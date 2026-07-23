# ===== test/ternary_opt.jl — bitwidth-specific FMA/FAA performance paths
#
# Gates for the ternary optimization layer:
#   1. ternary tables are definitionally bit-identical to the scalar path
#      (exhaustive over full ternary cubes, same-format and mixed-format)
#   2. the vmap ternary gather (eager band) matches the scalar path element-wise
#   3. the adaptive K=7 band builds only after TERNARY_BUILD_ELEMS and then serves
#      bit-identical results; the LRU byte budget holds and evicts
#   4. the sticky-head wide-spread escalation (StickyF) is differentially equal
#      to the MPFR reference (_bigfma/_bigsum3 + _finish) across every rounding
#      mode family including stochastic sub-grids, on random and adversarial
#      (exact-cancellation, tie-boundary, tiny-tail) operand triples
#   5. warm wide-spread scalar FMA/FAA allocate zero bytes
#   6. the threaded compute loop is element-identical to the sequential loop

using Test
using Random
using ByteFloats
using ByteFloats: apply_op, decode, rawvalue, get_table, _finish, _bigfma, _bigsum3,
    _faa_wide, _fma_wide, StickyF, TernaryKey, TERNARY_CACHE, TERNARY_USE,
    TERNARY_EAGER_BITS, TERNARY_ADAPTIVE_BITS, TERNARY_BUILD_ELEMS, TERNARY_CACHE_BYTES,
    THREAD_MIN_ELEMS, THREADED_KERNELS, empty_tables!

@testset "ternary tables ≡ scalar (exhaustive)" begin
    empty_tables!()
    ρs = (RNE_SatNone, RNE_SatFinite, ProjSpec(TowardZero(), SatNone()),
          ProjSpec(ToOdd(), SatFinite()))
    for T in (Binary3p1se, Binary4p2se), op in (:FMA, :FAA, :Clamp), ρ in ρs
        K = bitwidth(T)
        tbl = get_table(op, T, T, T, T, ρ)
        @test length(tbl) == 1 << (3K)
        for c1 in 0:(1<<K)-1, c2 in 0:(1<<K)-1, c3 in 0:(1<<K)-1
            idx = ((c1 << K | c2) << K) + c3 + 1
            r = apply_op(Val(op), T, ρ, 0,
                         decode(rawvalue(T, UInt8(c1))), decode(rawvalue(T, UInt8(c2))),
                         decode(rawvalue(T, UInt8(c3))))
            @test tbl[idx] == codepoint(r)
        end
    end
    # mixed formats: distinct result and operand formats
    let fr = Binary5p2ue, f1 = Binary4p2se, f2 = Binary3p1sf, f3 = Binary5p3se, ρ = RNE_SatNone
        tbl = get_table(:FMA, fr, f1, f2, f3, ρ)
        K2, K3 = bitwidth(f2), bitwidth(f3)
        @test length(tbl) == 1 << (bitwidth(f1) + K2 + K3)
        for c1 in 0:(1<<bitwidth(f1))-1, c2 in 0:(1<<K2)-1, c3 in 0:(1<<K3)-1
            idx = ((c1 << K2 | c2) << K3) + c3 + 1
            r = apply_op(Val(:FMA), fr, ρ, 0,
                         decode(rawvalue(f1, UInt8(c1))), decode(rawvalue(f2, UInt8(c2))),
                         decode(rawvalue(f3, UInt8(c3))))
            @test tbl[idx] == codepoint(r)
        end
    end
    # stochastic ρ never tabulable
    @test_throws ArgumentError get_table(:FMA, Binary4p2se, Binary4p2se, Binary4p2se,
                                         Binary4p2se, ProjSpec(StochasticA{4}(), SatNone()))
    empty_tables!()
end

@testset "vmap ternary gather ≡ scalar" begin
    empty_tables!()
    rng = Xoshiro(7)
    for T in (Binary4p2se, Binary6p3se), op in (:FMA, :FAA)
        K = bitwidth(T)
        A = [rawvalue(T, UInt8(rand(rng, 0:(1<<K)-1))) for _ in 1:500]
        B = [rawvalue(T, UInt8(rand(rng, 0:(1<<K)-1))) for _ in 1:500]
        C = [rawvalue(T, UInt8(rand(rng, 0:(1<<K)-1))) for _ in 1:500]
        D = vmap(op, T, RNE_SatNone, A, B, C)
        @test all(i -> D[i] === apply_op(Val(op), T, RNE_SatNone, 0,
                                         decode(A[i]), decode(B[i]), decode(C[i])),
                  eachindex(D))
    end
    # the eager band actually built tables (3·4 = 12 and 3·6 = 18 bits are both ≤ 18)
    @test !isempty(TERNARY_CACHE)
    # stochastic ternary vmap still runs (scalar path) and reproduces under a fixed rng
    let T = Binary4p2se, ρ = ProjSpec(StochasticC{8}(), SatNone())
        A = fill(T(1.5), 100); B = fill(T(1.25), 100); C = fill(T(0.125), 100)
        r1 = vmap(:FMA, T, ρ, A, B, C; rng=Xoshiro(3))
        r2 = vmap(:FMA, T, ρ, A, B, C; rng=Xoshiro(3))
        @test all(i -> r1[i] === r2[i], eachindex(r1))
    end
    empty_tables!()
end

@testset "adaptive K=7 band + LRU budget" begin
    empty_tables!()
    old_elems = TERNARY_BUILD_ELEMS[]
    TERNARY_BUILD_ELEMS[] = 2500
    try
        V7 = Binary7p3se                                  # 21 bits: adaptive band
        rng = Xoshiro(11)
        A7 = [rawvalue(V7, UInt8(rand(rng, 0:127))) for _ in 1:1000]
        FAA(V7, RNE_SatNone, A7, A7, A7)
        @test isempty(TERNARY_CACHE)                      # 1000 elems: below threshold
        FAA(V7, RNE_SatNone, A7, A7, A7)
        @test isempty(TERNARY_CACHE)                      # 2000: still below
        r = FAA(V7, RNE_SatNone, A7, A7, A7)              # 3000 ≥ 2500: builds
        @test length(TERNARY_CACHE) == 1
        @test table_bytes() == 1 << 21
        @test all(i -> r[i] === FAA(V7, RNE_SatNone, A7[i], A7[i], A7[i]), eachindex(r))
    finally
        TERNARY_BUILD_ELEMS[] = old_elems
    end
    # 24-bit signatures (all-K=8) never tabulate
    empty_tables!()
    W = Binary8p3se
    A8 = fill(W(1.5), 64)
    FMA(W, RNE_SatNone, A8, A8, A8)
    @test isempty(TERNARY_CACHE) && isempty(TERNARY_USE)
    # LRU: budget for three 4 KiB tables, insert six → three survive, budget holds
    old_bytes = TERNARY_CACHE_BYTES[]
    TERNARY_CACHE_BYTES[] = 3 * (1 << 12)
    try
        T = Binary4p2se
        for op in (Val(:FMA), Val(:FAA), Val(:Clamp)), S in (Binary4p2se, Binary4p3sf)
            vmap!(similar([T(1.0)], T, 4), op, T, RNE_SatNone,
                  fill(S(1.0), 4), fill(S(1.0), 4), fill(S(1.0), 4))
        end
        @test length(TERNARY_CACHE) == 3
        @test table_bytes() <= TERNARY_CACHE_BYTES[]
    finally
        TERNARY_CACHE_BYTES[] = old_bytes
    end
    empty_tables!()
end

@testset "sticky-head escalation ≡ MPFR reference" begin
    ρall = Any[RNE_SatNone, RNE_SatFinite, ProjSpec(TowardZero(), SatNone()),
               ProjSpec(TowardPositive(), SatFinite()), ProjSpec(TowardNegative(), SatNone()),
               ProjSpec(ToOdd(), SatFinite()), ProjSpec(NearestTiesToAway(), SatPropagate()),
               ProjSpec(StochasticA{4}(), SatNone()), ProjSpec(StochasticB{8}(), SatFinite()),
               ProjSpec(StochasticC{60}(), SatNone())]
    Rsweep(ρ) = (N = nrandbits(ρ); N == 0 ? (0,) : (0, 1, 1 << (N - 1), (1 << N) - 1))
    # random triples from the widest formats (the only ones that reach escalation)
    rng = Xoshiro(1)
    for W in (Binary8p1ue, Binary8p1se, Binary7p1ue, Binary8p2ue)
        vals = [decode(rawvalue(W, c)) for c in 0x00:UInt8((1 << bitwidth(W)) - 1)
                if isfinite(decode(rawvalue(W, c)))]
        for _ in 1:250
            x, y, z = rand(rng, vals), rand(rng, vals), rand(rng, vals)
            for ρ in ρall, R in Rsweep(ρ), (op, ref) in ((:FMA, _bigfma), (:FAA, _bigsum3))
                @test apply_op(Val(op), W, ρ, R, x, y, z) === _finish(W, ρ, R, ref(x, y, z))
            end
        end
    end
    # adversarial: exact cancellation, near-cancellation, tie-boundary heads, tiny tails
    W = Binary8p1se
    advs = Tuple{Float64,Float64,Float64}[]
    for e1 in (-100, -20, 0, 40, 120), e2 in (-127, -90, -30, 5)
        b = 2.0^e1; t = 2.0^e2
        append!(advs, [(b, -b, t), (b, -b, -t), (b, -b/2, t), (-b, b/2, t),
                       (b, t, -b), (t, b, -b), (b, b, t), (b, 3b, -t),
                       (1.5b, -b, t), (b, -0.75b, -t)])
    end
    for (x, y, z) in advs, ρ in ρall, R in Rsweep(ρ)
        for (op, ref) in ((:FMA, _bigfma), (:FAA, _bigsum3))
            @test apply_op(Val(op), W, ρ, R, x, y, z) === _finish(W, ρ, R, ref(x, y, z))
        end
    end
    # _faa_wide is exact-or-sticky, never silently wrong on Σ = 0
    @test _faa_wide(2.0^100, -(2.0^100), 0.0) == 0.0
    let s = _faa_wide(2.0^100, -(2.0^100), 2.0^-100)
        v = s isa StickyF ? Float64(s.v) : Float64(s)
        @test v == 2.0^-100
    end
end

@testset "wide-spread warm path: zero allocation" begin
    W = Binary8p1se
    a, b, c = W(2.0^60), W(2.0^-100), W(2.0^-100)
    fmaw(x, y, z) = FMA(W, RNE_SatNone, x, y, z)
    faaw(x, y, z) = FAA(W, RNE_SatNone, x, y, z)
    fmaw(a, b, c); faaw(a, b, c)                          # warm
    @test @allocated(fmaw(a, b, c)) == 0
    @test @allocated(faaw(a, b, c)) == 0
    @test Base.return_types(fmaw, NTuple{3,W}) == [W]
    @test Base.return_types(faaw, NTuple{3,W}) == [W]
end

@testset "threaded ternary compute ≡ sequential" begin
    W = Binary8p2se                                        # 24 bits: compute path
    rng = Xoshiro(2)
    A = [rawvalue(W, rand(rng, UInt8)) for _ in 1:50_000]
    B = shuffle(rng, A); C = shuffle(rng, A)
    old = THREAD_MIN_ELEMS[]
    try
        THREAD_MIN_ELEMS[] = 1 << 10
        Dt = FMA(W, RNE_SatNone, A, B, C)
        THREAD_MIN_ELEMS[] = typemax(Int)
        Ds = FMA(W, RNE_SatNone, A, B, C)
        @test all(i -> Dt[i] === Ds[i], eachindex(Dt))     # === : NaN codes compare too
    finally
        THREAD_MIN_ELEMS[] = old
    end
end
