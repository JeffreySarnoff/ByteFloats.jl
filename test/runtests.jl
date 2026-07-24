# ===== test/runtests.jl — consolidated exhaustive suite (design §9.1)
#
# Consolidates the seven development harnesses, verbatim in substance:
#   1. formats.jl vs an independent BigFloat transliteration of draft §4.7.2
#   2. projspec.jl construction/validation/printing + pipeline smoke
#   3. decode_encode.jl: encode round-trips, ordering three ways, Class, Next ops
#   4. ops_scalar.jl: registry, exhaustive Add/Subtract/Multiply vs references,
#      Annex A.3 division, stochastic R plumbing, Base register, Convert
#   5. oracle.jl: all unary ops exhaustive vs a 3072-bit run of the interval
#      protocol, monotonicity, draft pins, directed asymptotes, ternary exhaustive
#   6. tables.jl + kernels.jl: table ≡ scalar, cache semantics, vmap equivalence,
#      stochastic reproducibility, zero-allocation warm paths
#   7. blocks.jl + approx.jl + the Float128 revision plan's §7 gates
#      (carrier equivalence, differential Float128-vs-MPFR table builds,
#      width thresholds, CR-bracket soundness, envelope sanity, runtime switch): §5 composition vs a from-scratch reference,
#      ConvertToBlockMaxAbsFinite NOTEs, κ measurement/registration, conformance
#
# Assembled from the harnesses that each passed against the shipped sources
# (assertion counts in checkpoint.md); assembled file not re-executed at
# consolidation time per instruction — run with `Pkg.test("ByteFloats")`.

using Test
using Random
using ByteFloats
using Quadmath: Float128
using ByteFloats: project, project_interval, round_to_precision, encode, order_key,
    KIND_FIN, nan_code, posinf_code, neginf_code, signmask,
    apply_op, ωeval, BigExactF, EncloseF, StickyF, OP_REGISTRY, opinfo, _UNARY_OPS,
    get_table, blockdecode, _NAMED, Enclose128F, _USE_FLOAT128, _f128,
    _rtp_core, _rtp_f64, _extremal_SQ, _decode_compute, _decode_table, MaybeRNG,
    MaxFiniteOf, MinFiniteOf, MinPositiveOf, MaxSubnormalOf, MinNormalOf, formatname,
    nan_code, posinf_code, neginf_code, maxfinite_datum,
    get_table, _USE_FLOAT128, _f128, _UNARY_OPS, rawvalue, nan_code, TableKey,
    rawvalue, decode, _decode_compute, _decode_table, nan_code, Rounded, KIND_FIN,
    apply_op, MaybeRNG


const UN = collect(_UNARY_OPS)

# @testset "ByteFloats.jl" begin

# Independent reference decode straight from draft §4.7.2, in BigFloat
function refdecode(::Type{Binary{K,P,S,E}}, c::UInt8) where {K,P,S,E}
    S && c == UInt8(1 << (K-1)) && return :nan
    !S && c == UInt8((1<<K)-1) && return :nan
    if E
        S && c == UInt8((1<<(K-1))-1) && return :pinf
        S && c == UInt8((1<<K)-1) && return :ninf
        !S && c == UInt8((1<<K)-2) && return :pinf
    end
    neg = S && c >= UInt8(1<<(K-1)); m = neg ? c - UInt8(1<<(K-1)) : c
    T = Int(m) % (1 << (P-1)); Eb = Int(m) ÷ (1 << (P-1))
    B = S ? 1<<(K-P-1) : 1<<(K-P)
    X = Eb == 0 ? (0 + T*big(2.0)^(1-P)) * big(2.0)^(1-B) :
                  (1 + T*big(2.0)^(1-P)) * big(2.0)^(Eb-B)
    (neg ? -X : X)
end

nfmt = Ref(0); nchecked = Ref(0)
@testset "formats.jl vs draft" begin
    for K in 3:8, P in 1:K, S in (true,false), E in (true,false)
        S && P >= K && continue
        T = Binary{K,P,S,E}; nfmt[] += 1
        # decode agrees with independent reference on every code point + Float64 exactness
        for c in 0x00:UInt8((1<<K)-1)
            v = rawvalue(T, c); d = decode(v); r = refdecode(T, c)
            if r === :nan; @test isnan(d)
            elseif r === :pinf; @test d == Inf
            elseif r === :ninf; @test d == -Inf
            else
                @test big(d) == r            # exact carrier assertion
            end
            nchecked[] += 1
        end
        # extremal identities
        @test decode(MaxFiniteOf(T)) == maximum(filter(isfinite, [decode(rawvalue(T,c)) for c in 0x00:UInt8((1<<K)-1)]))
        @test decode(MinFiniteOf(T)) == minimum(filter(isfinite, [decode(rawvalue(T,c)) for c in 0x00:UInt8((1<<K)-1)]))
        @test decode(MinPositiveOf(T)) == minimum(filter(x -> isfinite(x) && x > 0, [decode(rawvalue(T,c)) for c in 0x00:UInt8((1<<K)-1)]))
        @test decode(MinNormalOf(T)) == big(2.0)^(1 - expbias(T))
        P > 1 && @test decode(MaxSubnormalOf(T)) == (big(2.0)^(1-P) * ((1<<(P-1))-1)) * big(2.0)^(1-expbias(T))
        # typemax/typemin/zero/one/eps
        @test decode(typemax(T)) == (E ? Inf : decode(MaxFiniteOf(T)))
        @test decode(typemin(T)) == (S ? (E ? -Inf : decode(MinFiniteOf(T))) : 0.0)
        @test iszero(zero(T)) && decode(one(T)) == 1.0
        @test decode(eps(T)) == 2.0^(1-P)
        # predicates on specials
        @test isnan(rawvalue(T, nan_code(T)))
        E && @test isinf(rawvalue(T, posinf_code(T)))
        (S && E) && @test decode(rawvalue(T, neginf_code(T))) == -Inf
        @test !signbit(zero(T)) && !signbit(rawvalue(T, nan_code(T)))
    end
end
println("formats verified: $(nfmt[]) formats, $(nchecked[]) code points, all exhaustive")

# Group M and the extremal queries also take a value; the value form must agree
# with the type form on every format and stay allocation-free.
@testset "Group M value forwarders" begin
    for K in 3:8, P in 1:K, S in (true,false), E in (true,false)
        S && P >= K && continue
        T = Binary{K,P,S,E}
        v = rawvalue(T, 0x01)
        for F in (BitwidthOf, PrecisionOf, SignednessOf, DomainOf, ExponentBiasOf,
                  ExponentBitwidthOf, TrailingSignificandBitwidthOf,
                  MaxFiniteOf, MinFiniteOf, MinPositiveOf, MaxSubnormalOf, MinNormalOf)
            @test F(v) === F(T)
        end
    end
    vb = Binary8p4se(1.5)
    fw(x) = BitwidthOf(x); fw(vb)
    @test @allocated(fw(vb)) == 0
    @test Base.return_types(fw, Tuple{Binary8p4se}) == [Int]
    mx(x) = MaxFiniteOf(x); mx(vb)
    @test @allocated(mx(vb)) == 0
    @test Base.return_types(mx, Tuple{Binary8p4se}) == [Binary8p4se]
end
println(Binary8p4se, "  ", Binary8p4se(2.0), "  ", formatname(Binary8p1uf))
@test Binary8p4se === Binary{8,4,true,true}
@test_throws ArgumentError Binary{9,4,true,true}(Val(:code), 0x00)
@test_throws ArgumentError Binary{8,8,true,true}(Val(:code), 0x00)
println("parameter validation OK")

# K=8 boundary: 0xff is a legitimate code point (−Inf for signed·extended)
@test decode(rawvalue(Binary8p4se, 0xff)) == -Inf

# ==========================================================================
# projspec.jl
# ==========================================================================

@testset "projspec" begin
    ρ = ProjSpec(NearestTiesToEven(), SatNone())
    @test ρ === RNE_SatNone && sizeof(ρ) == 0
    @test RoundOf(ρ) === NearestTiesToEven() && SatOf(ρ) === SatNone()
    @test !isstochastic(ρ) && nrandbits(ρ) == 0
    σ = ProjSpec(StochasticB{4}(), SatFinite())
    @test isstochastic(σ) && nrandbits(σ) == 4
    @test isstochastic(StochasticA{1}) && !isstochastic(ToOdd)
    @test_throws ArgumentError StochasticA{0}()
    @test_throws ArgumentError StochasticC{61}()
    @test_throws ArgumentError StochasticB{1.5}()
    @test projmode(RoundNearest) === NearestTiesToEven()
    @test projmode(RoundUp) === TowardPositive()
    @test projmode(RoundToZero) === TowardZero()
    @test projmode(ToOdd()) === ToOdd()
    @test default_projspec(Binary8p4se) === RNE_SatNone
    @test sprint(show, ρ) == "(NearestTiesToEven, SatNone)"
    @test sprint(show, σ) == "(StochasticB[4], SatFinite)"
    @test Base.issingletontype(typeof(RNE_SatFinite))
    T = Binary8p4se   # P=4 ⇒ ulp = 0.25 in [2,4)
    @test decode(project(T, RNE_SatNone, 2.1)) == 2.0
    @test decode(project(T, ProjSpec(TowardPositive(), SatNone()), 2.1)) == 2.25
    @test decode(project(T, RNE_SatNone, 2.125)) == 2.0          # exact tie → even
    @test decode(project(T, ProjSpec(NearestTiesToAway(), SatNone()), 2.125)) == 2.25
end

# ==========================================================================
# defaults.jl — session defaults and the projection/component coherence invariant
# ==========================================================================

@testset "defaults.jl" begin
    # initial values
    @test DefaultType() === Binary8p2se
    @test DefaultReturnType() === Binary8p2se
    @test DefaultAccumulatorType() === binary32 === Float32
    @test DefaultRoundingMode() === NearestTiesToEven()
    @test DefaultSaturationMode() === SatNone()
    @test DefaultProjection() === RNE_SatNone
    @test DefaultRNG() === Random.Xoshiro
    @test DefaultRbits() == 8

    coherent() = DefaultProjection() ===
                 ProjSpec(DefaultRoundingMode(), DefaultSaturationMode())
    @test coherent()

    # component setters update the projection
    @test DefaultRoundingMode!(TowardZero()) === TowardZero()
    @test DefaultProjection() === RTZ_SatNone
    @test coherent()
    @test DefaultSaturationMode!(SatFinite()) === SatFinite()
    @test DefaultProjection() === RTZ_SatFinite
    @test DefaultRoundingMode() === TowardZero()        # unchanged by the sat setter
    @test coherent()

    # type-argument convenience forms
    DefaultRoundingMode!(NearestTiesToAway)
    @test DefaultRoundingMode() === NearestTiesToAway()
    @test DefaultProjection() === RNA_SatFinite
    DefaultSaturationMode!(SatPropagate)
    @test DefaultProjection() === RNA_SatPropagate
    @test coherent()

    # direct projection setter updates both components
    DefaultProjection!(ProjSpec(StochasticA{8}(), SatNone()))
    @test DefaultRoundingMode() === StochasticA{8}()
    @test DefaultSaturationMode() === SatNone()
    @test coherent()
    DefaultProjection!(TowardPositive(), SatFinite)     # (mode, sat) convenience
    @test DefaultProjection() === RTP_SatFinite
    @test DefaultRoundingMode() === TowardPositive()
    @test DefaultSaturationMode() === SatFinite()
    @test coherent()

    # remaining setters, validation, and the RNG instance form
    @test DefaultType!(Binary5p3sf) === Binary5p3sf
    @test DefaultType() === Binary5p3sf
    @test_throws ArgumentError DefaultType!(Binary{8,8,true,true})   # invalid params
    @test DefaultReturnType!(Binary8p4se) === Binary8p4se
    @test DefaultReturnType() === Binary8p4se
    @test DefaultType() === Binary5p3sf                  # independent of the return type
    @test_throws ArgumentError DefaultReturnType!(Binary{8,8,true,true})
    @test DefaultAccumulatorType!(binary64) === Float64
    @test DefaultAccumulatorType() === binary64
    @test_throws MethodError DefaultAccumulatorType!(Int)   # not an AbstractFloat
    @test DefaultRbits!(16) == 16
    @test_throws ArgumentError DefaultRbits!(0)
    @test_throws ArgumentError DefaultRbits!(61)
    @test DefaultRbits() == 16                           # failed sets don't stick
    rng = Random.Xoshiro(42)
    @test DefaultRNG!(rng) === rng
    @test DefaultRNG() === rng

    # restore initial state for any later consumer
    DefaultType!(Binary8p2se)
    DefaultReturnType!(Binary8p2se)
    DefaultAccumulatorType!(binary32)
    DefaultProjection!(RNE_SatNone)
    DefaultRNG!(Random.Xoshiro)
    DefaultRbits!(8)
    @test coherent()

    # ---- consumption combinators: speculative fast path over a barrier
    mkval(T, x) = T(x)
    addρ(ρ, x, y) = Add(Binary8p4se, ρ, x, y)
    a4, b4 = Binary8p4se(1.5), Binary8p4se(0.25)

    # fast path (defaults at initial values): correct on every combinator
    @test with_default_type(mkval, 1.5) === Binary8p2se(1.5)
    @test with_default_returntype(mkval, 1.5) === Binary8p2se(1.5)
    @test with_default_accumulatortype(zero) === 0.0f0
    @test with_default_projection(addρ, a4, b4) === Add(Binary8p4se, RNE_SatNone, a4, b4)
    # Allocation contract (see defaults.jl): zero-alloc + concrete inference hold
    # when f's result type does not depend on the default — the projection
    # combinator's normal shape (caller fixes the formats, ρ steers rounding).
    # A result whose type IS the default (mkval above) boxes once at escape;
    # that box is irreducible for a runtime-chosen type, so it is not pinned.
    wdp(x, y) = with_default_projection(addρ, x, y); wdp(a4, b4)
    @test @allocated(wdp(a4, b4)) == 0
    @test Base.return_types(wdp, Tuple{Binary8p4se,Binary8p4se}) == [Binary8p4se]

    # slow path (defaults changed): same answers as passing the default explicitly
    DefaultType!(Binary6p3se)
    @test with_default_type(mkval, 1.5) === Binary6p3se(1.5)
    DefaultProjection!(RTZ_SatFinite)
    @test with_default_projection(addρ, a4, b4) === Add(Binary8p4se, RTZ_SatFinite, a4, b4)
    DefaultAccumulatorType!(binary64)
    @test with_default_accumulatortype(zero) === 0.0
    DefaultType!(Binary8p2se); DefaultProjection!(RNE_SatNone); DefaultAccumulatorType!(binary32)
    @test coherent()
end

# ==========================================================================
# decode_encode.jl
# ==========================================================================
allfmts = DataType[]
for K in 3:8, P in 1:K, S in (true,false), E in (true,false)
    (S && P >= K) || push!(allfmts, Binary{K,P,S,E})
end
npair = Ref(0)
@testset "decode_encode.jl §3" begin
for T in allfmts
    K = bitwidth(T); P = precision(T); B = expbias(T)
    codes = UInt8.(0:(1<<K)-1)
    vals  = [rawvalue(T, c) for c in codes]
    ds    = decode.(vals)

    # --- encode: round-trip every code through the canonical-form pipeline
    for (c, v, d) in zip(codes, vals, ds)
        if isfinite(d)
            r = round_to_precision(P, B, NearestTiesToEven(), d, 0, 0)   # exact datum ⇒ pure extraction
            @test r.kind == KIND_FIN
            @test encode(T, Int(r.sign), r.S, r.Q) == c
        end
        # ωProject of an exact datum is the identity on the code point (incl. NaN/±Inf)
        @test codepoint(project(T, RNE_SatNone, d)) == c
    end

    # --- order_key ⟺ TotalOrder ⟺ numeric-with-NaN-top, over all pairs
    for x in vals, y in vals
        npair[] += 1
        to = TotalOrder(x, y)
        @test to == (order_key(x) <= order_key(y))
        dx, dy = decode(x), decode(y)
        ref = isnan(dx) ? isnan(dy) : (isnan(dy) ? true : dx <= dy)
        @test to == ref
        # numeric comparisons: NaN unordered
        if isnan(dx) || isnan(dy)
            @test !(x == y) && !(x < y) && !(x <= y)
        else
            @test (x == y) == (dx == dy) && (x < y) == (dx < dy)
        end
    end

    # --- Class vs decode-derived reference
    for (v, d) in zip(vals, ds)
        want = isnan(d) ? ClassNaN : d == Inf ? ClassPosInf : d == -Inf ? ClassNegInf :
               d == 0 ? ClassZero :
               (abs(d) < 2.0^(1 - B) ? (d > 0 ? ClassPosSubnormal : ClassNegSubnormal) :
                                        (d > 0 ? ClassPosNormal   : ClassNegNormal))
        @test Class(v) == want
    end

    # --- Next ops vs sorted-enumeration reference:
    # NGT(x) = least datum > x, else NaN;  NLT(x) = greatest datum < x, else NaN;  NaN→NaN.
    finite_and_inf = sort([d for d in ds if !isnan(d)])
    for (v, d) in zip(vals, ds)
        g, l = NextGreaterThan(v), NextLessThan(v)
        if isnan(d)
            @test isnan(g) && isnan(l)
        else
            ups = filter(>(d), finite_and_inf)
            dns = filter(<(d), finite_and_inf)
            isempty(ups) ? (@test isnan(g)) : (@test decode(g) == first(ups))
            isempty(dns) ? (@test isnan(l)) : (@test decode(l) == last(dns))
        end
        @test nextfloat(v) === g && prevfloat(v) === l
    end
end
end
println("§3 verified over $(length(allfmts)) formats, $(npair[]) ordered pairs")

# ==========================================================================
# ops_scalar.jl
# ==========================================================================

# independent reference: exact BigFloat op → project
refbin(T, ρ, fop, x, y; R=0) = setprecision(BigFloat, 2400) do
    project(T, ρ, fop(BigFloat(decode(x)), BigFloat(decode(y))); R)
end

@testset "ops_scalar.jl §6" begin
    @test length(OP_REGISTRY) == 52                       # 31 unary(+Convert) + 18 binary + 3 ternary
    @test opinfo(:FMA).arity == 3 && opinfo(:Exp).arity == 1

    ρs = [RNE_SatNone, ProjSpec(NearestTiesToAway(), SatFinite()),
          ProjSpec(TowardPositive(), SatNone()), ProjSpec(TowardNegative(), SatPropagate()),
          ProjSpec(TowardZero(), SatNone()), ProjSpec(ToOdd(), SatNone())]

    # exhaustive Add / Subtract over K=5 (1024 pairs) × 6 ρ, vs BigFloat reference
    T = Binary5p2se
    vals = [rawvalue(T, UInt8(c)) for c in 0:31]
    for ρ in ρs, x in vals, y in vals
        dx, dy = decode(x), decode(y)
        want = if isnan(dx) || isnan(dy) || (isinf(dx) && isinf(dy) && dx != dy)
            rawvalue(T, nan_code(T))
        elseif isinf(dx) || isinf(dy)
            project(T, ρ, isinf(dx) ? dx : dy)
        else
            refbin(T, ρ, +, x, y)
        end
        @test codepoint(Add(T, ρ, x, y)) == codepoint(want)
        # Subtract vs reference (Inf algebra of x − y), plus the exact-negation identity,
        # which is a theorem only for finite y (projection may saturate Negate(±Inf)).
        wsub = if isnan(dx) || isnan(dy) || (isinf(dx) && isinf(dy) && dx == dy)
            rawvalue(T, nan_code(T))
        elseif isinf(dx) || isinf(dy)
            project(T, ρ, isinf(dx) ? dx : -dy)
        else
            refbin(T, ρ, -, x, y)
        end
        @test codepoint(Subtract(T, ρ, x, y)) == codepoint(wsub)
        if isfinite(dy)
            @test codepoint(Subtract(T, ρ, x, y)) == codepoint(Add(T, ρ, x, Negate(T, ρ, y)))
        end
    end

    # exhaustive Multiply on K=4 across mixed result format (Binary4p2se → Binary5p3se)
    S4 = Binary4p2se; R5 = Binary5p3se
    v4 = [rawvalue(S4, UInt8(c)) for c in 0:15]
    for ρ in ρs, x in v4, y in v4
        dx, dy = decode(x), decode(y)
        want = if isnan(dx) || isnan(dy) || (iszero(dx) && isinf(dy)) || (isinf(dx) && iszero(dy))
            rawvalue(R5, nan_code(R5))
        else
            refbin(R5, ρ, *, x, y)
        end
        @test codepoint(Multiply(R5, ρ, x, y)) == codepoint(want)
    end

    # Divide semantics (Annex A.3) + interval path
    z, o = zero(T), one(T)
    @test isnan(Divide(T, RNE_SatNone, o, z))                       # x/0 → NaN, all x
    @test isnan(Divide(T, RNE_SatNone, z, z))
    @test iszero(Divide(T, RNE_SatNone, o, rawvalue(T, posinf_code(T))))
    thr = T(3.0)
    @test decode(Divide(T, RNE_SatNone, o, thr)) ==
          decode(refbin(T, RNE_SatNone, /, o, thr))                 # 1/3 via enclosure

    # stochastic plumbing: explicit R deterministic; drawn R in range; rng honored
    σ = ProjSpec(StochasticA{2}(), SatNone())
    a, b = T(2.0), T(0.25)
    for r in 0:3
        @test Add(T, σ, a, b; R=r) === Add(T, σ, a, b; R=r)
    end
    @test_throws ArgumentError Add(T, σ, a, b; R=4)
    rng1, rng2 = Xoshiro(7), Xoshiro(7)
    @test Add(T, σ, a, b; rng=rng1) === Add(T, σ, a, b; rng=rng2)
    # stochastic R-sweep sanity: P=2 ⇒ ulp 1 in [2,4); 2+0.25 has ν=0.25,
    # StochasticA_{2,R} rounds away iff ⌊0.25·4⌋+R ≥ 4 ⇔ R=3: exactly 1 of 4 draws → 3.0
    @test [decode(Add(T, σ, a, b; R=r)) for r in 0:3] == [2.0, 2.0, 2.0, 3.0]

    # Base register ≡ spec register; mixed formats have no silent promotion
    @test codepoint(a + b) == codepoint(Add(T, RNE_SatNone, a, b))
    @test codepoint(exp(b)) == codepoint(Exp(T, RNE_SatNone, b))
    @test codepoint(-a) == codepoint(Negate(a))
    @test_throws Union{MethodError,ErrorException} a + one(Binary8p4se)

    # Convert: exact big Integer, Float32, BigFloat carrier, default constructor loop
    T8 = Binary8p4se
    n = Int64(2)^60 + 1
    @test decode(Convert(T8, RNE_SatNone, n)) ==
          setprecision(() -> decode(project(T8, RNE_SatNone, BigFloat(n))), BigFloat, 128)
    @test decode(Convert(T8, RNE_SatNone, Float32(2.1))) == decode(T8(Float64(Float32(2.1))))
    @test T8(2.0) === Convert(T8, RNE_SatNone, 2.0) && decode(one(T8)) == 1.0
    @test decode(eps(T8)) == 0.125

    # Tanh asymptote through the enclosure protocol (design §4.7):
    # tanh(maxfinite) = 1⁻ ⇒ RNE → 1, TowardNegative → greatest datum below 1
    mx = MaxFiniteOf(T8)
    @test decode(Tanh(T8, RNE_SatNone, mx)) == 1.0
    below1 = decode(NextLessThan(one(T8)))
    @test decode(Tanh(T8, ProjSpec(TowardNegative(), SatNone()), mx)) == below1
    @test decode(Tanh(T8, ProjSpec(TowardPositive(), SatNone()), Negate(mx))) == -below1
end
println("ops_scalar.jl verified")

# ==========================================================================
# oracle.jl
# ==========================================================================

# machinery cross-check: same defined result must emerge when the interval protocol
# starts at 3072 bits instead of 256 (catches escalation/precision bugs)
function highprec(T, ρ, res, R=0)
    res isa Float64 && return project(T, ρ, res; R)
    res isa Float128 && return project(T, ρ, res; R)
    res isa BigExactF && return project(T, ρ, res.f(); R)
    d, u = res.f(3072)
    isequal(d, u) && return project(T, ρ, d; R)
    cd = project(T, ρ, d; R, sticky=+1); cu = project(T, ρ, u; R, sticky=-1)
    codepoint(cd) == codepoint(cu) || error("unresolved at 3072 bits")
    cd
end

T8 = Binary8p3se
codes8 = [rawvalue(T8, UInt8(c)) for c in 0:255]
ρ4 = [RNE_SatNone, ProjSpec(TowardPositive(), SatNone()),
      ProjSpec(TowardNegative(), SatFinite()), ProjSpec(ToOdd(), SatNone())]

@testset "oracle.jl §5" begin
    # coverage: every registry op except Convert has ω-semantics
    for op in OP_REGISTRY
        op.name === :Convert && continue
        @test length(methods(ωeval, Tuple{Val{op.name}, Vararg{Float64, op.arity}})) >= 1
    end

    # exhaustive unary: all 30 ops × 256 codes × 4 ρ, vs the 3072-bit protocol
    @testset "unary exhaustive ×hp" begin
        for opn in UN, ρ in ρ4, v in codes8
            res = ωeval(Val(opn), decode(v))
            got = apply_op(Val(opn), T8, ρ, 0, decode(v))
            if codepoint(got) != codepoint(highprec(T8, ρ, res))
                println("T8=$T8,op=$opn ρ=$ρ v=$(decode(v)) got=$(decode(got)) want=$(decode(highprec(T8, ρ, res)))")
                error("exhaustive unary ×hp failed")
            end
            @test codepoint(got) == codepoint(highprec(T8, ρ, res))
        end
    end

    # monotone ops: RNE image must be nondecreasing over increasing finite inputs
    @testset "monotonicity" begin
        finv = sort(filter(v -> isfinite(decode(v)), codes8); by=decode)
        for opn in (:Exp, :Exp2, :ExpMinusOne, :Softplus, :Sinh, :Tanh, :ArcSinh,
                    :ArcTan, :ArcTanPi, :ArcSin, :ArcSinPi)
            prev = -Inf
            for v in finv
                d = decode(apply_op(Val(opn), T8, RNE_SatNone, 0, decode(v)))
                isnan(d) && continue
                @test d >= prev
                prev = d
            end
        end
        for opn in (:Log, :Log2, :LogOnePlus, :Sqrt, :ArcCosh)   # monotone on their domains
            prev = -Inf
            for v in finv
                dv = decode(v)
                d = decode(apply_op(Val(opn), T8, RNE_SatNone, 0, dv))
                isnan(d) && continue
                @test d >= prev
                prev = d
            end
        end
    end

    # semantic pins (draft rows / exact identities)
    @testset "pins" begin
        ρ0 = RNE_SatNone
        val(op, x, ρ=ρ0; R=0) = decode(apply_op(Val(op), T8, ρ, R, Float64(x)))
        val2(op, x, y, ρ=ρ0) = decode(apply_op(Val(op), T8, ρ, 0, Float64(x), Float64(y)))
        @test val(:Exp, 0) == 1 && val(:Log, 1) == 0
        @test val(:Log2, 4) == 2 && val(:Log2, 4, ProjSpec(TowardNegative(), SatNone())) == 2
        @test val(:Log, 0) == -Inf && isnan(val(:Log, -2)) && isnan(val(:Log, -Inf))
        @test val(:Sqrt, 4) == 2 && isnan(val(:Sqrt, -1)) && val(:Sqrt, Inf) == Inf
        @test val(:Recip, 4) == 0.25 && isnan(val(:Recip, 0)) && val(:Recip, Inf) == 0
        @test isnan(val(:RSqrt, 0)) && val(:RSqrt, 4) == 0.5 && val(:RSqrt, Inf) == 0
        @test val(:ExpMinusOne, -Inf) == -1 && val(:Softplus, -Inf) == 0
        @test val(:Tanh, Inf) == 1 && val(:ArcTanh, 1) == Inf && isnan(val(:ArcTanh, 1.5))
        @test isnan(val(:Sin, Inf)) && val(:Cos, 0) == 1
        @test val(:SinPi, 8) == 0 && val(:CosPi, 2.5) == 0 && val(:SinPi, 224) == 0  # exact mod-2, huge arg
        @test val(:TanPi, 0.5) == Inf && val(:TanPi, 1.5) == -Inf
        @test val(:TanPi, 0.25) == 1 && val(:TanPi, 0.75) == -1 && val(:TanPi, 5.25) == 1
        @test val(:TanPi, 0.25, ProjSpec(TowardNegative(), SatNone())) == 1   # exact even directed
        @test val2(:ArcTan2Pi, 3, 3) == 0.25 && val2(:ArcTan2Pi, 3, -3) == 0.75
        @test val2(:ArcTan2Pi, -3, 3) == -0.25 && val2(:ArcTan2Pi, -3, -3) == -0.75
        @test decode(apply_op(Val(:TanPi), Binary8p3sf, RNE_SatNone, 0, 0.5)) |> isnan  # Finite domain: ∞→NaN
        @test val(:ArcTanPi, 1) == 0.25 && val(:ArcTanPi, Inf) == 0.5
        @test val(:ArcCosPi, -1) == 1 && val(:ArcCosPi, 0) == 0.5
        @test val2(:Hypot, 3, 4) == 5 && val2(:Hypot, Inf, NaN) == Inf
        @test isnan(val2(:Divide, 1, 0)) && isnan(val2(:Divide, 0, 0)) && val2(:Divide, 1, Inf) == 0
        # ArcTan2 single-zero branch cuts: (0, x<0) → π-nearest; (y>0, 0) → π/2-nearest
        piT = decode(project(T8, ρ0, BigFloat(π)))
        @test val2(:ArcTan2, 0, -3) == piT
        @test val2(:ArcTan2, 0, 3) == 0 && val2(:ArcTan2, 0, 0) == 0
        halfpiT = decode(project(T8, ρ0, BigFloat(π) / 2))
        @test val2(:ArcTan2, 5, 0) == halfpiT
        @test val2(:ArcTan2Pi, 3, -Inf) == 1 && val2(:ArcTan2Pi, -3, 0) == -0.5
        # extremum family vectors (semantics, not formulas)
        @test isnan(val2(:Maximum, NaN, 3)) && val2(:MaximumNumber, NaN, 3) == 3
        @test val2(:MinimumMagnitude, -0.5, 4) == -0.5
        @test val2(:MaximumMagnitude, -4, 3) == -4
        @test val2(:MaximumMagnitudeNumber, NaN, -4) == -4
        @test val2(:MaximumFinite, Inf, 3) == 3          # finite preferred over ∞
        @test val2(:MaximumFinite, NaN, Inf) == Inf      # ∞ beats NaN
        @test isnan(val2(:MinimumFinite, NaN, NaN))
        @test val2(:CopySign, 3, -0.5) == -3 && isnan(val2(:CopySign, 3, NaN))
    end

    # directed asymptotes through the enclosure protocol
    @testset "asymptotes" begin
        mx = maxfinite_datum(T8); b1 = decode(NextLessThan(one(T8)))
        up = ProjSpec(TowardPositive(), SatNone()); dn = ProjSpec(TowardNegative(), SatNone())
        @test decode(apply_op(Val(:Tanh), T8, dn, 0, mx)) == b1
        @test decode(apply_op(Val(:Tanh), T8, RNE_SatNone, 0, mx)) == 1
        @test decode(apply_op(Val(:Exp), T8, up, 0, -mx)) == decode(MinPositiveOf(T8))
        @test decode(apply_op(Val(:Exp), T8, RNE_SatNone, 0, -mx)) == 0
        @test decode(apply_op(Val(:Exp), T8, ProjSpec(ToOdd(), SatNone()), 0, -mx)) == decode(MinPositiveOf(T8))
        @test decode(apply_op(Val(:ExpMinusOne), T8, up, 0, -mx)) == -b1
        @test decode(apply_op(Val(:Softplus), T8, up, 0, mx)) == Inf     # Mhi+ε rounds up, saturates
        @test decode(apply_op(Val(:Softplus), T8, RNE_SatNone, 0, mx)) == mx
        @test decode(apply_op(Val(:Softplus), T8, dn, 0, mx)) == mx
    end

    # exhaustive Divide over K=5 × ρ, vs high-precision protocol
    @testset "divide exhaustive" begin
        T5 = Binary5p2se
        v5 = [rawvalue(T5, UInt8(c)) for c in 0:31]
        for ρ in ρ4, x in v5, y in v5
            res = ωeval(Val(:Divide), decode(x), decode(y))
            @test codepoint(apply_op(Val(:Divide), T5, ρ, 0, decode(x), decode(y))) ==
                  codepoint(highprec(T5, ρ, res))
        end
    end

    # exhaustive ternary on K=4 vs exact BigFloat reference
    @testset "ternary exhaustive" begin
        T4 = Binary4p2se
        v4 = [decode(rawvalue(T4, UInt8(c))) for c in 0:15]
        for ρ in (RNE_SatNone, ProjSpec(TowardZero(), SatNone())), x in v4, y in v4, z in v4
            g_fma = decode(apply_op(Val(:FMA), T4, ρ, 0, x, y, z))
            w_fma = if isnan(x) || isnan(y) || isnan(z) ||
                       (iszero(x) && isinf(y)) || (isinf(x) && iszero(y))
                NaN
            else
                p = x * y
                if isinf(p) || isinf(z)
                    (isinf(p) && isinf(z) && p != z) ? NaN : (isinf(p) ? p : z)
                else
                    setprecision(() -> decode(project(T4, ρ, BigFloat(x) * BigFloat(y) + BigFloat(z))), BigFloat, 2400)
                end
            end
            @test isequal(g_fma, w_fma)
            g_faa = decode(apply_op(Val(:FAA), T4, ρ, 0, x, y, z))
            w_faa = if isnan(x) || isnan(y) || isnan(z)
                NaN
            elseif isinf(x) || isinf(y) || isinf(z)
                hp = any(==(Inf), (x, y, z)); hn = any(==(-Inf), (x, y, z))
                hp && hn ? NaN : (hp ? Inf : -Inf)
            else
                setprecision(() -> decode(project(T4, ρ, (BigFloat(x) + BigFloat(y)) + BigFloat(z))), BigFloat, 2400)
            end
            @test isequal(g_faa, w_faa)
            g_cl = decode(apply_op(Val(:Clamp), T4, ρ, 0, x, y, z))
            @test isequal(g_cl, (isnan(x) || isnan(y) || isnan(z)) ? NaN :
                                decode(project(T4, ρ, min(max(x, y), z))))
        end
    end

    # stochastic through the enclosure path: Exp under StochasticA{2}, full R-sweep,
    # vs an exact-ν decision computed independently at 3000 bits
    @testset "stochastic transcendental" begin
        σ = ProjSpec(StochasticA{2}(), SatNone())
        for v in codes8
            d = decode(v); (isfinite(d) && !iszero(d)) || continue
            for R in 0:3
                got = decode(Exp(T8, σ, v; R=R))
                res = ωeval(Val(:Exp), d)
                @test got == decode(highprec(T8, σ, res, R))
            end
        end
    end
end
println("oracle.jl verified")

# ==========================================================================
# tables+kernels
# ==========================================================================


T8 = Binary8p3se; T5 = Binary5p2se; S4 = Binary4p2se; R5 = Binary5p3se
ρup = ProjSpec(TowardPositive(), SatNone())

@testset "tables.jl + kernels.jl §7" begin
    empty_tables!()
    # --- table ≡ scalar path, exhaustively
    for (op, fr, f1, ρ) in ((:Exp, T8, T8, RNE_SatNone), (:Log, T8, T8, ρup),
                            (:Sqrt, T5, T5, RNE_SatNone), (:Convert, T5, T8, ρup))
        tbl = get_table(op, fr, f1, ρ)
        @test length(tbl) == 1 << bitwidth(f1)
        for c in 0:(1 << bitwidth(f1)) - 1
            want = op === :Convert ? codepoint(project(fr, ρ, decode(rawvalue(f1, UInt8(c))))) :
                                     codepoint(apply_op(Val(op), fr, ρ, 0, decode(rawvalue(f1, UInt8(c)))))
            @test tbl[c + 1] == want
        end
    end
    tbl2 = get_table(:Subtract, T5, T5, T5, RNE_SatNone)     # asymmetric op: catches index-order bugs
    @test length(tbl2) == 1 << 10
    for c1 in 0:31, c2 in 0:31
        want = codepoint(apply_op(Val(:Subtract), T5, RNE_SatNone, 0,
                                  decode(rawvalue(T5, UInt8(c1))), decode(rawvalue(T5, UInt8(c2)))))
        @test tbl2[(c1 << 5) + c2 + 1] == want
    end
    # mixed formats: Multiply Binary4p2se × Binary4p2se → Binary5p3se
    tblm = get_table(:Multiply, R5, S4, S4, RNE_SatNone)
    @test length(tblm) == 256

    # --- cache identity, byte accounting, reset
    @test get_table(:Exp, T8, T8, RNE_SatNone) === get_table(:Exp, T8, T8, RNE_SatNone)
    @test table_bytes() == 256 + 256 + 32 + 256 + 1024 + 256
    empty_tables!(); @test table_bytes() == 0
    @test_throws ArgumentError get_table(:Exp, T8, T8, ProjSpec(StochasticA{2}(), SatNone()))

    # --- vmap ≡ scalar map (unary, binary, mixed-format, asymmetric)
    codesA = [rawvalue(T5, UInt8(rand(0:31))) for _ in 1:4096]
    codesB = [rawvalue(T5, UInt8(rand(0:31))) for _ in 1:4096]
    out = Subtract(T5, RNE_SatNone, codesA, codesB)
    @test all(codepoint(out[i]) == codepoint(Subtract(T5, RNE_SatNone, codesA[i], codesB[i])) for i in eachindex(out))
    a4 = [rawvalue(S4, UInt8(rand(0:15))) for _ in 1:1000]
    b4 = [rawvalue(S4, UInt8(rand(0:15))) for _ in 1:1000]
    om = Multiply(R5, ρup, a4, b4)
    @test eltype(om) == R5
    @test all(codepoint(om[i]) == codepoint(Multiply(R5, ρup, a4[i], b4[i])) for i in eachindex(om))
    ou = Exp(T8, RNE_SatNone, [rawvalue(T8, UInt8(c)) for c in 0:255])
    @test all(codepoint(ou[c + 1]) == codepoint(Exp(T8, RNE_SatNone, rawvalue(T8, UInt8(c)))) for c in 0:255)
    oc = Convert(T5, ρup, [rawvalue(T8, UInt8(c)) for c in 0:255])
    @test all(codepoint(oc[c + 1]) == codepoint(Convert(T5, ρup, rawvalue(T8, UInt8(c)))) for c in 0:255)

    # --- ternary Shape B
    c4 = [rawvalue(S4, UInt8(rand(0:15))) for _ in 1:1000]
    of = FMA(S4, RNE_SatNone, a4, b4, c4)
    @test all(codepoint(of[i]) == codepoint(FMA(S4, RNE_SatNone, a4[i], b4[i], c4[i])) for i in eachindex(of))

    # --- stochastic arrays: reproducible under the same rng; matches a manual
    #     scalar loop consuming the identical draw sequence
    σ = ProjSpec(StochasticA{3}(), SatNone())
    o1 = Add(T5, σ, codesA, codesB; rng=Xoshiro(42))
    o2 = Add(T5, σ, codesA, codesB; rng=Xoshiro(42))
    @test all(codepoint.(o1) .== codepoint.(o2))
    rng = Xoshiro(42)
    o3 = [Add(T5, σ, codesA[i], codesB[i]; rng=rng) for i in eachindex(codesA)]
    @test all(codepoint.(o1) .== codepoint.(o3))
    @test any(codepoint.(o1) .!= codepoint.(Add(T5, σ, codesA, codesB; rng=Xoshiro(7))))  # R matters

    # --- warm-path allocation: gather loops allocate nothing beyond the output
    dest = similar(codesA); v = Val(:Add)
    vmap!(dest, v, T5, RNE_SatNone, codesA, codesB)
    @test (@allocated vmap!(dest, v, T5, RNE_SatNone, codesA, codesB)) == 0
    destu = similar(ou); vu = Val(:Exp); srcu = [rawvalue(T8, UInt8(rand(0:255))) for _ in 1:4096]
    du = similar(srcu)
    vmap!(du, vu, T8, RNE_SatNone, srcu)
    @test (@allocated vmap!(du, vu, T8, RNE_SatNone, srcu)) == 0

    # --- views and strides (AbstractArray contract)
    V = view(codesA, 100:2:900)
    ov = Exp(T5, RNE_SatNone, V)
    @test all(codepoint(ov[i]) == codepoint(Exp(T5, RNE_SatNone, V[i])) for i in eachindex(V))

    # --- throughput smoke (informational, not asserted)
    n = 1 << 18   # informational throughput smoke, kept small for CI
    bigA = [rawvalue(T8, UInt8(rand(0:255))) for _ in 1:n]
    bigB = [rawvalue(T8, UInt8(rand(0:255))) for _ in 1:n]
    bd = similar(bigA)
    vmap!(bd, Val(:Add), T8, RNE_SatNone, bigA, bigB)                      # warm + build 64 KiB table
    t = @elapsed vmap!(bd, Val(:Add), T8, RNE_SatNone, bigA, bigB)
    println("Shape-A 8×8 Add gather: ", round(3n / t / 1e9; digits=2), " GB/s effective (", n, " elems)")
    tu = @elapsed vmap!(du, Val(:Exp), T8, RNE_SatNone, srcu)
    println("Shape-A unary gather:   ", round(2 * length(srcu) / tu / 1e9; digits=2), " GB/s effective")
end
println("tables.jl + kernels.jl verified")

# ==========================================================================
# blocks.jl
# ==========================================================================

# ---------- independent reference: draft §5 composition, written from scratch ----------
# opv(prec, rd) must return the op's true value directed-rounded at prec (or ±Inf/NaN specials
# resolved to exact Float64 by the caller). Interval division and the sticky protocol below are
# deliberately a *separate* implementation from blocks.jl's.
function ref_elem(fr, ρ, R, opv, Sr::Float64)
    isnan(Sr) && return nan_code(fr)
    v0 = opv(64, RoundNearest)
    isnan(v0) && return nan_code(fr)
    iszero(Sr) && return codepoint(project(fr, ρ, 0.0; R))
    if isinf(Sr)
        sg = Float64(sign(v0)) * sign(Sr)
        return codepoint(project(fr, ρ, sg; R))
    end
    if isinf(v0)
        return codepoint(project(fr, ρ, Float64(sign(Sr)) * Float64(v0); R))
    end
    prec = 3072
    vd = opv(prec, RoundDown); vu = opv(prec, RoundUp)
    qd, qu = setprecision(BigFloat, prec) do
        if Sr > 0
            (setrounding(() -> vd / Sr, BigFloat, RoundDown), setrounding(() -> vu / Sr, BigFloat, RoundUp))
        else
            (setrounding(() -> vu / Sr, BigFloat, RoundDown), setrounding(() -> vd / Sr, BigFloat, RoundUp))
        end
    end
    isequal(qd, qu) && return codepoint(project(fr, ρ, qd; R))
    cd = project(fr, ρ, qd; R, sticky=+1); cu = project(fr, ρ, qu; R, sticky=-1)
    codepoint(cd) == codepoint(cu) || error("reference unresolved")
    codepoint(cd)
end
# decoded lane value with ωMultiply specials (independent expression)
function ref_lane(S::Float64, x::Float64)
    (isnan(S) || isnan(x)) && return NaN
    ((iszero(S) && isinf(x)) || (isinf(S) && iszero(x))) && return NaN
    p = S * x
    iszero(p) ? 0.0 : p              # ByteFloats single zero: reference must normalize −0.0 too
end

T5 = Binary5p2se; T8 = Binary8p3se; U1 = Binary8p1uf
rnd(T) = rawvalue(T, UInt8(rand(0:(1 << bitwidth(T)) - 1)))
ρ2 = (RNE_SatNone, ProjSpec(TowardPositive(), SatNone()))

@testset "blocks.jl §8" begin
    # --- blockdecode ≡ independent lane semantics, incl 0·∞ and NaN
    for _ in 1:200
        b = Block(rnd(T5), ntuple(_ -> rnd(T5), 3))
        X = blockdecode(b)
        for i in 1:3
            @test isequal(X[i], ref_lane(decode(b.s), decode(b.x[i])))
        end
    end
    binf = Block(rawvalue(T5, posinf_code(T5)), (zero(T5), one(T5)))
    @test isnan(blockdecode(binf)[1]) && blockdecode(binf)[2] == Inf

    # --- elementwise BlockOp vs the reference composition, random + special-heavy
    Random.seed!(1)
    for B in (1, 2, 3, 16, 31, 32), ρ in ρ2, trial in 1:(B <= 3 ? 60 : 12)
        FS = trial % 3 == 0 ? U1 : T5
        b1 = Block(rnd(FS), ntuple(_ -> rnd(T5), B))
        b2 = Block(rnd(FS), ntuple(_ -> rnd(T5), B))
        sr = rnd(FS)
        # BlockAdd (exact op values)
        got = BlockAdd(T8, ρ, b1, b2, sr)
        @test got.s === sr
        X1 = blockdecode(b1); X2 = blockdecode(b2)
        for i in 1:B
            x, y = X1[i], X2[i]
            opv = if isnan(x) || isnan(y) || (isinf(x) && isinf(y) && x != y)
                (p, rd) -> BigFloat(NaN)
            elseif isinf(x) || isinf(y)
                (p, rd) -> BigFloat(isinf(x) ? x : y)
            else
                (p, rd) -> setprecision(() -> BigFloat(x) + BigFloat(y), BigFloat, 2400)
            end
            @test codepoint(got.x[i]) == ref_elem(T8, ρ, 0, opv, decode(sr))
        end
        # BlockExp (enclosure op values)
        gexp = BlockExp(T8, ρ, b1, sr)
        for i in 1:B
            x = X1[i]
            opv = if isnan(x); (p, rd) -> BigFloat(NaN)
            elseif x == -Inf; (p, rd) -> BigFloat(0)
            elseif x == Inf;  (p, rd) -> BigFloat(Inf)
            else (p, rd) -> setprecision(() -> setrounding(() -> exp(BigFloat(x)), BigFloat, rd), BigFloat, p)
            end
            @test codepoint(gexp.x[i]) == ref_elem(T8, ρ, 0, opv, decode(sr))
        end
    end

    # --- BlockProject S-special rows (draft §5.1.2 NOTEs 1–2)
    nanT5 = rawvalue(T5, nan_code(T5)); infT5 = rawvalue(T5, posinf_code(T5))
    xs = (one(T5), Negate(one(T5)), zero(T5), nanT5)
    bz = ConvertToBlock(T5, T8, RNE_SatNone, xs, zero(T5))
    @test decode.(bz.x) === (0.0, 0.0, 0.0) .* 1 || all(i -> i == 4 ? isnan(decode(bz.x[i])) : decode(bz.x[i]) == 0.0, 1:4)
    bi = ConvertToBlock(T5, T8, RNE_SatNone, xs, infT5)
    @test decode(bi.x[1]) == 1.0 && decode(bi.x[2]) == -1.0 && decode(bi.x[3]) == 0.0 && isnan(decode(bi.x[4]))

    # --- ScaledOp ≡ B=1 BlockOp with unit result scale (draft §5.5)
    for _ in 1:100
        s1, x1, s2, x2 = rnd(T5), rnd(T5), rnd(T5), rnd(T5)
        r1 = ScaledAdd(T8, RNE_SatNone, s1, x1, s2, x2)
        r2 = BlockAdd(T8, RNE_SatNone, Block(s1, (x1,)), Block(s2, (x2,)), one(T5))
        @test codepoint(r1) == codepoint(r2.x[1])
        r3 = ScaledDivide(T8, RNE_SatNone, s1, x1, s2, x2)
        r4 = BlockDivide(T8, RNE_SatNone, Block(s1, (x1,)), Block(s2, (x2,)), one(T5))
        @test codepoint(r3) == codepoint(r4.x[1])
    end

    # --- reductions vs references
    Random.seed!(2)
    for B in (1, 3, 16, 32), trial in 1:40
        b = Block(rnd(T5), ntuple(_ -> rnd(T5), B))
        X = blockdecode(b)
        # ReduceAdd
        want = if any(isnan, X); NaN
        elseif any(==(Inf), X) && any(==(-Inf), X); NaN
        elseif any(isinf, X); X[findfirst(isinf, X)]
        else setprecision(() -> decode(project(T8, RNE_SatNone, sum(BigFloat, X; init=BigFloat(0)))), BigFloat, 3000)
        end
        @test isequal(decode(BlockReduceAdd(T8, RNE_SatNone, b)), want)
        # DotProduct incl >53-bit lane products
        by = Block(rnd(T5), ntuple(_ -> rnd(T5), B))
        Y = blockdecode(by)
        lanes = [ref_lane(1.0, 1.0)]; ok = true
        wdot = begin
            cls = [ (isnan(X[i]) || isnan(Y[i]) || (iszero(X[i]) && isinf(Y[i])) || (isinf(X[i]) && iszero(Y[i]))) ? NaN :
                    (isinf(X[i]) || isinf(Y[i])) ? sign(X[i]) * sign(Y[i]) * Inf : 1.0 for i in 1:B ]
            if any(isnan, cls); NaN
            elseif any(isinf, cls)
                (any(==(Inf), cls) && any(==(-Inf), cls)) ? NaN : cls[findfirst(isinf, cls)]
            else
                setprecision(BigFloat, 3000) do
                    decode(project(T8, RNE_SatNone, sum(BigFloat(X[i]) * BigFloat(Y[i]) for i in 1:B; init=BigFloat(0))))
                end
            end
        end
        @test isequal(decode(BlockDotProduct(T8, RNE_SatNone, b, by)), wdot)
    end
    # 64-bit lane-product stress: full-significand scales × full-significand elements
    smax = MaxFiniteOf(T5)
    bx = Block(smax, ntuple(_ -> MaxFiniteOf(T5), 4))
    by = Block(MinPositiveOf(T5), ntuple(_ -> MinPositiveOf(T5), 4))
    wexact = setprecision(BigFloat, 3000) do
        decode(project(T8, RNE_SatNone, 4 * BigFloat(decode(smax))^2 * BigFloat(decode(MinPositiveOf(T5)))^2))
    end
    @test decode(BlockDotProduct(T8, RNE_SatNone, bx, by)) == wexact
    # ReduceMultiply pins: 0 present with ∞ → NaN; sign of ∞ product; plain product
    b0i = Block(one(T5), (zero(T5), rawvalue(T5, posinf_code(T5))))
    @test isnan(decode(BlockReduceMultiply(T8, RNE_SatNone, b0i)))
    bni = Block(one(T5), (Negate(one(T5)), rawvalue(T5, posinf_code(T5))))
    @test decode(BlockReduceMultiply(T8, RNE_SatNone, bni)) == -Inf
    bp = Block(T5(2.0), (T5(3.0), T5(2.0)))
    @test decode(BlockReduceMultiply(T8, RNE_SatNone, bp)) == 24.0

    # --- ConvertToBlockMaxAbsFinite: the five draft NOTEs
    ρs_up = ProjSpec(TowardPositive(), SatNone())
    allnan = ntuple(_ -> nanT5, 3)
    r = ConvertToBlockMaxAbsFinite(T5, T8, RNE_SatNone, RNE_SatNone, allnan)
    @test isnan(decode(r.s)) && all(isnan ∘ decode, r.x)                       # NOTE 1
    allinf = ntuple(_ -> infT5, 3)
    r = ConvertToBlockMaxAbsFinite(T5, T8, RNE_SatNone, RNE_SatNone, allinf)
    @test decode(r.s) == Inf && all(v -> decode(v) == 1.0, r.x)                # NOTE 2 (SatNone: s=∞, elems ±1)
    mixed = (infT5, T5(2.0), Negate(infT5))
    r = ConvertToBlockMaxAbsFinite(T5, T8, RNE_SatNone, RNE_SatNone, mixed)
    @test decode(r.s) == 2.0                                                   # NOTE 3: ∞ doesn't set scale
    @test decode(r.x[1]) == Inf && decode(r.x[3]) == -Inf && decode(r.x[2]) == 1.0
    tiny = ntuple(_ -> MinPositiveOf(T5), 3)                                    # NOTE 4 shape: scale→0 ⇒ all zero
    rz = ConvertToBlockMaxAbsFinite(Binary3p1se, T8, ProjSpec(TowardZero(), SatFinite()), RNE_SatNone, tiny)
    if iszero(decode(rz.s)); @test all(v -> decode(v) == 0.0, rz.x); end
    r5 = ConvertToBlockMaxAbsFinite(U1, T8, ρs_up, RNE_SatFinite, (T5(3.0), T5(0.5), Negate(T5(2.0))))
    @test decode(r5.s) == 4.0                                                  # NOTE 5: TowardPositive P=1 scale
    @test decode(r5.x[1]) == 0.75 && decode(r5.x[2]) == 0.125 && decode(r5.x[3]) == -0.5

    # --- P=1 scale exactness: division by 2^k collapses to the exact path
    bU = Block(U1(4.0), (T5(3.0), T5(0.5)))
    g = BlockAdd(T8, RNE_SatNone, bU, Block(U1(1.0), (zero(T5), zero(T5))), U1(2.0))
    @test decode(g.x[1]) == 6.0 && decode(g.x[2]) == 1.0                       # (4·3+0)/2, (4·0.5+0)/2

    # --- stochastic block ops: reproducible per seeded rng
    σ = ProjSpec(StochasticA{2}(), SatNone())
    b1 = Block(T5(2.0), (T5(1.5), T5(3.0))); b2 = Block(one(T5), (T5(0.25), T5(0.25)))
    o1 = BlockAdd(T5, σ, b1, b2, one(T5); rng=Xoshiro(9))
    o2 = BlockAdd(T5, σ, b1, b2, one(T5); rng=Xoshiro(9))
    @test codepoint.(o1.x) == codepoint.(o2.x)

    # --- BlockVector SoA round trip and layout
    blocks = [Block(rnd(T5), ntuple(_ -> rnd(T5), 4)) for _ in 1:10]
    bv = BlockVector(blocks)
    @test length(bv) == 10 && all(codepoint.(bv[j].x) == codepoint.(blocks[j].x) && bv[j].s === blocks[j].s for j in 1:10)
    bv[3] = blocks[7]
    @test codepoint.(bv[3].x) == codepoint.(blocks[7].x)
    @test size(bv.elems) == (4, 10)
end
println("blocks.jl verified")

# ==========================================================================
# approx.jl
# ==========================================================================

T8 = Binary8p4se; T5 = Binary5p2se
@testset "approx.jl §9" begin
    # --- κ of the exact path is 0 (exhaustive, unary and binary)
    exact_exp(x) = Exp(T8, RNE_SatNone, x)
    @test measure_kappa(exact_exp, :Exp, T8, (T8,), RNE_SatNone) === (0.0, true)
    exact_add(x, y) = Add(T5, RNE_SatNone, x, y)
    @test measure_kappa(exact_add, :Add, T5, (T5, T5), RNE_SatNone) === (0.0, true)
    # ternary exhaustive at K=4 (2^12 inputs)
    T4 = Binary4p2se
    exact_fma(x, y, z) = FMA(T4, RNE_SatNone, x, y, z)
    @test measure_kappa(exact_fma, :FMA, T4, (T4, T4, T4), RNE_SatNone) === (0.0, true)

    # --- synthetic known-κ implementations: perturb the defined result by k steps
    step2(x) = (r = Exp(T8, RNE_SatNone, x);
                isfinite(decode(r)) ? NextGreaterThan(NextGreaterThan(r)) : r)
    κ2, exh = measure_kappa(step2, :Exp, T8, (T8,), RNE_SatNone)
    @test exh
    @test κ2 >= 2.0     # ≥: stepping can also cross into a region ≥2 keys away
    # NaN mismatch: return a number where the defined result is NaN
    denan(x) = (r = Log(T8, RNE_SatNone, x); isnan(decode(r)) ? zero(T8) : r)
    @test measure_kappa(denan, :Log, T8, (T8,), RNE_SatNone)[1] |> isnan
    # finite→Inf deviation is also NaN-κ
    toinf(x) = rawvalue(T8, posinf_code(T8))
    @test measure_kappa(toinf, :Abs, T8, (T8,), RNE_SatNone)[1] |> isnan

    # --- FTZ worked example (draft Annex): flush subnormal results, ties to zero
    ρf = RNE_SatFinite
    ftz = ftz_variant(:Exp, T8, T8, ρf)
    # behavior: every subnormal defined result maps to 0 or ±MinNormal, nearest, ties→0
    P = precision(T8); half = 1 << (P - 2); mn = decode(MinNormalOf(T8))
    for c in 0:255
        x = rawvalue(T8, UInt8(c))
        want = Exp(T8, ρf, x); got = ftz(x)
        if issubnormal(want)
            m = Int(codepoint(want) & ~signmask(T8))
            @test decode(got) == (m <= half ? 0.0 : mn)   # exp ≥ 0: no negative results
        else
            @test codepoint(got) == codepoint(want)
        end
    end
    κf, exhf = measure_kappa(ftz, :Exp, T8, (T8,), ρf)
    @test exhf && κf == 1 << (P - 2)                       # largest flushed-to-zero subnormal
    # P=1 target: no subnormals, FTZ is exact
    U1 = Binary8p1uf
    @test measure_kappa(ftz_variant(:Convert, U1, T8, ρf), :Convert, U1, (T8,), ρf) === (0.0, true)

    # --- registration semantics
    foreach(unregister_approx!, list_approx())
    impl = register_approx!(:exp_ftz_8p4, :Exp, T8, (T8,), ρf, ftz)   # κ auto-declared
    @test kappa(:exp_ftz_8p4) == κf && kappa_measured(impl) == κf && impl.exhaustive
    @test :exp_ftz_8p4 in list_approx()
    @test approx(:exp_ftz_8p4).fn === ftz
    @test_throws ArgumentError register_approx!(:exp_ftz_8p4, :Exp, T8, (T8,), ρf, ftz)  # duplicate
    @test_throws ArgumentError register_approx!(:lie, :Exp, T8, (T8,), ρf, ftz; κ=1)      # understated
    ok = register_approx!(:generous, :Exp, T8, (T8,), ρf, ftz; κ=10)                       # overstated OK
    @test kappa(ok) == 10 && kappa_measured(ok) == κf
    @test_throws ArgumentError register_approx!(:nanimpl, :Log, T8, (T8,), RNE_SatNone, denan)  # NaN needs κ=NaN
    nreg = register_approx!(:nanimpl, :Log, T8, (T8,), RNE_SatNone, denan; κ=NaN)
    @test isnan(kappa(nreg))
    @test_throws ArgumentError register_approx!(:badop, :Nope, T8, (T8,), ρf, ftz)
    @test_throws ArgumentError measure_kappa(ftz, :Exp, T8, (T8,), ProjSpec(StochasticA{2}(), SatNone()))
    unregister_approx!(:generous)
    @test !(:generous in list_approx())

    # --- conformance declaration reflects registry, cache, and approx state
    empty_tables!()
    get_table(:Exp, T8, T8, RNE_SatNone)
    get_table(:Add, T5, T5, T5, RNE_SatNone)
    c = conformance()
    @test length(c.formats) == 120 && :Binary8p4se in c.formats
    @test length(c.operations) == 52
    @test count(o -> o.arity == 3, c.operations) == 3
    @test length(c.cached_specializations) == 2
    @test :BlockDotProduct in c.block_surface && :ScaledFMA in c.block_surface && :BlockTanh in c.block_surface
    @test length(c.block_surface) == 51 * 2 + 6
    @test any(a -> a.name === :exp_ftz_8p4 && a.kappa == κf && a.exhaustive, c.approximate)
    d = conformance_dict(c)
    @test d["package"] == "ByteFloats.jl 0.1.0" && length(d["formats"]) == 120
    @test any(s -> s["op"] == "Add" && s["saturation"] == "SatNone", d["cached_specializations"])
    buf = IOBuffer(); conformance_report(buf, c); rep = String(take!(buf))
    @test occursin("κ verified exhaustively", rep) && occursin("Exp⟨", rep)
    @test occursin("Scalar operations (52)", rep)
end
println("approx.jl verified") 

@testset "Float128 revision plan §7" begin
    # --- carrier equivalence: project(Float64 d) ≡ project(Float128(d)),
    #     exhaustive over all datums of representative formats × ρ × R-sweep
    ρs = [RNE_SatNone, ProjSpec(TowardPositive(), SatNone()),
          ProjSpec(TowardNegative(), SatFinite()), ProjSpec(ToOdd(), SatNone()),
          ProjSpec(NearestTiesToAway(), SatPropagate())]
    for T in (Binary8p4se, Binary8p1uf, Binary5p2se, Binary3p1se), ρ in ρs
        for c in 0x00:UInt8((1 << bitwidth(T)) - 1)
            d = decode(rawvalue(T, c))
            @test codepoint(project(T, ρ, d)) == codepoint(project(T, ρ, Float128(d)))
        end
    end
    σ = ProjSpec(StochasticA{3}(), SatNone())
    for c in 0x00:0xff, R in (0, 3, 7)
        d = decode(rawvalue(Binary8p3se, c))
        @test codepoint(project(Binary8p3se, σ, d; R)) ==
              codepoint(project(Binary8p3se, σ, Float128(d); R))
    end
    # Float128 inputs must not double-round through Float64 at a target midpoint.
    x128 = Float128(0.53125) + ldexp(one(Float128), -101)
    @test Float64(x128) == 0.53125
    @test decode(Convert(Binary8p4se, RNE_SatNone, x128)) == 0.5625

    # The interval ladder must honor a non-power-of-two precision ceiling exactly.
    seen_precisions = Int[]
    unresolved = p -> (push!(seen_precisions, p); (BigFloat(0.7), BigFloat(0.8)))
    @test_throws ErrorException project_interval(Binary3p2se, RNE_SatNone, unresolved; maxprec=300)
    @test seen_precisions == [256, 300]
    @test_throws ArgumentError project_interval(Binary3p2se, RNE_SatNone, unresolved; maxprec=1)

    # sticky semantics on the Float128 carrier (the asymptote machinery)
    @test decode(project(Binary8p4se, ProjSpec(TowardNegative(), SatNone()),
                         Float128(1); sticky=-1)) == decode(NextLessThan(one(Binary8p4se)))
    @test decode(project(Binary8p4se, ProjSpec(TowardPositive(), SatNone()),
                         Float128(0); sticky=+1)) == decode(MinPositiveOf(Binary8p4se))

    # --- differential builds: Float128-first vs forced-MPFR, byte-identical
    function build_both(op, fr, args...)
        empty_tables!(); _USE_FLOAT128[] = true
        t1 = copy(get_table(op, fr, args...))
        empty_tables!(); _USE_FLOAT128[] = false
        t0 = copy(get_table(op, fr, args...))
        _USE_FLOAT128[] = true; empty_tables!()
        (t1, t0)
    end
    ρd = (RNE_SatNone, ProjSpec(TowardPositive(), SatNone()),
          ProjSpec(TowardNegative(), SatFinite()), ProjSpec(ToOdd(), SatNone()))
    for opn in _UNARY_OPS, ρ in ρd, T in (Binary8p3se, Binary8p4se)
        t1, t0 = build_both(opn, T, T, ρ)
        @test t1 == t0
    end
    for ρ in ρd
        for (fr, f1, f2) in ((Binary8p4se, Binary8p4se, Binary8p4se),
                             (Binary8p1uf, Binary8p1uf, Binary8p1uf))   # wide-spread stress
            t1, t0 = build_both(:Add, fr, f1, f2, ρ)
            @test t1 == t0
            t1, t0 = build_both(:Divide, fr, f1, f2, ρ)
            @test t1 == t0
        end
    end

    # --- width thresholds: escalation picks Float128 exactly within the band
    e40 = 2.0^40
    @test ωeval(Val(:Add), 2.0^100, 1.0) isa Float128
    @test ωeval(Val(:Add), 2.0^101, 1.0) isa BigExactF
    @test ωeval(Val(:Add), e40, 1.0) isa Float64              # within Float64: no escalation
    # and the Float128 band is *exact*: compare against the BigFloat truth
    for ΔE in (60, 80, 100)
        x, y = 2.0^ΔE, 1.0 + 0.5^7
        r = ωeval(Val(:Add), x, y)
        @test r isa Float128
        @test setprecision(() -> BigFloat(r) == BigFloat(x) + BigFloat(y), BigFloat, 300)
    end
    @test ωeval(Val(:FMA), 2.0^46, 2.0^46, 1.0) isa Float128       # ΔE(p,z) = 92
    let r = ωeval(Val(:FMA), 2.0^47, 2.0^46, 1.0)                  # 93 > 92: sticky head
        @test r isa StickyF{Float64} && r.v == 2.0^93 && r.sgn == 1
    end
    @test ωeval(Val(:FAA), 2.0^98, 1.0, 1.0) isa Float128
    let r = ωeval(Val(:FAA), 2.0^99, 1.0, 1.0)                     # 99 > 98: distilled, exact
        @test r isa Float128
        @test setprecision(() -> BigFloat(r) == BigFloat(2.0^99) + 2, BigFloat, 300)
    end
    let r = ωeval(Val(:FAA), 2.0^99, 2.0^-99, 2.0^-99)             # inexact tail → sticky
        @test r isa StickyF{Float128} && Float64(r.v) == 2.0^99 && r.sgn == 1
    end

    # --- CR enclosure soundness: inexact Divide now returns the Float64-CR fast
    #     enclosure — the eager estimate is the IEEE correctly-rounded quotient
    #     (within half an ulp of truth) and the MPFR ladder strictly brackets it
    for (x, y) in ((1.0, 3.0), (2.0^60, 3.0), (7.0, -11.0), (1.0, 2.0^-60 * 3))
        r = ωeval(Val(:Divide), x, y)
        @test r isa EncloseF
        t = setprecision(() -> BigFloat(x) / BigFloat(y), BigFloat, 500)
        @test r.yd == x / y                                   # the IEEE-CR quotient
        @test abs(t - BigFloat(r.yd)) <= 0.5 * eps(abs(r.yd)) # CR ⇒ ≤ half an ulp
        lo, hi = r.f(256)
        @test lo < t < hi                                     # ladder strictly brackets
    end
    s = ωeval(Val(:Sqrt), 2.0)
    @test s isa EncloseF && s.yd == sqrt(2.0)                 # IEEE-CR hardware sqrt
    ts = setprecision(() -> sqrt(BigFloat(2)), BigFloat, 500)
    @test abs(ts - BigFloat(s.yd)) <= 0.5 * eps(s.yd)         # CR ⇒ ≤ half an ulp
    slo, shi = s.f(256)
    @test slo < ts < shi                                      # ladder strictly brackets
    rs = ωeval(Val(:RSqrt), 2.0)
    trs = setprecision(() -> 1 / sqrt(BigFloat(2)), BigFloat, 500)
    @test rs isa EncloseF && rs.yd == 1.0 / sqrt(2.0)         # composed hardware-CR estimate
    @test abs(trs - BigFloat(rs.yd)) <= 1.5 * eps(rs.yd)      # ≤ ~1.5 ulp by composition
    rlo, rhi = rs.f(256)
    @test rlo < trs < rhi                                     # ladder strictly brackets

    # --- envelope sanity: libquadmath error ≪ the 2^-90 claim (assert < 2^-100)
    rng = Xoshiro(2026)
    for f in (exp, log, sin, tanh, atan, log1p, expm1, asinh)
        for _ in 1:200
            x = (f === log ? rand(rng) * 100 + 1e-6 :
                 f === log1p ? rand(rng) * 10 - 0.99 : randn(rng) * 3)
            y128 = f(Float128(x))
            ymp = setprecision(() -> f(BigFloat(x)), BigFloat, 256)
            iszero(ymp) && continue
            @test abs(BigFloat(y128) - ymp) <= abs(ymp) * big(2.0)^-100
        end
    end

    # --- runtime switch semantics: identical results either way
    a, b = Binary8p1uf(2.0^60), Binary8p1uf(2.0^-60)
    r_on = Add(Binary8p1uf, RNE_SatNone, a, b)
    _USE_FLOAT128[] = false
    r_off = Add(Binary8p1uf, RNE_SatNone, a, b)
    _USE_FLOAT128[] = true
    @test codepoint(r_on) == codepoint(r_off)
end

# ==========================================================================
# Bit-operations revision plan gates (Phase 0, K1–K5)
# ==========================================================================
modes_bo = [NearestTiesToEven(), NearestTiesToAway(), TowardPositive(), TowardNegative(), TowardZero(), ToOdd()]

@testset "bitops plan gates" begin
    # ---- K2: table ≡ compute ≡ bit-composed decode, all 13,296 code points
    for T in allfmts, c in 0x00:UInt8((1 << bitwidth(T)) - 1)
        v = rawvalue(T, c)
        @test isequal(decode(v), _decode_compute(v))
        @test isequal(_decode_table(T)[Int(c) + 1], _decode_compute(v))
    end
    # constant folding preserved (Group M identities still exact)
    @test ByteFloats.maxfinite_datum(Binary8p4se) == 224.0   # 1.75·2^7 (B = 8, Eb = 15)

    # ---- K4: _extremal_SQ ≡ Base.decompose of the extremal datum, all formats
    for T in allfmts
        S, Q = _extremal_SQ(T)
        (mS, mE, _) = Base.decompose(ByteFloats.maxfinite_datum(T))
        # normalize decompose's trailing zeros into the comparison
        @test S * big(2.0)^Q == mS * big(2.0)^mE
    end

    # ---- K3: bit rounding ≡ generic core, exhaustively over reachable exact values
    function eqrtp(P, B, μ, X, R, st)
        a = _rtp_f64(P, B, μ, X, R, st)
        b = _rtp_core(P, B, μ, X, R, st)
        a.kind == b.kind && a.sign == b.sign && a.S == b.S && a.Q == b.Q
    end
    Random.seed!(11)
    for T in (Binary8p4se, Binary8p1uf, Binary5p2se, Binary3p1se, Binary8p8uf)
        P = precision(T); B = ByteFloats.expbias(T)
        dats = [decode(rawvalue(T, UInt8(c))) for c in 0:(1 << bitwidth(T)) - 1]
        fin = filter(x -> isfinite(x) && !iszero(x), dats)
        # datums, pairwise sums/products (the reachable arithmetic values), sticky variants
        pool = Float64[]
        append!(pool, fin)
        for _ in 1:4000
            x, y = rand(fin), rand(fin)
            s = x + y; iszero(s) || !isfinite(s) || push!(pool, s)
            p = x * y; iszero(p) || !isfinite(p) || push!(pool, p)
        end
        for μ in modes_bo, X in pool, st in (-1, 0, 1)
            @test eqrtp(P, B, μ, X, 0, st)
        end
        # stochastic: full R-sweep at small N, boundary N ∈ {45, 60} sampled R
        for X in pool[1:200]
            for N in (1, 2), R in 0:(1 << N) - 1, SV in (StochasticA{N}(), StochasticB{N}(), StochasticC{N}())
                @test eqrtp(P, B, SV, X, R, 0)
            end
            for (N, Rs) in ((45, (0, 1, (1 << 45) - 1)), (60, (0, 7, (1 << 60) - 1)))
                for R in Rs, SV in (StochasticA{N}(), StochasticC{N}())
                    @test eqrtp(P, B, SV, X, R, 0)
                end
            end
        end
    end
    # subnormal-Float64 inputs route to the generic core and still project correctly
    @test decode(project(Binary8p4se, RNE_SatNone, 5.0e-324)) == 0.0
    @test decode(project(Binary8p4se, ProjSpec(TowardPositive(), SatNone()), 5.0e-324)) ==
          decode(MinPositiveOf(Binary8p4se))

    # ---- K1: comparisons and sorting
    for T in (Binary8p4se, Binary5p2se, Binary8p1uf)
        vals = [rawvalue(T, UInt8(c)) for c in 0:(1 << bitwidth(T)) - 1]
        for x in vals, y in vals
            dx, dy = decode(x), decode(y)
            refto = isnan(dx) ? isnan(dy) : (isnan(dy) ? true : dx <= dy)
            @test TotalOrder(x, y) == refto
            if isnan(dx) || isnan(dy)
                @test !(x == y) && !(x < y) && !(x <= y)
            else
                @test (x == y) == (dx == dy) && (x < y) == (dx < dy) && (x <= y) == (dx <= dy)
            end
        end
        A = [rawvalue(T, rand(UInt8) & UInt8((1 << bitwidth(T)) - 1)) for _ in 1:5000]
        s1 = sort(A)                                   # counting sort via defalg
        s2 = sort(A; alg=Base.Sort.DEFAULT_UNSTABLE)   # stock comparison sort
        @test codepoint.(s1) == codepoint.(s2)
        r1 = sort(A; rev=true)
        r2 = sort(A; alg=Base.Sort.DEFAULT_UNSTABLE, rev=true)
        @test codepoint.(r1) == codepoint.(r2)
        @test issorted(s1)
        sv = sort!(view(copy(A), 100:4000))
        @test issorted(sv)
        byp = sort(A; by=decode)                       # non-default ordering falls back safely
        @test decode.(byp[.!isnan.(decode.(byp))]) |> issorted
    end
    nanheavy = [Binary8p4se(1.0), rawvalue(Binary8p4se, nan_code(Binary8p4se)), Binary8p4se(0.5)]
    @test isnan(decode(sort(nanheavy)[end]))           # NaN last, Base convention
    @test isnan(decode(sort(nanheavy; rev=true)[1]))   # NaN first under rev

    # ---- Phase 0: dispatch/kwarg semantics preserved incl explicit rng
    σ = ProjSpec(StochasticA{2}(), SatNone())
    a, b = Binary5p2se(2.0), Binary5p2se(0.25)
    @test Add(Binary5p2se, σ, a, b; rng=Xoshiro(1)) === Add(Binary5p2se, σ, a, b; rng=Xoshiro(1))
    @test Add(Binary5p2se, σ, a, b; R=3) === Add(Binary5p2se, σ, a, b; R=3)
    @test Add(Binary5p2se, RNE_SatNone, a, b) === Add(Binary5p2se, RNE_SatNone, a, b)

    # ---- K5: packed storage
    for T in (Binary3p1se, Binary4p2se, Binary5p2se, Binary6p3se, Binary7p3se, Binary8p4se)
        for n in 0:66
            A = [rawvalue(T, UInt8(rand(0:(1 << bitwidth(T)) - 1))) for _ in 1:n]
            pv = PackedVector(A)
            @test length(pv) == n && codepoint.(collect(pv)) == codepoint.(A)
            for i in 1:n                                # setindex round trip
                pv[i] = A[n - i + 1]
            end
            @test codepoint.(collect(pv)) == codepoint.(reverse(A))
        end
        A = [rawvalue(T, UInt8(rand(0:(1 << bitwidth(T)) - 1))) for _ in 1:1000]
        pv = PackedVector(A)
        @test sizeof(pv.data) <= cld(1000 * bitwidth(T), 64) * 8 + 8
        o1 = ByteFloats.vmap(:Exp, T, RNE_SatNone, pv)       # kernel-through-packed ≡ bytes
        o2 = ByteFloats.vmap(:Exp, T, RNE_SatNone, A)
        @test codepoint.(o1) == codepoint.(o2)
    end
end

# ==========================================================================
# Specialization regressions (deterministic replacements for the resolved
# measurement flags — see checkpoint.md; concrete return types and zero
# warm-path allocation are the properties whose absence mimics slowness)
# ==========================================================================
@testset "specialization regressions" begin
    T = Binary8p4se
    a, b, c = T(1.5), T(0.25), T(2.0)
    σ = ProjSpec(StochasticA{3}(), SatNone())
    # concrete inferred return types at public entry points, kwarg paths included
    @test Base.return_types(project, Tuple{Type{T}, typeof(RNE_SatNone), Float64}) == [T]
    for f in (Add, Multiply, Divide, Exp, FMA)
        n = f === FMA ? 3 : (f in (Add, Multiply, Divide) ? 2 : 1)
        sig = Tuple{Type{T}, typeof(RNE_SatNone), ntuple(_ -> T, n)...}
        @test Base.return_types(f, sig) == [T]
    end
    @test Base.return_types((x, y) -> Add(T, σ, x, y; rng=Xoshiro(1), R=nothing), Tuple{T,T}) == [T]
    # zero allocation on warm public scalar paths (pure ρ) and the engine
    add2(x, y) = Add(T, RNE_SatNone, x, y); add2(a, b)
    prj(d) = project(T, RNE_SatNone, d); prj(2.3)
    fma3(x, y, z) = FMA(T, RNE_SatNone, x, y, z); fma3(a, b, c)
    cmp2(x, y) = x < y; cmp2(a, b)
    @test @allocated(add2(a, b)) == 0
    @test @allocated(prj(2.3)) == 0
    @test @allocated(fma3(a, b, c)) == 0
    @test @allocated(cmp2(a, b)) == 0
    # stochastic scalar with explicit R: deterministic and allocation-free
    addR(x, y) = Add(T, σ, x, y; R=5); addR(a, b)
    @test @allocated(addR(a, b)) == 0
    @test addR(a, b) === addR(a, b)
end

# ==========================================================================
# UInt8 code-point constructor
# ==========================================================================
@testset "UInt8 code-point constructor" begin
    # the motivating example, exactly
    a = Binary5p3sf(1.0)
    b = Binary5p3sf(0x08)
    @test a == b && a === b
    # exhaustive: T(c) ≡ rawvalue(T, c) over every format and every valid code
    for nm in keys(_NAMED)
        T = getfield(ByteFloats, nm)
        for c in 0x00:UInt8((1 << bitwidth(T)) - 1)
            @test T(c) === rawvalue(T, c)
        end
        # invalid codes throw for K < 8 (K = 8 accepts all of UInt8)
        if bitwidth(T) < 8
            @test_throws ArgumentError T(UInt8(1 << bitwidth(T)))
            @test_throws ArgumentError T(0xff)
        end
    end
    # semantic split: UInt8 is a code point; every other Integer is a numeric value
    @test Binary8p4se(0x02) === rawvalue(Binary8p4se, 0x02)
    @test decode(Binary8p4se(2)) == 2.0
    @test Binary8p4se(0x02) != Binary8p4se(2)
    @test decode(Binary8p4se(UInt16(2))) == 2.0        # only UInt8, not wider unsigned
    @test decode(Convert(Binary8p4se, RNE_SatNone, 0x02)) == 2.0   # Convert stays numeric
    # Rational disambiguation policy
    @test_throws ArgumentError Binary8p4se(1//2)
    # method-table hygiene and specialization
    @test isempty(Test.detect_ambiguities(ByteFloats))
    mk(c) = Binary5p3sf(c); mk(0x08)
    @test @allocated(mk(0x08)) == 0
    @test Base.return_types(mk, Tuple{UInt8}) == [Binary5p3sf]
end

# ==========================================================================
# Bitwidth-specific FMA/FAA performance paths (ternary tables, adaptive cache,
# sticky-head escalation, threaded compute)
# ==========================================================================
include("ternary_opt.jl")

# ==========================================================================
# Package hygiene (Aqua) and static error analysis (JET)
# ==========================================================================
include("quality.jl")

# end # @testset "ByteFloats.jl"
