# test_faa_64_32.jl — validation for Float64faa.jl and Float32faa.jl
#
# Oracles:
#   1. Exact value in high-precision BigFloat (4096 bits for Float64, 512 for
#      Float32 — both exceed the worst-case exact-sum width), with the
#      round-to-nearest / ties-to-even property checked directly against both
#      neighbours. Platform independent.
#   2. Native sequential (x + y) + z for special-value inputs with no NaN
#      (no rounding occurs, so fused == sequential bit-for-bit); the
#      deterministic NaN rule is checked directly against its specification.
#   3. Boldo-Melquiond round-to-odd construction from native hardware ops —
#      an independent correctly-rounded sum — mass-compared in a safe
#      exponent range; disagreements adjudicated by oracle 1.
#   4. Cross-format diagnostic: faa32(x,y,z) vs Float32(faa64(x,y,z)).
#      These are expected to deviate on genuine double-rounding cases (a
#      three-addend sum can land exactly on a Float32 tie point after the
#      Float64 rounding); every deviation must be adjudicated in faa32's
#      favor by oracle 1.
#
# Run:  julia test_faa_64_32.jl

using Test, Random, Printf

include("Float64faa.jl")
include("Float32faa.jl")
const faa64 = Float64FAA.faa
const faa32 = Float32FAA.faa

bits(x::Float64) = reinterpret(UInt64, x)
bits(x::Float32) = reinterpret(UInt32, x)
implicit(::Type{Float64}) = UInt64(1) << 52
implicit(::Type{Float32}) = UInt32(1) << 23
bigprec(::Type{Float64}) = 4096
bigprec(::Type{Float32}) = 512

"Round-to-nearest-even check of r against the exact x + y + z."
function nearest_ok(x::T, y::T, z::T, r::T) where {T<:Union{Float32,Float64}}
    (isfinite(x) && isfinite(y) && isfinite(z)) || return true
    setprecision(BigFloat, bigprec(T)) do
        exact = BigFloat(x) + BigFloat(y) + BigFloat(z)
        if !isfinite(r)
            isnan(r) && return false
            lim = BigFloat(floatmax(T)) +
                  (BigFloat(floatmax(T)) - BigFloat(prevfloat(floatmax(T)))) / 2
            return abs(exact) >= lim && (signbit(r) == (exact < 0))
        end
        if iszero(r)
            return iszero(exact) ? !signbit(r) || (signbit(x) && signbit(y) && signbit(z)) :
                   abs(exact) <= BigFloat(nextfloat(zero(T))) / 2
        end
        err = abs(exact - BigFloat(r))
        for nb in (prevfloat(r), nextfloat(r))
            isfinite(nb) || continue
            e2 = abs(exact - BigFloat(nb))
            e2 < err && return false
            if e2 == err
                trailing_zeros(bits(r) | implicit(T)) >=
                    trailing_zeros(bits(nb) | implicit(T)) || return false
            end
        end
        return true
    end
end

# --- oracle 3: Boldo-Melquiond via hardware ops ----------------------------
@inline function twosum(a::T, b::T) where {T}
    s = a + b
    ap = s - b; bp = s - ap
    s, (a - ap) + (b - bp)
end
@inline function roundodd_add(u::T, v::T) where {T}
    s, t = twosum(u, v)
    if t != zero(T) && (bits(s) & 1) == 0
        s = t > zero(T) ? nextfloat(s) : prevfloat(s)
    end
    s
end
function sum3_bm(x::T, y::T, z::T) where {T}
    uh, ul = twosum(x, y)
    th, tl = twosum(uh, z)
    th + roundodd_add(tl, ul)
end

rng = MersenneTwister(0x64_32)

expbits(::Type{Float64}) = 11; expbits(::Type{Float32}) = 8
fracbits(::Type{Float64}) = 52; fracbits(::Type{Float32}) = 23
uof(::Type{Float64}) = UInt64; uof(::Type{Float32}) = UInt32
function rand_finite(::Type{T}, emin::Int, emax::Int; possign=false) where {T}
    U = uof(T)
    ef = U(rand(rng, emin:emax))
    fr = rand(rng, U) & ((U(1) << fracbits(T)) - one(U))
    sn = (!possign && rand(rng, Bool)) ? U(1) << (8 * sizeof(U) - 1) : U(0)
    reinterpret(T, sn | (ef << fracbits(T)) | fr)
end
efmax(::Type{Float64}) = 2046; efmax(::Type{Float32}) = 254
function rand_nan(::Type{T}) where {T}
    U = uof(T)
    p = rand(rng, U) & ((U(1) << fracbits(T)) - 1)
    p == U(0) && (p = one(U))
    sn = rand(rng, Bool) ? U(1) << (8 * sizeof(U) - 1) : U(0)
    reinterpret(T, sn | (U(2^expbits(T) - 1) << fracbits(T)) | p)
end

function run_battery(::Type{T}, faaT) where {T}
    U = uof(T)
    EFM = efmax(T)
    check_big(x, y, z) = nearest_ok(x, y, z, faaT(x, y, z))
    p = fracbits(T) + 1                      # significand width

    @testset "$T: BigFloat oracle, broad random" begin
        ok = true
        for _ in 1:40_000
            ok &= check_big(rand_finite(T, 1, EFM), rand_finite(T, 1, EFM),
                            rand_finite(T, 1, EFM))
        end
        @test ok
    end

    @testset "$T: exponent-gap sweep (incl. carry region)" begin
        gset = T === Float64 ? (0,1,2,20,52,53,54,73,74,75,76,120,127,128,200,2000) :
                               (0,1,2,10,23,24,25,38,39,40,41,60,63,64,90,250)
        g2set = T === Float64 ? (0,1,53,74,75,128,200) : (0,1,24,39,40,64,90)
        ok = true
        for g1 in gset, g2 in g2set
            for _ in 1:80
                e1 = rand(rng, (EFM ÷ 2):(EFM ÷ 2 + 100))
                x = rand_finite(T, e1, e1)
                y = rand_finite(T, max(1, e1 - g1), max(1, e1 - g1))
                z = rand_finite(T, max(1, e1 - g1 - g2), max(1, e1 - g1 - g2))
                ok &= check_big(x, y, z)
                # all-positive variant hammers the accumulator-carry path
                ok &= check_big(abs(x), abs(y), abs(z))
            end
        end
        @test ok
    end

    @testset "$T: cancellation & return-c path" begin
        ok = true
        for _ in 1:30_000
            x = rand_finite(T, EFM ÷ 4, 3 * (EFM ÷ 4))
            z = rand_finite(T, 1, EFM)
            ok &= check_big(x, -x, z)
            ok &= bits(faaT(x, -x, z)) == bits(z)
            ok &= bits(faaT(z, x, -x)) == bits(z)
            y = rand_finite(T, EFM ÷ 4, 3 * (EFM ÷ 4))
            s = x + y
            isfinite(s) && (ok &= check_big(x, y, -s))
            w = nextfloat(x, rand(rng, -2:2))
            ok &= check_big(x, -w, rand_finite(T, 1, 30))
        end
        @test ok
    end

    @testset "$T: subnormal & overflow boundaries" begin
        ok = true
        for _ in 1:30_000
            ok &= check_big(rand_finite(T, 0, 40), rand_finite(T, 0, 40),
                            rand_finite(T, 0, 40))
            ok &= check_big(rand_finite(T, 0, 40), rand_finite(T, 1, EFM),
                            rand_finite(T, 0, 40))
            ok &= check_big(rand_finite(T, EFM - 40, EFM), rand_finite(T, EFM - 40, EFM),
                            rand_finite(T, EFM - 40, EFM))
        end
        tiny = nextfloat(zero(T))
        ok &= check_big(tiny, tiny, tiny)
        ok &= check_big(floatmin(T), -prevfloat(floatmin(T)), tiny)
        @test ok
        @test bits(faaT(floatmax(T), floatmax(T), -floatmax(T))) == bits(floatmax(T))
        @test isinf(faaT(floatmax(T), floatmax(T), zero(T)))
    end

    @testset "$T: ties & double rounding" begin
        one_ = one(T); u = T(2)^(-p); uu = T(2)^(-2p)
        @test faaT(one_, u, zero(T)) == one_                # exact tie -> even
        @test faaT(one_, u, uu) == nextfloat(one_)          # above the tie
        @test faaT(one_, u, -uu) == one_                    # below the tie
        @test (one_ + u) + uu == one_                       # sequential misses it
        ok = true
        for _ in 1:30_000
            e1 = rand(rng, (EFM ÷ 2):(EFM ÷ 2 + 60))
            fr = (rand(rng, uof(T)) & ((uof(T)(1) << (fracbits(T) ÷ 2)) - one(uof(T)))) <<
                 (fracbits(T) - fracbits(T) ÷ 2)
            x = reinterpret(T, (uof(T)(e1) << fracbits(T)) | fr)
            y = ldexp(one(T), e1 - (2^(expbits(T) - 1) - 1) - p)
            z = ldexp(T(rand(rng, (-1, 0, 1))),
                      e1 - (2^(expbits(T) - 1) - 1) - p - rand(rng, 1:2p))
            ok &= check_big(x, y, z)
            ok &= check_big(x, -y, z)
        end
        @test ok
    end

    @testset "$T: special values" begin
        U = uof(T)
        Q = U(1) << (fracbits(T) - 1)
        EXPM = U(2^expbits(T) - 1) << fracbits(T)
        SGN = U(1) << (8 * sizeof(U) - 1)
        mknan(pl; s=false, q=true) =
            reinterpret(T, (s ? SGN : U(0)) | EXPM | (q ? Q : U(0)) | U(pl))
        DEFAULT = reinterpret(T, SGN | EXPM | Q)
        # non-NaN specials: fused must equal native sequential bit-for-bit
        specials = [zero(T), -zero(T), one(T), -T(2.5), T(Inf), -T(Inf),
                    floatmin(T), -nextfloat(zero(T)), prevfloat(floatmin(T))]
        ok = true
        for x in specials, y in specials, z in specials
            ok &= bits(faaT(x, y, z)) == bits((x + y) + z)
        end
        @test ok
        # the deterministic NaN rule (checked against its specification)
        @test bits(faaT(mknan(0xF), mknan(0x1), one(T))) == bits(mknan(0xF))
        @test bits(faaT(mknan(0x1), one(T), mknan(0xF))) == bits(mknan(0xF))
        @test bits(faaT(mknan(0xF), one(T), mknan(0x1))) == bits(mknan(0xF))
        @test bits(faaT(mknan(0x5; s=true), one(T), mknan(0x5))) == bits(mknan(0x5; s=true))
        @test bits(faaT(mknan(0xA; q=false), one(T), one(T))) == bits(mknan(0xA))
        @test bits(faaT(mknan(0xA; q=false), one(T), mknan(0xC))) == bits(mknan(0xC))
        @test bits(faaT(mknan(0xA), one(T), mknan(0xC; q=false))) == bits(mknan(0xA))
        @test bits(faaT(T(Inf), one(T), -T(Inf))) == bits(DEFAULT)
        @test bits(faaT(T(Inf), -T(Inf), mknan(0xC; q=false))) == bits(DEFAULT)
        @test bits(faaT(T(Inf), one(T), mknan(0xC))) == bits(mknan(0xC))
        @test isinf(faaT(T(Inf), one(T), one(T)))
        # signed zeros
        @test bits(faaT(-zero(T), -zero(T), -zero(T))) == bits(-zero(T))
        @test bits(faaT(-zero(T), zero(T), -zero(T))) == bits(zero(T))
        @test bits(faaT(one(T), -one(T), -zero(T))) == bits(zero(T))
    end

    @testset "$T: Boldo-Melquiond cross-check (adjudicated)" begin
        lo, hi = T === Float64 ? (300, 1700) : (40, 210)   # BM-safe range
        n_bad = 0; n_dis = 0
        for _ in 1:10_000_000
            x = rand_finite(T, lo, hi); y = rand_finite(T, lo, hi)
            z = rand_finite(T, lo, hi)
            r = faaT(x, y, z)
            if bits(r) != bits(sum3_bm(x, y, z))
                n_dis += 1
                nearest_ok(x, y, z, r) || (n_bad += 1)
            end
        end
        @test n_bad == 0
        n_dis > 0 && @info "$T: BM deviated (adjudicated in faa's favor)" n_dis
    end

    @testset "$T: quality" begin
        a, b, c = T(1.5)^3, T(2.5)^-2, T(3.25)
        g(x, y, z) = faaT(x, y, z); g(a, b, c)
        @test @allocated(g(a, b, c)) == 0
        @test Base.return_types(faaT, Tuple{T,T,T}) == [T]
    end
end

@testset "faa: Float64 & Float32" begin
    run_battery(Float64, faa64)
    run_battery(Float32, faa32)

    @testset "cross-format: faa32 vs RN32(faa64) (adjudicated)" begin
        # NOT an identity: the 2p+2 innocuous-double-rounding theorem covers
        # single two-operand operations, but a THREE-addend exact sum can lie
        # arbitrarily close to a Float32 tie point, so rounding through
        # Float64 breaks ties the wrong way (~0.3% of random cases). This is
        # the very reason faa exists. The test requires that in every
        # deviation the BigFloat oracle rules faa32 correct.
        n_bad = 0; n_dis = 0
        for _ in 1:10_000_000
            x = rand_finite(Float32, 0, 254)
            y = rand_finite(Float32, 0, 254)
            z = rand_finite(Float32, 0, 254)
            r32 = faa32(x, y, z)
            v = Float32(faa64(Float64(x), Float64(y), Float64(z)))
            if bits(r32) != bits(v)
                n_dis += 1
                nearest_ok(x, y, z, r32) || (n_bad += 1)
            end
        end
        @test n_bad == 0
        n_dis > 0 && @info "cross-format deviations (adjudicated)" n_dis
    end
end