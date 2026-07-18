# Float128FAA test suite.
#
# There is no native "faaq" to compare against, so three oracles are used:
#   1. The exact value in 40000-bit BigFloat, checked for round-to-nearest /
#      ties-to-even directly against both Float128 neighbours (platform
#      independent; no reliance on any BigFloat -> Float128 conversion).
#   2. For special values (NaN/Inf/zero, where no rounding occurs), the
#      fused sum must be bit-identical to (x + y) + z evaluated with native
#      libquadmath addition (runs on any OS — addq works everywhere).
#   3. A Boldo-Melquiond round-to-odd construction of the correctly rounded
#      three-way sum built from native Float128 +/- (an independent fast
#      implementation), used for mass cross-checking in a safe exponent
#      range; any disagreement is adjudicated by the BigFloat oracle.

using Test
using Random
using Printf
using Quadmath
using Quadmath: Float128
using Float128FAA

const F = Float128
bits(x::F) = reinterpret(UInt128, x)
strict_same(a::F, b::F) = bits(a) == bits(b)

setprecision(BigFloat, 40000)

"Round-to-nearest-even check of r against the exact x + y + z (finite inputs)."
function nearest_ok(x::F, y::F, z::F, r::F)
    (isfinite(x) && isfinite(y) && isfinite(z)) || return true
    exact = BigFloat(x) + BigFloat(y) + BigFloat(z)
    if !isfinite(r)
        isnan(r) && return false
        lim = BigFloat(floatmax(F)) +
              (BigFloat(floatmax(F)) - BigFloat(prevfloat(floatmax(F)))) / 2
        return abs(exact) >= lim && (signbit(r) == (exact < 0))
    end
    if iszero(r)
        return iszero(exact) ? !signbit(r) || (signbit(x) && signbit(y) && signbit(z)) :
               abs(exact) <= BigFloat(nextfloat(F(0.0))) / 2   # rounded to zero
    end
    br = BigFloat(r)
    err = abs(exact - br)
    for nb in (prevfloat(r), nextfloat(r))
        isfinite(nb) || continue
        errn = abs(exact - BigFloat(nb))
        errn < err && return false                       # neighbour strictly closer
        if errn == err                                   # exact tie: must be even
            trailing_zeros(bits(r) | (UInt128(1) << 112)) >=
                trailing_zeros(bits(nb) | (UInt128(1) << 112)) || return false
        end
    end
    return true
end

# --- oracle 3: Boldo-Melquiond correctly-rounded sum3 from native ops -------
@inline function twosum(a::F, b::F)
    s = a + b
    ap = s - b; bp = s - ap
    s, (a - ap) + (b - bp)
end
@inline function roundodd_add(u::F, v::F)
    s, t = twosum(u, v)
    if t != zero(F) && (bits(s) & 1) == 0
        s = t > zero(F) ? nextfloat(s) : prevfloat(s)
    end
    s
end
function sum3_bm(x::F, y::F, z::F)
    uh, ul = twosum(x, y)
    th, tl = twosum(uh, z)
    th + roundodd_add(tl, ul)
end

rng = MersenneTwister(0xfa5eed)
rand_bits() = reinterpret(F, rand(rng, UInt128))
function rand_finite(emin::Int, emax::Int)
    ef = UInt128(rand(rng, emin:emax))
    fr = rand(rng, UInt128) & ((UInt128(1) << 112) - 1)
    sn = rand(rng, Bool) ? UInt128(1) << 127 : UInt128(0)
    reinterpret(F, sn | (ef << 112) | fr)
end
rand_normalish() = rand_finite(1, 32766)
rand_saferange() = rand_finite(2000, 30000)     # BM-safe: no under/overflow anywhere
rand_lowrange()  = rand_finite(0, 300)
rand_highrange() = rand_finite(32400, 32766)
function rand_nan()
    p = rand(rng, UInt128) & ((UInt128(1) << 112) - 1)
    p == 0 && (p = UInt128(1))
    sn = rand(rng, Bool) ? UInt128(1) << 127 : UInt128(0)
    reinterpret(F, sn | (UInt128(0x7fff) << 112) | p)
end

check_big(x, y, z) = nearest_ok(x, y, z, faa(x, y, z))

@testset "Float128FAA" begin

@testset "BigFloat oracle: broad random" begin
    ok = true
    for _ in 1:60_000
        ok &= check_big(rand_normalish(), rand_normalish(), rand_normalish())
    end
    @test ok
end

@testset "BigFloat oracle: exponent-gap sweep" begin
    # every alignment regime: gaps around 0, the 113/142 exactness borders,
    # the 143+ jam region, and huge gaps
    ok = true
    for g1 in (0, 1, 2, 50, 112, 113, 114, 141, 142, 143, 144, 200, 255, 256, 400, 30000),
        g2 in (0, 1, 113, 142, 143, 256, 400)
        for _ in 1:60
            e1 = rand(rng, 16000:16400)
            x = rand_finite(e1, e1)
            y = rand_finite(max(1, e1 - g1), max(1, e1 - g1))
            z = rand_finite(max(1, e1 - g1 - g2), max(1, e1 - g1 - g2))
            ok &= check_big(x, y, z)
        end
    end
    @test ok
end

@testset "cancellation & the return-c path" begin
    ok = true
    for _ in 1:40_000
        x = rand_saferange()
        z = rand_normalish()
        ok &= check_big(x, -x, z)                       # exact a+b cancel -> c
        ok &= strict_same(faa(x, -x, z), z)
        ok &= strict_same(faa(z, x, -x), z)             # any operand order
        y = rand_saferange()
        s = x + y
        ok &= check_big(x, y, -s)                       # residual of an add
        ok &= check_big(-s, x, y)
        w = nextfloat(x, rand(rng, -2:2))
        ok &= check_big(x, -w, rand_lowrange())         # near-cancel + tiny third
    end
    @test ok
end

@testset "subnormals & underflow boundary" begin
    ok = true
    for _ in 1:40_000
        ok &= check_big(rand_lowrange(), rand_lowrange(), rand_lowrange())
        ok &= check_big(rand_lowrange(), rand_normalish(), rand_lowrange())
    end
    tiny = nextfloat(F(0.0))
    ok &= check_big(tiny, tiny, tiny)
    ok &= check_big(floatmin(F), -prevfloat(floatmin(F)), tiny)
    ok &= check_big(tiny, tiny, -tiny)
    @test ok
end

@testset "overflow boundary" begin
    ok = true
    for _ in 1:30_000
        ok &= check_big(rand_highrange(), rand_highrange(), rand_highrange())
    end
    # fused stays finite where sequential addition would overflow
    @test strict_same(faa(floatmax(F), floatmax(F), -floatmax(F)), floatmax(F))
    @test isinf(faa(floatmax(F), floatmax(F), F(0.0)))
    ok &= check_big(floatmax(F), floatmax(F), -floatmax(F))
    ok &= check_big(floatmax(F), prevfloat(floatmax(F)), -floatmax(F))
    @test ok
end

@testset "ties & double rounding" begin
    # 1 + 2^-113 is an exact tie; the third addend breaks it invisibly to
    # sequential addition
    @test faa(F(1.0), F(2.0)^-113, F(0.0)) == F(1.0)                 # tie -> even
    @test faa(F(1.0), F(2.0)^-113, F(2.0)^-226) == nextfloat(F(1.0)) # above tie
    @test faa(F(1.0), F(2.0)^-113, -F(2.0)^-226) == F(1.0)           # below tie
    @test (F(1.0) + F(2.0)^-113) + F(2.0)^-226 == F(1.0)             # sequential misses it
    ok = true
    for _ in 1:40_000
        # x with a short significand, y half an ulp of x, z tiny of random sign
        e1 = rand(rng, 16300:16400)
        fr = (rand(rng, UInt128) & ((UInt128(1) << 56) - 1)) << 56
        x = reinterpret(F, (UInt128(e1) << 112) | fr)
        y = ldexp(F(1.0), e1 - 16383 - 113)
        z = ldexp(F(rand(rng, (-1, 0, 1))), e1 - 16383 - 113 - rand(rng, 1:200))
        ok &= check_big(x, y, z)
        ok &= check_big(x, -y, z)
    end
    @test ok
end

@testset "special values ≡ native (x + y) + z" begin
    # no rounding is involved for special inputs, so the fused result must be
    # bit-identical to sequential native addition (small finite magnitudes so
    # partial sums cannot overflow)
    EXPM = UInt128(0x7fff) << 112
    Q = UInt128(1) << 111
    mknan(p; s=false, q=true) =
        reinterpret(F, (s ? UInt128(1) << 127 : UInt128(0)) | EXPM |
                       (q ? Q : UInt128(0)) | UInt128(p))
    specials = [F(0.0), -F(0.0), F(1.0), -F(2.5), F(Inf), -F(Inf), F(NaN),
                mknan(0xA), mknan(0xB; s=true), mknan(0xC; q=false),
                mknan(0x1; q=false, s=true), floatmin(F), -nextfloat(F(0.0))]
    ok = true
    for x in specials, y in specials, z in specials
        ok &= strict_same(faa(x, y, z), (x + y) + z)
    end
    @test ok
    # fuzz: random specials mixed with modest finite values
    okf = true
    for _ in 1:200_000
        pick() = (r = rand(rng, 1:4);
                  r == 1 ? rand_nan() :
                  r == 2 ? rand_finite(16300, 16400) :
                  r == 3 ? reinterpret(F, rand(rng, (UInt128(0), UInt128(1) << 127,
                           EXPM, UInt128(1) << 127 | EXPM))) :
                  rand_finite(1, 32000))
        x, y, z = pick(), pick(), pick()
        (isnan(x) | isnan(y) | isnan(z) | isinf(x) | isinf(y) | isinf(z) |
         iszero(x) | iszero(y) | iszero(z)) || continue
        okf &= strict_same(faa(x, y, z), (x + y) + z)
    end
    @test okf
    # signed zeros
    @test bits(faa(-F(0.0), -F(0.0), -F(0.0))) == bits(-F(0.0))
    @test bits(faa(-F(0.0), F(0.0), -F(0.0)))  == bits(F(0.0))
    @test bits(faa(F(1.0), -F(1.0), -F(0.0)))  == bits(F(0.0))      # exact cancel -> +0
    @test bits(faa(F(1.0), F(2.0), -F(3.0)))   == bits(F(0.0))
end

@testset "Boldo-Melquiond cross-check (adjudicated)" begin
    n_faa_bad = 0
    n_disagree = 0
    for _ in 1:2_000_000
        x, y, z = rand_saferange(), rand_saferange(), rand_saferange()
        r = faa(x, y, z)
        b = sum3_bm(x, y, z)
        if !strict_same(r, b)
            n_disagree += 1
            nearest_ok(x, y, z, r) || (n_faa_bad += 1)
        end
    end
    @test n_faa_bad == 0
    # report (not a failure) how often the BM construction itself deviated
    n_disagree > 0 && @info "BM oracle deviated (adjudicated in faa's favor)" n_disagree
end

@testset "quality: allocations & inference" begin
    a, b, c = F(1.5)^7, F(2.5)^-3, F(3.25)
    f3(x, y, z) = faa(x, y, z); f3(a, b, c)
    @test @allocated(f3(a, b, c)) == 0
    @test Base.return_types(faa, Tuple{F,F,F}) == [F]
    @test faa128 === faa
end

end # outer testset
