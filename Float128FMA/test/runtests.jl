# Float128FMA test suite.
#
# Primary oracle (Linux/macOS): libquadmath's fmaq, compared with STRICT bit
# equality — including NaN sign/payload. Secondary, platform-independent
# oracle: the exact value in 40000-bit BigFloat, checked for the round-to-
# nearest / ties-to-even property directly against both Float128 neighbours
# (no reliance on any BigFloat -> Float128 conversion). On Windows only the
# BigFloat oracle runs; the NaN batteries then check the documented rule
# directly.

using Test
using Random
using Printf
using Quadmath
using Quadmath: Float128

# detect the native fmaq binding BEFORE loading our module
const HAVE_NATIVE = hasmethod(fma, Tuple{Float128,Float128,Float128})

using Float128FMA

const F = Float128
bits(x::F) = reinterpret(UInt128, x)
strict_same(a::F, b::F) = bits(a) == bits(b)

setprecision(BigFloat, 40000)

"Round-to-nearest-even check against the exact value (finite inputs)."
function nearest_ok(x::F, y::F, z::F, r::F)
    (isfinite(x) && isfinite(y) && isfinite(z)) || return true
    exact = BigFloat(x) * BigFloat(y) + BigFloat(z)
    if !isfinite(r)
        isnan(r) && return false
        lim = BigFloat(floatmax(F)) +
              (BigFloat(floatmax(F)) - BigFloat(prevfloat(floatmax(F)))) / 2
        return abs(exact) >= lim && (signbit(r) == (exact < 0))
    end
    if iszero(r) && iszero(exact)
        # exact zero: sign rules checked in the directed battery
        return true
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

const FRAC_MASK_T = (UInt128(1) << 112) - 1
rng = MersenneTwister(0x5eed)
rand_bits() = reinterpret(F, rand(rng, UInt128))
function rand_finite(emin::Int, emax::Int)
    ef = UInt128(rand(rng, emin:emax))
    fr = rand(rng, UInt128) & ((UInt128(1) << 112) - 1)
    sn = rand(rng, Bool) ? UInt128(1) << 127 : UInt128(0)
    reinterpret(F, sn | (ef << 112) | fr)
end
rand_normalish() = rand_finite(1, 32766)
rand_midrange()  = rand_finite(16000, 16700)
rand_lowrange()  = rand_finite(0, 300)
rand_highrange() = rand_finite(32400, 32766)
function rand_nan()
    payload = rand(rng, UInt128) & ((UInt128(1) << 112) - 1)
    payload == 0 && (payload = UInt128(1))               # keep it a NaN if quiet bit off
    sn = rand(rng, Bool) ? UInt128(1) << 127 : UInt128(0)
    reinterpret(F, sn | (UInt128(0x7fff) << 112) | payload)
end

check_native(x, y, z) = !HAVE_NATIVE || strict_same(fma128(x, y, z), fma(x, y, z))
check_big(x, y, z)    = nearest_ok(x, y, z, fma128(x, y, z))

@testset "Float128FMA" begin

@testset "random bit patterns (all classes)" begin
    ok = true
    for _ in 1:400_000
        ok &= check_native(rand_bits(), rand_bits(), rand_bits())
    end
    @test ok
end

@testset "normals: mid & broad range" begin
    ok = okb = true
    for i in 1:300_000
        a, b, c = rand_midrange(), rand_midrange(), rand_midrange()
        ok &= check_native(a, b, c); i <= 1500 && (okb &= check_big(a, b, c))
        a, b, c = rand_normalish(), rand_normalish(), rand_normalish()
        ok &= check_native(a, b, c); i <= 1500 && (okb &= check_big(a, b, c))
    end
    @test ok; @test okb
end

@testset "subnormals & underflow boundary" begin
    ok = okb = true
    for i in 1:200_000
        a, b, c = rand_lowrange(), rand_lowrange(), rand_lowrange()
        ok &= check_native(a, b, c); i <= 1500 && (okb &= check_big(a, b, c))
        ok &= check_native(rand_lowrange(), rand_midrange(), rand_lowrange())
    end
    # directed sweep: products landing within a few ulps of the subnormal
    # threshold, all guard/sticky flavours
    for e in -120:2:120, fr in (UInt128(0), UInt128(1), UInt128(3) << 50,
                                (UInt128(1) << 112) - 1)
        x = reinterpret(F, (UInt128(1) << 112) | fr)              # ef=1 * frac
        y = ldexp(F(1.0), e)                                      # product ~ 2^(-16382 + e)
        okb &= check_big(x, y, F(0.0))
        okb &= check_big(x, y, nextfloat(F(0.0)))
        okb &= check_big(x, y, -nextfloat(F(0.0)))
        ok  &= check_native(x, y, F(0.0))
        ok  &= check_native(x, y, nextfloat(F(0.0)))
    end
    @test ok; @test okb
end

@testset "overflow boundary" begin
    ok = okb = true
    for i in 1:150_000
        a, b, c = rand_highrange(), rand_highrange(), rand_highrange()
        ok &= check_native(a, b, c); i <= 1000 && (okb &= check_big(a, b, c))
        ok &= check_native(rand_highrange(), rand_midrange(), rand_highrange())
    end
    okb &= check_big(floatmax(F), F(1.0) + F(2.0)^-113, F(0.0))
    okb &= check_big(floatmax(F), nextfloat(F(1.0)), -floatmax(F))
    okb &= check_big(prevfloat(F(2.0)), prevfloat(F(2.0)), F(2.0)^-200)
    @test ok; @test okb
end

@testset "massive cancellation / division residuals" begin
    ok = okb = true
    for i in 1:200_000
        x = rand_midrange(); y = rand_midrange()
        z = -(x * y)
        isfinite(z) || continue
        ok &= check_native(x, y, z); i <= 1500 && (okb &= check_big(x, y, z))
        q = x / y
        isfinite(q) || continue
        ok &= check_native(q, y, -x)
        if BigFloat(q) * BigFloat(y) == BigFloat(x)
            ok &= iszero(fma128(q, y, -x))               # ByteFloats' witness
        end
    end
    @test ok; @test okb
end

@testset "constructed halfway ties at bit 113" begin
    ok = okb = true
    for i in 1:100_000
        fx = (rand(rng, UInt128) & ((UInt128(1) << 56) - 1)) << 56
        fy = (rand(rng, UInt128) & ((UInt128(1) << 56) - 1)) << 56 | UInt128(1) << 55
        x = reinterpret(F, (UInt128(16383) << 112) | fx)
        y = reinterpret(F, (UInt128(16383) << 112) | fy)
        z = rand(rng, Bool) ? F(0.0) :
            reinterpret(F, (UInt128(16383) << 112) | (rand(rng, UInt128) & FRAC_MASK_T))
        ok &= check_native(x, y, z); i <= 1500 && (okb &= check_big(x, y, z))
    end
    @test ok; @test okb
end

@testset "special-value grid" begin
    specials = [F(0.0), -F(0.0), F(1.0), -F(1.0), F(Inf), -F(Inf), F(NaN),
                floatmax(F), -floatmax(F), floatmin(F), -floatmin(F),
                nextfloat(F(0.0)), -nextfloat(F(0.0)), prevfloat(floatmin(F)),
                nextfloat(F(1.0)), prevfloat(F(1.0)), F(2.0)^-16000, F(2.0)^16000]
    ok = okb = true
    for x in specials, y in specials, z in specials
        ok &= check_native(x, y, z)
        okb &= check_big(x, y, z)
    end
    @test ok; @test okb
    # signed-zero results
    @test bits(fma128(F(0.0), F(1.0), F(0.0)))  == bits(F(0.0))
    @test bits(fma128(-F(0.0), F(1.0), -F(0.0))) == bits(-F(0.0))
    @test bits(fma128(-F(0.0), F(1.0), F(0.0)))  == bits(F(0.0))    # mixed -> +0 (RN)
    @test bits(fma128(F(1.0), F(1.0), -F(1.0)))  == bits(F(0.0))    # exact cancel -> +0
    # classic double-rounding trap, scaled across the exponent range
    for s in (-16000, -8000, 0, 8000, 16000)
        a = F(2.0)^60 + F(2.0)^-52
        b = F(2.0)^60 - F(2.0)^-52
        okd = check_big(a * F(2.0)^s, b, -F(2.0)^(120 + s))
        HAVE_NATIVE && (okd &= check_native(a * F(2.0)^s, b, -F(2.0)^(120 + s)))
        @test okd
    end
    @test fma128(F(2.0)^60 + F(2.0)^-52, F(2.0)^60 - F(2.0)^-52, -F(2.0)^120) == -F(2.0)^-104
end

@testset "NaN payloads & invalid operations (libquadmath semantics)" begin
    EXPM = UInt128(0x7fff) << 112
    Q = UInt128(1) << 111
    mknan(p; s=false, q=true) =
        reinterpret(F, (s ? UInt128(1) << 127 : UInt128(0)) | EXPM |
                       (q ? Q : UInt128(0)) | UInt128(p))
    DEFAULT = reinterpret(F, UInt128(1) << 127 | EXPM | Q)
    # documented rule, checked directly (works on Windows too)
    @test strict_same(fma128(mknan(0xF), mknan(0x1), F(1.0)), mknan(0xF))
    @test strict_same(fma128(mknan(0x1), mknan(0xF), F(1.0)), mknan(0xF))
    @test strict_same(fma128(mknan(0x1), F(1.0), mknan(0xF)), mknan(0xF))
    @test strict_same(fma128(mknan(0xF), F(1.0), mknan(0x1)), mknan(0xF))
    @test strict_same(fma128(mknan(0x5; s=true), F(1.0), mknan(0x5)), mknan(0x5; s=true))
    @test strict_same(fma128(mknan(0x5), F(1.0), mknan(0x5; s=true)), mknan(0x5))
    @test strict_same(fma128(mknan(0xA; q=false), F(1.0), F(1.0)), mknan(0xA))    # quieted
    @test strict_same(fma128(mknan(0xA; q=false), F(1.0), mknan(0xC)), mknan(0xC))
    @test strict_same(fma128(mknan(0xA), F(1.0), mknan(0xC; q=false)), mknan(0xA))
    @test strict_same(fma128(F(0.0), F(Inf), F(1.0)), DEFAULT)
    @test strict_same(fma128(F(0.0), F(Inf), mknan(0xC)), mknan(0xC))
    @test strict_same(fma128(F(0.0), F(Inf), mknan(0x0)), DEFAULT)                # tie -> product stage
    @test strict_same(fma128(F(0.0), F(Inf), mknan(0xC; q=false)), DEFAULT)       # frac C < Q
    @test strict_same(fma128(F(Inf), F(1.0), -F(Inf)), DEFAULT)
    @test strict_same(fma128(F(Inf), F(1.0), mknan(0xC)), mknan(0xC))
    # fuzz vs native
    if HAVE_NATIVE
        ok = true
        for _ in 1:200_000
            pick() = (r = rand(rng, 1:4);
                      r == 1 ? rand_nan() :
                      r == 2 ? rand_normalish() :
                      r == 3 ? reinterpret(F, rand(rng, (UInt128(0), UInt128(1) << 127,
                               EXPM, UInt128(1) << 127 | EXPM))) :
                      rand_bits())
            x, y, z = pick(), pick(), pick()
            ok &= strict_same(fma128(x, y, z), fma(x, y, z))
        end
        @test ok
    end
end

@testset "quality: allocations, inference, Base.fma wiring" begin
    a, b, c = F(1.5)^7, F(2.5)^-3, F(3.25)
    f3(x, y, z) = fma128(x, y, z); f3(a, b, c)
    @test @allocated(f3(a, b, c)) == 0
    @test Base.return_types(fma128, Tuple{F,F,F}) == [F]
    @test hasmethod(fma, Tuple{F,F,F})            # native (Unix) or installed (Windows)
end

end # outer testset
