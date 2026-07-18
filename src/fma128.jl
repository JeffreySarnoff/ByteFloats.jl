# Float128FMA.jl
#
# A pure-Julia, correctly-rounded (IEEE 754 round-to-nearest-even) fused
# multiply-add for Quadmath.Float128, intended for Windows, where Quadmath.jl
# does not define `fma(::Float128, ::Float128, ::Float128)` because calling
# libquadmath's `fmaq` corrupts floating-point state there
# (https://github.com/JuliaMath/Quadmath.jl/issues/31).
#
# Bit-for-bit compatible with libquadmath's `fmaq` — including NaN payload
# propagation and invalid-operation NaN generation (glibc soft-fp x86
# semantics) — so results on Windows are indistinguishable from Linux/macOS.
#
# Design: classic soft-float mulAdd on the raw binary128 bit patterns.
#   * decompose x, y, z into sign / unbiased exponent / 113-bit significand
#     (subnormals normalized);
#   * exact 113x113 -> 226-bit product via four 64x64 -> 128-bit multiplies;
#   * both product and addend placed in a 256-bit fixed-point accumulator
#     (a pair of UInt128) with the leading significand bit at bit 254;
#   * alignment shifts use "shift-right-jam" (sticky bit folded into the LSB).
#     Jamming is exact for round-to-nearest here: a jam can only occur when
#     the exponent gap exceeds ~30 (product shifted) or ~142 (addend shifted),
#     in which case catastrophic cancellation is impossible and the jam bit
#     stays >= 100 positions below the guard bit; conversely, whenever
#     cancellation can occur (|gap| <= 1) the alignment shift loses no bits
#     and the subtraction is exact;
#   * one round-to-nearest-even at 113 bits (guard bit 141, sticky 0..140),
#     with gradual underflow and overflow to +-Inf.
#
# NaN semantics (empirically matched against libquadmath, glibc soft-fp
# `_FP_CHOOSENAN` for x86): among competing NaNs the one with the larger raw
# 112-bit fraction field wins (quiet bit participates in the comparison);
# ties go to the product-stage NaN. Invalid operations (0 * Inf, Inf - Inf)
# generate the x86 "indefinite" NaN 0xffff8000...0, which competes by the
# same rule. The winner is returned with its quiet bit set, sign and payload
# preserved. (Floating-point exception *flags* are not modeled; Julia does
# not expose them for Float128.)
#
# No calls into libquadmath are made (only bit-level reinterpret), so the
# Windows FP-state corruption cannot occur. The code is platform-independent
# and allocation-free; it works on Julia >= 1.6 and is tested on Julia 1.12.
#
# Usage as a package:      pkg> dev path/to/Float128FMA ; using Float128FMA
# Usage as a single file:  include("Float128FMA.jl"); using .Float128FMA
#
# On Windows, loading the module installs `Base.fma` for Float128
# automatically. `fma128(x, y, z)` is always available on every OS and always
# uses this implementation; `Float128FMA.install!()` forces `Base.fma` to use
# it even where a native binding exists.

module Float128FMA

using Quadmath: Float128

export fma128

# ---------------------------------------------------------------------------
# binary128 layout constants
# ---------------------------------------------------------------------------
const SIGN_MASK = UInt128(1) << 127
const EXP_MASK  = UInt128(0x7fff) << 112          # exponent field
const FRAC_MASK = (UInt128(1) << 112) - 1          # stored fraction
const IMPLICIT  = UInt128(1) << 112                # implicit leading bit
const EXP_BIAS  = 16383
const QUIET_BIT = UInt128(1) << 111
# x86 "indefinite" NaN produced by invalid operations (sign set, payload 0)
const DEFAULT_NAN = SIGN_MASK | EXP_MASK | QUIET_BIT

# ---------------------------------------------------------------------------
# 256-bit helpers on (hi::UInt128, lo::UInt128) pairs
# ---------------------------------------------------------------------------

"""Exact 128x128 -> 256-bit product, returned as (hi, lo)."""
@inline function mul256(x::UInt128, y::UInt128)
    m64 = UInt128(typemax(UInt64))
    x0 = x & m64; x1 = x >> 64
    y0 = y & m64; y1 = y >> 64
    p00 = x0 * y0                     # each factor < 2^64: exact in UInt128
    p01 = x0 * y1
    p10 = x1 * y0
    p11 = x1 * y1
    mid = (p00 >> 64) + (p01 & m64) + (p10 & m64)      # < 3*2^64, no overflow
    lo  = (mid << 64) | (p00 & m64)
    hi  = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64)
    return hi, lo
end

@inline function add256(ah::UInt128, al::UInt128, bh::UInt128, bl::UInt128)
    lo = al + bl
    hi = ah + bh + (lo < al ? one(UInt128) : zero(UInt128))
    return hi, lo
end

"""a - b for 256-bit values, assuming a >= b."""
@inline function sub256(ah::UInt128, al::UInt128, bh::UInt128, bl::UInt128)
    lo = al - bl
    hi = ah - bh - (al < bl ? one(UInt128) : zero(UInt128))
    return hi, lo
end

"""Left shift by 0 <= s < 256 (callers never shift set bits past bit 255)."""
@inline function shl256(hi::UInt128, lo::UInt128, s::Int)
    if s == 0
        return hi, lo
    elseif s < 128
        return (hi << s) | (lo >> (128 - s)), lo << s
    else
        return lo << (s - 128), zero(UInt128)          # s == 128 gives (lo, 0)
    end
end

"""
Right shift by s >= 0 with "jamming": if any shifted-out bit is nonzero, the
result's least-significant bit is set (a sticky bit in the lowest position).
"""
@inline function shr256jam(hi::UInt128, lo::UInt128, s::Int)
    if s == 0
        return hi, lo
    elseif s < 128
        sticky = (lo << (128 - s)) != 0
        nlo = (lo >> s) | (hi << (128 - s))
        return hi >> s, nlo | (sticky ? one(UInt128) : zero(UInt128))
    elseif s < 256
        t = s - 128                                    # t == 0: hi << 128 == 0 in Julia
        sticky = (lo != 0) | ((hi << (128 - t)) != 0)
        return zero(UInt128), (hi >> t) | (sticky ? one(UInt128) : zero(UInt128))
    else
        sticky = (hi | lo) != 0
        return zero(UInt128), (sticky ? one(UInt128) : zero(UInt128))
    end
end

@inline lz256(hi::UInt128, lo::UInt128) =
    hi == 0 ? 128 + leading_zeros(lo) : leading_zeros(hi)

# ---------------------------------------------------------------------------
# decompose |x| (finite, nonzero) into (sig, e):
#   value == sig * 2^(e - 112),  sig in [2^112, 2^113)  (subnormals normalized)
# ---------------------------------------------------------------------------
@inline function decompose(a::UInt128)
    ef = Int((a >> 112) % UInt16)
    fr = a & FRAC_MASK
    if ef == 0                                          # subnormal (fr != 0)
        s = leading_zeros(fr) - 15                      # bring MSB to bit 112
        return fr << s, -16382 - s
    else
        return fr | IMPLICIT, ef - EXP_BIAS
    end
end

# ---------------------------------------------------------------------------
# NaN propagation, matching libquadmath (glibc soft-fp, x86 _FP_CHOOSENAN):
# among two NaNs the larger raw fraction field wins, ties to the first.
# ---------------------------------------------------------------------------
@inline choosenan(a::UInt128, b::UInt128) =
    (a & FRAC_MASK) >= (b & FRAC_MASK) ? a : b

# ---------------------------------------------------------------------------
# the fma itself
# ---------------------------------------------------------------------------

"""
    fma128(x::Float128, y::Float128, z::Float128) -> Float128

Correctly rounded (round-to-nearest, ties-to-even) `x*y + z` with a single
rounding, computed in pure Julia integer arithmetic; bit-for-bit compatible
with libquadmath's `fmaq`, including signed zeros, subnormals (gradual
underflow), infinities, NaN payload propagation, and invalid-operation NaN
generation.
"""
function fma128(x::Float128, y::Float128, z::Float128)
    ux = reinterpret(UInt128, x)
    uy = reinterpret(UInt128, y)
    uz = reinterpret(UInt128, z)
    ax = ux & ~SIGN_MASK; ay = uy & ~SIGN_MASK; az = uz & ~SIGN_MASK

    sp = ((ux ⊻ uy) & SIGN_MASK) != 0                  # sign of the product
    sz = (uz & SIGN_MASK) != 0

    # ---- specials: NaN propagation and invalid operations ---------------
    xn = ax > EXP_MASK; yn = ay > EXP_MASK; zn = az > EXP_MASK
    if xn | yn | zn ||
       ((ax == EXP_MASK) & (ay == 0)) || ((ay == EXP_MASK) & (ax == 0))
        # product-stage NaN (an actual NaN operand, or 0 * Inf -> indefinite)
        local t::UInt128
        have_t = true
        if xn & yn
            t = choosenan(ux, uy)
        elseif xn
            t = ux
        elseif yn
            t = uy
        elseif !zn                                     # 0 * Inf, z not NaN
            return reinterpret(Float128, DEFAULT_NAN)
        elseif ((ax == EXP_MASK) & (ay == 0)) || ((ay == EXP_MASK) & (ax == 0))
            t = DEFAULT_NAN                            # 0 * Inf competes with NaN z
        else
            have_t = false                             # only z is NaN
            t = zero(UInt128)
        end
        # soft-fp quiets the product-stage NaN when packing the intermediate
        # result, BEFORE it competes with z (whose fraction stays raw)
        t |= QUIET_BIT
        r = have_t ? (zn ? choosenan(t, uz) : t) : uz
        return reinterpret(Float128, r | QUIET_BIT)
    end
    if ax == EXP_MASK || ay == EXP_MASK                 # x or y infinite (no NaN)
        if az == EXP_MASK && sz != sp
            return reinterpret(Float128, DEFAULT_NAN)   # Inf - Inf
        end
        return reinterpret(Float128, (sp ? SIGN_MASK : zero(UInt128)) | EXP_MASK)
    end
    az == EXP_MASK && return z                          # finite*finite + Inf
    if ax == 0 || ay == 0                               # product is a zero
        if az == 0
            # (+-0) + (+-0): same signs keep the sign, else +0 (RN)
            return sp == sz ?
                reinterpret(Float128, sp ? SIGN_MASK : zero(UInt128)) :
                reinterpret(Float128, zero(UInt128))
        end
        return z
    end

    # ---- exact product in 256-bit fixed point ---------------------------
    sigx, ex = decompose(ax)
    sigy, ey = decompose(ay)
    ph, pl = mul256(sigx, sigy)          # value = P * 2^(ex + ey - 224), MSB at 224 or 225
    msb = (ph >> 97) != 0 ? 225 : 224    # bit 225 of P == bit 97 of ph
    Ep  = ex + ey + (msb - 224)          # value = (M/2^254) * 2^Ep after the shift below
    ph, pl = shl256(ph, pl, 254 - msb)   # normalize: leading bit at position 254

    local Mh::UInt128, Ml::UInt128
    local E::Int
    local sres::Bool

    if az == 0
        Mh, Ml, E, sres = ph, pl, Ep, sp
    else
        sigz, ez = decompose(az)
        zh, zl = shl256(zero(UInt128), sigz, 142)       # leading bit at 254
        d = Ep - ez
        if d >= 0
            zh, zl = shr256jam(zh, zl, min(d, 300))     # align z down to product
            E = Ep
        else
            ph, pl = shr256jam(ph, pl, min(-d, 300))    # align product down to z
            E = ez
        end
        if sp == sz                                     # effective addition
            Mh, Ml = add256(ph, pl, zh, zl)
            sres = sp
            if (Mh & SIGN_MASK) != 0                    # carry into bit 255
                Mh, Ml = shr256jam(Mh, Ml, 1)
                E += 1
            end
        else                                            # effective subtraction
            if ph > zh || (ph == zh && pl >= zl)
                Mh, Ml = sub256(ph, pl, zh, zl)
                sres = sp
            else
                Mh, Ml = sub256(zh, zl, ph, pl)
                sres = sz
            end
            if Mh == 0 && Ml == 0
                return reinterpret(Float128, zero(UInt128))   # exact cancel -> +0 (RN)
            end
            sh = lz256(Mh, Ml) - 1                      # restore leading bit to 254
            if sh > 0
                Mh, Ml = shl256(Mh, Ml, sh)
                E -= sh
            end
        end
    end

    # ---- round to nearest even at 113 bits ------------------------------
    be = E + EXP_BIAS                                   # tentative exponent field
    if be >= 0x7fff                                     # certain overflow
        return reinterpret(Float128,
            (sres ? SIGN_MASK : zero(UInt128)) | EXP_MASK)
    end
    if be <= 0                                          # subnormal range: pre-shift
        Mh, Ml = shr256jam(Mh, Ml, min(1 - be, 300))
        be = 0
    end

    sig    = Mh >> 14                                   # bits 142..254 -> 113-bit significand
    guard  = (Mh >> 13) & 1                             # bit 141
    sticky = ((Mh & ((UInt128(1) << 13) - 1)) | Ml) != 0  # bits 0..140
    if guard == 1 && (sticky || (sig & 1) == 1)
        sig += 1
    end

    local r::UInt128
    if be == 0
        r = sig            # if rounding carried to 2^112 this is exactly the min normal
    else
        if sig == (UInt128(1) << 113)                   # rounding carried out of 113 bits
            sig >>= 1
            be += 1
            be >= 0x7fff && return reinterpret(Float128,
                (sres ? SIGN_MASK : zero(UInt128)) | EXP_MASK)
        end
        r = (UInt128(be) << 112) | (sig & FRAC_MASK)
    end
    return reinterpret(Float128, (sres ? SIGN_MASK : zero(UInt128)) | r)
end

# ---------------------------------------------------------------------------
# wiring into Base
# ---------------------------------------------------------------------------
if Sys.iswindows()
    import Base: fma
    @eval Base.fma(x::Float128, y::Float128, z::Float128) = fma128(x, y, z)
end

end # module
