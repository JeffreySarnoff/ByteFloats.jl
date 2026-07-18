# Float128FAA.jl
#
# faa(x, y, z) — a "fused add-add" for Quadmath.Float128: the correctly
# rounded (IEEE 754 round-to-nearest, ties-to-even) value of x + y + z with a
# SINGLE rounding, in pure Julia. Companion to Float128FMA.jl; like it, this
# never calls into libquadmath (only bit-level reinterpret and integer
# arithmetic), so it works on Windows, where Float128 fma is unavailable
# (https://github.com/JuliaMath/Quadmath.jl/issues/31), and is identical on
# every platform. Julia >= 1.6; developed and validated on Julia 1.12.
#
# Design: soft-float three-way addition on the raw binary128 bit patterns.
#   * decompose the operands into sign / unbiased exponent / 113-bit
#     significand (subnormals normalized) and sort by magnitude |a|>=|b|>=|c|;
#   * accumulate in 256-bit fixed point (a pair of UInt128) anchored at the
#     largest operand, whose leading significand bit sits at bit 254;
#   * fold b and then c in with "shift-right-jam" alignment (any shifted-out
#     bit sets the LSB as a sticky);
#   * one round-to-nearest-even at 113 bits, with gradual underflow and
#     overflow to +-Inf.
#
# Why jamming is exact here (the three-addend hazard): two operands can
# cancel and leave the third as the whole answer. Analysis by cases, with
# d1 = e_a - e_b and d2 = e_a - e_c the alignment shifts:
#   * a jam requires a shift > 142, i.e. the jammed operand is < 2^112 in
#     accumulator units while the running sum stays >= 2^141 - 2^112, so the
#     jam bit ends >= 28 positions below the round/guard position after any
#     normalization shift — harmless sticky;
#   * catastrophic cancellation between a and b requires d1 <= 1, where the
#     alignment shift loses no bits and the subtraction is EXACT, a multiple
#     of 2^141; the single unrepairable case — a + b == 0 exactly with c far
#     below — is detected and returns c itself (trivially correctly rounded);
#   * cancellation between (a +- b) and c requires c's magnitude to be within
#     one bit of the partial sum, which forces d2 <= 114 — no jam — so that
#     subtraction is exact too.
#
# Special values follow IEEE 754 addition applied to the fused sum: only
# actual INPUT infinities count (finite operands whose partial sums would
# overflow do not create an infinity unless the fused result overflows);
# +Inf and -Inf together give the x86 "indefinite" NaN 0xffff8000...0; an
# exact zero result is +0 in round-to-nearest, and the sum of three zeros is
# -0 only when all three are -0. NaN propagation mirrors libquadmath's
# addition chain (x + y) + z (glibc soft-fp, x86 _FP_CHOOSENAN): among the
# stage-1 NaNs (x, y) the larger raw 112-bit fraction wins (ties to x), the
# survivor is quieted and then competes with z's raw fraction (ties to the
# stage-1 survivor); the winner is returned quieted, sign and payload
# preserved. For special-value inputs (where no rounding is involved) the
# result is therefore bit-identical to evaluating (x + y) + z with native
# libquadmath addition. Floating-point exception flags are not modeled.
#
# Usage as a package:      pkg> dev path/to/Float128FAA ; using Float128FAA
# Usage as a single file:  include("Float128FAA.jl"); using .Float128FAA
#
#     faa(Float128(1), Float128(2)^-113, Float128(2)^-226)   # nextfloat(1.0)
#     # sequential addition gives 1.0 — the fused sum sees past the tie
#
# `faa128` is an alias of `faa`. There is no `Base.faa` to extend; call the
# module's function directly.

module Float128FAA

using Quadmath: Float128

export faa, faa128

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

"""NaN choice, matching libquadmath: larger raw fraction wins, ties to first."""
@inline choosenan(a::UInt128, b::UInt128) =
    (a & FRAC_MASK) >= (b & FRAC_MASK) ? a : b

# ---------------------------------------------------------------------------
# the fused add-add
# ---------------------------------------------------------------------------

"""
    faa(x::Float128, y::Float128, z::Float128) -> Float128
    faa128(x, y, z)

Correctly rounded (round-to-nearest, ties-to-even) `x + y + z` with a single
rounding, computed in pure Julia integer arithmetic. Handles signed zeros,
subnormals (gradual underflow), infinities, and NaN payload propagation with
the same semantics as libquadmath addition applied to `(x + y) + z` — except
that, being fused, only actual input infinities produce an infinite result.
"""
function faa(x::Float128, y::Float128, z::Float128)
    ux = reinterpret(UInt128, x)
    uy = reinterpret(UInt128, y)
    uz = reinterpret(UInt128, z)
    ax = ux & ~SIGN_MASK; ay = uy & ~SIGN_MASK; az = uz & ~SIGN_MASK
    sx = (ux & SIGN_MASK) != 0
    sy = (uy & SIGN_MASK) != 0
    sz = (uz & SIGN_MASK) != 0

    xinf = ax == EXP_MASK; yinf = ay == EXP_MASK; zinf = az == EXP_MASK
    xn = ax > EXP_MASK; yn = ay > EXP_MASK; zn = az > EXP_MASK

    # ---- NaN propagation, mirroring (x + y) + z through soft-fp ----------
    if xn | yn | zn
        local t::UInt128
        have_t = true
        if xn & yn
            t = choosenan(ux, uy)
        elseif xn
            t = ux
        elseif yn
            t = uy
        elseif xinf & yinf & (sx != sy)
            t = DEFAULT_NAN                            # stage 1: Inf - Inf
        else
            have_t = false                             # only z is NaN
            t = zero(UInt128)
        end
        # soft-fp quiets the stage-1 NaN when packing the intermediate
        # result, BEFORE it competes with z (whose fraction stays raw)
        t |= QUIET_BIT
        r = have_t ? (zn ? choosenan(t, uz) : t) : uz
        return reinterpret(Float128, r | QUIET_BIT)
    end

    # ---- infinities (no NaN): opposite input infinities are invalid ------
    if xinf | yinf | zinf
        hasp = (xinf & !sx) | (yinf & !sy) | (zinf & !sz)
        hasn = (xinf & sx) | (yinf & sy) | (zinf & sz)
        (hasp & hasn) && return reinterpret(Float128, DEFAULT_NAN)
        return reinterpret(Float128, (hasn ? SIGN_MASK : zero(UInt128)) | EXP_MASK)
    end

    # ---- zeros ----------------------------------------------------------
    if ax == 0 && ay == 0 && az == 0
        # -0 only when every addend is -0 (IEEE addition, RN)
        return reinterpret(Float128, (sx & sy & sz) ? SIGN_MASK : zero(UInt128))
    end

    # ---- sort by magnitude: (m1,s1,u1) >= (m2,s2,u2) >= (m3,s3,u3) -------
    # (for finite values the |bits| order IS the magnitude order; zeros sink)
    m1, sg1, u1 = ax, sx, ux
    m2, sg2, u2 = ay, sy, uy
    m3, sg3, u3 = az, sz, uz
    if m2 > m1
        m1, sg1, u1, m2, sg2, u2 = m2, sg2, u2, m1, sg1, u1
    end
    if m3 > m1
        m1, sg1, u1, m3, sg3, u3 = m3, sg3, u3, m1, sg1, u1
    end
    if m3 > m2
        m2, sg2, u2, m3, sg3, u3 = m3, sg3, u3, m2, sg2, u2
    end
    m2 == 0 && return reinterpret(Float128, u1)         # single nonzero addend

    # ---- accumulate in 256-bit fixed point, anchored at m1 ---------------
    sig1, e1 = decompose(m1)
    Mh, Ml = shl256(zero(UInt128), sig1, 142)           # leading bit at 254
    E = e1
    sres = sg1

    sig2, e2 = decompose(m2)
    Bh, Bl = shl256(zero(UInt128), sig2, 142)
    Bh, Bl = shr256jam(Bh, Bl, min(e1 - e2, 300))
    if sg2 == sres
        Mh, Ml = add256(Mh, Ml, Bh, Bl)
        if (Mh & SIGN_MASK) != 0                        # carry into bit 255:
            Mh, Ml = shr256jam(Mh, Ml, 1)               # renormalize NOW so that
            E += 1                                      # adding c cannot overflow
        end
    else
        Mh, Ml = sub256(Mh, Ml, Bh, Bl)                 # |a| >= |b|: no borrow
        if (Mh | Ml) == 0
            # a + b cancels exactly; the answer is c itself (or +0)
            return m3 == 0 ? reinterpret(Float128, zero(UInt128)) :
                             reinterpret(Float128, u3)
        end
    end

    if m3 != 0
        sig3, e3 = decompose(m3)
        Ch, Cl = shl256(zero(UInt128), sig3, 142)
        Ch, Cl = shr256jam(Ch, Cl, min(E - e3, 300))    # align against CURRENT anchor
        if sg3 == sres
            Mh, Ml = add256(Mh, Ml, Ch, Cl)
        elseif (Mh > Ch) || (Mh == Ch && Ml >= Cl)
            Mh, Ml = sub256(Mh, Ml, Ch, Cl)
            (Mh | Ml) == 0 && return reinterpret(Float128, zero(UInt128))  # exact cancel -> +0
        else
            Mh, Ml = sub256(Ch, Cl, Mh, Ml)             # sign flips to c's
            sres = sg3
        end
    end

    if (Mh & SIGN_MASK) != 0                            # carry into bit 255
        Mh, Ml = shr256jam(Mh, Ml, 1)
        E += 1
    end
    sh = lz256(Mh, Ml) - 1                              # restore leading bit to 254
    if sh > 0
        Mh, Ml = shl256(Mh, Ml, sh)
        E -= sh
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

const faa128 = faa

end # module
