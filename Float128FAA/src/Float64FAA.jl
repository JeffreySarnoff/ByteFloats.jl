# Float64faa.jl
#
# faa(x, y, z) for Float64 — the correctly rounded (IEEE 754 round-to-nearest,
# ties-to-even) value of x + y + z with a SINGLE rounding: the additive
# analogue of fma. Pure Julia bit-level integer arithmetic, no dependencies,
# allocation-free, identical on every platform. Same algorithm as
# Float128FAA.jl, instantiated for binary64; because the significand is
# 53 bits, a single UInt128 serves as the whole fixed-point accumulator.
#
# Design: decompose the operands into sign / unbiased exponent / 53-bit
# significand (subnormals normalized); sort by magnitude |a| >= |b| >= |c|;
# accumulate in a UInt128 anchored at the largest operand with its leading
# significand bit at position 126 (bit 127 is carry headroom); fold b and c
# in with shift-right-jam alignment; round once to nearest-even at 53 bits
# (guard bit 73, sticky bits 0..72), with gradual underflow and overflow.
#
# The three-addend hazards (same as the binary128 analysis, rescaled):
#   * accumulator overflow: three same-sign operands can exceed 2^128, so the
#     bit-127 carry is renormalized after folding in b, BEFORE c is aligned
#     against the updated anchor;
#   * jam vs cancellation: a jam requires an alignment shift >= 75, making
#     the jammed operand < 2^52 while every reachable partial sum stays
#     >= 2^73 - 2^52, so jam bits end below bit 54 < guard bit 73 after any
#     normalization; whenever deep cancellation is possible (gap <= 1) the
#     alignment loses no bits and the subtraction is exact;
#   * a + b == 0 exactly: the answer is c itself, returned directly.
#
# Special values follow IEEE 754 addition applied to the fused sum: only
# actual INPUT infinities produce an infinite result (finite operands whose
# partial sums would overflow do not: faa(floatmax, floatmax, -floatmax) ==
# floatmax); +Inf and -Inf together give the x86 indefinite NaN
# 0xfff8000000000000; an exact zero result is +0, and three zeros give -0
# only when all three are -0. NaN payload propagation is DETERMINISTIC and
# uses the same rule as the Float128 family (glibc soft-fp x86
# _FP_CHOOSENAN): among the stage-1 NaNs (x, y) the larger raw fraction
# field wins (ties to x); the survivor is quieted and competes with z's raw
# fraction (ties to the survivor); the winner is returned quieted, sign and
# payload preserved. Hardware NaN propagation is platform-specific (x86 and
# ARM differ, and compilers may commute addition operands), so for NaN
# payloads the fused result is well-defined here rather than
# hardware-matching; for every non-NaN special input it is bit-identical to
# native (x + y) + z. Exception flags are not modeled.
#
# Usage:  include("Float64faa.jl"); using .Float64FAA
#         faa(1.0, 2.0^-53, 2.0^-106)   # nextfloat(1.0); sequential gives 1.0
# `faa64` is an alias (use it if Float32FAA is loaded alongside, since both
# modules export a function named `faa`).

module Float64FAA

export faa, faa64

const SIGN_MASK = 0x8000000000000000                 # UInt64
const EXP_MASK  = UInt64(0x7ff) << 52
const FRAC_MASK = (UInt64(1) << 52) - 1
const IMPLICIT  = UInt64(1) << 52
const EXP_BIAS  = 1023
const QUIET_BIT = UInt64(1) << 51
const DEFAULT_NAN = SIGN_MASK | EXP_MASK | QUIET_BIT # 0xfff8000000000000

const ANCHOR_SHIFT = 74                              # sig bits at 74..126
const GUARD_BIT    = 73
const STICKY_MASK  = (UInt128(1) << 73) - 1
const CARRY_BIT    = UInt128(1) << 127

"""Right shift with jamming: any shifted-out bit sets the result's LSB."""
@inline function shrjam(m::UInt128, s::Int)
    s == 0 && return m
    s >= 128 && return UInt128(m != 0)
    (m >> s) | UInt128((m << (128 - s)) != 0)
end

"""(sig, e) with value == sig * 2^(e - 52), sig in [2^52, 2^53)."""
@inline function decompose(a::UInt64)
    ef = Int(a >> 52)                                 # sign already stripped
    fr = a & FRAC_MASK
    if ef == 0                                        # subnormal (fr != 0)
        s = leading_zeros(fr) - 11                    # bring MSB to bit 52
        return fr << s, -1022 - s
    else
        return fr | IMPLICIT, ef - EXP_BIAS
    end
end

"""NaN choice: larger raw fraction field wins, ties to the first argument."""
@inline choosenan(a::UInt64, b::UInt64) =
    (a & FRAC_MASK) >= (b & FRAC_MASK) ? a : b

"""
    faa(x::Float64, y::Float64, z::Float64) -> Float64
    faa64(x, y, z)

Correctly rounded (round-to-nearest, ties-to-even) `x + y + z` with a single
rounding.
"""
function faa(x::Float64, y::Float64, z::Float64)
    ux = reinterpret(UInt64, x)
    uy = reinterpret(UInt64, y)
    uz = reinterpret(UInt64, z)
    ax = ux & ~SIGN_MASK; ay = uy & ~SIGN_MASK; az = uz & ~SIGN_MASK
    sx = (ux & SIGN_MASK) != 0
    sy = (uy & SIGN_MASK) != 0
    sz = (uz & SIGN_MASK) != 0

    xinf = ax == EXP_MASK; yinf = ay == EXP_MASK; zinf = az == EXP_MASK
    xn = ax > EXP_MASK; yn = ay > EXP_MASK; zn = az > EXP_MASK

    # ---- NaN propagation (deterministic; mirrors (x + y) + z staging) ----
    if xn | yn | zn
        local t::UInt64
        have_t = true
        if xn & yn
            t = choosenan(ux, uy)
        elseif xn
            t = ux
        elseif yn
            t = uy
        elseif xinf & yinf & (sx != sy)
            t = DEFAULT_NAN                           # stage 1: Inf - Inf
        else
            have_t = false                            # only z is NaN
            t = zero(UInt64)
        end
        t |= QUIET_BIT                                # stage-1 NaN is quieted
        r = have_t ? (zn ? choosenan(t, uz) : t) : uz
        return reinterpret(Float64, r | QUIET_BIT)
    end

    # ---- infinities (no NaN): opposite input infinities are invalid ------
    if xinf | yinf | zinf
        hasp = (xinf & !sx) | (yinf & !sy) | (zinf & !sz)
        hasn = (xinf & sx) | (yinf & sy) | (zinf & sz)
        (hasp & hasn) && return reinterpret(Float64, DEFAULT_NAN)
        return reinterpret(Float64, (hasn ? SIGN_MASK : zero(UInt64)) | EXP_MASK)
    end

    # ---- zeros ----------------------------------------------------------
    if ax == 0 && ay == 0 && az == 0
        return reinterpret(Float64, (sx & sy & sz) ? SIGN_MASK : zero(UInt64))
    end

    # ---- sort by magnitude (|bits| order IS magnitude order; zeros sink) -
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
    m2 == 0 && return reinterpret(Float64, u1)        # single nonzero addend

    # ---- accumulate in a UInt128, anchored at m1 -------------------------
    sig1, e1 = decompose(m1)
    M = UInt128(sig1) << ANCHOR_SHIFT                 # leading bit at 126
    E = e1
    sres = sg1

    sig2, e2 = decompose(m2)
    B = shrjam(UInt128(sig2) << ANCHOR_SHIFT, min(e1 - e2, 200))
    if sg2 == sres
        M += B
        if (M & CARRY_BIT) != 0                       # carry into bit 127:
            M = shrjam(M, 1)                          # renormalize NOW so that
            E += 1                                    # adding c cannot overflow
        end
    else
        M -= B                                        # |a| >= |b|: no borrow
        if M == 0
            # a + b cancels exactly; the answer is c itself (or +0)
            return m3 == 0 ? 0.0 : reinterpret(Float64, u3)
        end
    end

    if m3 != 0
        sig3, e3 = decompose(m3)
        C = shrjam(UInt128(sig3) << ANCHOR_SHIFT, min(E - e3, 200))
        if sg3 == sres
            M += C
        elseif M >= C
            M -= C
            M == 0 && return 0.0                      # exact cancel -> +0
        else
            M = C - M                                 # sign flips to c's
            sres = sg3
        end
    end

    if (M & CARRY_BIT) != 0                           # final carry into bit 127
        M = shrjam(M, 1)
        E += 1
    end
    sh = leading_zeros(M) - 1                         # restore leading bit to 126
    if sh > 0
        M <<= sh
        E -= sh
    end

    # ---- round to nearest even at 53 bits --------------------------------
    be = E + EXP_BIAS                                 # tentative exponent field
    if be >= 0x7ff                                    # certain overflow
        return reinterpret(Float64, (sres ? SIGN_MASK : zero(UInt64)) | EXP_MASK)
    end
    if be <= 0                                        # subnormal range: pre-shift
        M = shrjam(M, min(1 - be, 200))
        be = 0
    end

    sig    = (M >> ANCHOR_SHIFT) % UInt64             # bits 74..126 -> 53 bits
    guard  = (M >> GUARD_BIT) & 1
    sticky = (M & STICKY_MASK) != 0
    if guard == 1 && (sticky || (sig & 1) == 1)
        sig += 1
    end

    local r::UInt64
    if be == 0
        r = sig           # a carry to 2^52 is exactly the min normal
    else
        if sig == (UInt64(1) << 53)                   # rounding carried out of 53 bits
            sig >>= 1
            be += 1
            be >= 0x7ff && return reinterpret(Float64,
                (sres ? SIGN_MASK : zero(UInt64)) | EXP_MASK)
        end
        r = (UInt64(be) << 52) | (sig & FRAC_MASK)
    end
    return reinterpret(Float64, (sres ? SIGN_MASK : zero(UInt64)) | r)
end

const faa64 = faa

end # module