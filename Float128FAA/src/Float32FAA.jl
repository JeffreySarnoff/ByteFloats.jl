# Float32faa.jl
#
# faa(x, y, z) for Float32 — the correctly rounded (IEEE 754 round-to-nearest,
# ties-to-even) value of x + y + z with a SINGLE rounding: the additive
# analogue of fma. Pure Julia bit-level integer arithmetic, no dependencies,
# allocation-free, identical on every platform. Same algorithm as
# Float128FAA.jl / Float64faa.jl, instantiated for binary32; the 24-bit
# significand lets a single UInt64 serve as the fixed-point accumulator.
#
# Design: decompose into sign / unbiased exponent / 24-bit significand
# (subnormals normalized); sort by magnitude |a| >= |b| >= |c|; accumulate in
# a UInt64 anchored at the largest operand with its leading significand bit
# at position 62 (bit 63 is carry headroom); fold b and c in with
# shift-right-jam alignment; round once to nearest-even at 24 bits (guard
# bit 38, sticky bits 0..37), with gradual underflow and overflow.
#
# Hazards, rescaled from the binary128 analysis: the bit-63 carry is
# renormalized after folding in b, before c is aligned against the updated
# anchor; a jam requires a shift >= 40, making the jammed operand < 2^23
# while every reachable partial sum stays >= 2^38 - 2^23, so jam bits end
# below bit 26 < guard bit 38 after any normalization; deep cancellation
# (gap <= 1) is always exact; a + b == 0 exactly returns c directly.
#
# Special values follow IEEE 754 addition applied to the fused sum: only
# actual INPUT infinities produce an infinite result
# (faa(floatmax, floatmax, -floatmax) == floatmax); +Inf and -Inf together
# give the x86 indefinite NaN 0xffc00000; an exact zero result is +0, and
# three zeros give -0 only when all three are -0. NaN payload propagation is
# DETERMINISTIC, using the same rule as the Float128/Float64 companions
# (larger raw fraction wins, ties to the earlier stage, survivor quieted) —
# hardware NaN propagation is platform-specific, so payloads are well-defined
# here rather than hardware-matching; every non-NaN special input is
# bit-identical to native (x + y) + z. Exception flags are not modeled.
#
# Usage:  include("Float32faa.jl"); using .Float32FAA
#         faa(1.0f0, 2.0f0^-24, 2.0f0^-48)   # nextfloat(1f0); sequential: 1f0
# `faa32` is an alias (use it if Float64FAA is loaded alongside, since both
# modules export a function named `faa`).

module Float32FAA

export faa, faa32

const SIGN_MASK = 0x80000000                          # UInt32
const EXP_MASK  = UInt32(0xff) << 23
const FRAC_MASK = (UInt32(1) << 23) - one(UInt32)   # NB: -1 would promote to Int64
const IMPLICIT  = UInt32(1) << 23
const EXP_BIAS  = 127
const QUIET_BIT = UInt32(1) << 22
const DEFAULT_NAN = SIGN_MASK | EXP_MASK | QUIET_BIT  # 0xffc00000

const ANCHOR_SHIFT = 39                               # sig bits at 39..62
const GUARD_BIT    = 38
const STICKY_MASK  = (UInt64(1) << 38) - 1
const CARRY_BIT    = UInt64(1) << 63

"""Right shift with jamming: any shifted-out bit sets the result's LSB."""
@inline function shrjam(m::UInt64, s::Int)
    s == 0 && return m
    s >= 64 && return UInt64(m != 0)
    (m >> s) | UInt64((m << (64 - s)) != 0)
end

"""(sig, e) with value == sig * 2^(e - 23), sig in [2^23, 2^24)."""
@inline function decompose(a::UInt32)
    ef = Int(a >> 23)                                 # sign already stripped
    fr = a & FRAC_MASK
    if ef == 0                                        # subnormal (fr != 0)
        s = leading_zeros(fr) - 8                     # bring MSB to bit 23
        return fr << s, -126 - s
    else
        return fr | IMPLICIT, ef - EXP_BIAS
    end
end

"""NaN choice: larger raw fraction field wins, ties to the first argument."""
@inline choosenan(a::UInt32, b::UInt32) =
    (a & FRAC_MASK) >= (b & FRAC_MASK) ? a : b

"""
    faa(x::Float32, y::Float32, z::Float32) -> Float32
    faa32(x, y, z)

Correctly rounded (round-to-nearest, ties-to-even) `x + y + z` with a single
rounding.
"""
function faa(x::Float32, y::Float32, z::Float32)
    ux = reinterpret(UInt32, x)
    uy = reinterpret(UInt32, y)
    uz = reinterpret(UInt32, z)
    ax = ux & ~SIGN_MASK; ay = uy & ~SIGN_MASK; az = uz & ~SIGN_MASK
    sx = (ux & SIGN_MASK) != 0
    sy = (uy & SIGN_MASK) != 0
    sz = (uz & SIGN_MASK) != 0

    xinf = ax == EXP_MASK; yinf = ay == EXP_MASK; zinf = az == EXP_MASK
    xn = ax > EXP_MASK; yn = ay > EXP_MASK; zn = az > EXP_MASK

    # ---- NaN propagation (deterministic; mirrors (x + y) + z staging) ----
    if xn | yn | zn
        local t::UInt32
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
            t = zero(UInt32)
        end
        t |= QUIET_BIT                                # stage-1 NaN is quieted
        r = have_t ? (zn ? choosenan(t, uz) : t) : uz
        return reinterpret(Float32, r | QUIET_BIT)
    end

    # ---- infinities (no NaN): opposite input infinities are invalid ------
    if xinf | yinf | zinf
        hasp = (xinf & !sx) | (yinf & !sy) | (zinf & !sz)
        hasn = (xinf & sx) | (yinf & sy) | (zinf & sz)
        (hasp & hasn) && return reinterpret(Float32, DEFAULT_NAN)
        return reinterpret(Float32, (hasn ? SIGN_MASK : zero(UInt32)) | EXP_MASK)
    end

    # ---- zeros ----------------------------------------------------------
    if ax == 0 && ay == 0 && az == 0
        return reinterpret(Float32, (sx & sy & sz) ? SIGN_MASK : zero(UInt32))
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
    m2 == 0 && return reinterpret(Float32, u1)        # single nonzero addend

    # ---- accumulate in a UInt64, anchored at m1 --------------------------
    sig1, e1 = decompose(m1)
    M = UInt64(sig1) << ANCHOR_SHIFT                  # leading bit at 62
    E = e1
    sres = sg1

    sig2, e2 = decompose(m2)
    B = shrjam(UInt64(sig2) << ANCHOR_SHIFT, min(e1 - e2, 100))
    if sg2 == sres
        M += B
        if (M & CARRY_BIT) != 0                       # carry into bit 63:
            M = shrjam(M, 1)                          # renormalize NOW so that
            E += 1                                    # adding c cannot overflow
        end
    else
        M -= B                                        # |a| >= |b|: no borrow
        if M == 0
            # a + b cancels exactly; the answer is c itself (or +0)
            return m3 == 0 ? 0.0f0 : reinterpret(Float32, u3)
        end
    end

    if m3 != 0
        sig3, e3 = decompose(m3)
        C = shrjam(UInt64(sig3) << ANCHOR_SHIFT, min(E - e3, 100))
        if sg3 == sres
            M += C
        elseif M >= C
            M -= C
            M == 0 && return 0.0f0                    # exact cancel -> +0
        else
            M = C - M                                 # sign flips to c's
            sres = sg3
        end
    end

    if (M & CARRY_BIT) != 0                           # final carry into bit 63
        M = shrjam(M, 1)
        E += 1
    end
    sh = leading_zeros(M) - 1                         # restore leading bit to 62
    if sh > 0
        M <<= sh
        E -= sh
    end

    # ---- round to nearest even at 24 bits --------------------------------
    be = E + EXP_BIAS                                 # tentative exponent field
    if be >= 0xff                                     # certain overflow
        return reinterpret(Float32, (sres ? SIGN_MASK : zero(UInt32)) | EXP_MASK)
    end
    if be <= 0                                        # subnormal range: pre-shift
        M = shrjam(M, min(1 - be, 100))
        be = 0
    end

    sig    = (M >> ANCHOR_SHIFT) % UInt32             # bits 39..62 -> 24 bits
    guard  = (M >> GUARD_BIT) & 1
    sticky = (M & STICKY_MASK) != 0
    if guard == 1 && (sticky || (sig & 1) == 1)
        sig += 1
    end

    local r::UInt32
    if be == 0
        r = sig           # a carry to 2^23 is exactly the min normal
    else
        if sig == (UInt32(1) << 24)                   # rounding carried out of 24 bits
            sig >>= 1
            be += 1
            be >= 0xff && return reinterpret(Float32,
                (sres ? SIGN_MASK : zero(UInt32)) | EXP_MASK)
        end
        r = (UInt32(be) << 23) | (sig & FRAC_MASK)
    end
    return reinterpret(Float32, (sres ? SIGN_MASK : zero(UInt32)) | r)
end

const faa32 = faa

end # module
