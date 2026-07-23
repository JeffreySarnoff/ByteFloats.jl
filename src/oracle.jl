# ===== oracle.jl — the ω-semantics catalog (design §9.2, architecture §5,
#                   Float128 revision plan Sites B/C/D)
#
# Rigorous evaluation of every operation's defined result. Never on the hot path.
#
# Result protocol (consumed by ops_scalar.jl's _finish):
#   Float64      exact result: draft-tabled specials and exactly representable finites
#   Float128     exact result by width analysis (Class R: sums of decoded datums whose
#                required significand width fits 113 bits — see the _DE_* thresholds)
#   BigExactF    f() → BigFloat, exact at 2200 bits (the wide-spread tail)
#   Enclose128F  correctly-rounded Float128 bracket (Class R: IEEE-CR sqrt/rsqrt;
#                Divide/Recip now resolve via EncloseF's yd — see below) with the
#                MPFR closure as grid-straddle fallback
#   EncloseF     f(prec) → MPFR directed (lo, hi); optional fq Float128 pre-filter
#                (Class E: |truth − fq()| ≤ |fq()|·2^-90, ≥2^18 slack over published
#                libquadmath bounds, discharged by the differential-build tests);
#                optional yd eager Float64 estimate (faithful libm or IEEE-CR
#                quotient: |truth − yd| ≤ |yd|·2^-45) tried before fq
#
# Special-value rows are written out explicitly even where the substrate coincides:
# the rows are the spec, not an optimization. Interpretations beyond the readable
# draft text are marked [interp] and tracked in checkpoint.md.

const _BIGP = 2200   # exact-arithmetic precision: > full Float64 exponent span + slack

# ---- error-free transforms (Class-1 arithmetic, design §5.1)
@inline function _twosum(a::Float64, b::Float64)
    s = a + b
    bb = s - a
    (s, (a - (s - bb)) + (b - bb))
end

# ---- Float128 exactness-by-width thresholds (plan Site B, Class R).
# Decoded datums carry ≤8-bit significands (Float64-exact products ≤17); a sum is
# exactly representable in Float128 (113-bit) when operand bits + exponent spread
# + carry fit. Escalation is only reached when the Float64 residual ≠ 0 (spread > 44),
# so Float128 covers the whole escalation band except the P ∈ {1,2} extreme tail.
const _DE_ADD = 100    # 9 + ΔE ≤ 113, margin 4
const _DE_FMA = 92     # 17-bit exact product + 8-bit addend: 18 + ΔE ≤ 113, margin 3
const _DE_FAA = 98     # three 8-bit terms: 11 + span ≤ 113, margin 4
@inline _expdiff(a::Float64, b::Float64) = abs(exponent(a) - exponent(b))
@inline function _span3(x::Float64, y::Float64, z::Float64)
    lo = typemax(Int); hi = typemin(Int)
    for v in (x, y, z)
        iszero(v) && continue
        e = exponent(v)
        lo = min(lo, e); hi = max(hi, e)
    end
    hi - lo
end

# Float128 twin of _twosum (IEEE-CR Float128 add/sub, so the Knuth transform is exact)
@inline function _twosum128(a::Float128, b::Float128)
    s = a + b
    bb = s - a
    (s, (a - (s - bb)) + (b - bb))
end

# ---- wide-spread sticky escalations (non-allocating replacements for the
# BigFloat tail; see StickyF's soundness note in ops_scalar.jl).
#
# FMA, ΔE(p, z) > _DE_FMA = 92: the exact 17-bit product p and the ≤8-bit addend z
# are bit-disjoint by > 75 binades, so the larger term is the head and the smaller
# contributes only its sign. Threshold-grid bound (head h, tail w): h is a multiple
# of 2^(e_h−16) while |w| < |h|·2^-92, far below both the finest stochastic
# sub-grid unit (≥ 2^(e_h−7−61)) and h's distance to any off-grid threshold
# (≥ 2^(e_h−16)).
@inline function _fma_wide(p::Float64, z::Float64)
    abs(p) >= abs(z) ? StickyF(p, signbit(z) ? -1 : 1) :
                       StickyF(z, signbit(p) ? -1 : 1)
end

# FAA, span > _DE_FAA = 98 (cancellation among x, y, z possible): distill the three
# exact Float128 terms into (v1, v2, v3) with v1+v2+v3 invariant (each 2sum pass is
# an exact transform) until v1 = fl128(Σ) up to its own lsb and the residual
# v2 + v3 has a determinable sign strictly below lsb(v1) = 2^(e₁−112):
#   · 2sum guarantees |v2| ≤ ulp(v1)/2 = 2^(e₁−113);
#   · accepting requires |v3| < |v2| (sign = sign(v2), |v2+v3| < 2^(e₁−112))
#     or v3 = 0 (sign = sign(v2)) or v2 = 0 with |v3| < 2^(e₁−112) (sign = sign(v3)).
# Priest-style distillation on 3 sorted terms converges in ≤ 3 sweeps; 6 is margin,
# with the (unreachable in testing) MPFR fallback preserving semantics regardless.
function _faa_wide(x::Float64, y::Float64, z::Float64)
    v1, v2, v3 = Float128(x), Float128(y), Float128(z)
    for _ in 1:6
        # sort descending by magnitude (convergence, not exactness, needs this)
        if abs(v2) < abs(v3); v2, v3 = v3, v2; end
        if abs(v1) < abs(v2); v1, v2 = v2, v1; end
        if abs(v2) < abs(v3); v2, v3 = v3, v2; end
        t, tt = _twosum128(v2, v3)
        s, e  = _twosum128(v1, t)
        v1, v2, v3 = s, e, tt
        if iszero(v2) && iszero(v3)
            return v1                                        # exact (covers Σ = 0)
        end
        iszero(v1) && continue                               # total hiding in v2+v3; resort
        if iszero(v3)
            return StickyF(v1, signbit(v2) ? -1 : 1)
        elseif iszero(v2)
            abs(v3) < ldexp(one(Float128), Base.exponent(v1) - 112) &&
                return StickyF(v1, signbit(v3) ? -1 : 1)
        elseif abs(v3) < abs(v2)
            return StickyF(v1, signbit(v2) ? -1 : 1)
        end
    end
    _bigsum3(x, y, z)
end

_bigsum2(x::Float64, y::Float64) =
    BigExactF(() -> setprecision(() -> BigFloat(x) + BigFloat(y), BigFloat, _BIGP))
_bigfma(x::Float64, y::Float64, z::Float64) =
    BigExactF(() -> setprecision(() -> BigFloat(x) * BigFloat(y) + BigFloat(z), BigFloat, _BIGP))
_bigsum3(x::Float64, y::Float64, z::Float64) =
    BigExactF(() -> setprecision(() -> (BigFloat(x) + BigFloat(y)) + BigFloat(z), BigFloat, _BIGP))

# ---- MPFR directed-enclosure closures (the rigorous ladder; unchanged semantics)
_mpfr1(f::F, x::Float64) where {F} = prec -> setprecision(BigFloat, prec) do
    (setrounding(() -> f(BigFloat(x)), BigFloat, RoundDown),
     setrounding(() -> f(BigFloat(x)), BigFloat, RoundUp))
end
_mpfr2(f::F, x::Float64, y::Float64) where {F} = prec -> setprecision(BigFloat, prec) do
    (setrounding(() -> f(BigFloat(x), BigFloat(y)), BigFloat, RoundDown),
     setrounding(() -> f(BigFloat(x), BigFloat(y)), BigFloat, RoundUp))
end
# enclosure builders: MPFR ladder (f) + optional Float128 pre-filter (fq) +
# optional eager Float64 estimate (yd); see EncloseF's docstring for the
# three-stage resolution each of these feeds
_encl1(f::F, x::Float64; fq=nothing, yd=NaN) where {F} = EncloseF(_mpfr1(f, x), fq, yd)
_encl2(f::F, x::Float64, y::Float64; fq=nothing, yd=NaN) where {F} = EncloseF(_mpfr2(f, x, y), fq, yd)

# Quadmath omits acos(::Float128); π/2 − asin is exact well within the 2^-90 envelope
# (asin faithful + π half-ulp + one subtraction: ≲ 2^-108 relative).
@inline _acos128(x::Float128) = ldexp(Float128(π), -1) - asin(x)

# π · num / 2^k with exact power-of-two scaling; sign-aware directed bounds.
function _encl_piscale(num::Int, k::Int)
    fq = () -> ldexp(Float128(π) * num, -k)
    # eager Float64 estimate: π half-ulp + one multiply + exact ldexp ≤ ~1.5 ulp
    yd = ldexp(Float64(π) * num, -k)
    EncloseF(prec -> setprecision(BigFloat, prec) do
        lo = setrounding(BigFloat, RoundDown) do
            p = num > 0 ? BigFloat(π) : setrounding(() -> BigFloat(π), BigFloat, RoundUp)
            ldexp(p * num, -k)
        end
        hi = setrounding(BigFloat, RoundUp) do
            p = num > 0 ? BigFloat(π) : setrounding(() -> BigFloat(π), BigFloat, RoundDown)
            ldexp(p * num, -k)
        end
        (lo, hi)
    end, fq, yd)
end

# g(x…)/π for the arc-Pi variants: MPFR ladder with sign-aware denominator bounds;
# gq is the Float128 counterpart of g for the pre-filter (envelope covers the
# composed faithful-op + π-constant + CR-division errors with ≥2^16 slack).
function _encl_divpi(g::F, gq::G, xs::Float64...) where {F,G}
    fq = () -> gq(map(Float128, xs)...) / Float128(π)
    # eager Float64 estimate: faithful g + π half-ulp + CR division ≤ ~2 ulp; the
    # division by π ≈ 3.14 is well-conditioned, so no amplification anywhere
    yd = g(xs...) / Float64(π)
    EncloseF(prec -> setprecision(BigFloat, prec) do
        lo = setrounding(BigFloat, RoundDown) do
            n = g(map(BigFloat, xs)...)
            p = n >= 0 ? setrounding(() -> BigFloat(π), BigFloat, RoundUp) : BigFloat(π)
            n / p
        end
        hi = setrounding(BigFloat, RoundUp) do
            n = g(map(BigFloat, xs)...)
            p = n >= 0 ? setrounding(() -> BigFloat(π), BigFloat, RoundDown) : BigFloat(π)
            n / p
        end
        (lo, hi)
    end, fq, yd)
end

# f(π·r) for the *Pi trig family on the exactly reduced r ∈ (0,2) with the
# extremum/pole/Niven points already peeled: the tiny m-interval contains no interior
# extremum (an ≤8-bit-significand r is never within 2^-9 of one), so the 4-combo
# min/max of directed endpoint evaluations encloses. Pre-filter: one Float128 eval
# (π-constant half-ulp + CR product half-ulp + trig envelope ≪ 2^-90).
# The eager Float64 estimate is supplied by the CALLER as the Base.*pi native
# (sinpi/cospi/tanpi): those are relative-faithful in the *result* domain even
# adjacent to zeros/poles. Evaluating f(Float64(π)·r) instead would amplify the
# π-rounding argument error by |πr|/|f| near the zeros of sin — measured at only
# ~2× envelope slack for the closest format-reachable r — so it is NOT used.
function _encl_pitrig(f::F, r::Float64; yd::Float64=NaN) where {F}
    fq = () -> f(Float128(π) * Float128(r))
    EncloseF(prec -> setprecision(BigFloat, prec) do
        ml = setrounding(() -> BigFloat(π) * r, BigFloat, RoundDown)
        mh = setrounding(() -> BigFloat(π) * r, BigFloat, RoundUp)
        a = setrounding(() -> (f(ml), f(mh)), BigFloat, RoundDown)
        b = setrounding(() -> (f(ml), f(mh)), BigFloat, RoundUp)
        (min(a[1], a[2]), max(b[1], b[2]))
    end, fq, yd)
end

# exact mod-2 reduction for the *Pi family: Float64 `rem` is exact by IEEE
@inline function _mod2(x::Float64)
    r = rem(x, 2.0)
    r < 0.0 ? r + 2.0 : r
end

# ============================================================================
# Group "exact selection": Abs, Negate, CopySign, extremum family, Clamp
# ============================================================================
ωeval(::Val{:Abs}, x::Float64) = isnan(x) ? NaN : abs(x)
ωeval(::Val{:Negate}, x::Float64) = isnan(x) ? NaN : (iszero(x) ? 0.0 : -x)
function ωeval(::Val{:CopySign}, x::Float64, y::Float64)
    (isnan(x) | isnan(y)) && return NaN                     # [interp]: NaN sign source → NaN
    v = copysign(x, y)
    iszero(v) ? 0.0 : v                                     # single zero
end

ωeval(::Val{:Maximum}, x::Float64, y::Float64) = (isnan(x) | isnan(y)) ? NaN : max(x, y)
ωeval(::Val{:Minimum}, x::Float64, y::Float64) = (isnan(x) | isnan(y)) ? NaN : min(x, y)
ωeval(::Val{:MaximumNumber}, x::Float64, y::Float64) =
    isnan(x) ? (isnan(y) ? NaN : y) : (isnan(y) ? x : max(x, y))
ωeval(::Val{:MinimumNumber}, x::Float64, y::Float64) =
    isnan(x) ? (isnan(y) ? NaN : y) : (isnan(y) ? x : min(x, y))
function ωeval(::Val{:MaximumMagnitude}, x::Float64, y::Float64)
    (isnan(x) | isnan(y)) && return NaN
    abs(x) > abs(y) ? x : abs(y) > abs(x) ? y : max(x, y)
end
function ωeval(::Val{:MinimumMagnitude}, x::Float64, y::Float64)
    (isnan(x) | isnan(y)) && return NaN
    abs(x) < abs(y) ? x : abs(y) < abs(x) ? y : min(x, y)
end
function ωeval(::Val{:MaximumMagnitudeNumber}, x::Float64, y::Float64)
    isnan(x) && return isnan(y) ? NaN : y
    isnan(y) && return x
    abs(x) > abs(y) ? x : abs(y) > abs(x) ? y : max(x, y)
end
function ωeval(::Val{:MinimumMagnitudeNumber}, x::Float64, y::Float64)
    isnan(x) && return isnan(y) ? NaN : y
    isnan(y) && return x
    abs(x) < abs(y) ? x : abs(y) < abs(x) ? y : min(x, y)
end
# Finite variants (§4.11.3): prefer finite operands; then infinities beat NaN.
# These are the reduction semantics ConvertToBlockMaxAbsFinite's NaN seed relies on.
function ωeval(::Val{:MaximumFinite}, x::Float64, y::Float64)
    fx, fy = isfinite(x), isfinite(y)
    fx & fy && return max(x, y)
    fx && return x
    fy && return y
    isnan(x) && return isnan(y) ? NaN : y
    isnan(y) && return x
    max(x, y)
end
function ωeval(::Val{:MinimumFinite}, x::Float64, y::Float64)
    fx, fy = isfinite(x), isfinite(y)
    fx & fy && return min(x, y)
    fx && return x
    fy && return y
    isnan(x) && return isnan(y) ? NaN : y
    isnan(y) && return x
    min(x, y)
end
function ωeval(::Val{:Clamp}, x::Float64, lo::Float64, hi::Float64)
    (isnan(x) | isnan(lo) | isnan(hi)) && return NaN
    min(max(x, lo), hi)
end

# ============================================================================
# Group "exact arithmetic": Add, Subtract, Multiply, FMA, FAA
# ============================================================================
function ωeval(::Val{:Add}, x::Float64, y::Float64)
    (isnan(x) | isnan(y)) && return NaN
    if isinf(x) || isinf(y)
        (isinf(x) && isinf(y) && x != y) && return NaN      # ∞ + (−∞) → NaN
        return isinf(x) ? x : y
    end
    s, e = _twosum(x, y)
    e == 0.0 && return iszero(s) ? 0.0 : s
    (_f128() && _expdiff(x, y) <= _DE_ADD) && return Float128(x) + Float128(y)   # exact by width
    _bigsum2(x, y)
end
ωeval(::Val{:Subtract}, x::Float64, y::Float64) =
    ωeval(Val(:Add), x, isnan(y) ? y : (iszero(y) ? 0.0 : -y))
function ωeval(::Val{:Multiply}, x::Float64, y::Float64)
    (isnan(x) | isnan(y)) && return NaN
    ((iszero(x) && isinf(y)) || (isinf(x) && iszero(y))) && return NaN   # 0·∞ → NaN
    p = x * y
    iszero(p) ? 0.0 : p                # exact: ≤16-bit significands, exponents in range
end
function ωeval(::Val{:FMA}, x::Float64, y::Float64, z::Float64)
    (isnan(x) | isnan(y) | isnan(z)) && return NaN
    ((iszero(x) && isinf(y)) || (isinf(x) && iszero(y))) && return NaN
    p = x * y                                                # exact when finite
    if isinf(p) || isinf(z)
        (isinf(p) && isinf(z) && p != z) && return NaN
        return isinf(p) ? p : z
    end
    s, e = _twosum(p, z)
    e == 0.0 && return iszero(s) ? 0.0 : s
    if _f128()                                               # e ≠ 0 ⇒ p, z both nonzero
        _expdiff(p, z) > _DE_FMA && return _fma_wide(p, z)   # sticky head, no allocation
        return Float128(p) + Float128(z)                     # p exact ⇒ sum exact by width
    end
    _bigfma(x, y, z)
end
function ωeval(::Val{:FAA}, x::Float64, y::Float64, z::Float64)
    (isnan(x) | isnan(y) | isnan(z)) && return NaN
    if isinf(x) || isinf(y) || isinf(z)
        hasp = (x == Inf) | (y == Inf) | (z == Inf)
        hasn = (x == -Inf) | (y == -Inf) | (z == -Inf)
        (hasp & hasn) && return NaN
        return hasp ? Inf : -Inf
    end
    s1, e1 = _twosum(x, y)
    s2, e2 = _twosum(s1, z)
    (e1 == 0.0 && e2 == 0.0) && return iszero(s2) ? 0.0 : s2
    if _f128()
        _span3(x, y, z) <= _DE_FAA &&
            return (Float128(x) + Float128(y)) + Float128(z) # every partial exact by width
        return _faa_wide(x, y, z)                            # sticky head, no allocation
    end
    _bigsum3(x, y, z)
end

# ============================================================================
# Group "quotient": Divide, Recip, Sqrt, RSqrt   (Annex A.3 zero/∞ semantics)
# All four rest on IEEE correct rounding (Class R, plan Site C), split two ways:
#   Divide/Recip — the *Float64* quotient is CR (≤ half an ulp), so it serves as
#     EncloseF's eager yd estimate directly; no Float128 arithmetic on this path.
#   Sqrt/RSqrt   — Float128 sqrt/inv are CR, so an inexact nearest-CR result q
#     brackets the truth in (prevfloat(q), nextfloat(q)) → Enclose128F.
# ============================================================================
function ωeval(::Val{:Divide}, x::Float64, y::Float64)
    (isnan(x) | isnan(y)) && return NaN
    iszero(y) && return NaN                                  # x/0 → NaN, all x (A.3)
    (isinf(x) && isinf(y)) && return NaN
    isinf(y) && return 0.0                                   # x/±∞ → 0 (single zero)
    isinf(x) && return y > 0 ? x : -x
    q = x / y
    (isfinite(q) && fma(q, y, -x) == 0.0) && return (iszero(q) ? 0.0 : q)   # exact quotient
    # IEEE-CR Float64 quotient as eager estimate; degenerate q ⇒ yd = NaN ⇒ the
    # Float128 filter / ladder decide. Single call site keeps the union narrow;
    # the former 113-bit-exact branch is dead (a dyadic x/y is exact at ≤ p_x ≤ 53
    # bits whenever it is finite and representable).
    yd = (isfinite(q) && !iszero(q)) ? q : NaN
    _encl2(/, x, y; fq=() -> Float128(x) / Float128(y), yd)
end
function ωeval(::Val{:Recip}, x::Float64)
    isnan(x) && return NaN
    iszero(x) && return NaN                                  # Recip(0) → NaN (A.3)
    isinf(x) && return 0.0                                   # Recip(±∞) → 0 (A.3)
    q = 1.0 / x
    (isfinite(q) && fma(q, x, -1.0) == 0.0) && return q      # exact ⇔ x a power of two
    # IEEE-CR Float64 quotient (≤ half an ulp) is the ideal eager estimate; a
    # degenerate q (over/underflow) sends yd = NaN so _finish skips straight to
    # the Float128 filter / MPFR ladder. Single call site keeps the return union
    # to {Float64, EncloseF} — the former Float128/Enclose128F branch is dead for
    # Float64 inputs (1/x exact at 113 bits ⇔ x a power of two ⇔ exact at 53).
    yd = (isfinite(q) && !iszero(q)) ? q : NaN
    _encl1(inv, x; fq=() -> inv(Float128(x)), yd)
end
function ωeval(::Val{:Sqrt}, x::Float64)
    isnan(x) && return NaN
    x < 0.0 && return NaN
    iszero(x) && return 0.0
    isinf(x) && return Inf
    s = sqrt(x)
    fma(s, s, -x) == 0.0 && return s                          # exact square root
    # IEEE-CR hardware sqrt (≤ half an ulp, unconditional) is the ideal eager
    # estimate; s is always finite/nonzero/normal here (x positive finite), so no
    # degenerate guard is needed. Single call site keeps the return union to
    # {Float64, EncloseF} — the former Float128/Enclose128F branch is dead for
    # Float64 inputs: a dyadic √x has ≤ ⌈53/2⌉ = 27 significand bits, so any
    # 113-bit-exact root is already 53-bit-exact and caught by the fma test above.
    _encl1(sqrt, x; fq=() -> sqrt(Float128(x)), yd=s)
end
_mpfr_rsqrt(x::Float64) = prec -> setprecision(BigFloat, prec) do
    lo = setrounding(BigFloat, RoundDown) do
        s = setrounding(() -> sqrt(BigFloat(x)), BigFloat, RoundUp)
        1 / s
    end
    hi = setrounding(BigFloat, RoundUp) do
        s = setrounding(() -> sqrt(BigFloat(x)), BigFloat, RoundDown)
        1 / s
    end
    (lo, hi)
end
function ωeval(::Val{:RSqrt}, x::Float64)
    isnan(x) && return NaN
    x < 0.0 && return NaN
    iszero(x) && return NaN                                   # 1/√0 = 1/0 → NaN
    isinf(x) && return 0.0
    s = sqrt(x)
    if fma(s, s, -x) == 0.0                                   # √x exact in Float64
        r = 1.0 / s
        fma(r, s, -1.0) == 0.0 && return r                    # 1/√x exact too ⇔ x = 4^k
    end
    # Composition of two IEEE-CR ops (hardware sqrt, then divide, each ≤ half an
    # ulp) gives |truth − yd| ≤ ~1.5 ulp ≈ 2⁻⁵¹·⁴ — ≥ 2⁶ slack under the 2⁻⁴⁵
    # envelope. The fq stage composes the same two CR ops at 113 bits (≪ 2⁻⁹⁰).
    # The former Float128 branch is dead for Float64 inputs (1/√x dyadic ⇔ x a
    # power of 4, caught above); deleting it narrows the union to {Float64, EncloseF}.
    EncloseF(_mpfr_rsqrt(x), () -> inv(sqrt(Float128(x))), 1.0 / s)
end

# ============================================================================
# Group "exponential / logarithmic"
# ============================================================================
function ωeval(::Val{:Exp}, x::Float64)
    isnan(x) && return NaN
    isinf(x) && return x > 0.0 ? Inf : 0.0
    iszero(x) && return 1.0
    _encl1(exp, x; fq=() -> exp(Float128(x)), yd=exp(x))
end
function ωeval(::Val{:Exp2}, x::Float64)
    isnan(x) && return NaN
    isinf(x) && return x > 0.0 ? Inf : 0.0
    iszero(x) && return 1.0
    _encl1(exp2, x; fq=() -> exp2(Float128(x)), yd=exp2(x))
end
function ωeval(::Val{:ExpMinusOne}, x::Float64)
    isnan(x) && return NaN
    x == -Inf && return -1.0
    x == Inf && return Inf
    iszero(x) && return 0.0
    _encl1(expm1, x; fq=() -> expm1(Float128(x)), yd=expm1(x))
end
function ωeval(::Val{:Log}, x::Float64)          # draft §4.11 Log rows
    isnan(x) && return NaN
    x == -Inf && return NaN
    x == Inf && return Inf
    x < 0.0 && return NaN
    iszero(x) && return -Inf
    x == 1.0 && return 0.0
    _encl1(log, x; fq=() -> log(Float128(x)), yd=log(x))
end
function ωeval(::Val{:Log2}, x::Float64)
    isnan(x) && return NaN
    x == -Inf && return NaN
    x == Inf && return Inf
    x < 0.0 && return NaN
    iszero(x) && return -Inf
    # exact powers of two: the fq envelope would straddle the exact integer grid
    # value and defeat the filter; MPFR resolves them via lo == hi, but the dyadic
    # screen here is cheaper and exact (Class R).
    xe = exponent(x)
    x == ldexp(1.0, xe) && return Float64(xe)
    _encl1(log2, x; fq=() -> log2(Float128(x)), yd=log2(x))
end
function ωeval(::Val{:LogOnePlus}, x::Float64)
    isnan(x) && return NaN
    x < -1.0 && return NaN
    x == -1.0 && return -Inf
    x == Inf && return Inf
    iszero(x) && return 0.0
    _encl1(log1p, x; fq=() -> log1p(Float128(x)), yd=log1p(x))
end
function ωeval(::Val{:Softplus}, x::Float64)
    isnan(x) && return NaN
    x == -Inf && return 0.0
    x == Inf && return Inf
    fq = x > 0.0 ? (() -> (q = Float128(x); q + log1p(exp(-q)))) :
                   (() -> log1p(exp(Float128(x))))
    yd = x > 0.0 ? x + log1p(exp(-x)) : log1p(exp(x))
    EncloseF(prec -> setprecision(BigFloat, prec) do
        # monotone-↑ composition: whole-expression directed rounding is an enclosure;
        # the x > 0 form keeps it tight (and exact-side correct) for large x
        f = x > 0.0 ? (b -> b + log1p(exp(-b))) : (b -> log1p(exp(b)))
        (setrounding(() -> f(BigFloat(x)), BigFloat, RoundDown),
         setrounding(() -> f(BigFloat(x)), BigFloat, RoundUp))
    end, fq, yd)
end

# ============================================================================
# Group "trigonometric / hyperbolic"
# ============================================================================
for (name, bf) in ((:Sin, :sin), (:Cos, :cos), (:Tan, :tan))
    zval = name === :Cos ? 1.0 : 0.0
    @eval function ωeval(::Val{$(QuoteNode(name))}, x::Float64)
        isnan(x) && return NaN
        isinf(x) && return NaN                               # trig of ±∞ → NaN
        iszero(x) && return $zval
        # huge-argument reduction accuracy in libquadmath is undocumented; keep the
        # pre-filter inside a conservative window and let MPFR own the rest (§9.2:
        # "undocumented envelope ⇒ fall back")
        abs(x) <= 1.0e15 ? _encl1($bf, x; fq=() -> $bf(Float128(x)), yd=$bf(x)) : _encl1($bf, x)
    end
end
function ωeval(::Val{:ArcSin}, x::Float64)
    isnan(x) && return NaN
    abs(x) > 1.0 && return NaN
    iszero(x) && return 0.0
    x == 1.0 && return _encl_piscale(1, 1)                    #  π/2
    x == -1.0 && return _encl_piscale(-1, 1)                  # −π/2
    _encl1(asin, x; fq=() -> asin(Float128(x)), yd=asin(x))
end
function ωeval(::Val{:ArcCos}, x::Float64)
    isnan(x) && return NaN
    abs(x) > 1.0 && return NaN
    x == 1.0 && return 0.0
    x == -1.0 && return _encl_piscale(1, 0)                   # π
    iszero(x) && return _encl_piscale(1, 1)                   # π/2
    _encl1(acos, x; fq=() -> _acos128(Float128(x)), yd=acos(x))
end
function ωeval(::Val{:ArcTan}, x::Float64)
    isnan(x) && return NaN
    x == Inf && return _encl_piscale(1, 1)
    x == -Inf && return _encl_piscale(-1, 1)
    iszero(x) && return 0.0
    _encl1(atan, x; fq=() -> atan(Float128(x)), yd=atan(x))
end
for (name, bf, inf2) in ((:Sinh, :sinh, :same), (:Cosh, :cosh, :pos), (:Tanh, :tanh, :one))
    @eval function ωeval(::Val{$(QuoteNode(name))}, x::Float64)
        isnan(x) && return NaN
        if isinf(x)
            $(inf2 === :same ? :(return x) :
              inf2 === :pos  ? :(return Inf) :
                               :(return x > 0 ? 1.0 : -1.0))
        end
        iszero(x) && return $(name === :Cosh ? 1.0 : 0.0)
        _encl1($bf, x; fq=() -> $bf(Float128(x)), yd=$bf(x))
    end
end
function ωeval(::Val{:ArcSinh}, x::Float64)
    isnan(x) && return NaN
    isinf(x) && return x
    iszero(x) && return 0.0
    _encl1(asinh, x; fq=() -> asinh(Float128(x)), yd=asinh(x))
end
function ωeval(::Val{:ArcCosh}, x::Float64)
    isnan(x) && return NaN
    x < 1.0 && return NaN
    x == 1.0 && return 0.0
    x == Inf && return Inf
    _encl1(acosh, x; fq=() -> acosh(Float128(x)), yd=acosh(x))
end
function ωeval(::Val{:ArcTanh}, x::Float64)
    isnan(x) && return NaN
    abs(x) > 1.0 && return NaN
    x == 1.0 && return Inf
    x == -1.0 && return -Inf
    iszero(x) && return 0.0
    _encl1(atanh, x; fq=() -> atanh(Float128(x)), yd=atanh(x))
end

# ============================================================================
# Group "π-scaled trigonometric" — exact mod-2 reduction, then enclosure
# ============================================================================
function ωeval(::Val{:SinPi}, x::Float64)
    isnan(x) && return NaN
    isinf(x) && return NaN
    r = _mod2(x)
    r == 0.0 && return 0.0
    r == 0.5 && return 1.0
    r == 1.0 && return 0.0
    r == 1.5 && return -1.0
    _encl_pitrig(sin, r; yd=sinpi(r))
end
function ωeval(::Val{:CosPi}, x::Float64)
    isnan(x) && return NaN
    isinf(x) && return NaN
    r = _mod2(x)
    r == 0.0 && return 1.0
    r == 0.5 && return 0.0
    r == 1.0 && return -1.0
    r == 1.5 && return 0.0
    _encl_pitrig(cos, r; yd=cospi(r))
end
function ωeval(::Val{:TanPi}, x::Float64)
    isnan(x) && return NaN
    isinf(x) && return NaN
    r = _mod2(x)
    r == 0.0 && return 0.0
    r == 0.5 && return Inf                                    # 754-2019 parity convention
    r == 1.0 && return 0.0
    r == 1.5 && return -Inf
    # Niven: for dyadic r the only remaining rational (hence dyadic, hence
    # grid-hittable) values of tan(πr) are ±1 at the quarter-integers — peel them,
    # or the enclosure protocol would chase an interior grid point forever.
    r == 0.25 && return 1.0
    r == 0.75 && return -1.0
    r == 1.25 && return 1.0
    r == 1.75 && return -1.0
    _encl_pitrig(tan, r; yd=tanpi(r))
end
function ωeval(::Val{:ArcSinPi}, x::Float64)
    isnan(x) && return NaN
    abs(x) > 1.0 && return NaN
    iszero(x) && return 0.0
    x == 1.0 && return 0.5
    x == -1.0 && return -0.5
    _encl_divpi(asin, asin, x)
end
function ωeval(::Val{:ArcCosPi}, x::Float64)
    isnan(x) && return NaN
    abs(x) > 1.0 && return NaN
    x == 1.0 && return 0.0
    x == -1.0 && return 1.0
    iszero(x) && return 0.5
    _encl_divpi(acos, _acos128, x)
end
function ωeval(::Val{:ArcTanPi}, x::Float64)
    isnan(x) && return NaN
    x == Inf && return 0.5
    x == -Inf && return -0.5
    iszero(x) && return 0.0
    x == 1.0 && return 0.25
    x == -1.0 && return -0.25
    _encl_divpi(atan, atan, x)
end

# ============================================================================
# Hypot, ArcTan2, ArcTan2Pi   (operand order: (y, x), as IEEE atan2)
# ============================================================================
function ωeval(::Val{:Hypot}, x::Float64, y::Float64)
    (isinf(x) || isinf(y)) && return Inf                     # ∞ dominates even NaN
    (isnan(x) | isnan(y)) && return NaN
    iszero(y) && return abs(x)                               # exact
    iszero(x) && return abs(y)
    _encl2(hypot, x, y; fq=() -> hypot(Float128(x), Float128(y)), yd=hypot(x, y))
end
function ωeval(::Val{:ArcTan2}, y::Float64, x::Float64)      # [interp]: single-zero branch cuts
    (isnan(y) | isnan(x)) && return NaN
    if iszero(y)
        (x > 0.0 || iszero(x)) && return 0.0                 # atan2(0, x≥0) = 0
        return _encl_piscale(1, 0)                           # atan2(0, x<0 or −∞) = π
    end
    if isinf(y)
        x == Inf && return _encl_piscale(y > 0 ? 1 : -1, 2)      # ±π/4
        x == -Inf && return _encl_piscale(y > 0 ? 3 : -3, 2)     # ±3π/4
        return _encl_piscale(y > 0 ? 1 : -1, 1)                  # ±π/2
    end
    x == Inf && return 0.0
    x == -Inf && return _encl_piscale(y > 0 ? 1 : -1, 0)         # ±π
    iszero(x) && return _encl_piscale(y > 0 ? 1 : -1, 1)         # ±π/2
    _encl2(atan, y, x; fq=() -> atan(Float128(y), Float128(x)), yd=atan(y, x))  # MPFR atan2, correct per mode
end
function ωeval(::Val{:ArcTan2Pi}, y::Float64, x::Float64)
    (isnan(y) | isnan(x)) && return NaN
    if iszero(y)
        (x > 0.0 || iszero(x)) && return 0.0
        return 1.0
    end
    if isinf(y)
        x == Inf && return y > 0 ? 0.25 : -0.25
        x == -Inf && return y > 0 ? 0.75 : -0.75
        return y > 0 ? 0.5 : -0.5
    end
    x == Inf && return 0.0
    x == -Inf && return y > 0 ? 1.0 : -1.0
    iszero(x) && return y > 0 ? 0.5 : -0.5
    # diagonals give exactly ±1/4, ±3/4 — dyadic, so peel (Niven-termination):
    y == x && return y > 0 ? 0.25 : -0.75
    y == -x && return y > 0 ? 0.75 : -0.25
    _encl_divpi(atan, atan, y, x)
end
