# ===== project.jl — the projection engine (design §5)
# One generic pipeline: ωRoundToPrecision → ωSaturate → ωEncode, shared by the
# Float64 exact path, the BigFloat exact path, and the interval oracle.
#
# `sticky ∈ {-1,0,+1}` is symbolic information: the true value equals the carried
# value plus sticky·ε for an infinitesimal ε > 0. It lets a correctly-rounded
# (MPFR) endpoint stand in for an irrational true value without losing the
# direction information that directed/tie/stochastic rounding needs (design §5.1/§9.2).

# ---- Rounded value in canonical integer form (design §5.2)
const KIND_FIN = 0x00; const KIND_NAN = 0x01; const KIND_PINF = 0x02; const KIND_NINF = 0x03
"""Result of ωRoundToPrecision in canonical integer form: a kind tag plus
`sign · S · 2^Q` for finite results (design §5.2). The saturation stage consumes
this without ever re-touching the carrier value."""
struct Rounded
    kind::UInt8
    sign::Int8       # ±1 (finite, nonzero); +1 for zero
    S::Int64         # significand, 0 ≤ S ≤ 2^P (2^P = carry, normalized below)
    Q::Int64         # value = sign · S · 2^Q
end
const HUGEQ = Int64(1) << 40   # "finite but beyond every format" exponent sentinel

@inline function _floorexp(a)      # ⌊log₂ a⌋ for a > 0 finite (Float64 or BigFloat)
    Int64(Base.exponent(a))
end

"""
    round_to_precision(P, B, μ, X, R, sticky) -> Rounded

ωRoundToPrecision (draft §4.7.4), exact for exact carriers; `R` supplies the
stochastic random bits (0 ≤ R < 2^N), ignored for pure modes.
"""
# Float64 carrier (bitops plan K3): mask-extracted guard/round/sticky. Specials,
# zeros, sticky-zeros, and (Convert-only-reachable) subnormal Float64 inputs bail
# to the generic core; the exhaustive equivalence gate pins bit ≡ generic.
function round_to_precision(P::Int, B::Int, μ::RoundingMode3109, X::Float64, R::Int, sticky::Int)
    (isnan(X) | isinf(X) | iszero(X)) && return _rtp_core(P, B, μ, X, R, sticky)
    ((reinterpret(UInt64, X) >> 52) & 0x7ff) == 0x000 && return _rtp_core(P, B, μ, X, R, sticky)
    _rtp_f64(P, B, μ, X, R, sticky)
end

const _HALF128 = UInt128(1) << 127
# sticky-aware fixed-point comparisons: true fraction = νfix·2^-128 (+ lost bits) + νs·ε
@inline _bgt(ν::UInt128, lost::Bool, νs::Int, c::UInt128) = ν > c || (ν == c && (lost | (νs > 0)))
@inline _beq(ν::UInt128, lost::Bool, νs::Int, c::UInt128) = ν == c && !lost && νs == 0
@inline _bge(ν::UInt128, lost::Bool, νs::Int, c::UInt128) = ν > c || (ν == c && (lost | (νs >= 0)))
@inline function _bfloorscaled(ν::UInt128, lost::Bool, νs::Int, N::Int)
    k = Int64(ν >> (128 - N))
    (νs < 0 && !lost && (ν << N) == 0 && k > 0) && return k - 1   # exact grid hit, true below
    k
end
@inline function _brnite(ν::UInt128, lost::Bool, νs::Int, N::Int)
    k = Int64(ν >> (128 - N))
    rbit = (ν >> (127 - N)) & 0x1
    rbit == 0 && return k                                # fr < ½ regardless of ε
    low = ν << (N + 1)
    (low != 0 || lost || νs > 0) && return k + 1         # fr > ½
    νs < 0 && return k                                   # fr = ½ − ε
    isodd(k) ? k + 1 : k                                 # exact tie → even
end

# ROUND-AWAY PREDICATES, fixed-point family: these are the bit-path twins of the
# generic `_roundaway` family below — one predicate per mode, identical semantics
# on (fraction, sticky) evidence, differing only in carrier ((ν::UInt128, lost)
# here vs exact float ν there). Any semantic edit must land in BOTH families;
# the exhaustive bit ≡ generic equivalence gate in the suite enforces this.
@inline _rab(::TowardZero, ν, lost, νs, Sfl, Q, B, P, sign, R) = false
@inline _rab(::TowardPositive, ν, lost, νs, Sfl, Q, B, P, sign, R) =
    _bgt(ν, lost, νs, UInt128(0)) && sign > 0
@inline _rab(::TowardNegative, ν, lost, νs, Sfl, Q, B, P, sign, R) =
    _bgt(ν, lost, νs, UInt128(0)) && sign < 0
@inline _rab(::NearestTiesToAway, ν, lost, νs, Sfl, Q, B, P, sign, R) = _bge(ν, lost, νs, _HALF128)
@inline _rab(::NearestTiesToEven, ν, lost, νs, Sfl, Q, B, P, sign, R) =
    _bgt(ν, lost, νs, _HALF128) || (_beq(ν, lost, νs, _HALF128) && !_codeiseven(Sfl, Q, B, P))
@inline _rab(::ToOdd, ν, lost, νs, Sfl, Q, B, P, sign, R) =
    _bgt(ν, lost, νs, UInt128(0)) && _codeiseven(Sfl, Q, B, P)
@inline _rab(::StochasticA{N}, ν, lost, νs, Sfl, Q, B, P, sign, R) where {N} =
    _bfloorscaled(ν, lost, νs, N) + R >= Int64(1) << N
@inline _rab(::StochasticB{N}, ν, lost, νs, Sfl, Q, B, P, sign, R) where {N} =
    _bfloorscaled(ν, lost, νs, N + 1) + (2R + 1) >= Int64(1) << (N + 1)
@inline _rab(::StochasticC{N}, ν, lost, νs, Sfl, Q, B, P, sign, R) where {N} =
    _brnite(ν, lost, νs, N) + R >= Int64(1) << N

function _rtp_f64(P::Int, B::Int, μ::RoundingMode3109, X::Float64, R::Int, sticky::Int)
    sign = Int8(X > 0.0 ? 1 : -1)
    u = reinterpret(UInt64, abs(X))
    e = Int(u >> 52) - 1023                                  # normal input guaranteed
    m = (u & ((UInt64(1) << 52) - 1)) | (UInt64(1) << 52)    # 53-bit significand
    Q = Int64(max(e, 1 - B) - P + 1)
    d = e - Int(Q)                                           # units-bit position; d ≤ P−1 ≤ 7
    local Sfl::Int64, νfix::UInt128
    lost = false
    if d >= 0
        t = 52 - d                                           # t ∈ [45, 52]
        Sfl = Int64(m >> t)
        νfix = UInt128(m & ((UInt64(1) << t) - 1)) << (128 - t)
    else
        Sfl = 0
        t = 52 - d                                           # > 52
        if t <= 128
            νfix = UInt128(m) << (128 - t)
        else
            sh = t - 128
            if sh >= 64
                νfix = UInt128(0); lost = true               # m ≠ 0 here
            else
                νfix = UInt128(m >> sh)
                lost = (m & ((UInt64(1) << sh) - 1)) != 0
            end
        end
    end
    νs = sticky == 0 ? 0 : sticky * Int(sign)
    # step-down for "true value just below an exact dyadic": mirrors the identical
    # block in _rtp_core (binade edge borrows from Q; otherwise decrement S)
    if νs < 0 && νfix == UInt128(0) && !lost                 # true just below the dyadic
        if Sfl == Int64(1) << (P - 1) && Q > Int64(2 - B - P)
            Q -= 1
            Sfl = (Int64(1) << P) - 1
        else
            Sfl -= 1                                          # Sfl ≥ 1 here (X normal ⇒ d ≥ 0 path)
        end
        νfix = typemax(UInt128); lost = false; νs = -1        # ν = 1⁻ in fixed point
    end
    away = _rab(μ, νfix, lost, νs, Sfl, Q, B, P, sign, R)
    S = away ? Sfl + 1 : Sfl
    if S == Int64(1) << P
        S = Int64(1) << (P - 1); Q += 1
    end
    S == 0 && return Rounded(KIND_FIN, Int8(1), 0, 0)
    Rounded(KIND_FIN, sign, S, Q)
end

# Float128 carrier (Float128 revision plan, Site A): every operation _rtp_core uses
# (exponent, exact power-of-two ldexp scaling, floor, exact fraction subtraction,
# comparisons against dyadic constants) is exact on Float128, and the ν-exactness
# argument transfers verbatim (S̃ < 2^P after exact scaling). ν granularity on a
# Float128 carrier is ≥ 2^(P−113), so every stochastic grid up to N = 60 stays
# aligned — no mode restriction. No precision ceremony needed.
round_to_precision(P::Int, B::Int, μ::RoundingMode3109, X::Float128, R::Int, sticky::Int) =
    _rtp_core(P, B, μ, X, R, sticky)
# BigFloat carriers may exceed the ambient MPFR default precision; ldexp/floor/ν
# arithmetic must run at (at least) the operand's precision or bits are silently
# truncated. Function barrier sets it for the whole finite path.
round_to_precision(P::Int, B::Int, μ::RoundingMode3109, X::BigFloat, R::Int, sticky::Int) =
    setprecision(() -> _rtp_core(P, B, μ, X, R, sticky), BigFloat, Base.precision(X) + 8)

function _rtp_core(P::Int, B::Int, μ::RoundingMode3109, X, R::Int, sticky::Int)
    if isnan(X)
        return Rounded(KIND_NAN, Int8(1), 0, 0)
    elseif isinf(X)
        if X > 0
            sticky < 0 && return Rounded(KIND_FIN, Int8(1), (Int64(1) << P) - 1, HUGEQ)
            return Rounded(KIND_PINF, Int8(1), 0, 0)
        else
            sticky > 0 && return Rounded(KIND_FIN, Int8(-1), (Int64(1) << P) - 1, HUGEQ)
            return Rounded(KIND_NINF, Int8(-1), 0, 0)
        end
    end
    local sign::Int8, Sfl::Int64, Q::Int64, ν, νs::Int
    if iszero(X)
        sticky == 0 && return Rounded(KIND_FIN, Int8(1), 0, 0)
        sign = Int8(sticky > 0 ? 1 : -1)
        Q = Int64(2 - B - P)
        Sfl = 0
        ν = zero(float(typeof(X))); νs = 1                  # |true| ∈ (0, ε)
    else
        sign = Int8(X > 0 ? 1 : -1)
        a = abs(X)
        e = _floorexp(a)
        Q = max(e, Int64(1 - B)) - P + 1
        S̃ = ldexp(a, -Int(Q))                               # exact power-of-two scaling
        f = floor(S̃)
        Sfl = Int64(f)
        ν = S̃ - f                                           # exact (see design §5.2)
        νs = sticky == 0 ? 0 : (sticky * Int(sign))          # sticky in |value| direction
        if νs < 0 && iszero(ν)                               # true magnitude just below S̃·2^Q
            if Sfl == Int64(1) << (P - 1) && Q > Int64(2 - B - P)   # binade edge
                Q -= 1
                Sfl = (Int64(1) << P) - 1
            else
                Sfl -= 1                                     # Sfl ≥ 1 here (X ≠ 0)
            end
            ν = one(ν); νs = -1                              # ν = 1⁻
        end
    end
    away = _roundaway(μ, ν, νs, Sfl, Q, B, P, sign, R)
    S = away ? Sfl + 1 : Sfl
    if S == Int64(1) << P                                    # next-binade carry (NOTE 4)
        S = Int64(1) << (P - 1); Q += 1
    end
    S == 0 && return Rounded(KIND_FIN, Int8(1), 0, 0)
    Rounded(KIND_FIN, sign, S, Q)
end

# sticky-aware ν comparisons: true fraction = ν + νs·ε
@inline _νgt(ν, νs, c) = ν > c || (ν == c && νs > 0)
@inline _νeq(ν, νs, c) = ν == c && νs == 0
@inline _νge(ν, νs, c) = ν > c || (ν == c && νs >= 0)
@inline function _νfloorscaled(ν, νs, N)                     # ⌊(ν+νs·ε)·2^N⌋
    t = ldexp(ν, N); ft = floor(t)
    (νs < 0 && t == ft && ft > 0) && return Int64(ft) - 1
    Int64(ft)
end
@inline function _νrnite(ν, νs, N)                            # RNITE((ν+νs·ε)·2^N)
    t = ldexp(ν, N); ft = floor(t); fr = t - ft
    if fr > 0.5 || (fr == 0.5 && (νs > 0 || (νs == 0 && isodd(Int64(ft)))))
        return Int64(ft) + 1
    end
    Int64(ft)
end

@inline function _codeiseven(Sfl::Int64, Q::Int64, B::Int, P::Int)
    P > 1 ? iseven(Sfl) : (Sfl == 0 || iseven(Q + B))
end

# ROUND-AWAY PREDICATES, generic-carrier family (exact float ν): twins of the
# fixed-point `_rab` family above — see the note there; edits must land in both.
@inline _roundaway(::TowardZero, ν, νs, Sfl, Q, B, P, sign, R) = false
@inline _roundaway(::TowardPositive, ν, νs, Sfl, Q, B, P, sign, R) = _νgt(ν, νs, 0) && sign > 0
@inline _roundaway(::TowardNegative, ν, νs, Sfl, Q, B, P, sign, R) = _νgt(ν, νs, 0) && sign < 0
@inline _roundaway(::NearestTiesToAway, ν, νs, Sfl, Q, B, P, sign, R) = _νge(ν, νs, 0.5)
@inline _roundaway(::NearestTiesToEven, ν, νs, Sfl, Q, B, P, sign, R) =
    _νgt(ν, νs, 0.5) || (_νeq(ν, νs, 0.5) && !_codeiseven(Sfl, Q, B, P))
@inline _roundaway(::ToOdd, ν, νs, Sfl, Q, B, P, sign, R) =
    _νgt(ν, νs, 0) && _codeiseven(Sfl, Q, B, P)
@inline _roundaway(::StochasticA{N}, ν, νs, Sfl, Q, B, P, sign, R) where {N} =
    _νfloorscaled(ν, νs, N) + R >= Int64(1) << N
@inline _roundaway(::StochasticB{N}, ν, νs, Sfl, Q, B, P, sign, R) where {N} =
    _νfloorscaled(ν, νs, N + 1) + (2R + 1) >= Int64(1) << (N + 1)
@inline _roundaway(::StochasticC{N}, ν, νs, Sfl, Q, B, P, sign, R) where {N} =
    _νrnite(ν, νs, N) + R >= Int64(1) << N

# ---- dyadic comparison for saturation (design §5.3)
@inline _nbits(S::Int64) = 64 - leading_zeros(UInt64(S))
# compare s·S·2^Q against extremal datum m (as Float64, exact); |result| may exceed Float64 range,
# so compare in (sign, msb-position, aligned significand) space.
@inline function _cmp_rounded_datum(sign::Int8, S::Int64, Q::Int64, m::Float64)
    S == 0 && return m == 0.0 ? 0 : (m > 0 ? -1 : 1)
    m == 0.0 && return Int(sign)
    sm = m > 0 ? 1 : -1
    Int(sign) != sm && return Int(sign) < sm ? -1 : 1
    am = abs(m)
    (mS, mE, _) = Base.decompose(am)                      # am = mS · 2^mE exactly
    mS64 = Int64(mS); mE64 = Int64(mE)
    p1 = Q + _nbits(S); p2 = mE64 + _nbits(mS64)          # msb positions
    if p1 != p2
        c = p1 < p2 ? -1 : 1
    else
        Δ = Q - mE64                                       # |Δ| ≤ 64 when positions equal
        a1 = Δ >= 0 ? (Int128(S) << Δ) : Int128(S)
        a2 = Δ >= 0 ? Int128(mS64) : (Int128(mS64) << (-Δ))
        c = a1 < a2 ? -1 : (a1 > a2 ? 1 : 0)
    end
    Int(sign) > 0 ? c : -c
end

# Extremal magnitude in canonical (S, Q) integer form — a pure function of the type
# parameters, so it constant-folds (bitops plan K4). Verified by enumeration against
# Base.decompose(maxfinite_datum(T)) in the gate suite.
@inline function _extremal_SQ(::Type{T}) where {K,P,SGN,EXT,T<:Binary{K,P,SGN,EXT}}
    c = codepoint(MaxFiniteOf(T))
    Tt = c & UInt8((1 << (P - 1)) - 1)
    Eb = Int(c >> (P - 1))
    B = expbias(T)
    S = Eb == 0 ? Int64(Tt) : Int64(Tt) + (Int64(1) << (P - 1))
    Q = Int64((Eb == 0 ? 1 : Eb) - B - (P - 1))
    (S, Q)
end

# ---- ωSaturate (draft §4.7.5): the pattern rows, specialized by (Sat, Round, Σ, Δ).
# Range tests are two-integer comparisons against folded type constants: rounded
# canonical forms and the extremal form share one (Q, then S) lexicographic order
# (subnormals and the lowest normal binade share Q = 2−B−P), signed formats are
# sign–magnitude symmetric (|Mlo| = Mhi), unsigned underflow is simply sign < 0,
# and the HUGEQ sentinel exceeds every Q_hi. Returns :asis/:mhi/:mlo/:pinf/:ninf/:nan.
function saturate(::Type{T}, ρ::ProjSpec{RM,SM}, r::Rounded) where
        {K,P,SGN,EXT,T<:Binary{K,P,SGN,EXT},RM,SM}
    r.kind == KIND_NAN && return :nan
    over = false; under = false
    if r.kind == KIND_FIN
        r.S == 0 && return :asis                                     # zero always in range
        Shi, Qhi = _extremal_SQ(T)
        overmag = (r.Q > Qhi) | ((r.Q == Qhi) & (r.S > Shi))
        over = overmag & (r.sign > 0)
        under = SGN ? (overmag & (r.sign < 0)) : (r.sign < 0)
        (!over & !under) && return :asis                             # row 2
    end
    sat = SM(); isP = r.kind == KIND_PINF; isN = r.kind == KIND_NINF
    if sat isa SatFinite
        isP && return :mhi
        isN && return :mlo
        return under ? :mlo : :mhi
    elseif sat isa SatPropagate
        isP && return EXT ? :pinf : :mhi
        isN && return (SGN && EXT) ? :ninf : :mlo
        return under ? :mlo : :mhi
    else # SatNone
        rm = RM()
        if r.kind == KIND_FIN
            if over && (rm isa TowardZero || rm isa TowardNegative)
                return :mhi
            elseif under && (rm isa TowardZero || rm isa TowardPositive)
                return :mlo
            end
            if EXT
                over && return :pinf
                if under
                    return SGN ? :ninf : :nan
                end
            end
            return :nan                                              # Finite catch-all
        else
            if EXT
                isP && return :pinf
                isN && return SGN ? :ninf : :nan
            end
            return :nan                                              # Finite catch-all
        end
    end
end

# ---- ωProject (draft §4.7.3): full pipeline to a code point
"""
    project(T, ρ, X; R=0, sticky=0) -> T

Project a closed-extended-real carrier value into format `T` under ρ.
`X::Float64` must be *exact* (the Float64 fast path); `X::BigFloat` likewise
(or an interval endpoint with `sticky ∈ {-1,+1}`). `R` is the stochastic draw.
"""
function project(::Type{T}, ρ::ProjSpec{RM,SM}, X; R::Int=0, sticky::Int=0) where
        {K,P,SGN,EXT,T<:Binary{K,P,SGN,EXT},RM,SM}
    isnan(X) && return rawvalue(T, nan_code(T))                       # ωProject row 1
    B = expbias(T)
    r = round_to_precision(P, B, RM(), X, R, sticky)
    a = saturate(T, ρ, r)
    a === :asis && return rawvalue(T, encode(T, Int(r.sign), r.S, r.Q))
    a === :mhi  && return MaxFiniteOf(T)
    a === :mlo  && return MinFiniteOf(T)
    a === :pinf && return rawvalue(T, posinf_code(T))
    a === :ninf && return rawvalue(T, neginf_code(T))
    return rawvalue(T, nan_code(T))
end

"""
    project_interval(T, ρ, f; R=0, maxprec=4096) -> T

Project the true value of `f`, where `f(prec)` returns an enclosure
`(d::BigFloat, u::BigFloat)` with the true value in `[d, u]` (endpoints from
MPFR directed rounding). Rigor (design §9.2): if d == u the value is exact —
project it. Otherwise the true value lies in the *open* interval (d, u)
(directed correct rounding of a non-representable value moves strictly);
projection is monotone in the value for every mode at fixed R, so if
project(d, sticky=+1) == project(u, sticky=-1) that common code is the answer;
otherwise a projection grid point sits inside the interval and precision escalates.
"""
function project_interval(::Type{T}, ρ::ProjSpec, f; R::Int=0, maxprec::Int=4096) where {T<:Binary}
    prec = 256
    while true
        d, u = f(prec)
        if isequal(d, u)
            return project(T, ρ, d; R, sticky=0)
        end
        cd = project(T, ρ, d; R, sticky=+1)
        cu = project(T, ρ, u; R, sticky=-1)
        codepoint(cd) == codepoint(cu) && return cd
        prec >= maxprec && error("project_interval: unresolved at $maxprec bits (op enclosure too wide or true value pathologically near a grid point)")
        prec *= 2
    end
end
