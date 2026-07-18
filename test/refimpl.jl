# Independent reference implementation of ωRoundToPrecision/ωSaturate/ωEncode using
# Rational{BigInt} arithmetic — shares no code with the engine under test.
using P3109
using P3109: nan_code, posinf_code, neginf_code, rawvalue, codepoint

const RQ = Rational{BigInt}

function ilog2rat(x::RQ)
    @assert x > 0
    e = ndigits(numerator(x); base=2) - ndigits(denominator(x); base=2)
    p = big(2)^abs(e); v = e >= 0 ? RQ(p) : RQ(1, p)
    while v > x
        e -= 1; v //= 2
    end
    while 2v <= x
        e += 1; v *= 2
    end
    e
end
pow2r(q::Integer) = q >= 0 ? RQ(big(2)^q) : RQ(1, big(2)^(-q))

# margin: for approximate inputs assert ν is far from every decision constant
function refround(P::Int, B::Int, mode, R::Int, X::RQ; approx::Bool=false)
    X == 0 && return (0, big(0), 0)
    sgn = X > 0 ? 1 : -1
    a = abs(X)
    e = ilog2rat(a)
    Q = max(e, 1 - B) - P + 1
    S̃ = a / pow2r(Q)
    fl = floor(BigInt, S̃)
    ν = S̃ - fl
    if approx
        for c in (RQ(0), RQ(1, 2), RQ(1))
            @assert ν == c || abs(ν - c) > RQ(1, big(2)^2500) "ν too close to a decision boundary"
        end
    end
    codeiseven = P > 1 ? iseven(fl) : (fl == 0 || iseven(Q + B))
    away = if mode isa TowardZero
        false
    elseif mode isa TowardPositive
        ν > 0 && sgn > 0
    elseif mode isa TowardNegative
        ν > 0 && sgn < 0
    elseif mode isa NearestTiesToAway
        ν >= RQ(1, 2)
    elseif mode isa NearestTiesToEven
        ν > RQ(1, 2) || (ν == RQ(1, 2) && !codeiseven)
    elseif mode isa ToOdd
        ν > 0 && codeiseven
    elseif mode isa StochasticA
        N = P3109.nrandbits(mode)
        floor(BigInt, ν * big(2)^N) + R >= big(2)^N
    elseif mode isa StochasticB
        N = P3109.nrandbits(mode)
        floor(BigInt, ν * big(2)^(N + 1)) + (2R + 1) >= big(2)^(N + 1)
    elseif mode isa StochasticC
        N = P3109.nrandbits(mode)
        t = ν * big(2)^N; ft = floor(BigInt, t); fr = t - ft
        rnite = (fr < RQ(1, 2) || (fr == RQ(1, 2) && iseven(ft))) ? ft : ft + 1
        rnite + R >= big(2)^N
    else
        error("mode $mode")
    end
    S = fl + (away ? 1 : 0)
    (sgn, S, Q)
end

datum_rq(v) = RQ(decode(v))   # decode is exact in Float64 (asserted elsewhere)

function refproject(::Type{T}, ρ::ProjSpec, X::Union{RQ,Float64}; R::Int=0, approx::Bool=false) where {T<:Binary}
    mode = P3109.roundingmode(ρ); sat = P3109.saturationmode(ρ)
    K = bitwidth(T); P = precision(T); B = P3109.expbias(T)
    SGN = P3109.issigned(T); EXT = P3109.isextended(T)
    nanc = nan_code(T)
    X isa Float64 && isnan(X) && return nanc
    if X isa Float64 && isinf(X)
        Z = X
    else
        Xr = X isa Float64 ? RQ(X) : X
        sgn, S, Q = refround(P, B, mode, R, Xr; approx)
        Z = S == 0 ? RQ(0) : sgn * RQ(S) * pow2r(Q)
        Z = Z::RQ
    end
    mhi = datum_rq(MaxFiniteOf(T)); mlo = datum_rq(MinFiniteOf(T))
    # ωSaturate rows, in draft order
    out = if Z isa RQ && mlo <= Z <= mhi
        Z
    elseif sat isa SatFinite
        Z isa Float64 ? (Z > 0 ? mhi : mlo) : (Z < mlo ? mlo : mhi)
    elseif sat isa SatPropagate
        if Z isa Float64            # ±Inf
            Z > 0 ? (EXT ? :pinf : mhi) : ((SGN && EXT) ? :ninf : mlo)
        else
            Z < mlo ? mlo : mhi
        end
    else # SatNone
        if Z isa RQ && Z > mhi && (mode isa TowardZero || mode isa TowardNegative)
            mhi
        elseif Z isa RQ && Z < mlo && (mode isa TowardZero || mode isa TowardPositive)
            mlo
        elseif EXT && ((Z isa Float64 && Z > 0) || (Z isa RQ && Z > mhi))
            :pinf
        elseif EXT && ((Z isa Float64 && Z < 0) || (Z isa RQ && Z < mlo))
            SGN ? :ninf : :nan
        else
            :nan
        end
    end
    out === :pinf && return posinf_code(T)
    out === :ninf && return neginf_code(T)
    out === :nan && return nanc
    # encode by exhaustive search over the datum set (independent of ωEncode)
    for c in 0x00:UInt8((1 << K) - 1)
        c == nanc && continue
        v = rawvalue(T, c)
        isfinite(v) || continue
        datum_rq(v) == out && return c
    end
    error("refproject: value $out not in datum set of $(T)")
end

# reference ω semantics for Add/Multiply/Divide on datums (exact rational)
function refop2(name::Symbol, x, y)::Union{RQ,Float64}   # Float64 = special (NaN/±Inf)
    dx, dy = decode(x), decode(y)
    (isnan(dx) || isnan(dy)) && return NaN
    if name === :Add
        if isinf(dx) || isinf(dy)
            (isinf(dx) && isinf(dy) && sign(dx) != sign(dy)) && return NaN
            return isinf(dx) ? dx : dy
        end
        return RQ(dx) + RQ(dy)
    elseif name === :Multiply
        if isinf(dx) || isinf(dy)
            (dx == 0 || dy == 0) && return NaN
            return sign(dx) * sign(dy) > 0 ? Inf : -Inf
        end
        return RQ(dx) * RQ(dy)
    elseif name === :Divide
        dy == 0 && return NaN
        if isinf(dx)
            isinf(dy) && return NaN
            return sign(dx) * sign(dy) > 0 ? Inf : -Inf
        end
        isinf(dy) && return 0.0
        dx == 0 && return 0.0
        return RQ(dx) / RQ(dy)
    end
    error("refop2 $name")
end
