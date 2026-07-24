# ===== formats.jl — the Binary{K,P,SGN,EXT} type, Group M, naming, Base API (design §2, §3.1)

"""
    Binary{K,P,SGN,EXT} <: AbstractFloat

A P3109 floating-point value: a code point in the format with bitwidth `K ∈ 3:8`,
precision `P`, signedness `SGN::Bool` (`true` = Signed), domain `EXT::Bool`
(`true` = Extended, i.e. the datum set includes infinities).

The code point occupies the low `K` bits of the payload byte; the high `8-K` bits
are maintained zero as a representation invariant. Construct raw code points with
`T(c::UInt8)` (validated code-point construction), `rawvalue(T, c)` (unchecked
kernel route), or `Binary{...}(Val(:code), x::UInt8)`; construct from numeric
values with `T(x::Real)` (default projection spec) or `Convert`. UInt8 is
the one argument type meaning *code point*; all other Reals mean *value*.
"""
struct Binary{K,P,SGN,EXT} <: AbstractFloat
    x::UInt8
    @inline function Binary{K,P,SGN,EXT}(::Val{:code}, x::UInt8) where {K,P,SGN,EXT}
        @boundscheck begin
            checkformat(K, P, SGN, EXT)
            Int(x) < (1 << K) || throw(ArgumentError("code point $x out of range for K=$K"))
        end
        new{K,P,SGN,EXT}(x)
    end
end

function checkformat(K, P, SGN, EXT)
    (K isa Int && P isa Int && SGN isa Bool && EXT isa Bool) ||
        throw(ArgumentError("Binary parameters must be (K::Int, P::Int, SGN::Bool, EXT::Bool)"))
    3 <= K <= 8 || throw(ArgumentError("bitwidth K=$K outside supported range 3:8"))
    if SGN
        0 < P < K || throw(ArgumentError("signed format requires 0 < P < K (got P=$P, K=$K)"))
    else
        0 < P <= K || throw(ArgumentError("unsigned format requires 0 < P ≤ K (got P=$P, K=$K)"))
    end
    nothing
end

"""Unsafe raw constructor used by kernels after invariants are established."""
@inline rawvalue(::Type{Binary{K,P,SGN,EXT}}, x::UInt8) where {K,P,SGN,EXT} =
    @inbounds Binary{K,P,SGN,EXT}(Val(:code), x)

"""
    (T::Type{<:Binary})(c::UInt8) -> T

Construct a value from its **code point** (validated: for K < 8, `c` must be
`< 2^K` or an `ArgumentError` is thrown; the range check folds away at K = 8 and
for constant codes). `Binary5p3sf(0x08) == Binary5p3sf(1.0)`.

`UInt8` is the *only* argument type with code-point semantics — every other
`Real` (including other `Integer`s) constructs by projecting the numeric value,
and `Convert` is numeric for all integers: `Binary8p4se(0x02)` is code point 2
(the second-smallest subnormal), while `Binary8p4se(2)` is the number 2.0.
`rawvalue(T, c)` remains the unchecked kernel-internal route.
"""
@inline (::Type{Binary{K,P,SGN,EXT}})(c::UInt8) where {K,P,SGN,EXT} =
    Binary{K,P,SGN,EXT}(Val(:code), c)

@inline Base.codepoint(v::Binary) = v.x   # extends Base.codepoint (Char); avoids export clash

# ---- Group M (meta) operations: pure functions of the type parameters (design §2.3)
"Format bitwidth K (3–8)."
bitwidth(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = K
Base.precision(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = P
"Whether the format is Signed (has a sign bit and negative datums)."
issigned(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = S
"Whether the format's domain is Extended (datum set includes infinities)."
isextended(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = E
"Exponent bias: 2^(K−P−1) signed, 2^(K−P) unsigned."
expbias(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = S ? (1 << (K - P - 1)) : (1 << (K - P))
"Width of the exponent field in bits: (K − signbit) − (P − 1)."
expbitwidth(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = (S ? K - 1 : K) - (P - 1)
"Trailing-significand width P − 1 (the stored fraction bits)."
trailingsigbits(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = P - 1

const BitwidthOf = bitwidth
const PrecisionOf = precision
const SignednessOf = issigned
const DomainOf = isextended
const ExponentBiasOf = expbias
const ExponentBitwidthOf = expbitwidth
const TrailingSignificandBitwidthOf = trailingsigbits

# Special code points (literals after constant folding)
@inline nan_code(::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    S ? UInt8(1 << (K - 1)) : UInt8((1 << K) - 1)
@inline posinf_code(::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    S ? UInt8((1 << (K - 1)) - 1) : UInt8((1 << K) - 2)   # meaningful only when E
@inline neginf_code(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = UInt8((1 << K) - 1)  # signed+E only
@inline signmask(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = UInt8(1 << (K - 1))

# Extremal *code points* (draft Group M returns format values).
# Largest finite magnitude code: the greatest code below the NaN/Inf slots.
#   signed·extended  : +Inf at 2^(K-1)-1        → maxfinite = 2^(K-1)-2
#   signed·finite    : NaN  at 2^(K-1)          → maxfinite = 2^(K-1)-1
#   unsigned·extended: NaN 2^K-1, +Inf 2^K-2    → maxfinite = 2^K-3
#   unsigned·finite  : NaN  at 2^K-1            → maxfinite = 2^K-2
@inline function MaxFiniteOf(T::Type{Binary{K,P,S,E}}) where {K,P,S,E}
    c = S ? (E ? UInt8((1 << (K - 1)) - 2) : UInt8((1 << (K - 1)) - 1)) :
            (E ? UInt8((1 << K) - 3)       : UInt8((1 << K) - 2))
    rawvalue(T, c)
end
@inline function MinFiniteOf(T::Type{Binary{K,P,S,E}}) where {K,P,S,E}
    S || return rawvalue(T, 0x00)                              # unsigned: 0
    rawvalue(T, codepoint(MaxFiniteOf(T)) | signmask(T))       # most negative finite
end
@inline MinPositiveOf(T::Type{<:Binary}) = rawvalue(T, 0x01)
@inline MaxSubnormalOf(T::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    rawvalue(T, UInt8((1 << (P - 1)) - 1))          # P=1 formats have no subnormals ⇒ code 0
@inline MinNormalOf(T::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    rawvalue(T, UInt8(1 << (P - 1)))

# Datum-valued companions (design §2.3): Float64 is the universal exact carrier
maxfinite_datum(T::Type{<:Binary}) = decode(MaxFiniteOf(T))
minfinite_datum(T::Type{<:Binary}) = decode(MinFiniteOf(T))

# ---- Value-argument forwarders for Group M and the extremal queries.
# The answer is a pure function of the type parameters, so a value carries it:
# `BitwidthOf(x)` ≡ `BitwidthOf(typeof(x))`. Restricted to `Binary` — a `where {T}`
# signature would claim every type in the language. Constant-folds to the same
# literal as the type-argument form.
for f in (:bitwidth, :issigned, :isextended, :expbias, :expbitwidth, :trailingsigbits,
          :MaxFiniteOf, :MinFiniteOf, :MinPositiveOf, :MaxSubnormalOf, :MinNormalOf)
    @eval @inline $f(x::Binary) = $f(typeof(x))
end
# `PrecisionOf` aliases `Base.precision`, so its forwarder must be spelled out.
@inline Base.precision(x::Binary) = precision(typeof(x))

# ---- Naming grid, draft §3.2: BinaryKpP + s|u + e|f
const _NAMED = Dict{Symbol,DataType}()
for K in 3:8, P in 1:K, S in (true, false), E in (true, false)
    S && P >= K && continue
    name = Symbol("Binary", K, "p", P, S ? "s" : "u", E ? "e" : "f")
    T = Binary{K,P,S,E}
    @eval const $name = $T
    _NAMED[name] = T
end
"""`formatname(T)` — the draft §3.2 name of a format type."""
formatname(::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    Symbol("Binary", K, "p", P, S ? "s" : "u", E ? "e" : "f")

# Print fully-instantiated formats by their draft name; anything else (UnionAlls,
# TypeVar-parameterized types met during stacktrace printing) defers to Base —
# a parametric `::Type{Binary{K,P,S,E}}` method here can be handed unbound
# static parameters by the printing machinery and crash (found by test).
function Base.show(io::IO, T::Type{<:Binary})
    if T isa DataType && length(T.parameters) == 4 &&
       T.parameters[1] isa Int && T.parameters[2] isa Int &&
       T.parameters[3] isa Bool && T.parameters[4] isa Bool
        print(io, formatname(T))
    else
        invoke(show, Tuple{IO,Type}, io, T)
    end
end
function Base.show(io::IO, v::Binary)
    T = typeof(v)
    print(io, formatname(T), "(")
    d = decode(v)
    isnan(d) ? print(io, "NaN") : print(io, d)
    print(io, " ≡ 0x", string(codepoint(v); base=16, pad=2), ")")
end

# ---- Base numeric API on the type (defined via Group M / decode; see also ops_scalar.jl)
Base.zero(T::Type{<:Binary}) = rawvalue(T, 0x00)
Base.zero(::T) where {T<:Binary} = zero(T)
Base.iszero(v::Binary) = codepoint(v) == 0x00
Base.typemax(T::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    E ? rawvalue(T, posinf_code(T)) : MaxFiniteOf(T)
Base.typemin(T::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    S ? (E ? rawvalue(T, neginf_code(T)) : MinFiniteOf(T)) : zero(T)
Base.floatmax(T::Type{<:Binary}) = MaxFiniteOf(T)
Base.floatmin(T::Type{<:Binary}) = MinNormalOf(T)
Base.eps(T::Type{Binary{K,P,S,E}}) where {K,P,S,E} = T(2.0^(1 - P))

Base.isnan(v::Binary) = codepoint(v) == nan_code(typeof(v))
function Base.isinf(v::Binary{K,P,S,E}) where {K,P,S,E}
    E || return false
    c = codepoint(v)
    c == posinf_code(typeof(v)) && return true
    (S && c == neginf_code(typeof(v))) && return true
    false
end
Base.isfinite(v::Binary) = !isnan(v) & !isinf(v)
function Base.signbit(v::Binary{K,P,S,E}) where {K,P,S,E}
    S || return false
    isnan(v) && return false
    codepoint(v) >= signmask(typeof(v)) && codepoint(v) != nan_code(typeof(v))
end
function issubnormal_3109(v::Binary{K,P,S,E}) where {K,P,S,E}
    (isnan(v) | isinf(v) | iszero(v)) && return false
    m = S ? (codepoint(v) & ~signmask(typeof(v))) : codepoint(v)
    m < UInt8(1 << (P - 1))
end
Base.issubnormal(v::Binary) = issubnormal_3109(v)

# `one` must be representable; it always is for these formats (biased exp fits by construction).
Base.one(T::Type{<:Binary}) = T(1.0)
Base.one(::T) where {T<:Binary} = one(T)

# Promotion (design §2.4): Binary ⋄ external float promotes to Float64 (the exact carrier);
# no automatic promotion between distinct Binary formats.
Base.promote_rule(::Type{<:Binary}, ::Type{Float64}) = Float64
Base.promote_rule(::Type{<:Binary}, ::Type{Float32}) = Float64
Base.promote_rule(::Type{<:Binary}, ::Type{Float16}) = Float64
Base.promote_rule(::Type{<:Binary}, ::Type{<:Integer}) = Float64
Base.Float64(v::Binary) = decode(v)
Base.Float32(v::Binary) = Float32(decode(v))   # exact: all K≤8 datums fit Float32
(::Type{T})(v::Binary) where {T<:Binary} = _convert_default(T, v)
(::Type{T})(x::Real) where {T<:Binary} = _convert_default(T, x)
# Disambiguates against Base's (::Type{T})(::Rational) where T<:AbstractFloat
# (found by Test.detect_ambiguities). Consistent with the Convert policy: inputs
# must be exactly projectable; Rationals are rejected rather than double-rounded.
(::Type{T})(x::Rational{S}) where {S,T<:Binary} =
    throw(ArgumentError("cannot exactly project a Rational; convert explicitly, e.g. $(T)(Float64(x)), and own the double rounding"))
Base.convert(::Type{T}, v::Binary) where {T<:Binary} = _convert_default(T, v)
Base.convert(::Type{T}, x::Real) where {T<:Binary} = _convert_default(T, x)
# _convert_default is defined in ops_scalar.jl (needs the projection engine).
