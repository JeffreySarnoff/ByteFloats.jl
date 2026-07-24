# ===== juliacompat.jl — the Base register (design §10.2)
#
# Base-function veneers for every draft operation that has a Base counterpart.
# Every method here is exactly one same-format spec-register call under the
# session default projection; there is no third semantics.
#
# The mapping is a declarative table partitioned over the op lists in
# ops_scalar.jl (_UNARY_OPS / _BINARY_OPS / _TERNARY_OPS): each op is either
# mapped below or listed in _NO_BASE_COUNTERPART with the reason class. The
# suite asserts that partition is exhaustive, so adding an op to a list forces
# a decision here — the same non-divergence mechanism as the registry.
#
# Scope decisions:
#   · Same-format operands only. Promotion between distinct Binary formats is
#     deliberately absent (design §2.4) — mixing formats is an explicit
#     `Convert`, not a silent widening. Binary ⋄ Float64/32/16/Integer promotes
#     to Float64 through the rules in formats.jl, as before.
#   · The extremum veneers map Base's IEEE semantics: `min`/`max` propagate NaN,
#     which is exactly draft Minimum/Maximum. The NaN-ignoring (Number),
#     magnitude, and finite families have no Base spelling — call them by their
#     draft names.

# ---- unary: plain Base name ⇒ draft op
const _BASE_UNARY = (
    :Abs => :abs,   :Recip => :inv,  :Sqrt => :sqrt,
    :Exp => :exp,   :Exp2 => :exp2,  :ExpMinusOne => :expm1,
    :Log => :log,   :Log2 => :log2,  :LogOnePlus => :log1p,
    :Sin => :sin,   :Cos => :cos,    :Tan => :tan,
    :ArcSin => :asin, :ArcCos => :acos, :ArcTan => :atan,
    :Sinh => :sinh, :Cosh => :cosh,  :Tanh => :tanh,
    :ArcSinh => :asinh, :ArcCosh => :acosh, :ArcTanh => :atanh,
    :SinPi => :sinpi, :CosPi => :cospi, :TanPi => :tanpi,
)
for (op, bf) in _BASE_UNARY
    @eval Base.$bf(x::Binary) = $op(x)
end

# ---- operators and irregular spellings (arity or argument order differs)
Base.:-(x::Binary) = Negate(x)                                  # unary minus
Base.:+(x::T, y::T) where {T<:Binary} = Add(x, y)
Base.:-(x::T, y::T) where {T<:Binary} = Subtract(x, y)
Base.:*(x::T, y::T) where {T<:Binary} = Multiply(x, y)
Base.:/(x::T, y::T) where {T<:Binary} = Divide(x, y)
Base.atan(y::T, x::T) where {T<:Binary} = ArcTan2(y, x)         # Base's (y, x) order
const _BASE_OPERATOR = (:Negate, :Add, :Subtract, :Multiply, :Divide, :ArcTan2)

# ---- binary: plain Base name ⇒ draft op, same format
const _BASE_BINARY = (
    :CopySign => :copysign, :Hypot => :hypot,
    :Maximum => :max, :Minimum => :min,
)
for (op, bf) in _BASE_BINARY
    @eval Base.$bf(x::T, y::T) where {T<:Binary} = $op(x, y)
end

# ---- ternary
Base.fma(x::T, y::T, z::T) where {T<:Binary} = FMA(x, y, z)
Base.muladd(x::T, y::T, z::T) where {T<:Binary} = FMA(x, y, z)  # one rounding, like fma
Base.clamp(x::T, lo::T, hi::T) where {T<:Binary} = Clamp(x, lo, hi)
const _BASE_TERNARY = (:FMA, :Clamp)

# ---- composites: Base functions with no single draft op, defined componentwise
# from draft ops so each component carries its defined result
Base.sincos(x::Binary) = (Sin(x), Cos(x))
Base.sincospi(x::Binary) = (SinPi(x), CosPi(x))
Base.minmax(x::T, y::T) where {T<:Binary} = (Minimum(x, y), Maximum(x, y))

# ---- draft ops with no Base counterpart — call these by their draft names.
# Kept as one auditable set; the suite asserts mapped ∪ this == every op list.
const _NO_BASE_COUNTERPART = (
    :RSqrt, :Softplus,                                   # no Base function
    :ArcSinPi, :ArcCosPi, :ArcTanPi, :ArcTan2Pi,         # Base has sinpi-family only
    :MaximumNumber, :MinimumNumber,                      # Base has no NaN-ignoring pair
    :MaximumMagnitude, :MinimumMagnitude,
    :MaximumMagnitudeNumber, :MinimumMagnitudeNumber,
    :MinimumFinite, :MaximumFinite,
    :FAA,                                                # no Base fused add-add
)
