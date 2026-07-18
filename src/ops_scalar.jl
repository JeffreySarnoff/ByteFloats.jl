# ===== ops_scalar.jl — operation registry and scalar API (design §7.1/§10.2, architecture §6)
#
# The registry is the single source from which the package generates the
# spec-register functions (draft names, draft parameter order), the same-format
# convenience methods, and the Base-register veneers. Block/scaled variants,
# table enumeration, exhaustive test sets, and conformance() are generated from
# the same rows (blocks.jl, tables.jl, approx.jl).
#
# Evaluation contract (architecture §5.1): oracle.jl supplies
#     ωeval(::Val{op}, xs::Float64...) ::Union{Float64, BigExactF, EncloseF}
# and this file supplies `apply_op`, which finishes each result kind through the
# projection engine. `Float64` results are *exact* by construction; `BigExactF.f()`
# returns an exact BigFloat; `EncloseF.f(prec)` returns a directed (lo, hi) enclosure.

using Random: AbstractRNG, default_rng

# Runtime switch (Float128 revision plan §5): every Float128 path fronts a complete
# MPFR path with identical semantics, so disabling costs speed, never correctness.
# Set ENV["P3109_FLOAT128"] = "disable" before loading (read in __init__).
const _USE_FLOAT128 = Ref(true)
@inline _f128() = _USE_FLOAT128[]

"""Deferred exact result: `f()` returns the value as an exact `BigFloat`."""
struct BigExactF{F}
    f::F
end
"""
Deferred enclosure: `f(prec)` returns `(lo, hi)::NTuple{2,BigFloat}` bracketing the
true value (MPFR directed rounding — the rigorous ladder). `fq`, when present, is a
zero-argument Float128 estimator used as a pre-filter (plan Site D, Class E): the
true value is taken to lie within `|y|·2^-90` of `y = fq()` — a bound ≥ 2^18 slacker
than any published libquadmath transcendental error, discharged empirically by the
differential-build tests. Filter disagreement always falls through to `f`.
"""
struct EncloseF{F,G}
    f::F
    fq::G
end
EncloseF(f) = EncloseF(f, nothing)
"""
Correctly-rounded Float128 bracket (plan Site C, Class R): IEEE mandates correct
rounding for Float128 `/`, `sqrt`, `fma`, so a nearest-CR result of a value known
inexact brackets the truth in the *open* interval `(lo, hi) = (prevfloat(q), nextfloat(q))`
with no envelope assumption. `f` is the MPFR fallback for (near-impossible at P ≤ 8)
grid-straddling brackets.
"""
struct Enclose128F{F}
    lo::Float128
    hi::Float128
    f::F
end

const _F128_RELEXP = -90        # envelope exponent: E = |y|·2^-90

@inline _finish(::Type{fr}, ρ::ProjSpec, R::Int, v::Float64) where {fr<:Binary} =
    project(fr, ρ, v; R)
@inline _finish(::Type{fr}, ρ::ProjSpec, R::Int, v::Float128) where {fr<:Binary} =
    project(fr, ρ, v; R)
@inline _finish(::Type{fr}, ρ::ProjSpec, R::Int, b::BigExactF) where {fr<:Binary} =
    project(fr, ρ, b.f(); R)
function _finish(::Type{fr}, ρ::ProjSpec, R::Int, e::EncloseF) where {fr<:Binary}
    if e.fq !== nothing && _f128()
        y = e.fq()::Float128
        if isfinite(y) && !iszero(y)
            E = ldexp(abs(y), _F128_RELEXP)
            d = y - E
            u = y + E
            cd = project(fr, ρ, d; R, sticky=+1)
            cu = project(fr, ρ, u; R, sticky=-1)
            codepoint(cd) == codepoint(cu) && return cd
        end
        # non-finite/zero estimate or filter disagreement: rigorous ladder decides
    end
    project_interval(fr, ρ, e.f; R)
end
function _finish(::Type{fr}, ρ::ProjSpec, R::Int, z::Enclose128F) where {fr<:Binary}
    isequal(z.lo, z.hi) && return project(fr, ρ, z.lo; R)
    cd = project(fr, ρ, z.lo; R, sticky=+1)
    cu = project(fr, ρ, z.hi; R, sticky=-1)
    codepoint(cd) == codepoint(cu) && return cd
    project_interval(fr, ρ, z.f; R)
end

@inline function apply_op(op::Val, ::Type{fr}, ρ::ProjSpec, R::Int, xs::Float64...) where {fr<:Binary}
    res = ωeval(op, xs...)
    # bitops plan Phase 0(a): explicit fast split — Class-1/selection results are
    # Float64 for every ordinary input; keep the widened union off the hot path.
    # Justification: like-for-like measurement (both variants under identical
    # harness conditions) showed the split alone recovering 399 → 269 ns/elem;
    # see checkpoint.md "Resolution of the two flagged measurements".
    res isa Float64 && return project(fr, ρ, res; R)
    _finish_slow(fr, ρ, R, res)
end
@noinline _finish_slow(::Type{fr}, ρ::ProjSpec, R::Int, res) where {fr<:Binary} =
    _finish(fr, ρ, R, res)

# ---- stochastic draw plumbing (design §5.5)
# bitops plan Phase 0(b): the rng default is `nothing`, resolved to the task-local
# default only when a stochastic draw is actually taken.
# Justification is SEMANTIC only — pure-ρ calls never touch RNG state, and array
# kernels resolve the rng once per call instead of per element. The original
# performance justification is withdrawn: controlled A/B (checkpoint.md,
# "Resolution of the two flagged measurements") showed the previous eager
# `default_rng()` kwarg default cost ≈ nothing in specialized code (25.4 vs
# 26.5 ns/elem); the 1,347 ns reading that motivated it was dynamic keyword
# dispatch through a non-const global in the measurement harness, not this code.
const MaybeRNG = Union{Nothing,AbstractRNG}
@inline function _drawR(ρ::ProjSpec, rng::MaybeRNG, R::Union{Nothing,Int})
    isstochastic(ρ) || return 0
    N = nrandbits(ρ)
    if R === nothing
        r = rng === nothing ? default_rng() : rng
        return Int(rand(r, UInt64) & ((UInt64(1) << N) - 1))
    end
    0 <= R < (1 << N) || throw(ArgumentError("explicit R=$R outside 0:$(2^N - 1) for N=$N random bits"))
    return R
end

# ---- the operation registry (draft §5.4 substitution lists)
struct OpInfo
    name::Symbol
    arity::Int
    group::Symbol   # :A arithmetic, :B elementary, :C extremum/misc, :conv
end
const OP_REGISTRY = OpInfo[]
register_op!(name::Symbol, arity::Int, group::Symbol) = push!(OP_REGISTRY, OpInfo(name, arity, group))
opinfo(name::Symbol) = OP_REGISTRY[findfirst(o -> o.name === name, OP_REGISTRY)]

const _UNARY_OPS = (:Abs, :Negate, :Sqrt, :RSqrt, :Recip, :Exp, :Log, :ExpMinusOne,
    :LogOnePlus, :Exp2, :Log2, :Sin, :Cos, :Tan, :ArcSin, :ArcCos, :ArcTan,
    :Sinh, :Cosh, :Tanh, :ArcSinh, :ArcCosh, :ArcTanh,
    :SinPi, :CosPi, :TanPi, :ArcSinPi, :ArcCosPi, :ArcTanPi, :Softplus)
const _BINARY_OPS = (:CopySign, :Add, :Subtract, :Multiply, :Divide, :Hypot,
    :ArcTan2, :ArcTan2Pi, :Maximum, :Minimum, :MaximumNumber, :MinimumNumber,
    :MaximumMagnitude, :MinimumMagnitude, :MaximumMagnitudeNumber,
    :MinimumMagnitudeNumber, :MinimumFinite, :MaximumFinite)
const _TERNARY_OPS = (:FMA, :FAA, :Clamp)

for n in _UNARY_OPS;   register_op!(n, 1, n in (:Abs, :Negate) ? :A : :B); end
for n in _BINARY_OPS;  register_op!(n, 2, n in (:Add, :Subtract, :Multiply, :Divide, :CopySign) ? :A : :C); end
for n in _TERNARY_OPS; register_op!(n, 3, :A); end
register_op!(:Convert, 1, :conv)

# ---- generated spec register + same-format convenience methods
# Spec form follows the draft's parameterization order: Op(f_r, ρ, operands...).
for op in OP_REGISTRY
    op.name === :Convert && continue
    name = op.name; V = Val{name}
    if op.arity == 1
        @eval begin
            @inline function $name(fr::Type{<:Binary}, ρ::ProjSpec, x::Binary;
                                   rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
                apply_op($V(), fr, ρ, _drawR(ρ, rng, R), decode(x))
            end
            @inline $name(x::T; kw...) where {T<:Binary} = $name(T, default_projspec(T), x; kw...)
        end
    elseif op.arity == 2
        @eval begin
            @inline function $name(fr::Type{<:Binary}, ρ::ProjSpec, x::Binary, y::Binary;
                                   rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
                apply_op($V(), fr, ρ, _drawR(ρ, rng, R), decode(x), decode(y))
            end
            @inline $name(x::T, y::T; kw...) where {T<:Binary} = $name(T, default_projspec(T), x, y; kw...)
        end
    else
        @eval begin
            @inline function $name(fr::Type{<:Binary}, ρ::ProjSpec, x::Binary, y::Binary, z::Binary;
                                   rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
                apply_op($V(), fr, ρ, _drawR(ρ, rng, R), decode(x), decode(y), decode(z))
            end
            @inline $name(x::T, y::T, z::T; kw...) where {T<:Binary} = $name(T, default_projspec(T), x, y, z; kw...)
        end
    end
end

# ---- Convert (draft §4.9): the one op accepting external operands
"""
    Convert(fr, ρ, x) -> fr

Draft §4.9 Convert⟨f_x, f_r, ρ⟩. Accepts `Binary`, IEEE floats (widened exactly to
the Float64 carrier), `Integer` (exact via a sufficiently wide BigFloat), and
`BigFloat` (projected directly; the caller warrants the value is exact).
"""
@inline function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::Binary;
                         rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
    project(fr, ρ, decode(x); R=_drawR(ρ, rng, R))
end
@inline function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::Union{Float64,Float32,Float16};
                         rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
    project(fr, ρ, Float64(x); R=_drawR(ρ, rng, R))   # exact widening
end
function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::Integer;
                 rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
    b = BigFloat(x; precision=max(64, ndigits(x, base=2) + 8))       # exact
    project(fr, ρ, b; R=_drawR(ρ, rng, R))
end
function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::BigFloat;
                 rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
    project(fr, ρ, x; R=_drawR(ρ, rng, R))
end
Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::AbstractFloat; kw...) = Convert(fr, ρ, Float64(x); kw...)

# closes the constructor loop declared in formats.jl
@inline _convert_default(::Type{T}, v::T) where {T<:Binary} = v
@inline _convert_default(::Type{T}, x) where {T<:Binary} = Convert(T, default_projspec(T), x)

# ---- Base register (design §10.2): same-format, default-ρ veneers only.
# Every method is exactly one spec-register call; there is no third semantics.
Base.:+(x::T, y::T) where {T<:Binary} = Add(x, y)
Base.:-(x::T, y::T) where {T<:Binary} = Subtract(x, y)
Base.:*(x::T, y::T) where {T<:Binary} = Multiply(x, y)
Base.:/(x::T, y::T) where {T<:Binary} = Divide(x, y)
Base.:-(x::T) where {T<:Binary} = Negate(x)
Base.abs(x::Binary) = Abs(x)
Base.copysign(x::T, y::T) where {T<:Binary} = CopySign(x, y)
Base.inv(x::Binary) = Recip(x)
Base.sqrt(x::Binary) = Sqrt(x)
Base.fma(x::T, y::T, z::T) where {T<:Binary} = FMA(x, y, z)
Base.muladd(x::T, y::T, z::T) where {T<:Binary} = FMA(x, y, z)
Base.min(x::T, y::T) where {T<:Binary} = Minimum(x, y)
Base.max(x::T, y::T) where {T<:Binary} = Maximum(x, y)
Base.clamp(x::T, lo::T, hi::T) where {T<:Binary} = Clamp(x, lo, hi)
Base.hypot(x::T, y::T) where {T<:Binary} = Hypot(x, y)
Base.atan(y::T, x::T) where {T<:Binary} = ArcTan2(y, x)
for (bf, op) in ((:exp, :Exp), (:exp2, :Exp2), (:expm1, :ExpMinusOne),
                 (:log, :Log), (:log2, :Log2), (:log1p, :LogOnePlus),
                 (:sin, :Sin), (:cos, :Cos), (:tan, :Tan),
                 (:asin, :ArcSin), (:acos, :ArcCos), (:atan, :ArcTan),
                 (:sinh, :Sinh), (:cosh, :Cosh), (:tanh, :Tanh),
                 (:asinh, :ArcSinh), (:acosh, :ArcCosh), (:atanh, :ArcTanh),
                 (:sinpi, :SinPi), (:cospi, :CosPi), (:tanpi, :TanPi))
    @eval Base.$bf(x::Binary) = $op(x)
end
