# ===== ops_scalar.jl — operation registry and scalar API (design §7.1/§10.2, architecture §6)
#
# The registry is the single source from which the package generates the
# spec-register functions (draft names, draft parameter order), the same-format
# convenience methods, and the Base-register veneers. Block/scaled variants,
# table enumeration, exhaustive test sets, and conformance() are generated from
# the same rows (blocks.jl, tables.jl, approx.jl).
#
# Evaluation contract (architecture §5.1): oracle.jl supplies
#     ωeval(::Val{op}, xs::Float64...)
#         ::Union{Float64, Float128, BigExactF, EncloseF, Enclose128F}
# and this file supplies `apply_op`, which finishes each result kind through the
# projection engine. `Float64`/`Float128` results are *exact* by construction;
# `BigExactF.f()` returns an exact BigFloat; `EncloseF` resolves through up to
# three stages (eager Float64 estimate → Float128 pre-filter → MPFR ladder);
# `Enclose128F` carries a correctly-rounded Float128 bracket with an MPFR fallback.

using Random: AbstractRNG, default_rng

# Runtime switch (Float128 revision plan §5): every Float128 path fronts a complete
# MPFR path with identical semantics, so disabling costs speed, never correctness.
# Set ENV["ByteFloats_Float128"] = "disable" before loading (read in __init__).
const _USE_FLOAT128 = Ref(true)
@inline _f128() = _USE_FLOAT128[]

"""Deferred exact result: `f()` returns the value as an exact `BigFloat`."""
struct BigExactF{F}
    f::F
end
"""
Sticky-head exact result (wide-spread tail, Float128 revision follow-on): the true
value is `v + sgn·ε` for an infinitesimal `ε > 0` — `v` carries every bit the
projection can consume and `sgn ∈ {-1,+1}` the direction of the neglected tail.
Sound because `v` is either exactly on a rounding threshold of the target grid
(where `sticky` decides, for every mode including stochastic sub-grids) or
strictly farther from the nearest threshold than the tail magnitude; the emitting
sites in oracle.jl discharge that bound (operand significands ≤ 17 bits, target
grids ≤ P−1+N ≤ 67 fractional bits, spreads > the _DE_* thresholds). Non-allocating:
replaces the BigFloat escalation for FMA/FAA."""
struct StickyF{T<:Union{Float64,Float128}}
    v::T
    sgn::Int
end
"""
Deferred enclosure, resolved by `_finish` through up to three stages, cheapest first:

1. `yd` — an *eager* Float64 estimate (NaN ⇒ absent) whose libm evaluation is
   faithful (≤ 1 ulp); the true value is taken to lie within `|yd|·2^-45`
   (`_F64_RELEXP`), ≥ 2^7 slacker than any faithful-libm error.
2. `fq` — a zero-argument Float128 estimator (plan Site D, Class E; `nothing` ⇒
   absent): true value within `|y|·2^-90` of `y = fq()` (`_F128_RELEXP`), a bound
   ≥ 2^18 slacker than any published libquadmath transcendental error, discharged
   empirically by the differential-build tests.
3. `f(prec)` — the rigorous MPFR ladder: directed `(lo, hi)::NTuple{2,BigFloat}`
   strictly bracketing the true value.

Each estimate stage returns only when the two-sided sticky projection of its
envelope agrees on a single code point; any disagreement (or an absent/degenerate
estimate) falls through to the next stage. Stage 3 always decides.
"""
struct EncloseF{F,G}
    f::F
    fq::G
    yd::Float64
end
EncloseF(f) = EncloseF(f, nothing, NaN)
EncloseF(f, fq) = EncloseF(f, fq, NaN)
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
const _F64_RELEXP  = -45        # Float64 pre-filter envelope: E = |y|·2^-45
# The Float64 stage is sound for ops whose Float64 libm evaluation is faithful
# (error ≤ 1 ulp ≈ 2^-52 relative on normal results): 2^-45 gives ≥ 2^7 slack.
# Estimates that are non-finite, zero, or too near the subnormal range (where
# the relative-error model breaks) skip the stage; disagreement at the sticky
# gate falls through to the Float128 filter and the rigorous ladder unchanged.
const _F64_MINNORMISH = 6.7e-290   # ≈ 2^-960: comfortably clear of subnormals

@inline _finish(::Type{fr}, ρ::ProjSpec, R::Int, v::Float64) where {fr<:Binary} =
    project(fr, ρ, v; R)
@inline _finish(::Type{fr}, ρ::ProjSpec, R::Int, v::Float128) where {fr<:Binary} =
    project(fr, ρ, v; R)
@inline _finish(::Type{fr}, ρ::ProjSpec, R::Int, b::BigExactF) where {fr<:Binary} =
    project(fr, ρ, b.f(); R)
@inline _finish(::Type{fr}, ρ::ProjSpec, R::Int, s::StickyF) where {fr<:Binary} =
    project(fr, ρ, s.v; R, sticky=s.sgn)
function _finish(::Type{fr}, ρ::ProjSpec, R::Int, e::EncloseF) where {fr<:Binary}
    yd = e.yd
    if yd == yd && isfinite(yd) && abs(yd) >= _F64_MINNORMISH   # yd==yd: not NaN
        Ed = ldexp(abs(yd), _F64_RELEXP)
        dd = project(fr, ρ, yd - Ed; R, sticky=+1)
        du = project(fr, ρ, yd + Ed; R, sticky=-1)
        codepoint(dd) == codepoint(du) && return dd
    end
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

"""apply_op(Val(op), fr, ρ, R, xs...) — evaluate `op`'s ω-semantics on decoded
Float64 operands and project the result into `fr` under ρ (with random bits R)."""
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
"""An explicit stochastic draw `R`, or `nothing` to draw one from the rng."""
const MaybeR = Union{Nothing,Int}

"""Resolve a caller-supplied rng, falling back to the task-local default. Call
sites hoist this out of loops so array kernels resolve once per call, not per
element."""
@inline _resolve_rng(rng::MaybeRNG) = rng === nothing ? default_rng() : rng
"""The rng an operation under ρ will actually draw from — `nothing` for pure ρ,
which must never touch RNG state."""
@inline _rng_for(ρ::ProjSpec, rng::MaybeRNG) = isstochastic(ρ) ? _resolve_rng(rng) : nothing
@inline function _drawR(ρ::ProjSpec, rng::MaybeRNG, R::MaybeR)
    isstochastic(ρ) || return 0
    N = nrandbits(ρ)
    if R === nothing
        r = _resolve_rng(rng)
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
# One shape, generated at every arity: the three hand-written branches were the
# same two methods with one, two, or three operands spelled out. Arity now comes
# from the registry row, so a change to the calling convention — the keyword set,
# the draw, the decode step — cannot land unevenly across arities.
for op in OP_REGISTRY
    op.name === :Convert && continue
    name = op.name; V = Val{name}
    xs = [Symbol(:x, i) for i in 1:op.arity]
    spec_args = [:($x::Binary) for x in xs]               # spec form: any formats
    same_args = [:($x::T) for x in xs]                    # convenience form: one format
    decoded = [:(decode($x)) for x in xs]
    @eval begin
        @inline function $name(fr::Type{<:Binary}, ρ::ProjSpec, $(spec_args...);
                               rng::MaybeRNG=nothing, R::MaybeR=nothing)
            apply_op($V(), fr, ρ, _drawR(ρ, rng, R), $(decoded...))
        end
        @inline $name($(same_args...); kw...) where {T<:Binary} =
            $name(T, default_projspec(T), $(xs...); kw...)
    end
end

# ---- Convert (draft §4.9): the one op accepting external operands
"""
    Convert(fr, ρ, x) -> fr

Draft §4.9 Convert⟨f_x, f_r, ρ⟩. Accepts `Binary`, IEEE binary16/32/64
floats (widened exactly to the Float64 carrier), `Float128` (projected directly),
`Integer` (exact via a sufficiently wide BigFloat), and `BigFloat` (projected
directly; the caller warrants the value is exact).
"""
@inline function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::Binary;
                         rng::MaybeRNG=nothing, R::MaybeR=nothing)
    project(fr, ρ, decode(x); R=_drawR(ρ, rng, R))
end
@inline function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::Union{Float64,Float32,Float16};
                         rng::MaybeRNG=nothing, R::MaybeR=nothing)
    project(fr, ρ, Float64(x); R=_drawR(ρ, rng, R))   # exact widening
end
@inline function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::Float128;
                         rng::MaybeRNG=nothing, R::MaybeR=nothing)
    project(fr, ρ, x; R=_drawR(ρ, rng, R))            # preserve all 113 significand bits
end
function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::Integer;
                 rng::MaybeRNG=nothing, R::MaybeR=nothing)
    b = BigFloat(x; precision=max(64, ndigits(x, base=2) + 8))       # exact
    project(fr, ρ, b; R=_drawR(ρ, rng, R))
end
function Convert(fr::Type{<:Binary}, ρ::ProjSpec, x::BigFloat;
                 rng::MaybeRNG=nothing, R::MaybeR=nothing)
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
