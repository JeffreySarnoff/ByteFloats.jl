# ===== defaults.jl — session-wide mutable defaults
#
# Settable defaults behind `Ref`s, read with `DefaultX()` and written with
# `DefaultX!(v)`. The projection default and its two components are kept
# coherent in both directions:
#   - `DefaultRoundingMode!` / `DefaultSaturationMode!` rebuild `DefaultProjection`
#     from the new component and the other one's current value;
#   - `DefaultProjection!` decomposes its argument back into both components.
# The invariant `DefaultProjection() === ProjSpec(DefaultRoundingMode(),
# DefaultSaturationMode())` therefore holds after any setter.
#
# These are *session* defaults for interactive use and library entry points that
# opt in to them. They are deliberately not consulted by the hot paths: the
# convenience methods (`a + b`, `Exp(x)`, …) continue to use the static
# `default_projspec`, so their specialization and tabling behavior cannot be
# changed — or broken — by a global. Plain `Ref` reads/writes, not atomic;
# set defaults from one task, not concurrently.
#
# Consumption discipline (the performant world-age-free pattern): a consumer
# never works directly on a `DefaultX()` result — the dispatch-steering reads
# are abstractly typed, and computing on them un-specializes the caller.
# Instead it goes through the `with_default_*` combinators below, which layer
# a speculative fast path over a function barrier:
#   1. read the Ref ONCE into a local (two reads could straddle a concurrent
#      set — the coupled setters keep the stored ProjSpec coherent, so guard
#      the pair, never the two component Refs separately);
#   2. `===`-test against the shipped initial value (`_GUARD_*`): on hit the
#      call is statically compiled against a constant — zero dispatch, zero
#      allocation, inference sees the concrete type;
#   3. on miss, cross the `@noinline` `_default_barrier`: one dynamic dispatch
#      + one boxed return, everything inside fully specialized.
# `DefaultRNG`/`DefaultRbits` need no combinator: their reads don't steer
# dispatch (`Ref{Int}` is concrete; an rng is passed through as a value).

# Speculation guards: the shipped initial values, `===`-tested on every
# combinator call. Changing an initial value means changing its guard — they
# must stay identical, so both come from these constants.
const _GUARD_TYPE       = Binary8p2se
const _GUARD_RETURN     = Binary8p2se
const _GUARD_ACCUM      = binary32
const _GUARD_PROJECTION = ProjSpec(NearestTiesToEven(), SatNone())

const _DEFAULT_TYPE       = Ref{Type{<:Binary}}(_GUARD_TYPE)
const _DEFAULT_RETURN     = Ref{Type{<:Binary}}(_GUARD_RETURN)
const _DEFAULT_ACCUM      = Ref{Type{<:AbstractFloat}}(_GUARD_ACCUM)
const _DEFAULT_ROUNDING   = Ref{RoundingMode3109}(roundingmode(_GUARD_PROJECTION))
const _DEFAULT_SATURATION = Ref{SaturationMode}(saturationmode(_GUARD_PROJECTION))
const _DEFAULT_PROJECTION = Ref{ProjSpec}(_GUARD_PROJECTION)
const _DEFAULT_RNG        = Ref{Union{Type{<:AbstractRNG},AbstractRNG}}(Random.Xoshiro)
const _DEFAULT_RBITS      = Ref{Int}(8)

"""
    DefaultType() -> Type{<:Binary}
    DefaultType!(T::Type{<:Binary}) -> T

The session's default format. Initialized to `Binary8p2se`.
"""
DefaultType() = _DEFAULT_TYPE[]
function DefaultType!(T::Type{Binary{K,P,S,E}}) where {K,P,S,E}
    checkformat(K, P, S, E)     # Binary{9,…} is a legal type object; only the params are validated
    _DEFAULT_TYPE[] = T
    T
end

"""
    DefaultReturnType() -> Type{<:Binary}
    DefaultReturnType!(T::Type{<:Binary}) -> T

The session's default *result* format — the format an operation projects into
when the caller does not name one. Initialized to `Binary8p2se`. Independent of
[`DefaultType`](@ref), which is the default *operand* format.
"""
DefaultReturnType() = _DEFAULT_RETURN[]
function DefaultReturnType!(T::Type{Binary{K,P,S,E}}) where {K,P,S,E}
    checkformat(K, P, S, E)
    _DEFAULT_RETURN[] = T
    T
end

"""
    DefaultAccumulatorType() -> Type{<:AbstractFloat}
    DefaultAccumulatorType!(T::Type{<:AbstractFloat}) -> T

The session's default accumulator carrier for wide-precision work (reductions,
dot products). Initialized to `binary32` (the exported IEEE 754 alias for
`Float32`). Any `AbstractFloat` type is accepted (`binary64`, `Float128`,
`BigFloat`, …).
"""
DefaultAccumulatorType() = _DEFAULT_ACCUM[]
DefaultAccumulatorType!(T::Type{<:AbstractFloat}) = (_DEFAULT_ACCUM[] = T; T)

"""
    DefaultRoundingMode() -> RoundingMode3109
    DefaultRoundingMode!(m) -> m

The session's default rounding mode. Initialized to `NearestTiesToEven()`.
Setting it rebuilds [`DefaultProjection`](@ref) from the new mode and the
current [`DefaultSaturationMode`](@ref). Accepts an instance or a
fully-parameterized type (`NearestTiesToAway`, `StochasticA{8}`, …).
"""
DefaultRoundingMode() = _DEFAULT_ROUNDING[]
function DefaultRoundingMode!(m::RoundingMode3109)
    _DEFAULT_ROUNDING[] = m
    _DEFAULT_PROJECTION[] = ProjSpec(m, _DEFAULT_SATURATION[])
    m
end
DefaultRoundingMode!(M::Type{<:RoundingMode3109}) = DefaultRoundingMode!(M())

"""
    DefaultSaturationMode() -> SaturationMode
    DefaultSaturationMode!(s) -> s

The session's default saturation mode. Initialized to `SatNone()`.
Setting it rebuilds [`DefaultProjection`](@ref) from the current
[`DefaultRoundingMode`](@ref) and the new mode. Accepts an instance or a type.
"""
DefaultSaturationMode() = _DEFAULT_SATURATION[]
function DefaultSaturationMode!(s::SaturationMode)
    _DEFAULT_SATURATION[] = s
    _DEFAULT_PROJECTION[] = ProjSpec(_DEFAULT_ROUNDING[], s)
    s
end
DefaultSaturationMode!(S::Type{<:SaturationMode}) = DefaultSaturationMode!(S())

"""
    DefaultProjection() -> ProjSpec
    DefaultProjection!(ρ::ProjSpec) -> ρ
    DefaultProjection!(m, s) -> ProjSpec

The session's default projection specification. Initialized to
`(DefaultRoundingMode(), DefaultSaturationMode())` = `RNE_SatNone`.
Setting it directly decomposes ρ into [`DefaultRoundingMode`](@ref) and
[`DefaultSaturationMode`](@ref), so the three stay coherent in both directions.
"""
DefaultProjection() = _DEFAULT_PROJECTION[]
function DefaultProjection!(ρ::ProjSpec)
    _DEFAULT_PROJECTION[] = ρ
    _DEFAULT_ROUNDING[] = roundingmode(ρ)
    _DEFAULT_SATURATION[] = saturationmode(ρ)
    ρ
end
DefaultProjection!(m, s) = DefaultProjection!(ProjSpec(projmode(m), s isa Type ? s() : s))

"""
    DefaultRNG() -> Type{<:AbstractRNG} | AbstractRNG
    DefaultRNG!(rng) -> rng

The session's default random-number generator for stochastic rounding.
Initialized to the `Xoshiro` *type* (a fresh generator per use); may be set to
an `AbstractRNG` instance for a reproducible stream.
"""
DefaultRNG() = _DEFAULT_RNG[]
DefaultRNG!(rng::Union{Type{<:AbstractRNG},AbstractRNG}) = (_DEFAULT_RNG[] = rng; rng)

"""
    DefaultRbits() -> Int
    DefaultRbits!(n::Int) -> n

The session's default random-bit budget N for the stochastic rounding families
(`StochasticA{N}`, `StochasticB{N}`, `StochasticC{N}`). Initialized to 8;
must satisfy 1 ≤ N ≤ 60.
"""
DefaultRbits() = _DEFAULT_RBITS[]
DefaultRbits!(n::Int) = (_check_nrandbits(n); _DEFAULT_RBITS[] = n; n)

# ---------------------------------------------------------------------------
# Consumption combinators: speculative fast path over a function barrier
# ---------------------------------------------------------------------------

# The slow path. `@noinline` so speculation failure cannot bloat or
# de-specialize the caller; `where {F}` forces specialization on the closure,
# `where {T}` / dispatch on ρ recovers the concrete type from the abstract Ref
# read. Cost: one dynamic dispatch on entry + one boxed return.
@noinline _default_barrier(f::F, ::Type{T}, args...) where {F,T} = f(T, args...)
@noinline _default_barrier(f::F, ρ::ProjSpec, args...) where {F} = f(ρ, args...)

"""
    with_default_type(f, args...)          -> f(DefaultType(), args...)
    with_default_returntype(f, args...)    -> f(DefaultReturnType(), args...)
    with_default_accumulatortype(f, args...) -> f(DefaultAccumulatorType(), args...)
    with_default_projection(f, args...)    -> f(DefaultProjection(), args...)

Call `f` with the named session default as its first argument — the supported
way to *consume* a default. While the default still holds its initial value
(the overwhelmingly common case) the call is statically compiled against that
constant: zero dispatch, zero allocation, concrete inferred result. After the
default is changed, the call crosses a function barrier instead: one dynamic
dispatch plus one boxed return, with everything inside fully specialized.

```julia-repl
julia> with_default_type((T, x) -> T(x), 1.5)
Binary8p2se(1.5 ≡ 0x2e)

julia> with_default_projection((ρ, x, y) -> Add(Binary8p4se, ρ, x, y),
                               Binary8p4se(1.5), Binary8p4se(0.25))
Binary8p4se(1.75 ≡ 0x46)
```
"""
@inline function with_default_type(f::F, args...) where {F}
    T = DefaultType()                       # single Ref read
    T === _GUARD_TYPE && return f(_GUARD_TYPE, args...)
    _default_barrier(f, T, args...)
end

@inline function with_default_returntype(f::F, args...) where {F}
    T = DefaultReturnType()
    T === _GUARD_RETURN && return f(_GUARD_RETURN, args...)
    _default_barrier(f, T, args...)
end

@inline function with_default_accumulatortype(f::F, args...) where {F}
    T = DefaultAccumulatorType()
    T === _GUARD_ACCUM && return f(_GUARD_ACCUM, args...)
    _default_barrier(f, T, args...)
end

@inline function with_default_projection(f::F, args...) where {F}
    ρ = DefaultProjection()                 # guard the coherent pair, not the component Refs
    ρ === _GUARD_PROJECTION && return f(_GUARD_PROJECTION, args...)
    _default_barrier(f, ρ, args...)
end

@doc (@doc with_default_type) with_default_returntype
@doc (@doc with_default_type) with_default_accumulatortype
@doc (@doc with_default_type) with_default_projection
