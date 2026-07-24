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
#      call is statically compiled against a constant — no dynamic dispatch;
#   3. on miss, cross the `@noinline` `_default_barrier`: one dynamic dispatch,
#      everything inside fully specialized.
# What the guard can and cannot buy: the combinator's return type is the UNION
# over both branches. When `f`'s result type does not depend on the default —
# the projection combinator's normal shape, where the caller fixes the formats
# and ρ steers only the rounding — both branches infer the same concrete type
# and the whole call is zero-allocation (pinned in the suite). When the result's
# type IS the default (`with_default_type` as a constructor), the slow branch
# infers `Any`, so the value boxes once at escape even on the fast path: the
# irreducible cost of a runtime-chosen type. Only `@eval` method redefinition
# (world age) could remove that box; these combinators deliberately don't.
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
const _DEFAULT_RBITS      = Ref{Int}(DEFAULT_RBITS)   # one source for the budget "8"

"""
    DefaultType() -> Type{<:Binary}
    DefaultType!(T::Type{<:Binary}) -> T

The session's default format. Initialized to `Binary8p2se`.
"""
# Shared by the two format-valued setters: `Binary{9,…}` is a perfectly legal
# type object, so it is the *parameters* that must be validated, not the type.
function _set_format_default!(ref, ::Type{Binary{K,P,S,E}}) where {K,P,S,E}
    checkformat(K, P, S, E)
    ref[] = Binary{K,P,S,E}
end

DefaultType() = _DEFAULT_TYPE[]
DefaultType!(T::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    _set_format_default!(_DEFAULT_TYPE, T)

"""
    DefaultReturnType() -> Type{<:Binary}
    DefaultReturnType!(T::Type{<:Binary}) -> T

The session's default *result* format — the format an operation projects into
when the caller does not name one. Initialized to `Binary8p2se`. Independent of
[`DefaultType`](@ref), which is the default *operand* format.
"""
DefaultReturnType() = _DEFAULT_RETURN[]
DefaultReturnType!(T::Type{Binary{K,P,S,E}}) where {K,P,S,E} =
    _set_format_default!(_DEFAULT_RETURN, T)

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
# read, so `f` runs fully specialized inside. Cost: one dynamic dispatch on
# entry (plus a boxed return only when `f`'s result type is the default —
# see the allocation contract below).
@noinline _default_barrier(f::F, ::Type{T}, args...) where {F,T} = f(T, args...)
@noinline _default_barrier(f::F, ρ::ProjSpec, args...) where {F} = f(ρ, args...)

# The shared speculation. `current` is the *already-read* default: each
# combinator reads its Ref exactly once, at the call below, so the read-once
# discipline lives in one place rather than in four copies. `guard` carries its
# value in its type (a `Type{…}` or a singleton `ProjSpec`), so the hit branch
# is compiled against a constant.
@inline function _with_default(f::F, current, guard, args...) where {F}
    current === guard && return f(guard, args...)
    _default_barrier(f, current, args...)
end

"""
    with_default_type(f, args...)          -> f(DefaultType(), args...)
    with_default_returntype(f, args...)    -> f(DefaultReturnType(), args...)
    with_default_accumulatortype(f, args...) -> f(DefaultAccumulatorType(), args...)
    with_default_projection(f, args...)    -> f(DefaultProjection(), args...)

Call `f` with the named session default as its first argument — the supported
way to *consume* a default. While the default still holds its initial value
(the overwhelmingly common case) the call is statically compiled against that
constant: no dynamic dispatch, `f` fully specialized. After the default is
changed, the call crosses a function barrier instead: one dynamic dispatch,
everything inside fully specialized.

Allocation contract: the combinator's return type unions both paths. When `f`'s
result type does not depend on the default — e.g. `with_default_projection`
with the formats fixed by the caller — the union is concrete and the call is
**zero-allocation with a concretely inferred result**. When the result's type
*is* the default (`with_default_type((T, x) -> T(x), …)`), the value is computed
on the specialized path but boxes once at escape — the irreducible cost of a
runtime-chosen type.

```julia-repl
julia> with_default_type((T, x) -> T(x), 1.5)
Binary8p2se(1.5 ≡ 0x41)

julia> with_default_projection((ρ, x, y) -> Add(Binary8p4se, ρ, x, y),
                               Binary8p4se(1.5), Binary8p4se(0.25))
Binary8p4se(1.75 ≡ 0x46)
```
"""
@inline with_default_type(f::F, args...) where {F} =
    _with_default(f, DefaultType(), _GUARD_TYPE, args...)

@inline with_default_returntype(f::F, args...) where {F} =
    _with_default(f, DefaultReturnType(), _GUARD_RETURN, args...)

@inline with_default_accumulatortype(f::F, args...) where {F} =
    _with_default(f, DefaultAccumulatorType(), _GUARD_ACCUM, args...)

# Reads the coherent pair, never the two component Refs separately: a concurrent
# set could otherwise hand the guard a torn (rounding, saturation) combination.
@inline with_default_projection(f::F, args...) where {F} =
    _with_default(f, DefaultProjection(), _GUARD_PROJECTION, args...)

@doc (@doc with_default_type) with_default_returntype
@doc (@doc with_default_type) with_default_accumulatortype
@doc (@doc with_default_type) with_default_projection
