# ===== rand.jl — Random-API integration: rand / randn for Binary formats
#
# rand(T) — the format's answer to "uniform on [0, 1)": draw the 53-bit Float64
# uniform and project it onto T's grid under the `projection` keyword. The
# default is RTZ_SatNone — floor on [0, 1) (≡ TowardNegative there) — chosen so
# each code point receives exactly the real measure of its floor interval and
# the result is provably < 1, preserving Julia's documented rand contract. A
# nearest-mode projection is available through the keyword, but it can return
# 1.0 (mass near the top rounds up), which is why it is not the default.
#
# randn(T) — a standard-normal Float64 draw, projected under `projection`,
# default RNE_SatFinite: round-to-nearest with tail draws beyond MaxFiniteOf(T)
# clamped to the extremal finite datum, so randn never returns ±Inf (extended
# formats) or NaN (finite formats). Signed formats only — an unsigned format
# cannot represent half the mass, and folding it onto zero would be a silent lie.
#
# The `projection` keyword lives on the scalar `::Type` methods — the only
# places Random's dispatch can carry it (the generic array/dims wrappers do not
# forward keywords, and the sampler hook is invoked without them). A stochastic
# projection draws its R from the SAME rng as the uniform/normal draw (threaded
# into `Convert`), so seeded streams stay reproducible. For arrays under a
# non-default projection, draw scalars: `[rand(rng, T; projection=ρ) for _ in 1:n]`.
#
# Wiring: the sampler hook serves every derived form — rand(T, dims), rand!(A),
# randn(T, dims), randn!(A), … — at the default semantics; the rng-less forms
# use Julia's task-local default generator, a Xoshiro. Binary <: AbstractFloat,
# so rand routes through the float sampler protocol
# (SamplerTrivial{CloseOpen01{T}}, the hook the IEEE types implement — not
# SamplerType, which never fires for floats). Every method produces its value
# through `Convert`, i.e. the projection engine — the single write path into a
# code point.

Random.rand(rng::AbstractRNG, ::Random.SamplerTrivial{Random.CloseOpen01{T}}) where {T<:Binary} =
    Convert(T, RTZ_SatNone, rand(rng, Float64))

# rand([rng], T; projection=RTZ_SatNone) / randn([rng], T; projection=RNE_SatFinite):
# the scalar ::Type forms take any ProjSpec via the keyword. The defaults keep
# the Julia contracts (floor ⇒ always < 1; nearest+SatFinite ⇒ never ±Inf/NaN);
# a stochastic projection draws its random bits from the same rng. Documented in
# the User Guide/Examples rather than a docstring: attaching documentation to
# the Base `rand`/`randn` bindings drags Base's whole doc group (and its
# unresolvable @refs) into this package's reference page.
@inline Random.rand(rng::AbstractRNG, ::Type{T}; projection::ProjSpec=RTZ_SatNone) where {T<:Binary} =
    Convert(T, projection, rand(rng, Float64); rng)::T
@inline Random.rand(::Type{T}; projection::ProjSpec=RTZ_SatNone) where {T<:Binary} =
    rand(default_rng(), T; projection)::T

function Random.randn(rng::AbstractRNG, ::Type{T}; projection::ProjSpec=RNE_SatFinite) where {T<:Binary}
    issigned(T) || throw(ArgumentError(
        "randn requires a signed format; $(formatname(T)) cannot represent negative draws"))
    Convert(T, projection, randn(rng); rng)::T
end
@inline Random.randn(::Type{T}; projection::ProjSpec=RNE_SatFinite) where {T<:Binary} =
    randn(default_rng(), T; projection)::T
