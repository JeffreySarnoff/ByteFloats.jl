# ===== rand.jl — Random-API integration: rand / randn for Binary formats
#
# rand(T) — the format's answer to "uniform on [0, 1)": draw the 53-bit Float64
# uniform and floor-project it onto T's grid (TowardZero ≡ TowardNegative on
# [0, 1)). Each code point receives exactly the real measure of its floor
# interval — P(code c) = width of [d_c, d_next) — so the CDF agrees with the
# real uniform at every datum boundary. Results are always in [0, 1): the
# projection never rounds up and the Float64 draw is < 1. The Float64 grid
# (2^-53) is far finer than any K ≤ 8 format's, so the discretization it adds
# is negligible against the format's own.
#
# randn(T) — a standard-normal Float64 draw, projected round-to-nearest with
# SatFinite: a tail draw beyond MaxFiniteOf(T) clamps to the extremal finite
# datum, so randn never returns ±Inf (extended formats) or NaN (finite formats).
# Signed formats only — an unsigned format cannot represent half the mass, and
# folding it onto zero would be a silent lie.
#
# Wiring: one scalar method through each of Random's extension points. Every
# derived form — rand(T), rand(T, dims), rand!(A), randn(T, dims), randn!(A),
# and the rng-taking variants — then works through Random's generic machinery,
# and the rng-less forms use Julia's task-local default generator, a Xoshiro.
# Both methods produce values through `Convert`, i.e. the projection engine —
# the single write path into a code point.

# Binary <: AbstractFloat, so Random routes `rand(rng, T)` through the float
# sampler protocol (SamplerTrivial{CloseOpen01{T}}), the same hook the IEEE
# float types implement — not SamplerType, which never fires for floats.
Random.rand(rng::AbstractRNG, ::Random.SamplerTrivial{Random.CloseOpen01{T}}) where {T<:Binary} =
    Convert(T, RTZ_SatNone, rand(rng, Float64))

function Random.randn(rng::AbstractRNG, ::Type{T}) where {T<:Binary}
    issigned(T) || throw(ArgumentError(
        "randn requires a signed format; $(formatname(T)) cannot represent negative draws"))
    Convert(T, RNE_SatFinite, randn(rng))
end
