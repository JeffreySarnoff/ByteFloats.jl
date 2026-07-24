# ===== test/quality.jl — package-hygiene and static-analysis gates
#
# Aqua: the community package-quality checks (method ambiguities, undefined
# exports, Project.toml/test-extras agreement, stale deps, compat bounds, type
# piracy, unbound type parameters, persistent tasks).
#
# JET: abstract-interpretation error analysis, run two ways —
#   1. whole-package (`report_package`), which analyzes every method against its
#      *declared* signature;
#   2. concrete entry points (`@test_call`), which analyze the calls users
#      actually make, with the format types statically known.
# Both are needed: (1) catches unreachable-method bugs in code no test exercises,
# (2) is the analysis that matches this package's specialization doctrine — and
# it is the one that can see through correlated type parameters.

using Test
using Random
using ByteFloats
using ByteFloats: project
using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(ByteFloats)
end

# `_vmap_packed` builds `Vector{fr}` from its `::Type{fr}` argument and views it
# per tile. Analyzed as a generic method, JET widens the correlation between the
# type parameter and the container's element type — it sees `SubArray{Any}` where
# every specialization has `SubArray{fr}` — and reports a method error that the
# concrete-call gate below shows does not exist. Nothing else may report.
const _JET_KNOWN_GENERIC_WIDENING = (:_vmap_packed,)

_report_method(r) = let li = r.vst[end].linfo
    li isa Core.MethodInstance && li.def isa Method ? li.def.name : nothing
end

@testset "JET (whole package)" begin
    result = JET.report_package(ByteFloats; target_modules = (ByteFloats,),
                                toplevel_logger = nothing)
    reports = filter(r -> _report_method(r) ∉ _JET_KNOWN_GENERIC_WIDENING,
                     JET.get_reports(result))
    @test isempty(reports)
    isempty(reports) || foreach(display, reports)
end

@testset "JET (concrete entry points)" begin
    T = Binary8p4se
    S = Binary8p3se
    a, b = T(1.5), T(0.25)
    A = T.(randn(64)); B = T.(randn(64)); D = similar(A)
    bx = Block(one(S), (a, b, a, b))
    pv = PackedVector(A)
    σ = ProjSpec(StochasticA{8}(), SatNone())

    JET.@test_call Add(T, RNE_SatNone, a, b)
    JET.@test_call FMA(T, RNE_SatNone, a, b, a)
    JET.@test_call Exp(T, RNE_SatNone, a)
    JET.@test_call Convert(S, RNE_SatNone, a)
    JET.@test_call a + b
    JET.@test_call project(T, RNE_SatNone, 1.6)
    JET.@test_call vmap!(D, Val(:Exp), T, RNE_SatNone, A)
    JET.@test_call vmap!(D, Val(:Add), T, RNE_SatNone, A, B)
    JET.@test_call vmap!(D, Val(:FMA), T, RNE_SatNone, A, B, A)
    JET.@test_call BlockDotProduct(T, RNE_SatNone, bx, bx)
    JET.@test_call BlockAdd(T, RNE_SatNone, bx, bx, one(S))
    JET.@test_call ScaledAdd(T, RNE_SatNone, one(S), a, one(S), b)
    # the method the package-wide pass cannot verify, verified where it is called
    JET.@test_call vmap(:Exp, T, RNE_SatNone, pv)
    JET.@test_call vmap(:Exp, T, σ, pv; rng = MersenneTwister(1))
    # session-default combinators: the fast path must be statically clean
    JET.@test_call with_default_type((F, x) -> F(x), 1.5)
    JET.@test_call with_default_projection((ρ, x, y) -> Add(T, ρ, x, y), a, b)
end
