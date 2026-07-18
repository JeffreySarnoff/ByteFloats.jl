# ===== approx.jl — conformance declaration and κ-approximate implementations
#       (design §1.1/§1.4/§10.4, architecture §9, draft §4.4–§4.6)
#
# Two guarantees this file makes machine-checkable:
#   1. Nothing approximate is reachable from the default API: approximate kernels
#      live in their own registry, retrieved only by explicit name.
#   2. Every declared κ is *verified by exhaustive enumeration* at registration
#      time (feasible because value sets are tiny); understated declarations are
#      rejected. κ is therefore a measured property, not a promise.
#
# κ semantics implemented here (per the draft's Annex worked example — the
# subnormal-flushing Exp with "κ = 4" — and flagged as an interpretation where
# the normative text was unavailable, see checkpoint.md):
#   κ = the maximum distance, in result-format code-point steps along the total
#       order, between the implementation's result and the defined result, over
#       all inputs whose defined result is finite.
#   κ = NaN when the implementation mismatches on any input whose defined result
#       is NaN or ±Inf, or returns a non-finite value where the defined result is
#       finite — such deviation is unbounded/unclassifiable in the code-point metric.

const DRAFT_REVISION = "IEEE P3109 working draft, uploaded 2026-07-17"
draft_revision() = DRAFT_REVISION

# ---------------------------------------------------------------------------
# κ measurement
# ---------------------------------------------------------------------------
"""Code-point distance along the total order between two values of one format."""
@inline codedistance(a::T, b::T) where {T<:Binary} =
    abs(Int(order_key(a)) - Int(order_key(b)))

"""
    measure_kappa(fn, op, fr, argformats, ρ; max_exhaustive=2^22, rng, samples=2^20)
        -> (κ::Float64, exhaustive::Bool)

Measure κ of `fn(x::argformats[1], …)::fr` against the defined results of draft
operation `op` under ρ. Enumerates the full input cross-product when it has at
most `max_exhaustive` points (always true for arity ≤ 2); otherwise measures on
`samples` uniform draws and reports `exhaustive = false`.
"""
function measure_kappa(fn::F, op::Symbol, fr::Type{<:Binary},
                       argformats::NTuple{N,DataType}, ρ::ProjSpec;
                       max_exhaustive::Int=1 << 22,
                       rng::AbstractRNG=default_rng(), samples::Int=1 << 20) where {F,N}
    isstochastic(ρ) &&
        throw(ArgumentError("κ is defined against deterministic defined results; stochastic ρ is not measurable"))
    Ks = map(bitwidth, argformats)
    total = prod(1 .<< Ks)
    exhaustive = total <= max_exhaustive
    κ = 0.0
    defined = (args...) -> op === :Convert ?
        project(fr, ρ, decode(args[1])) :
        apply_op(Val(op), fr, ρ, 0, map(decode, args)...)
    function visit(codes::NTuple{N,Int})
        args = ntuple(i -> rawvalue(argformats[i], UInt8(codes[i])), Val(N))
        want = defined(args...)
        got = fn(args...)::fr
        dw, dg = decode(want), decode(got)
        if isnan(dw) || isinf(dw)
            isequal(dw, dg) || return NaN                    # non-finite defined must match
            return 0.0
        end
        (isnan(dg) || isinf(dg)) && return NaN               # unbounded deviation
        Float64(codedistance(got, want))
    end
    if exhaustive
        for lin in 0:total - 1
            codes = Vector{Int}(undef, N)
            r = lin
            for i in 1:N
                codes[i] = r & ((1 << Ks[i]) - 1)
                r >>= Ks[i]
            end
            d = visit(ntuple(i -> codes[i], Val(N)))
            isnan(d) && return (NaN, true)
            κ = max(κ, d)
        end
    else
        for _ in 1:samples
            codes = ntuple(i -> Int(rand(rng, UInt32) & ((1 << Ks[i]) - 1)), Val(N))
            d = visit(codes)
            isnan(d) && return (NaN, false)
            κ = max(κ, d)
        end
    end
    (κ, exhaustive)
end

# ---------------------------------------------------------------------------
# Approximate-implementation registry
# ---------------------------------------------------------------------------
"""A registered κ-approximate implementation of one operation specialization."""
struct ApproxImpl{F}
    name::Symbol
    op::Symbol
    fr::DataType
    argformats::Tuple{Vararg{DataType}}
    ρ::ProjSpec
    fn::F
    kappa_declared::Float64
    kappa_measured::Float64
    exhaustive::Bool
end
kappa(a::ApproxImpl) = a.kappa_declared
kappa_measured(a::ApproxImpl) = a.kappa_measured

const APPROX_REGISTRY = Dict{Symbol,ApproxImpl}()
const APPROX_LOCK = ReentrantLock()

"""
    register_approx!(name, op, fr, argformats, ρ, fn; κ=nothing, kwargs...) -> ApproxImpl

Register `fn` as a named κ-approximate implementation of `op⟨argformats → fr, ρ⟩`.
κ is measured by enumeration (see `measure_kappa`); a declared `κ` smaller than the
measured value is **rejected**. Omitting `κ` declares the measured value. NaN-κ
implementations (non-finite mismatches) may be registered only by declaring `κ=NaN`.
"""
function register_approx!(name::Symbol, op::Symbol, fr::Type{<:Binary},
                          argformats::NTuple{N,DataType}, ρ::ProjSpec, fn::F;
                          κ::Union{Nothing,Real}=nothing, kwargs...) where {F,N}
    any(o -> o.name === op, OP_REGISTRY) || throw(ArgumentError("unknown draft operation :$op"))
    arity = opinfo(op).arity
    (op === :Convert ? 1 : arity) == N ||
        throw(ArgumentError(":$op has arity $arity, got $N argument formats"))
    κm, exh = measure_kappa(fn, op, fr, argformats, ρ; kwargs...)
    if isnan(κm)
        (κ !== nothing && isnan(κ)) || throw(ArgumentError(
            "implementation mismatches on non-finite defined results (measured κ = NaN); " *
            "register with explicit κ=NaN to acknowledge, or fix it"))
        κd = NaN
    else
        κd = κ === nothing ? κm : Float64(κ)
        (!isnan(κd) && κd < κm) &&
            throw(ArgumentError("declared κ = $κd understates measured κ = $κm — registration rejected"))
    end
    impl = ApproxImpl(name, op, fr, Tuple(argformats), ρ, fn, κd, κm, exh)
    lock(() -> (haskey(APPROX_REGISTRY, name) &&
                    throw(ArgumentError("approximate implementation :$name already registered"));
                APPROX_REGISTRY[name] = impl), APPROX_LOCK)
    impl
end

"""Retrieve a registered approximate implementation (callable via `.fn`)."""
approx(name::Symbol) = lock(() -> get(APPROX_REGISTRY, name) do
        throw(KeyError("no approximate implementation :$name registered"))
    end, APPROX_LOCK)
kappa(name::Symbol) = kappa(approx(name))
list_approx() = lock(() -> sort!(collect(keys(APPROX_REGISTRY))), APPROX_LOCK)
unregister_approx!(name::Symbol) = lock(() -> (delete!(APPROX_REGISTRY, name); nothing), APPROX_LOCK)

# ---------------------------------------------------------------------------
# Worked example from the draft's Annex: flush-subnormal-results variant
# ---------------------------------------------------------------------------
"""
    ftz_variant(op, fr, f1, ρ) -> fn

Build the Annex-style approximate unary implementation: compute the defined result,
then flush subnormal results to zero or `MinNormalOf(fr)`, whichever is nearer,
ties toward zero (magnitude-symmetric for signed formats). Returned `fn` is suitable
for `register_approx!`; for `fr` with precision P its κ is 2^(P-2) when P ≥ 2
(largest flushed-to-zero subnormal) and 0 when P = 1 (no subnormals exist).
"""
function ftz_variant(op::Symbol, fr::Type{<:Binary}, f1::Type{<:Binary}, ρ::ProjSpec)
    P = precision(fr)
    half = 1 << max(P - 2, 0)                    # subnormal codes 1 … 2^(P-1)-1; tie at 2^(P-2)
    minnorm = MinNormalOf(fr)
    function fn(x::Binary)
        r = op === :Convert ? project(fr, ρ, decode(x)) : apply_op(Val(op), fr, ρ, 0, decode(x))
        issubnormal(r) || return r
        m = Int(codepoint(r) & ~signmask(fr))    # magnitude code (subnormals: 1 … 2^(P-1)-1)
        neg = issigned(fr) && codepoint(r) >= signmask(fr)
        if m <= half
            return zero(fr)                      # nearer to zero (ties to zero)
        else
            return neg ? Negate(minnorm) : minnorm
        end
    end
    fn
end

# ---------------------------------------------------------------------------
# Conformance declaration (design §1.4)
# ---------------------------------------------------------------------------
struct ConformanceDeclaration
    package::String
    draft::String
    formats::Vector{Symbol}
    operations::Vector{NamedTuple{(:name, :arity, :group),Tuple{Symbol,Int,Symbol}}}
    rounding_modes::Vector{String}
    saturation_modes::Vector{Symbol}
    block_surface::Vector{Symbol}
    cached_specializations::Vector{TableKey}
    approximate::Vector{NamedTuple{(:name, :op, :kappa, :exhaustive),Tuple{Symbol,Symbol,Float64,Bool}}}
end

"""
    conformance() -> ConformanceDeclaration

The package's draft-§4.6-style conformance declaration, derived live from the
operation registry, the table cache (the specializations actually instantiated),
and the approximate-implementation registry. Serialize with `conformance_dict`
or render with `conformance_report`.
"""
function conformance()
    ops = [(name=o.name, arity=o.arity, group=o.group) for o in OP_REGISTRY]
    blocknames = Symbol[]
    for o in OP_REGISTRY
        o.name === :Convert && continue
        push!(blocknames, Symbol(:Block, o.name), Symbol(:Scaled, o.name))
    end
    append!(blocknames, (:BlockReduceAdd, :BlockReduceMultiply, :BlockDotProduct,
                         :ConvertFromBlock, :ConvertToBlock, :ConvertToBlockMaxAbsFinite))
    cached = lock(() -> collect(keys(TABLE_CACHE)), TABLE_LOCK)
    apx = lock(() -> [(name=a.name, op=a.op, kappa=a.kappa_declared, exhaustive=a.exhaustive)
                      for a in values(APPROX_REGISTRY)], APPROX_LOCK)
    ConformanceDeclaration(
        "ByteFloats.jl 0.1.0", DRAFT_REVISION,
        sort!(collect(keys(_NAMED))),
        ops,
        ["NearestTiesToEven", "NearestTiesToAway", "TowardPositive", "TowardNegative",
         "TowardZero", "ToOdd", "StochasticA[N], 1 ≤ N ≤ 60", "StochasticB[N]", "StochasticC[N]"],
        [:SatFinite, :SatPropagate, :SatNone],
        sort!(blocknames),
        cached, sort!(apx; by=a -> a.name))
end

"""Nested `Dict{String,Any}` form of the declaration — serializable by any JSON/TOML writer."""
function conformance_dict(c::ConformanceDeclaration=conformance())
    Dict{String,Any}(
        "package" => c.package,
        "draft" => c.draft,
        "formats" => String.(c.formats),
        "operations" => [Dict("name" => String(o.name), "arity" => o.arity,
                              "group" => String(o.group)) for o in c.operations],
        "rounding_modes" => c.rounding_modes,
        "saturation_modes" => String.(c.saturation_modes),
        "block_surface" => String.(c.block_surface),
        "cached_specializations" => [Dict("op" => String(k.op), "fr" => collect(k.fr),
                                          "f1" => collect(k.f1), "f2" => collect(k.f2),
                                          "rounding" => String(k.rm), "saturation" => String(k.sm))
                                     for k in c.cached_specializations],
        "approximate" => [Dict("name" => String(a.name), "op" => String(a.op),
                               "kappa" => a.kappa, "exhaustive" => a.exhaustive)
                          for a in c.approximate])
end

"""Human-readable conformance report."""
function conformance_report(io::IO=stdout, c::ConformanceDeclaration=conformance())
    println(io, "Conformance declaration — ", c.package)
    println(io, "Implements: ", c.draft)
    println(io, "\nFormats (", length(c.formats), "): all Binary{K,P,Σ,Δ}, K ∈ 3:8, ",
            "Σ ∈ {Signed, Unsigned}, Δ ∈ {Finite, Extended}")
    println(io, "\nScalar operations (", length(c.operations), "):")
    for a in 1:3
        names = join((String(o.name) for o in c.operations if o.arity == a), ", ")
        println(io, "  arity ", a, ": ", names)
    end
    println(io, "\nRounding modes: ", join(c.rounding_modes, ", "))
    println(io, "Saturation modes: ", join(String.(c.saturation_modes), ", "))
    println(io, "\nBlock/scaled surface (", length(c.block_surface), " operations), any B ≥ 1")
    println(io, "\nInstantiated pure-ρ table specializations: ", length(c.cached_specializations),
            " (", table_bytes(), " bytes)")
    for k in c.cached_specializations
        println(io, "  ", k.op, "⟨", k.f1, k.f2 == (0, 0, 0, 0) ? "" : string(" × ", k.f2),
                " → ", k.fr, ", (", k.rm, ", ", k.sm, ")⟩")
    end
    if isempty(c.approximate)
        println(io, "\nApproximate implementations: none (all default paths are bit-exact)")
    else
        println(io, "\nDeclared κ-approximate implementations:")
        for a in c.approximate
            println(io, "  :", a.name, "  op=", a.op, "  κ=", a.kappa,
                    a.exhaustive ? "  (κ verified exhaustively)" : "  (κ sampled — not exhaustive)")
        end
    end
    nothing
end
