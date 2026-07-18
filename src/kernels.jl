# ===== kernels.jl — array kernels (design §8, architecture §7)
#
# Two hot-loop shapes, written once and instantiated by the registry:
#   Shape A (gather):  out[i] = tbl[key(a[i]…)]   — pure ρ, arity ≤ 2
#   Shape B (compute): per-element scalar path    — ternary ops and stochastic ρ
# The table getter is @noinline and called ONCE per array call, hoisted out of
# the loop; loop bodies index a local Memory{UInt8} with no dict lookups, locks,
# or global loads per element. Threading and the K=8×8 compute-vs-table crossover
# are the Phase-2 items tracked in checkpoint.md; v1 loops are sequential.

"""
    vmap!(dest, Val(op), fr, ρ, A[, B[, C]]; rng) -> dest
    vmap(op, fr, ρ, A...; rng)

Elementwise draft operation over arrays of `Binary` values, projecting into `fr`
under ρ. Pure ρ with arity ≤ 2 runs the Shape-A gather; ternary and stochastic
specializations run the scalar path per element (each element drawing its own R).
"""
function vmap! end

# ---- Shape A: unary gather
function vmap!(dest::AbstractArray{FR}, ::Val{op}, ::Type{FR}, ρ::ProjSpec,
               A::AbstractArray{F1}) where {op,FR<:Binary,F1<:Binary}
    axes(dest) == axes(A) || throw(DimensionMismatch("dest and A must share axes"))
    if isstochastic(ρ)
        return _vmap_scalar!(dest, Val(op), FR, ρ, A)
    end
    tbl = get_table(op, FR, F1, ρ)                       # hoisted; @noinline
    @inbounds for i in eachindex(dest, A)
        dest[i] = rawvalue(FR, tbl[Int(codepoint(A[i])) + 1])
    end
    dest
end

# ---- Shape A: binary gather, index = (c1 << K2) | c2
function vmap!(dest::AbstractArray{FR}, ::Val{op}, ::Type{FR}, ρ::ProjSpec,
               A::AbstractArray{F1}, B::AbstractArray{F2}) where {op,FR<:Binary,F1<:Binary,F2<:Binary}
    axes(dest) == axes(A) == axes(B) || throw(DimensionMismatch("dest, A, B must share axes"))
    if isstochastic(ρ)
        return _vmap_scalar!(dest, Val(op), FR, ρ, A, B)
    end
    tbl = get_table(op, FR, F1, F2, ρ)
    K2 = bitwidth(F2)
    @inbounds for i in eachindex(dest, A, B)
        dest[i] = rawvalue(FR, tbl[(Int(codepoint(A[i])) << K2) + Int(codepoint(B[i])) + 1])
    end
    dest
end

# ---- Shape B: ternary (never tabulable at 2^(3K)) and the stochastic fallback
function vmap!(dest::AbstractArray{FR}, ::Val{op}, ::Type{FR}, ρ::ProjSpec,
               A::AbstractArray{<:Binary}, B::AbstractArray{<:Binary}, C::AbstractArray{<:Binary};
               rng::MaybeRNG=nothing) where {op,FR<:Binary}
    axes(dest) == axes(A) == axes(B) == axes(C) || throw(DimensionMismatch("operand axes must match"))
    rr = isstochastic(ρ) ? (rng === nothing ? default_rng() : rng) : nothing
    @inbounds for i in eachindex(dest, A, B, C)
        R = _drawR(ρ, rr, nothing)
        dest[i] = apply_op(Val(op), FR, ρ, R, decode(A[i]), decode(B[i]), decode(C[i]))
    end
    dest
end
function _vmap_scalar!(dest, ::Val{op}, ::Type{FR}, ρ::ProjSpec, A;
                       rng::MaybeRNG=nothing) where {op,FR<:Binary}
    rr = isstochastic(ρ) ? (rng === nothing ? default_rng() : rng) : nothing  # hoisted resolve
    @inbounds for i in eachindex(dest, A)
        dest[i] = apply_op(Val(op), FR, ρ, _drawR(ρ, rr, nothing), decode(A[i]))
    end
    dest
end
function _vmap_scalar!(dest, ::Val{op}, ::Type{FR}, ρ::ProjSpec, A, B;
                       rng::MaybeRNG=nothing) where {op,FR<:Binary}
    rr = isstochastic(ρ) ? (rng === nothing ? default_rng() : rng) : nothing
    @inbounds for i in eachindex(dest, A, B)
        dest[i] = apply_op(Val(op), FR, ρ, _drawR(ρ, rr, nothing), decode(A[i]), decode(B[i]))
    end
    dest
end
# stochastic entry points that thread the caller's rng through the Shape-A dispatchers
function vmap!(dest::AbstractArray{FR}, v::Val, ::Type{FR}, ρ::ProjSpec,
               A::AbstractArray{<:Binary}, B::AbstractArray{<:Binary}, rng::MaybeRNG) where {FR<:Binary}
    isstochastic(ρ) ? _vmap_scalar!(dest, v, FR, ρ, A, B; rng) : vmap!(dest, v, FR, ρ, A, B)
end
function vmap!(dest::AbstractArray{FR}, v::Val, ::Type{FR}, ρ::ProjSpec,
               A::AbstractArray{<:Binary}, rng::MaybeRNG) where {FR<:Binary}
    isstochastic(ρ) ? _vmap_scalar!(dest, v, FR, ρ, A; rng) : vmap!(dest, v, FR, ρ, A)
end

@inline function vmap(op::Symbol, fr::Type{<:Binary}, ρ::ProjSpec, As::AbstractArray...;
                      rng::MaybeRNG=nothing)
    dest = similar(first(As), fr)
    isstochastic(ρ) ? vmap!(dest, Val(op), fr, ρ, As..., rng) : vmap!(dest, Val(op), fr, ρ, As...)
end

# ---- registry-generated array surface for the spec register:
#      Op(fr, ρ, A::AbstractArray...) mirrors the scalar signature
for op in OP_REGISTRY
    op.name === :Convert && continue
    name = op.name
    if op.arity == 1
        @eval $name(fr::Type{<:Binary}, ρ::ProjSpec, A::AbstractArray{<:Binary};
                    rng::MaybeRNG=nothing) = vmap($(QuoteNode(name)), fr, ρ, A; rng)
    elseif op.arity == 2
        @eval $name(fr::Type{<:Binary}, ρ::ProjSpec, A::AbstractArray{<:Binary},
                    B::AbstractArray{<:Binary}; rng::MaybeRNG=nothing) =
            vmap($(QuoteNode(name)), fr, ρ, A, B; rng)
    else
        @eval $name(fr::Type{<:Binary}, ρ::ProjSpec, A::AbstractArray{<:Binary},
                    B::AbstractArray{<:Binary}, C::AbstractArray{<:Binary};
                    rng::MaybeRNG=nothing) = vmap($(QuoteNode(name)), fr, ρ, A, B, C; rng)
    end
end
Convert(fr::Type{<:Binary}, ρ::ProjSpec, A::AbstractArray{<:Binary}) =
    vmap(:Convert, fr, ρ, A)
