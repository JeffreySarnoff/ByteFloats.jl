# ===== kernels.jl — array kernels (design §8, architecture §7)
#
# Two hot-loop shapes, written once and instantiated by the registry:
#   Shape A (gather):  out[i] = tbl[key(a[i]…)]   — pure ρ, arity ≤ 2
#   Shape B (compute): per-element scalar path    — untabled ternary + stochastic ρ
# The table getter is @noinline and called ONCE per array call, hoisted out of
# the loop; loop bodies index a local Memory{UInt8} with no dict lookups, locks,
# or global loads per element. Ternary ops ride Shape A when the bitwidth policy
# grants a table (eager ≤ 2^18 entries, adaptive ≤ 2^21 — tables.jl) and an
# optionally threaded Shape B above that (the K=8 band).

"""
    vmap!(dest, Val(op), fr, ρ, A[, B[, C]]; rng) -> dest
    vmap(op, fr, ρ, A...; rng)

Elementwise draft operation over arrays of `Binary` values, projecting into `fr`
under ρ. Pure ρ runs the Shape-A gather whenever a table exists or the ternary
bitwidth policy grants one; stochastic specializations (and untabled ternary
signatures) run the scalar path per element (each element drawing its own R).
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

# ---- ternary: Shape-A gather where the policy grants a table (eager K≤6 band,
# adaptive K=7 band — see tables.jl), Shape-B compute otherwise. The pure-ρ
# compute loop is deterministic per element, so it threads for long arrays;
# stochastic ρ stays sequential (single rng stream, reproducible draws).
"""Minimum element count before the pure-ρ ternary compute loop threads."""
const THREAD_MIN_ELEMS = Ref(1 << 15)
"""Master switch for threaded ternary compute loops."""
const THREADED_KERNELS = Ref(true)

function vmap!(dest::AbstractArray{FR}, ::Val{op}, ::Type{FR}, ρ::ProjSpec,
               A::AbstractArray{F1}, B::AbstractArray{F2}, C::AbstractArray{F3};
               rng::MaybeRNG=nothing) where {op,FR<:Binary,F1<:Binary,F2<:Binary,F3<:Binary}
    axes(dest) == axes(A) == axes(B) == axes(C) || throw(DimensionMismatch("operand axes must match"))
    if !isstochastic(ρ)
        tbl = _ternary_table_for(op, FR, F1, F2, F3, ρ, length(dest))  # hoisted; @noinline
        if tbl !== nothing
            K2, K3 = bitwidth(F2), bitwidth(F3)
            @inbounds for i in eachindex(dest, A, B, C)
                idx = ((Int(codepoint(A[i])) << K2 | Int(codepoint(B[i]))) << K3) +
                      Int(codepoint(C[i])) + 1
                dest[i] = rawvalue(FR, tbl[idx])
            end
            return dest
        end
        inds = eachindex(dest, A, B, C)
        if THREADED_KERNELS[] && Threads.nthreads() > 1 &&
           length(inds) >= THREAD_MIN_ELEMS[] && inds isa AbstractUnitRange
            Threads.@threads for i in inds
                @inbounds dest[i] = apply_op(Val(op), FR, ρ, 0,
                                             decode(A[i]), decode(B[i]), decode(C[i]))
            end
            return dest
        end
        @inbounds for i in inds
            dest[i] = apply_op(Val(op), FR, ρ, 0, decode(A[i]), decode(B[i]), decode(C[i]))
        end
        return dest
    end
    rr = rng === nothing ? default_rng() : rng
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
function vmap!(dest::AbstractArray{FR}, v::Val, ::Type{FR}, ρ::ProjSpec,
               A::AbstractArray{<:Binary}, B::AbstractArray{<:Binary},
               C::AbstractArray{<:Binary}, rng::MaybeRNG) where {FR<:Binary}
    vmap!(dest, v, FR, ρ, A, B, C; rng)      # the ternary method handles both ρ kinds
end

@inline function vmap(op::Symbol, fr::Type{<:Binary}, ρ::ProjSpec, As::AbstractArray...;
                      rng::MaybeRNG=nothing)
    dest = similar(first(As), fr)
    # Stochastic ρ appends rng as a trailing *positional* argument, selecting the
    # rng-threading vmap! methods above; pure ρ takes the plain Shape-A/B methods.
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
# Convert has no ω-semantics (registry group :conv) so the loop above skips it;
# its array form still rides the Shape-A gather — the :Convert table is built by
# _scalar_code's bare-projection branch. Pure ρ only (stochastic Convert is per-R).
Convert(fr::Type{<:Binary}, ρ::ProjSpec, A::AbstractArray{<:Binary}) =
    vmap(:Convert, fr, ρ, A)
