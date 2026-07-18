# ===== tables.jl — table lifecycle for pure-ρ specializations (design §7.3/§8.3, architecture §7)
#
# Any *pure* (non-stochastic) specialization of a unary or binary operation is a
# total function on ≤ 2^K / 2^(K1+K2) code points: a ≤256-byte / ≤64 KiB table.
# The builder walks every input through the scalar path (oracle-backed), so a
# table entry IS the defined result — no residual correctness reasoning at use
# sites. Stochastic ρ is never tabulable (the result is a distribution over R);
# builders reject it loudly.
#
# Cache: value-keyed Dict under a ReentrantLock; kernels fetch via a @noinline
# getter once per array call (hoisted), then index a local Memory{UInt8}.

"""Value key for the table cache: op + format parameters + ρ parameters (all values,
so the cache itself is type-stable and non-allocating to query)."""
struct TableKey
    op::Symbol
    fr::NTuple{4,Int}          # result format (K, P, signed, extended)
    f1::NTuple{4,Int}          # first operand format
    f2::NTuple{4,Int}          # second operand format; (0,0,0,0) for unary
    rm::Symbol                 # rounding-mode type name
    sm::Symbol                 # saturation-mode type name
end

_fkey(::Type{Binary{K,P,S,E}}) where {K,P,S,E} = (K, P, Int(S), Int(E))
_rmname(ρ::ProjSpec) = nameof(typeof(roundingmode(ρ)))
_smname(ρ::ProjSpec) = nameof(typeof(saturationmode(ρ)))

const TABLE_CACHE = Dict{TableKey,Memory{UInt8}}()
const TABLE_LOCK = ReentrantLock()

"""Total bytes currently held by the table cache."""
table_bytes() = lock(() -> sum(length, values(TABLE_CACHE); init=0), TABLE_LOCK)
"""Drop every cached table (they rebuild lazily on next use)."""
empty_tables!() = lock(() -> (empty!(TABLE_CACHE); nothing), TABLE_LOCK)

@inline function _scalar_code(op::Val{name}, fr::Type{<:Binary}, ρ::ProjSpec, xs::Float64...) where {name}
    name === :Convert ? codepoint(project(fr, ρ, xs[1])) :
                        codepoint(apply_op(op, fr, ρ, 0, xs...))
end

function _build_unary(op::Symbol, fr::Type{<:Binary}, f1::Type{<:Binary}, ρ::ProjSpec)
    K1 = bitwidth(f1)
    tbl = Memory{UInt8}(undef, 1 << K1)
    V = Val(op)
    for c in 0:(1 << K1) - 1
        tbl[c + 1] = _scalar_code(V, fr, ρ, decode(rawvalue(f1, UInt8(c))))
    end
    tbl
end

function _build_binary(op::Symbol, fr::Type{<:Binary}, f1::Type{<:Binary}, f2::Type{<:Binary}, ρ::ProjSpec)
    K1, K2 = bitwidth(f1), bitwidth(f2)
    tbl = Memory{UInt8}(undef, 1 << (K1 + K2))
    V = Val(op)
    for c1 in 0:(1 << K1) - 1
        x1 = decode(rawvalue(f1, UInt8(c1)))
        base = c1 << K2
        for c2 in 0:(1 << K2) - 1
            tbl[base + c2 + 1] = _scalar_code(V, fr, ρ, x1, decode(rawvalue(f2, UInt8(c2))))
        end
    end
    tbl
end

# Double-checked pattern: probe under lock, build OUTSIDE the lock (builds may run
# MPFR escalations), insert under lock; a racing duplicate build is benign and rare.
@noinline function get_table(op::Symbol, fr::Type{<:Binary}, f1::Type{<:Binary}, ρ::ProjSpec)::Memory{UInt8}
    isstochastic(ρ) && throw(ArgumentError("stochastic ρ $ρ is not tabulable (design §5.4)"))
    key = TableKey(op, _fkey(fr), _fkey(f1), (0, 0, 0, 0), _rmname(ρ), _smname(ρ))
    t = lock(() -> get(TABLE_CACHE, key, nothing), TABLE_LOCK)
    t !== nothing && return t
    built = _build_unary(op, fr, f1, ρ)
    lock(() -> get!(TABLE_CACHE, key, built), TABLE_LOCK)
end
@noinline function get_table(op::Symbol, fr::Type{<:Binary}, f1::Type{<:Binary}, f2::Type{<:Binary}, ρ::ProjSpec)::Memory{UInt8}
    isstochastic(ρ) && throw(ArgumentError("stochastic ρ $ρ is not tabulable (design §5.4)"))
    key = TableKey(op, _fkey(fr), _fkey(f1), _fkey(f2), _rmname(ρ), _smname(ρ))
    t = lock(() -> get(TABLE_CACHE, key, nothing), TABLE_LOCK)
    t !== nothing && return t
    built = _build_binary(op, fr, f1, f2, ρ)
    lock(() -> get!(TABLE_CACHE, key, built), TABLE_LOCK)
end
