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

# ---- ternary tables (bitwidth-specific FMA/FAA/Clamp optimization) -----------
# A pure-ρ ternary specialization is a total function on 2^(K1+K2+K3) code
# points. That is 512 B at K=3 and 256 KiB at K=6 — cheap, cache-friendly, built
# eagerly on first array use — but 2 MiB at K=7 (built only for demonstrably hot
# signatures, LRU-evicted under a byte budget) and 16 MiB at K=8 (never built;
# the compute kernel with the sticky-head escalation serves that band).

"""Value key for a ternary table: op + result + three operand formats + ρ."""
struct TernaryKey
    op::Symbol
    fr::NTuple{4,Int}
    f1::NTuple{4,Int}
    f2::NTuple{4,Int}
    f3::NTuple{4,Int}
    rm::Symbol
    sm::Symbol
end
mutable struct TernaryEntry
    const tbl::Memory{UInt8}
    tick::Int                  # LRU stamp (monotone; larger = more recent)
end

const TERNARY_CACHE = Dict{TernaryKey,TernaryEntry}()
const TERNARY_USE   = Dict{TernaryKey,Int}()   # cumulative elements seen (adaptive band)
const TERNARY_TICK  = Ref(0)

"""K1+K2+K3 up to which a ternary table is built eagerly on first array call
(default 18 bits = 256 KiB — covers every all-K≤6 signature)."""
const TERNARY_EAGER_BITS = Ref(18)
"""K1+K2+K3 up to which a ternary table may be built adaptively once the
signature has processed `TERNARY_BUILD_ELEMS[]` elements (default 21 bits = 2 MiB
— the K=7 band). Above this, ternary ops always run the compute kernel."""
const TERNARY_ADAPTIVE_BITS = Ref(21)
"""Cumulative element count at which an adaptive-band signature earns its table."""
const TERNARY_BUILD_ELEMS = Ref(2_000_000)
"""Byte budget for the ternary cache; least-recently-used tables evict first."""
const TERNARY_CACHE_BYTES = Ref(32 * 1024 * 1024)

_ternary_bytes_locked() = sum(e -> length(e.tbl), values(TERNARY_CACHE); init=0)

"""Total bytes currently held by the table caches (unary/binary + ternary)."""
table_bytes() = lock(() -> sum(length, values(TABLE_CACHE); init=0) + _ternary_bytes_locked(),
                     TABLE_LOCK)
"""Drop every cached table and adaptive counter (they rebuild lazily on next use)."""
empty_tables!() = lock(TABLE_LOCK) do
    empty!(TABLE_CACHE); empty!(TERNARY_CACHE); empty!(TERNARY_USE)
    nothing
end

# One table entry = one trip through the scalar path. Convert is the sole op with
# no ω-semantics (registry group :conv): it is a bare projection of the decoded
# operand, so it bypasses apply_op. R = 0 is safe: stochastic ρ never reaches here.
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
    X2 = _decode_table(f2)
    for c1 in 0:(1 << K1) - 1
        x1 = decode(rawvalue(f1, UInt8(c1)))
        base = c1 << K2
        for c2 in 0:(1 << K2) - 1
            tbl[base + c2 + 1] = _scalar_code(V, fr, ρ, x1, X2[c2 + 1])
        end
    end
    tbl
end

"""
    get_table(op, fr, f1, [f2,] ρ) -> Memory{UInt8}

Fetch (building and caching on first use) the complete result table for the pure-ρ
specialization `op⟨fr; f1[, f2]⟩ under ρ`: entry `c + 1` (unary) or
`(c1 << K2) + c2 + 1` (binary) holds the result code point for those operand
code points. Throws for stochastic ρ. `@noinline` by design: kernels call this
once per array operation and index the returned table in their hot loop.
"""
function get_table end

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

# ---- ternary build / fetch / policy ------------------------------------------

function _build_ternary(op::Symbol, fr::Type{<:Binary}, f1::Type{<:Binary},
                        f2::Type{<:Binary}, f3::Type{<:Binary}, ρ::ProjSpec)
    K1, K2, K3 = bitwidth(f1), bitwidth(f2), bitwidth(f3)
    tbl = Memory{UInt8}(undef, 1 << (K1 + K2 + K3))
    V = Val(op)
    X2, X3 = _decode_table(f2), _decode_table(f3)
    for c1 in 0:(1 << K1) - 1
        x1 = decode(rawvalue(f1, UInt8(c1)))
        for c2 in 0:(1 << K2) - 1
            x2 = X2[c2 + 1]
            base = ((c1 << K2) | c2) << K3
            for c3 in 0:(1 << K3) - 1
                tbl[base + c3 + 1] = _scalar_code(V, fr, ρ, x1, x2, X3[c3 + 1])
            end
        end
    end
    tbl
end

_tkey(op::Symbol, fr, f1, f2, f3, ρ::ProjSpec) =
    TernaryKey(op, _fkey(fr), _fkey(f1), _fkey(f2), _fkey(f3), _rmname(ρ), _smname(ρ))

# Insert under the byte budget, evicting least-recently-used ternary tables first.
# The new table is always inserted (even if alone it exceeds the budget: the caller
# earned it; the budget then simply holds this one table).
function _ternary_insert!(key::TernaryKey, tbl::Memory{UInt8})
    lock(TABLE_LOCK) do
        e = get(TERNARY_CACHE, key, nothing)
        e !== nothing && return e.tbl                        # racing duplicate build
        budget = TERNARY_CACHE_BYTES[]
        while !isempty(TERNARY_CACHE) && _ternary_bytes_locked() + length(tbl) > budget
            victim = argmin(k -> TERNARY_CACHE[k].tick, keys(TERNARY_CACHE))
            delete!(TERNARY_CACHE, victim)
        end
        TERNARY_CACHE[key] = TernaryEntry(tbl, TERNARY_TICK[] += 1)
        tbl
    end
end

"""
    get_table(op, fr, f1, f2, f3, ρ) -> Memory{UInt8}

Ternary form: entry `((c1 << K2 | c2) << K3) + c3 + 1` holds the result code
point. Builds unconditionally (policy lives in `_ternary_table_for`); pure ρ only.
"""
@noinline function get_table(op::Symbol, fr::Type{<:Binary}, f1::Type{<:Binary},
                             f2::Type{<:Binary}, f3::Type{<:Binary}, ρ::ProjSpec)::Memory{UInt8}
    isstochastic(ρ) && throw(ArgumentError("stochastic ρ $ρ is not tabulable (design §5.4)"))
    key = _tkey(op, fr, f1, f2, f3, ρ)
    t = lock(TABLE_LOCK) do
        e = get(TERNARY_CACHE, key, nothing)
        e === nothing ? nothing : (e.tick = (TERNARY_TICK[] += 1); e.tbl)
    end
    t !== nothing && return t
    _ternary_insert!(key, _build_ternary(op, fr, f1, f2, f3, ρ))   # build outside the lock
end

"""
    _ternary_table_for(op, fr, f1, f2, f3, ρ, nelems) -> Union{Nothing,Memory{UInt8}}

The kernel-facing policy gate, called once per array operation with the call's
element count. Eager band (Σ Kᵢ ≤ `TERNARY_EAGER_BITS[]`): fetch/build now.
Adaptive band (≤ `TERNARY_ADAPTIVE_BITS[]`): return the cached table if present,
otherwise accumulate `nelems` against the signature and build only once it has
earned `TERNARY_BUILD_ELEMS[]`. Beyond the adaptive band (the K=8 territory):
always `nothing` — the compute kernel is the right tradeoff there.
"""
@noinline function _ternary_table_for(op::Symbol, fr::Type{<:Binary}, f1::Type{<:Binary},
                                      f2::Type{<:Binary}, f3::Type{<:Binary}, ρ::ProjSpec,
                                      nelems::Int)::Union{Nothing,Memory{UInt8}}
    isstochastic(ρ) && return nothing
    SB = bitwidth(f1) + bitwidth(f2) + bitwidth(f3)
    SB <= TERNARY_EAGER_BITS[] && return get_table(op, fr, f1, f2, f3, ρ)
    SB <= TERNARY_ADAPTIVE_BITS[] || return nothing
    key = _tkey(op, fr, f1, f2, f3, ρ)
    hit = lock(TABLE_LOCK) do
        e = get(TERNARY_CACHE, key, nothing)
        e === nothing ? nothing : (e.tick = (TERNARY_TICK[] += 1); e.tbl)
    end
    hit !== nothing && return hit
    n = lock(() -> (TERNARY_USE[key] = get(TERNARY_USE, key, 0) + nelems), TABLE_LOCK)
    n >= TERNARY_BUILD_ELEMS[] || return nothing
    _ternary_insert!(key, _build_ternary(op, fr, f1, f2, f3, ρ))
end
