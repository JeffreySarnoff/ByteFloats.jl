# ===== packed.jl — sub-byte packed storage (design §6.4, bitops plan K5)
#
# Compute unpacked, store packed, convert at boundaries. `PackedVector{F}` packs
# code points at K = bitwidth(F) bits with 64-bit word granularity; elementwise
# kernels run through unpack → byte kernels → (re)pack — never in-place packed
# arithmetic (a deliberate scope limit). Portable shift/mask splices only; a
# BMI2 pdep/pext tile variant is a tracked follow-up, not a dependency.

"""
    PackedVector{F<:Binary} <: AbstractVector{F}

Bit-packed storage of `F` code points at `bitwidth(F)` bits per element.
Construct with `PackedVector(A::AbstractVector{F})`; recover bytes with
`Vector(pv)` or `collect(pv)`. Memory: ⌈n·K/64⌉ words vs n bytes unpacked.
"""
struct PackedVector{F<:Binary} <: AbstractVector{F}
    data::Vector{UInt64}
    n::Int
end

# Element i at K bits/element begins at bit p = (i-1)K: word w = p ÷ 64 (1-based)
# at bit offset off = p mod 64. Elements with off + K > 64 splice across w, w+1.
@inline _wordpos(K::Int, i::Int) = (((i - 1) * K) >> 6 + 1, ((i - 1) * K) & 63)

function PackedVector(A::AbstractVector{F}) where {F<:Binary}
    K = bitwidth(F)
    n = length(A)
    words = zeros(UInt64, cld(n * K, 64))
    @inbounds for (i, v) in enumerate(A)
        w, off = _wordpos(K, i)
        c = UInt64(codepoint(v))
        words[w] |= c << off
        if _crosses_word(off, K)                          # cross-word splice
            words[w + 1] |= c >> (64 - off)
        end
    end
    PackedVector{F}(words, n)
end

Base.size(pv::PackedVector) = (pv.n,)
Base.@propagate_inbounds function Base.getindex(pv::PackedVector{F}, i::Int) where {F}
    @boundscheck checkbounds(pv, i)
    K = bitwidth(F)
    mask = _codemask(K)
    w, off = _wordpos(K, i)
    c = @inbounds pv.data[w] >> off
    if _crosses_word(off, K)
        c |= @inbounds(pv.data[w + 1]) << (64 - off)
    end
    rawvalue(F, UInt8(c & mask))
end
Base.@propagate_inbounds function Base.setindex!(pv::PackedVector{F}, v::F, i::Int) where {F}
    @boundscheck checkbounds(pv, i)
    K = bitwidth(F)
    mask = _codemask(K)
    w, off = _wordpos(K, i)
    c = UInt64(codepoint(v))
    @inbounds pv.data[w] = (pv.data[w] & ~(mask << off)) | (c << off)
    if _crosses_word(off, K)
        hi = K - (64 - off)                               # bits spilling into word w+1
        @inbounds pv.data[w + 1] = (pv.data[w + 1] & ~((UInt64(1) << hi) - 1)) | (c >> (64 - off))
    end
    pv
end

# Two facts every packed access needs: the low-K code mask, and whether an element
# starting at bit `off` spills into the next word.
@inline _codemask(K::Int) = UInt64((1 << K) - 1)
@inline _crosses_word(off::Int, K::Int) = off + K > 64

const _PACK_TILE = 256

"""Unpack `pv[first:first+len-1]` into `buf` (byte scratch); tile-granular kernel entry."""
function unpack_tile!(buf::AbstractVector{F}, pv::PackedVector{F}, first::Int, len::Int) where {F}
    @inbounds for j in 1:len
        buf[j] = pv[first + j - 1]
    end
    buf
end

"""
    vmap(op, fr, ρ, pv::PackedVector; rng=nothing) -> Vector{fr}

Elementwise `op` over packed storage: unpack a `$( _PACK_TILE)`-element tile into a
byte scratch, run the ordinary byte kernel (`vmap!`), emit, repeat. Never computes
on packed words directly — the deliberate compute-unpacked/store-packed boundary.
"""
function vmap(op::Symbol, fr::Type{<:Binary}, ρ::ProjSpec, pv::PackedVector{F};
              rng::MaybeRNG=nothing) where {F<:Binary}
    _vmap_packed(op, fr, ρ, pv, rng)
end

# Function barrier (module docstring's specialization rule): with the result
# format as a type parameter the tile loop specializes, so `out` is a concrete
# Vector and the views reach the `vmap!` methods statically. JET's *package*
# analysis still reports this method — it widens correlated type parameters, so
# it sees `SubArray{Any}` where a specialized call has `SubArray{fr}`; the
# concrete-call gate in test/quality.jl is what verifies the real path.
function _vmap_packed(op::Symbol, ::Type{fr}, ρ::ProjSpec, pv::PackedVector{F},
                      rng::MaybeRNG) where {fr<:Binary,F<:Binary}
    out = Vector{fr}(undef, pv.n)
    isempty(pv) && return out
    buf = Vector{F}(undef, _PACK_TILE)
    i = 1
    while i <= pv.n
        len = min(_PACK_TILE, pv.n - i + 1)
        unpack_tile!(buf, pv, i, len)
        bv = view(buf, 1:len)
        dest = view(out, i:i + len - 1)
        isstochastic(ρ) ? vmap!(dest, Val(op), fr, ρ, bv, rng) : vmap!(dest, Val(op), fr, ρ, bv)
        i += len
    end
    out
end
