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

function PackedVector(A::AbstractVector{F}) where {F<:Binary}
    K = bitwidth(F)
    n = length(A)
    words = zeros(UInt64, max(1, cld(n * K, 64)))
    @inbounds for (i, v) in enumerate(A)
        p = (i - 1) * K
        w = (p >> 6) + 1
        off = p & 63
        c = UInt64(codepoint(v))
        words[w] |= c << off
        if off + K > 64                                   # cross-word splice
            words[w + 1] |= c >> (64 - off)
        end
    end
    PackedVector{F}(words, n)
end

Base.size(pv::PackedVector) = (pv.n,)
Base.@propagate_inbounds function Base.getindex(pv::PackedVector{F}, i::Int) where {F}
    @boundscheck checkbounds(pv, i)
    K = bitwidth(F)
    mask = UInt64((1 << K) - 1)
    p = (i - 1) * K
    w = (p >> 6) + 1
    off = p & 63
    c = @inbounds pv.data[w] >> off
    if off + K > 64
        c |= @inbounds(pv.data[w + 1]) << (64 - off)
    end
    rawvalue(F, UInt8(c & mask))
end
Base.@propagate_inbounds function Base.setindex!(pv::PackedVector{F}, v::F, i::Int) where {F}
    @boundscheck checkbounds(pv, i)
    K = bitwidth(F)
    mask = UInt64((1 << K) - 1)
    p = (i - 1) * K
    w = (p >> 6) + 1
    off = p & 63
    c = UInt64(codepoint(v))
    @inbounds pv.data[w] = (pv.data[w] & ~(mask << off)) | (c << off)
    if off + K > 64
        hi = K - (64 - off)
        @inbounds pv.data[w + 1] = (pv.data[w + 1] & ~((UInt64(1) << hi) - 1)) | (c >> (64 - off))
    end
    pv
end

const _PACK_TILE = 256

"""Unpack `pv[first:first+len-1]` into `buf` (byte scratch); tile-granular kernel entry."""
function unpack_tile!(buf::AbstractVector{F}, pv::PackedVector{F}, first::Int, len::Int) where {F}
    @inbounds for j in 1:len
        buf[j] = pv[first + j - 1]
    end
    buf
end

# elementwise kernels through packed storage: tile-unpack → byte Shape-A/B → emit bytes
function vmap(op::Symbol, fr::Type{<:Binary}, ρ::ProjSpec, pv::PackedVector{F};
              rng::MaybeRNG=nothing) where {F<:Binary}
    out = Vector{fr}(undef, pv.n)
    buf = Vector{F}(undef, _PACK_TILE)
    seg = Vector{fr}(undef, _PACK_TILE)
    i = 1
    while i <= pv.n
        len = min(_PACK_TILE, pv.n - i + 1)
        unpack_tile!(buf, pv, i, len)
        bv = view(buf, 1:len)
        sv = view(seg, 1:len)
        isstochastic(ρ) ? vmap!(sv, Val(op), fr, ρ, bv, rng) : vmap!(sv, Val(op), fr, ρ, bv)
        copyto!(out, i, seg, 1, len)
        i += len
    end
    out
end
