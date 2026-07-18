# ===== decode_encode.jl — ωDecode / ωEncode / ordering keys / Next ops (design §3)

"""
    decode(v::Binary) -> Float64

ωDecode (draft §4.7.2). Float64 is the universal exact carrier for K ≤ 8 datums
(≤ 8-bit significands, |exponent| ≲ 2^7); exactness is asserted by exhaustive test.

Implemented as a constant-tuple lookup (bitops plan K2): the per-format table is
generated once from the computational decode below, so the two are correct by
construction and asserted equivalent exhaustively; constant inputs still fold.
"""
@inline function _decode_compute(v::Binary{K,P,SGN,EXT})::Float64 where {K,P,SGN,EXT}
    c = codepoint(v)
    if SGN
        c == UInt8(1 << (K - 1)) && return NaN
        if EXT
            c == UInt8((1 << (K - 1)) - 1) && return Inf
            c == UInt8((1 << K) - 1) && return -Inf
        end
    else
        c == UInt8((1 << K) - 1) && return NaN
        EXT && c == UInt8((1 << K) - 2) && return Inf
    end
    neg = SGN && (c >= UInt8(1 << (K - 1)))
    m = neg ? c - UInt8(1 << (K - 1)) : c
    T = m & UInt8((1 << (P - 1)) - 1)
    Eb = Int(m >> (P - 1))
    B = SGN ? (1 << (K - P - 1)) : (1 << (K - P))
    sig = Eb == 0 ? Int(T) : Int(T) + (1 << (P - 1))
    e = (Eb == 0 ? 1 : Eb) - B + (1 - P)
    # bit assembly (bitops plan K2): every datum exponent is deep inside Float64's
    # normal range (|e + nb − 1| ≤ ~260), so no subnormal/overflow cases exist and
    # ldexp's generality is pure overhead.
    sig == 0 && return 0.0
    nb = 64 - leading_zeros(UInt64(sig))
    mant = (UInt64(sig) << (53 - nb)) & ((UInt64(1) << 52) - 1)
    bits = (UInt64(e + nb - 1 + 1023) << 52) | mant
    neg && (bits |= UInt64(1) << 63)
    return reinterpret(Float64, bits)
end

@generated function _decode_table(::Type{Binary{K,P,S,E}}) where {K,P,S,E}
    t = ntuple(i -> _decode_compute(rawvalue(Binary{K,P,S,E}, UInt8(i - 1))), 1 << K)
    :($t)
end

@inline decode(v::Binary{K,P,SGN,EXT}) where {K,P,SGN,EXT} =
    @inbounds _decode_table(Binary{K,P,SGN,EXT})[Int(codepoint(v)) + 1]



"""
    encode(T, sign, S, Q) -> UInt8   (private; design §3.3)

ωEncode from canonical integer form: value = sign · S · 2^Q, with S ∈ 0:2^P
(S == 2^P is the next-binade carry, draft §4.7.4 NOTE 4). Precondition: the value
is in the datum set of `T` (guaranteed by RoundToPrecision ∘ Saturate).
"""
@inline function encode(::Type{Binary{K,P,SGN,EXT}}, sign::Int, S::Int64, Q::Int64) where {K,P,SGN,EXT}
    S == 0 && return 0x00
    if S == (Int64(1) << P)                    # carry into next binade
        S = Int64(1) << (P - 1); Q += 1
    end
    B = SGN ? (1 << (K - P - 1)) : (1 << (K - P))
    local c::UInt8
    if S < (Int64(1) << (P - 1))               # subnormal: Q must equal 2-B-P
        c = UInt8(S)
    else
        E = Int(Q) + P - 1
        Eb = E + B
        c = UInt8((S & ((Int64(1) << (P - 1)) - 1)) + (Int64(Eb) << (P - 1)))
    end
    (SGN && sign < 0) && (c |= UInt8(1 << (K - 1)))
    c
end

# ---- Total-order key (design §3.1): sign–magnitude → monotone unsigned key.
# NaN (at the −0 slot for signed formats / top code for unsigned) sorts ABOVE +Inf
# [interpretation; draft §4.12.1 text unavailable in upload — see checkpoint].
@inline function order_key(v::Binary{K,P,SGN,EXT}) where {K,P,SGN,EXT}
    c = codepoint(v)
    isnan(v) && return typemax(UInt16)
    if SGN
        neg = c >= UInt8(1 << (K - 1))
        return neg ? UInt16((1 << (K - 1))) - UInt16(c - UInt8(1 << (K - 1))) :
                     UInt16(1 << (K - 1)) + UInt16(c) + UInt16(1)
    else
        return UInt16(c) + UInt16(1)
    end
end

"""TotalOrder⟨fx,fy⟩ (draft §4.12.1): x ≤ y in the total order (single NaN largest).
Same-format comparisons run on the integer order key (bitops plan K1) — proven
equivalent to the decode order exhaustively over all pairs; cross-format keys are
not comparable, so mixed formats keep the decode path."""
@inline TotalOrder(x::T, y::T) where {T<:Binary} = order_key(x) <= order_key(y)
function TotalOrder(x::Binary, y::Binary)
    dx, dy = decode(x), decode(y)
    isnan(dx) && return isnan(dy) ? true : false
    isnan(dy) && return true
    dx <= dy
end
# key strict-< reproduces the old TotalOrder-derived isless exactly, including
# NaN-last (key(NaN) = typemax): isless(x, NaN) = true, isless(NaN, NaN) = false.
Base.isless(x::T, y::T) where {T<:Binary} = order_key(x) < order_key(y)

# Numeric comparisons (NaN unordered; keys are order-isomorphic to datums off NaN)
Base.:(==)(x::T, y::T) where {T<:Binary} = (isnan(x) | isnan(y)) ? false : order_key(x) == order_key(y)
Base.:(<)(x::T, y::T) where {T<:Binary}  = (isnan(x) | isnan(y)) ? false : order_key(x) < order_key(y)
Base.:(<=)(x::T, y::T) where {T<:Binary} = (isnan(x) | isnan(y)) ? false : order_key(x) <= order_key(y)

# ---- counting sort over the key space (bitops plan K1): ≤ 2^K + 1 distinct keys,
# equal keys ⇒ identical code points, so stability is moot; O(n) one-pass counts.
struct CodeCountingSort <: Base.Sort.Algorithm end
Base.Sort.defalg(::AbstractArray{<:Binary}) = CodeCountingSort()
# any ordering we don't specialize falls back to the stock algorithm
Base.sort!(v::AbstractVector{T}, lo::Int, hi::Int, ::CodeCountingSort,
           o::Base.Order.Ordering) where {T<:Binary} =
    sort!(v, lo, hi, Base.Sort.DEFAULT_UNSTABLE, o)
function Base.sort!(v::AbstractVector{T}, lo::Int, hi::Int, ::CodeCountingSort,
                    o::Union{Base.Order.ForwardOrdering,
                             Base.Order.ReverseOrdering{Base.Order.ForwardOrdering}}) where {T<:Binary}
    K = bitwidth(T)
    nk = (1 << K) + 1                              # keys 1..2^K plus NaN sentinel bucket
    counts = zeros(Int, nk + 1)
    key2code = Vector{UInt8}(undef, nk + 1)
    for c in 0x00:UInt8((1 << K) - 1)              # key ↔ code inversion, 2^K iterations
        k = order_key(rawvalue(T, c))
        b = k == typemax(UInt16) ? nk + 1 : Int(k)
        key2code[b] = c
    end
    @inbounds for i in lo:hi
        k = order_key(v[i])
        counts[k == typemax(UInt16) ? nk + 1 : Int(k)] += 1
    end
    rev = o isa Base.Order.ReverseOrdering
    i = rev ? hi : lo
    step = rev ? -1 : 1
    @inbounds for b in 1:nk + 1
        c = counts[b]
        c == 0 && continue
        # ascending buckets emitted backward under Reverse puts the NaN bucket
        # (largest key) at the front — exactly Base's rev=true isless semantics
        val = rawvalue(T, key2code[b])
        for _ in 1:c
            v[i] = val
            i += step
        end
    end
    v
end

# ---- Class (draft §4.13.1)
@enum FPClass::UInt8 ClassNaN ClassNegInf ClassNegNormal ClassNegSubnormal ClassZero ClassPosSubnormal ClassPosNormal ClassPosInf
function Class(v::Binary)
    isnan(v) && return ClassNaN
    d = decode(v)
    d == Inf && return ClassPosInf
    d == -Inf && return ClassNegInf
    iszero(v) && return ClassZero
    sub = issubnormal_p3109(v)
    if d > 0
        return sub ? ClassPosSubnormal : ClassPosNormal
    else
        return sub ? ClassNegSubnormal : ClassNegNormal
    end
end

# ---- NextGreaterThan / NextLessThan (draft §4.16): ±1 steps on magnitude code points
function NextGreaterThan(v::T) where {K,P,SGN,EXT,T<:Binary{K,P,SGN,EXT}}
    isnan(v) && return v
    c = codepoint(v)
    maxfin = codepoint(MaxFiniteOf(T))
    if EXT
        c == posinf_code(T) && return rawvalue(T, nan_code(T))         # Inf → NaN
        SGN && c == neginf_code(T) && return MinFiniteOf(T)            # -Inf → MinFinite
    end
    if !EXT && c == maxfin
        return rawvalue(T, nan_code(T))                                # Finite: MaxFinite → NaN
    end
    if SGN && signbit(v)
        c == (signmask(T) | 0x01) && return zero(T)                    # SmallestNegative → 0
        return rawvalue(T, c - 0x01)
    end
    rawvalue(T, c + 0x01)
end
function NextLessThan(v::T) where {K,P,SGN,EXT,T<:Binary{K,P,SGN,EXT}}
    isnan(v) && return v
    c = codepoint(v)
    if !SGN
        c == 0x00 && return rawvalue(T, nan_code(T))
        EXT && c == posinf_code(T) && return MaxFiniteOf(T)
        (!EXT && c == codepoint(MinFiniteOf(T))) && return rawvalue(T, nan_code(T))
        return rawvalue(T, c - 0x01)
    end
    if EXT
        c == neginf_code(T) && return rawvalue(T, nan_code(T))          # -Inf → NaN
        c == codepoint(MinFiniteOf(T)) && return rawvalue(T, neginf_code(T))
        c == posinf_code(T) && return MaxFiniteOf(T)
    else
        c == codepoint(MinFiniteOf(T)) && return rawvalue(T, nan_code(T))
    end
    c == 0x00 && return rawvalue(T, signmask(T) | 0x01)                 # 0 → SmallestNegative
    signbit(v) ? rawvalue(T, c + 0x01) : rawvalue(T, c - 0x01)
end
Base.nextfloat(v::Binary) = NextGreaterThan(v)
Base.prevfloat(v::Binary) = NextLessThan(v)
