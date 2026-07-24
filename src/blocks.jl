# ===== blocks.jl — blocks, scaled operations, reductions (design §6, architecture §8)
#
# A block is (scale, elements): draft §5. The elementwise §5.4 schema is
# implemented ONCE (decode all operand blocks → ω-op lanewise → BlockProject
# against the supplied result scale), and Block*/Scaled* variants for every
# registry operation are generated from it — the non-divergence mechanism.
#
# Correctness of the division-by-result-scale composition, per result kind:
#   Float64 exact V   : one correctly rounded directed division IS the enclosure
#                       (collapses to exact when V/S is representable, e.g. every
#                       P=1 power-of-two scale — the draft's own NOTE — for free)
#   Float128 exact V  : exact BigFloat conversion (≥128 bits), then one directed division
#   BigExactF         : evaluated exactly (2200 bits); dividing at ≥ that precision
#                       makes the conversion exact, so again one directed division
#   EncloseF (d,u)    : interval-divided by the exact scale, sign-aware
#   Enclose128F       : its MPFR fallback `f` is an EncloseF-shaped ladder — reuse it
# Reductions use an exact wide-precision accumulator (precision chosen from the
# operand structure ⇒ provably exact); the register-resident superaccumulator is
# the tracked Phase-3 optimization (checkpoint.md).

"""
    Block{B, FS<:Binary, FE<:Binary}

A P3109 block `(s, [x₁ … x_B])`: scale factor in format `FS`, `B ≥ 1` elements in
format `FE` (draft §5). `isbits`; blocks live in registers/stack.
"""
struct Block{B,FS<:Binary,FE<:Binary}
    s::FS
    x::NTuple{B,FE}
    # The element type is spelled out on the first lane rather than as
    # `NTuple{B,FE}`: an empty tuple would leave FE unbound (Aqua's unbound-type-
    # parameter gate). B ≥ 1 is therefore a signature property, and the empty
    # case falls to the explicit method below.
    function Block(s::FS, x::Tuple{FE,Vararg{FE,Bm1}}) where {Bm1,FS<:Binary,FE<:Binary}
        new{Bm1 + 1,FS,FE}(s, x)
    end
end
Block(::Binary, ::Tuple{}) = throw(ArgumentError("block size B must be ≥ 1"))
Block(s::Binary, xs::Binary...) = Block(s, xs)
blocksize(::Block{B}) where {B} = B
scaleformat(::Block{B,FS,FE}) where {B,FS,FE} = FS
elemformat(::Block{B,FS,FE}) where {B,FS,FE} = FE

"""ωBlockDecode (draft §5.1.1): per lane ωMultiply(decode(s), decode(xᵢ)) — exact Float64."""
@inline function blockdecode(b::Block{B}) where {B}
    S = decode(b.s)
    ntuple(i -> ωeval(Val(:Multiply), S, decode(b.x[i]))::Float64, Val(B))
end

# ---- ωBlockProject element pipeline (draft §5.1.2)
_res_sign(v::Float64) = sign(v)
_res_sign(v::Float128) = Float64(sign(v))
_res_sign(b::BigExactF) = Float64(sign(b.f()))
function _res_sign(e::EncloseF)
    d, u = e.f(256)
    s = sign(d)
    iszero(s) ? sign(u) : s
end
function _res_sign(z::Enclose128F)
    s = Float64(sign(z.lo))
    iszero(s) ? Float64(sign(z.hi)) : s
end
_res_isnan(v::Float64) = isnan(v)
_res_isnan(v::Float128) = isnan(v)
_res_isnan(::BigExactF) = false      # exact-arithmetic escalations are finite by construction
_res_isnan(::EncloseF) = false
_res_isnan(::Enclose128F) = false

# enclosure of (value of res)/S for finite nonzero exact S
function _encl_div_scale(res::Float64, S::Float64)
    EncloseF(prec -> setprecision(BigFloat, prec) do
        (setrounding(() -> BigFloat(res) / S, BigFloat, RoundDown),
         setrounding(() -> BigFloat(res) / S, BigFloat, RoundUp))
    end)
end
function _encl_div_scale(res::BigExactF, S::Float64)
    EncloseF(prec -> begin
        v = res.f()
        p = max(prec, Base.precision(v) + 8)                 # exact conversion ⇒ one rounded op
        setprecision(BigFloat, p) do
            (setrounding(() -> BigFloat(v) / S, BigFloat, RoundDown),
             setrounding(() -> BigFloat(v) / S, BigFloat, RoundUp))
        end
    end)
end
function _encl_div_scale(res::EncloseF, S::Float64)
    EncloseF(prec -> begin
        d, u = res.f(prec)
        setprecision(BigFloat, prec) do
            if S > 0
                (setrounding(() -> d / S, BigFloat, RoundDown),
                 setrounding(() -> u / S, BigFloat, RoundUp))
            else
                (setrounding(() -> u / S, BigFloat, RoundDown),
                 setrounding(() -> d / S, BigFloat, RoundUp))
            end
        end
    end)
end
function _encl_div_scale(res::Float128, S::Float64)
    EncloseF(prec -> begin
        p = max(prec, 128)                                    # exact conversion ⇒ one rounded op
        setprecision(BigFloat, p) do
            b = BigFloat(res)
            (setrounding(() -> b / S, BigFloat, RoundDown),
             setrounding(() -> b / S, BigFloat, RoundUp))
        end
    end)
end
_encl_div_scale(res::Enclose128F, S::Float64) = _encl_div_scale(EncloseF(res.f), S)

"""One element of ωBlockProject: the draft's S-special rows, then ωDivide + ωProject.

Fast-path cascade after the special rows, ordered cheapest-first by result kind:
exact Float64 quotient → CR Float128 quotient (exact or one-ulp bracket) →
CR-divided Enclose128F bracket → fq-composition envelope (2^-89) — each resolved
by the two-sided sticky gate; any miss falls through to the rigorous MPFR
interval at the bottom."""
function _bp_element(fr::Type{<:Binary}, ρ::ProjSpec, R::Int, res, Sdat::Float64)
    (isnan(Sdat) || _res_isnan(res)) && return rawvalue(fr, nan_code(fr))
    iszero(Sdat) && return project(fr, ρ, 0.0; R)
    if isinf(Sdat)
        return project(fr, ρ, _res_sign(res) * sign(Sdat); R)     # sgn(Xᵢ)·sgn(S) ∈ {−1,0,1}
    end
    if res isa Float64
        isinf(res) && return project(fr, ρ, sign(Sdat) * res; R)  # ωDivide(±∞, finite)
        iszero(res) && return project(fr, ρ, 0.0; R)
        q = res / Sdat
        (isfinite(q) && fma(q, Sdat, -res) == 0.0) && return project(fr, ρ, q; R)
        if _f128()                                                 # plan Site E, Class R:
            qr, qs = Float128(res), Float128(Sdat)                 # CR Float128 division
            q128 = qr / qs
            fma(q128, qs, -qr) == 0 && return project(fr, ρ, q128; R)
            return _finish(fr, ρ, R,
                           Enclose128F(prevfloat(q128), nextfloat(q128),
                                       _encl_div_scale(res, Sdat).f))
        end
    elseif res isa Float128 && _f128()
        q = res / Float128(Sdat)                              # CR division of an exact value
        fma(q, Float128(Sdat), -res) == 0 && return project(fr, ρ, q; R)
        return _finish(fr, ρ, R,
                       Enclose128F(prevfloat(q), nextfloat(q), _encl_div_scale(res, Sdat).f))
    elseif res isa Enclose128F && _f128()
        l128, h128 = Float128(Sdat) > 0 ? (res.lo, res.hi) : (res.hi, res.lo)
        lo = prevfloat(l128 / Float128(Sdat))                 # CR half-ulp each ⇒ one outward step
        hi = nextfloat(h128 / Float128(Sdat))
        cd = project(fr, ρ, lo; R, sticky=+1)
        cu = project(fr, ρ, hi; R, sticky=-1)
        codepoint(cd) == codepoint(cu) && return cd
    elseif res isa EncloseF && res.fq !== nothing && _f128()
        # Class E composition: op envelope 2^-90 + CR division half-ulp ⇒ 2^-89 covers
        y = res.fq()::Float128
        if isfinite(y) && !iszero(y)
            q = y / Float128(Sdat)
            E = ldexp(abs(q), _F128_RELEXP + 1)
            cd = project(fr, ρ, q - E; R, sticky=+1)
            cu = project(fr, ρ, q + E; R, sticky=-1)
            codepoint(cd) == codepoint(cu) && return cd
        end
    end
    project_interval(fr, ρ, _encl_div_scale(res, Sdat).f; R)
end

"""ωBlockProject (draft §5.1.2) over a tuple of result kinds, scale supplied as a value."""
function blockproject(fr::Type{<:Binary}, ρ::ProjSpec, sr::Binary, Z::NTuple{B,Any};
                      rng::MaybeRNG=nothing) where {B}
    Sdat = decode(sr)
    rr = isstochastic(ρ) ? (rng === nothing ? default_rng() : rng) : nothing
    elems = ntuple(i -> _bp_element(fr, ρ, _drawR(ρ, rr, nothing), Z[i], Sdat), Val(B))
    Block(sr, elems)
end

# ---- generated elementwise BlockOp / ScaledOp surface (draft §5.4 / §5.5)
for op in OP_REGISTRY
    op.name === :Convert && continue
    name = op.name; bname = Symbol(:Block, name); sname = Symbol(:Scaled, name); V = Val{name}
    if op.arity == 1
        @eval begin
            function $bname(fr::Type{<:Binary}, ρ::ProjSpec, b1::Block{B}, sr::Binary;
                            rng::MaybeRNG=nothing) where {B}
                X1 = blockdecode(b1)
                Z = ntuple(i -> ωeval($V(), X1[i]), Val(B))
                blockproject(fr, ρ, sr, Z; rng)
            end
            function $sname(fr::Type{<:Binary}, ρ::ProjSpec, s1::Binary, x1::Binary;
                            rng::MaybeRNG=nothing)
                X = ωeval(Val(:Multiply), decode(s1), decode(x1))::Float64
                _bp_element(fr, ρ, _drawR(ρ, rng, nothing), ωeval($V(), X), 1.0)
            end
        end
    elseif op.arity == 2
        @eval begin
            function $bname(fr::Type{<:Binary}, ρ::ProjSpec, b1::Block{B}, b2::Block{B}, sr::Binary;
                            rng::MaybeRNG=nothing) where {B}
                X1 = blockdecode(b1); X2 = blockdecode(b2)
                Z = ntuple(i -> ωeval($V(), X1[i], X2[i]), Val(B))
                blockproject(fr, ρ, sr, Z; rng)
            end
            function $sname(fr::Type{<:Binary}, ρ::ProjSpec, s1::Binary, x1::Binary,
                            s2::Binary, x2::Binary; rng::MaybeRNG=nothing)
                Xa = ωeval(Val(:Multiply), decode(s1), decode(x1))::Float64
                Xb = ωeval(Val(:Multiply), decode(s2), decode(x2))::Float64
                _bp_element(fr, ρ, _drawR(ρ, rng, nothing), ωeval($V(), Xa, Xb), 1.0)
            end
        end
    else
        @eval begin
            function $bname(fr::Type{<:Binary}, ρ::ProjSpec, b1::Block{B}, b2::Block{B},
                            b3::Block{B}, sr::Binary; rng::MaybeRNG=nothing) where {B}
                X1 = blockdecode(b1); X2 = blockdecode(b2); X3 = blockdecode(b3)
                Z = ntuple(i -> ωeval($V(), X1[i], X2[i], X3[i]), Val(B))
                blockproject(fr, ρ, sr, Z; rng)
            end
            function $sname(fr::Type{<:Binary}, ρ::ProjSpec, s1::Binary, x1::Binary, s2::Binary,
                            x2::Binary, s3::Binary, x3::Binary; rng::MaybeRNG=nothing)
                Xa = ωeval(Val(:Multiply), decode(s1), decode(x1))::Float64
                Xb = ωeval(Val(:Multiply), decode(s2), decode(x2))::Float64
                Xc = ωeval(Val(:Multiply), decode(s3), decode(x3))::Float64
                _bp_element(fr, ρ, _drawR(ρ, rng, nothing), ωeval($V(), Xa, Xb, Xc), 1.0)
            end
        end
    end
end

# ---- reductions (draft §5.3): specials by the fold algebra, then an exact sum/product
const _REDPREC = 2400   # ≥ Float64 span (~2100) + 64-bit terms + log₂B slack: exact for any B here

# Span filter (plan Site E, Class R): decoded lanes carry ≤17-bit significands; a
# B-term sum is exactly representable in Float128 when
#   17 + (exponent span) + ⌈log₂B⌉ + 1 ≤ 113.
# One integer pass decides; the BigFloat accumulator remains for the wide tail.
@inline function _expspan(X)
    lo = typemax(Int); hi = typemin(Int)
    for v in X
        iszero(v) && continue
        e = exponent(v)
        lo = min(lo, e); hi = max(hi, e)
    end
    lo > hi ? 0 : hi - lo
end
@inline _log2ceil(B::Int) = 8 * sizeof(Int) - leading_zeros(B - 1 > 0 ? B - 1 : 1)

function _reduce_add_datum(X)
    any(isnan, X) && return NaN
    hasp = any(==(Inf), X); hasn = any(==(-Inf), X)
    (hasp & hasn) && return NaN
    hasp && return Inf
    hasn && return -Inf
    all(iszero, X) && return 0.0
    B = length(X)
    if _f128() && _expspan(X) + _log2ceil(B) <= 92
        acc = Float128(0)
        for v in X
            acc += Float128(v)                               # every partial exact by width
        end
        return acc
    end
    BigExactF(() -> setprecision(() -> sum(BigFloat, X; init=BigFloat(0)), BigFloat, _REDPREC))
end

"""BlockReduceAdd (draft §5.3.1): project(reduce(ωAdd, [0, X…]))."""
function BlockReduceAdd(fr::Type{<:Binary}, ρ::ProjSpec, b::Block;
                        rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
    _finish(fr, ρ, _drawR(ρ, rng, R), _reduce_add_datum(blockdecode(b)))
end

"""BlockReduceMultiply (draft §5.3.1): project(reduce(ωMultiply, [1, X…]))."""
function BlockReduceMultiply(fr::Type{<:Binary}, ρ::ProjSpec, b::Block{B};
                             rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing) where {B}
    X = blockdecode(b)
    res = if any(isnan, X) || (any(iszero, X) && any(isinf, X))
        NaN                                                     # 0·∞ arises in the fold → NaN
    elseif any(isinf, X)
        s = prod(sign, X)
        s * Inf
    elseif any(iszero, X)
        0.0
    elseif _f128() && 16B + 8 <= 112                            # exact product by width (B ≤ 6)
        acc = Float128(1)
        for v in X
            acc *= Float128(v)
        end
        acc
    else
        BigExactF(() -> setprecision(() -> prod(BigFloat, X; init=BigFloat(1)), BigFloat, 16B + 128))
    end
    _finish(fr, ρ, _drawR(ρ, rng, R), res)
end

"""BlockDotProduct (draft §5.3.2). Lane products can carry 64 significant bits, so the
products and their sum are formed in the exact accumulator, with the ∞/NaN fold algebra
resolved on the Float64 classifications first."""
function BlockDotProduct(fr::Type{<:Binary}, ρ::ProjSpec, bx::Block{B}, by::Block{B};
                         rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing) where {B}
    X = blockdecode(bx); Y = blockdecode(by)
    pcls = ntuple(Val(B)) do i
        x, y = X[i], Y[i]
        (isnan(x) | isnan(y)) && return NaN
        ((iszero(x) && isinf(y)) || (isinf(x) && iszero(y))) && return NaN
        (isinf(x) || isinf(y)) && return sign(x) * sign(y) * Inf
        1.0                                                     # finite lane marker
    end
    res = if any(isnan, pcls)
        NaN
    elseif any(isinf, pcls)
        hasp = any(==(Inf), pcls); hasn = any(==(-Inf), pcls)
        (hasp & hasn) ? NaN : (hasp ? Inf : -Inf)
    else
        # product-span filter (plan Site E): lane products carry ≤33-bit significands;
        # the sum is exact in Float128 when 33 + span + ⌈log₂B⌉ + 2 ≤ 113
        spanlo = typemax(Int); spanhi = typemin(Int)
        for i in 1:B
            (iszero(X[i]) || iszero(Y[i])) && continue
            e = exponent(X[i]) + exponent(Y[i])
            spanlo = min(spanlo, e); spanhi = max(spanhi, e)
        end
        span = spanlo > spanhi ? 0 : (spanhi - spanlo) + 1
        if _f128() && span + _log2ceil(B) <= 76
            acc = Float128(0)
            for i in 1:B
                acc += Float128(X[i]) * Float128(Y[i])          # exact products, exact sum
            end
            acc
        else
            BigExactF(() -> setprecision(BigFloat, _REDPREC) do
                acc = BigFloat(0)
                for i in 1:B
                    acc += BigFloat(X[i]) * BigFloat(Y[i])      # exact products, exact sum
                end
                acc
            end)
        end
    end
    _finish(fr, ρ, _drawR(ρ, rng, R), res)
end

# ---- conversion family (draft §5.2)
"""ConvertFromBlock (§5.2.1): blockdecode, then project each element (no scale division)."""
function ConvertFromBlock(fr::Type{<:Binary}, ρ::ProjSpec, b::Block{B};
                          rng::MaybeRNG=nothing) where {B}
    X = blockdecode(b)
    ntuple(i -> project(fr, ρ, X[i]; R=_drawR(ρ, rng, nothing)), Val(B))
end

"""ConvertToBlock (§5.2.2): decode elements, BlockProject against the supplied scale."""
function ConvertToBlock(fs::Type{<:Binary}, fr::Type{<:Binary}, ρ::ProjSpec,
                        xs::NTuple{B,<:Binary}, s::Binary;
                        rng::MaybeRNG=nothing) where {B}
    s isa fs || throw(ArgumentError("scale operand must be in format $fs"))
    blockproject(fr, ρ, s, map(decode, xs); rng)
end

"""ConvertToBlockMaxAbsFinite (§5.2.3): S = reduce(ωMaximumFinite, [NaN, |X₁|…]),
project S under ρs, then BlockProject the elements under ρ."""
function ConvertToBlockMaxAbsFinite(fs::Type{<:Binary}, fr::Type{<:Binary},
                                    ρs::ProjSpec, ρ::ProjSpec, xs::NTuple{B,<:Binary};
                                    rng::MaybeRNG=nothing) where {B}
    X = map(decode, xs)
    M = map(x -> ωeval(Val(:Abs), x)::Float64, X)
    S = foldl((a, m) -> ωeval(Val(:MaximumFinite), a, m)::Float64, M; init=NaN)
    s = project(fs, ρs, S; R=_drawR(ρs, rng, nothing))
    blockproject(fr, ρ, s, X; rng)
end

# ---- SoA array-of-blocks container (design §6.3): scales and elements in planes
"""
    BlockVector{B,FS,FE} <: AbstractVector{Block{B,FS,FE}}

Structure-of-arrays storage: `scales::Vector{FS}`, `elems::Matrix{FE}` (B × n,
column-major so each block's elements are one contiguous column).
"""
struct BlockVector{B,FS<:Binary,FE<:Binary} <: AbstractVector{Block{B,FS,FE}}
    scales::Vector{FS}
    elems::Matrix{FE}
    function BlockVector{B}(scales::Vector{FS}, elems::Matrix{FE}) where {B,FS<:Binary,FE<:Binary}
        size(elems, 1) == B || throw(DimensionMismatch("elems must be $(B)×n"))
        size(elems, 2) == length(scales) || throw(DimensionMismatch("one scale per block"))
        new{B,FS,FE}(scales, elems)
    end
end
BlockVector(blocks::AbstractVector{Block{B,FS,FE}}) where {B,FS,FE} =
    BlockVector{B}([b.s for b in blocks],
                   FE[blocks[j].x[i] for i in 1:B, j in eachindex(blocks)])
Base.size(bv::BlockVector) = (length(bv.scales),)
Base.@propagate_inbounds function Base.getindex(bv::BlockVector{B,FS,FE}, j::Int) where {B,FS,FE}
    Block(bv.scales[j], ntuple(i -> bv.elems[i, j], Val(B)))
end
Base.@propagate_inbounds function Base.setindex!(bv::BlockVector{B}, b::Block{B}, j::Int) where {B}
    bv.scales[j] = b.s
    for i in 1:B
        bv.elems[i, j] = b.x[i]
    end
    bv
end
