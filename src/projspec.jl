# ===== projspec.jl — projection specifications as types (design §4, architecture §2)
#
# A projection specification ρ = (rounding mode, saturation mode) is a zero-size
# value whose *type* carries both modes, so kernels specialize on ρ exactly as
# they specialize on formats. The stochastic modes carry the random-bit budget N
# in the type: a kernel's randomness consumption is a compile-time fact, and
# pure-vs-stochastic dispatch is static (pure ⇒ tabulable, stochastic ⇒ never).

abstract type RoundingMode3109 end

"Round to nearest; ties to the even code point (IEEE default; draft §4.7.1)."
struct NearestTiesToEven <: RoundingMode3109 end
"Round to nearest; ties away from zero (draft §4.7.1)."
struct NearestTiesToAway <: RoundingMode3109 end
"Directed rounding toward +∞ (draft §4.7.2)."
struct TowardPositive    <: RoundingMode3109 end
"Directed rounding toward −∞ (draft §4.7.2)."
struct TowardNegative    <: RoundingMode3109 end
"Directed rounding toward zero — truncation (draft §4.7.2)."
struct TowardZero        <: RoundingMode3109 end
"Round inexact results to the nearest *odd* code point (draft §4.7.3; sticky-friendly)."
struct ToOdd             <: RoundingMode3109 end

# Stochastic variants (draft §4.7.4): R random bits with 0 ≤ R < 2^N are supplied
# per projected element. N is capped at 60 so every predicate stays in Int64
# arithmetic (StochasticB shifts by N+1; see project.jl).
@inline function _check_nrandbits(N)
    (N isa Int && 1 <= N <= 60) ||
        throw(ArgumentError("stochastic rounding requires an Int random-bit budget N with 1 ≤ N ≤ 60 (got $N)"))
    nothing
end
"""Stochastic rounding, variant A of draft §4.7.4: RoundAway ⟺ ⌊ν·2^N⌋ + R ≥ 2^N."""
struct StochasticA{N} <: RoundingMode3109
    StochasticA{N}() where {N} = (_check_nrandbits(N); new{N}())
end
"""Stochastic rounding, variant B: RoundAway ⟺ ⌊ν·2^(N+1)⌋ + (2R+1) ≥ 2^(N+1)."""
struct StochasticB{N} <: RoundingMode3109
    StochasticB{N}() where {N} = (_check_nrandbits(N); new{N}())
end
"""Stochastic rounding, variant C: RoundAway ⟺ RNITE(ν·2^N) + R ≥ 2^N."""
struct StochasticC{N} <: RoundingMode3109
    StochasticC{N}() where {N} = (_check_nrandbits(N); new{N}())
end

const RoundingModes_Ties = Union{NearestTiesToEven,NearestTiesToAway}
const RoundingModes_Directed = Union{TowardPositive,TowardNegative,TowardZero,ToOdd}
const RoundingModes_Stochastic = Union{StochasticA,StochasticB,StochasticC}

const RoundingModes_Deterministic = Union{RoundingModes_Ties,RoundingModes_Directed}
const RoundingModes = Union{RoundingModes_Deterministic,RoundingModes_Stochastic}

abstract type SaturationMode end

struct SatFinite    <: SaturationMode end   # clamp everything to the finite range
struct SatPropagate <: SaturationMode end   # keep representable infinities, clamp the rest
struct SatNone      <: SaturationMode end   # draft's direction/signedness/domain-governed rows

const SaturationModes_NonSaturating = Union{SatNone}
const SaturationModes_Saturating = Union{SatFinite,SatPropagate}

const SaturationModes = Union{SatFinite,SatPropagate,SatNone}

"""
    ProjSpec{R<:RoundingMode3109, S<:SaturationMode}

A projection specification ρ = (rounding mode, saturation mode), draft §4.2.
Zero-size; construct as `ProjSpec(NearestTiesToEven(), SatNone())` or via the
exported constants. Kernels specialize on its type.
"""
struct ProjSpec{R<:RoundingMode3109,S<:SaturationMode} end

ProjSpec(::R, ::S) where {R<:RoundingMode3109,S<:SaturationMode} = ProjSpec{R,S}()

roundingmode(::ProjSpec{R,S}) where {R,S} = R()

saturationmode(::ProjSpec{R,S}) where {R,S} = S()

"""Draft-named accessor: RoundOf(ρ) — the rounding mode of ρ (§4.2)."""
const RoundOf = roundingmode
"""Draft-named accessor: SatOf(ρ) — the saturation mode of ρ (§4.2)."""
const SatOf = saturationmode

# ---- queries
isstochastic(::Type{<:RoundingMode3109}) = false
isstochastic(::Type{<:RoundingModes_Stochastic}) = true
isstochastic(m::RoundingMode3109) = isstochastic(typeof(m))
isstochastic(::ProjSpec{R,S}) where {R,S} = isstochastic(R)

"""Number of random bits N consumed per projected element; 0 for pure modes."""
nrandbits(::Type{<:RoundingMode3109}) = 0
nrandbits(::Type{StochasticA{N}}) where {N} = N
nrandbits(::Type{StochasticB{N}}) where {N} = N
nrandbits(::Type{StochasticC{N}}) where {N} = N
nrandbits(m::RoundingMode3109) = nrandbits(typeof(m))
nrandbits(::ProjSpec{R,S}) where {R,S} = nrandbits(R)

# ---- predefined projections
"""(NearestTiesToEven, SatFinite)."""
const RNE_SatFinite    = ProjSpec{NearestTiesToEven, SatFinite}()
"""(NearestTiesToEven, SatPropagate)."""
const RNE_SatPropagate = ProjSpec{NearestTiesToEven, SatPropagate}()
"""(NearestTiesToEven, SatNone) — the package-wide default ρ."""
const RNE_SatNone      = ProjSpec{NearestTiesToEven, SatNone}()

"""(NearestTiesToAway, SatFinite)."""
const RNA_SatFinite    = ProjSpec{NearestTiesToAway, SatFinite}()
"""(NearestTiesToAway, SatPropagate)."""
const RNA_SatPropagate = ProjSpec{NearestTiesToAway, SatPropagate}()
"""(NearestTiesToAway, SatNone)."""
const RNA_SatNone      = ProjSpec{NearestTiesToAway, SatNone}()

"""(TowardPositive, SatFinite)."""
const RTP_SatFinite    = ProjSpec{TowardPositive, SatFinite}()
"""(TowardPositive, SatPropagate)."""
const RTP_SatPropagate = ProjSpec{TowardPositive, SatPropagate}()
"""(TowardPositive, SatNone)."""
const RTP_SatNone      = ProjSpec{TowardPositive, SatNone}()

"""(TowardNegative, SatFinite)."""
const RTN_SatFinite    = ProjSpec{TowardNegative, SatFinite}()
"""(TowardNegative, SatPropagate)."""
const RTN_SatPropagate = ProjSpec{TowardNegative, SatPropagate}()
"""(TowardNegative, SatNone)."""
const RTN_SatNone      = ProjSpec{TowardNegative, SatNone}()

"""(TowardZero, SatFinite)."""
const RTZ_SatFinite    = ProjSpec{TowardZero, SatFinite}()
"""(TowardZero, SatPropagate)."""
const RTZ_SatPropagate = ProjSpec{TowardZero, SatPropagate}()
"""(TowardZero, SatNone)."""
const RTZ_SatNone      = ProjSpec{TowardZero, SatNone}()

"""(ToOdd, SatFinite)."""
const RTO_SatFinite    = ProjSpec{ToOdd, SatFinite}()
"""(ToOdd, SatPropagate)."""
const RTO_SatPropagate = ProjSpec{ToOdd, SatPropagate}()
"""(ToOdd, SatNone)."""
const RTO_SatNone      = ProjSpec{ToOdd, SatNone}()

"""Default random-bit budget used by no-argument stochastic projection constructors."""
const DEFAULT_RBITS = 8

const RSA_SatFiniteType    = ProjSpec{StochasticA{N}, SatFinite} where {N}
const RSA_SatPropagateType = ProjSpec{StochasticA{N}, SatPropagate} where {N}
const RSA_SatNoneType      = ProjSpec{StochasticA{N}, SatNone} where {N}

RSA_SatFinite() = RSA_SatFinite(Val(DEFAULT_RBITS))
RSA_SatFinite(::Val{N}) where {N} = (_check_nrandbits(N); RSA_SatFiniteType{N}())
RSA_SatFinite(N::Int) = RSA_SatFinite(Val(N))

RSA_SatPropagate() = RSA_SatPropagate(Val(DEFAULT_RBITS))
RSA_SatPropagate(::Val{N}) where {N} = (_check_nrandbits(N); RSA_SatPropagateType{N}())
RSA_SatPropagate(N::Int) = RSA_SatPropagate(Val(N))

RSA_SatNone() = RSA_SatNone(Val(DEFAULT_RBITS))
RSA_SatNone(::Val{N}) where {N} = (_check_nrandbits(N); RSA_SatNoneType{N}())
RSA_SatNone(N::Int) = RSA_SatNone(Val(N))

const RSB_SatFiniteType    = ProjSpec{StochasticB{N}, SatFinite} where {N}
const RSB_SatPropagateType = ProjSpec{StochasticB{N}, SatPropagate} where {N}
const RSB_SatNoneType      = ProjSpec{StochasticB{N}, SatNone} where {N}

RSB_SatFinite() = RSB_SatFinite(Val(DEFAULT_RBITS))
RSB_SatFinite(::Val{N}) where {N} = (_check_nrandbits(N); RSB_SatFiniteType{N}())
RSB_SatFinite(N::Int) = RSB_SatFinite(Val(N))

RSB_SatPropagate() = RSB_SatPropagate(Val(DEFAULT_RBITS))
RSB_SatPropagate(::Val{N}) where {N} = (_check_nrandbits(N); RSB_SatPropagateType{N}())
RSB_SatPropagate(N::Int) = RSB_SatPropagate(Val(N))

RSB_SatNone() = RSB_SatNone(Val(DEFAULT_RBITS))
RSB_SatNone(::Val{N}) where {N} = (_check_nrandbits(N); RSB_SatNoneType{N}())
RSB_SatNone(N::Int) = RSB_SatNone(Val(N))

const RSC_SatFiniteType    = ProjSpec{StochasticC{N}, SatFinite} where {N}
const RSC_SatPropagateType = ProjSpec{StochasticC{N}, SatPropagate} where {N}
const RSC_SatNoneType      = ProjSpec{StochasticC{N}, SatNone} where {N}

RSC_SatFinite() = RSC_SatFinite(Val(DEFAULT_RBITS))
RSC_SatFinite(::Val{N}) where {N} = (_check_nrandbits(N); RSC_SatFiniteType{N}())
RSC_SatFinite(N::Int) = RSC_SatFinite(Val(N))

RSC_SatPropagate() = RSC_SatPropagate(Val(DEFAULT_RBITS))
RSC_SatPropagate(::Val{N}) where {N} = (_check_nrandbits(N); RSC_SatPropagateType{N}())
RSC_SatPropagate(N::Int) = RSC_SatPropagate(Val(N))

RSC_SatNone() = RSC_SatNone(Val(DEFAULT_RBITS))
RSC_SatNone(::Val{N}) where {N} = (_check_nrandbits(N); RSC_SatNoneType{N}())
RSC_SatNone(N::Int) = RSC_SatNone(Val(N))

default_projspec(::Type{<:Binary}) = RNE_SatNone

# ---- Base.RoundingMode compatibility shim (API boundaries only)
"""
Translate a `Base.RoundingMode` to its `RoundingMode3109` counterpart (identity on the latter).
Used only at API boundaries.
Internals always carry RoundingMode3109.
"""
projmode(::RoundingMode{:Nearest})         = NearestTiesToEven()
projmode(::RoundingMode{:NearestTiesAway}) = NearestTiesToAway()
projmode(::RoundingMode{:Up})              = TowardPositive()
projmode(::RoundingMode{:Down})            = TowardNegative()
projmode(::RoundingMode{:ToZero})          = TowardZero()
projmode(m::RoundingMode3109)              = m

# ---- draft-style printing: "(NearestTiesToEven, SatNone)"
_modename(m) = String(nameof(typeof(m)))
_modename(::StochasticA{N}) where {N} = "StochasticA[$N]"
_modename(::StochasticB{N}) where {N} = "StochasticB[$N]"
_modename(::StochasticC{N}) where {N} = "StochasticC[$N]"
Base.show(io::IO, ρ::ProjSpec) =
    print(io, "(", _modename(roundingmode(ρ)), ", ", _modename(saturationmode(ρ)), ")")
