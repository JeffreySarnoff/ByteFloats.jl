# ===== projspec.jl — projection specifications as types (design §4, architecture §2)
#
# A projection specification ρ = (rounding mode, saturation mode) is a zero-size
# value whose *type* carries both modes, so kernels specialize on ρ exactly as
# they specialize on formats. The stochastic modes carry the random-bit budget N
# in the type: a kernel's randomness consumption is a compile-time fact, and
# pure-vs-stochastic dispatch is static (pure ⇒ tabulable, stochastic ⇒ never).

abstract type RoundingMode3109 end

struct NearestTiesToEven <: RoundingMode3109 end
struct NearestTiesToAway <: RoundingMode3109 end
struct TowardPositive    <: RoundingMode3109 end
struct TowardNegative    <: RoundingMode3109 end
struct TowardZero        <: RoundingMode3109 end
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
const AnyStochastic = Union{StochasticA,StochasticB,StochasticC}

abstract type SaturationMode end
struct SatFinite    <: SaturationMode end   # clamp everything to the finite range
struct SatPropagate <: SaturationMode end   # keep representable infinities, clamp the rest
struct SatNone      <: SaturationMode end   # draft's direction/signedness/domain-governed rows

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
isstochastic(::Type{<:AnyStochastic}) = true
isstochastic(m::RoundingMode3109) = isstochastic(typeof(m))
isstochastic(::ProjSpec{R,S}) where {R,S} = isstochastic(R)

"""Number of random bits N consumed per projected element; 0 for pure modes."""
nrandbits(::Type{<:RoundingMode3109}) = 0
nrandbits(::Type{StochasticA{N}}) where {N} = N
nrandbits(::Type{StochasticB{N}}) where {N} = N
nrandbits(::Type{StochasticC{N}}) where {N} = N
nrandbits(m::RoundingMode3109) = nrandbits(typeof(m))
nrandbits(::ProjSpec{R,S}) where {R,S} = nrandbits(R)

# ---- defaults (design §10.2)
"""(NearestTiesToEven, SatNone) — the package-wide default ρ."""
const RNE_SatNone   = ProjSpec{NearestTiesToEven,SatNone}()
const RNE_SatFinite = ProjSpec{NearestTiesToEven,SatFinite}()
default_projspec(::Type{<:Binary}) = RNE_SatNone

# ---- Base.RoundingMode compatibility shim (API boundaries only)
projmode(::RoundingMode{:Nearest})         = NearestTiesToEven()
projmode(::RoundingMode{:NearestTiesAway}) = NearestTiesToAway()
projmode(::RoundingMode{:Up})              = TowardPositive()
projmode(::RoundingMode{:Down})            = TowardNegative()
projmode(::RoundingMode{:ToZero})          = TowardZero()
projmode(m::RoundingMode3109)             = m

# ---- draft-style printing: "(NearestTiesToEven, SatNone)"
_modename(m) = String(nameof(typeof(m)))
_modename(::StochasticA{N}) where {N} = "StochasticA[$N]"
_modename(::StochasticB{N}) where {N} = "StochasticB[$N]"
_modename(::StochasticC{N}) where {N} = "StochasticC[$N]"
Base.show(io::IO, ρ::ProjSpec) =
    print(io, "(", _modename(roundingmode(ρ)), ", ", _modename(saturationmode(ρ)), ")")
