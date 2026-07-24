# ===== ByteFloats.jl — module root (design §10.1, architecture §11)
"""
    ByteFloats

A conforming, performance-oriented Julia implementation of the IEEE P3109 draft
standard for arithmetic formats for machine learning (bitwidths 3–8).

Bit-exact defined results on every default path; the projection engine
(`RoundToPrecision → Saturate → Encode`) is the single write path into a code
point; approximate fast paths exist only behind the explicit κ registry
(`register_approx!` / `approx`), never substituted silently. See
`conformance()` for the live declaration and `ByteFloats.draft_revision()` for
the implemented draft.

# Performance note
Pass format types through `const` bindings, type parameters, or function
arguments. All exported entry points fully specialize when the format type is
statically known (measured: a complete scalar `Add` ≈ 26 ns, `project` ≈ 13 ns,
zero allocations). Calling through a **non-`const` global** format type forces
dynamic dispatch on every call — measured at ~1 µs per scalar keyword call —
which is Julia semantics, not package overhead; a single function barrier
(`f(::Type{T}, args...) where {T}`) restores full speed. The test suite pins
the specialization properties (concrete return types, zero warm-path
allocation) as deterministic regressions.
"""
module ByteFloats

const binary16 = Float16
const binary32 = Float32
const binary64 = Float64

using Random: Random, AbstractRNG, default_rng
using PrecompileTools: @setup_workload, @compile_workload
using Quadmath: Float128

include("fma128.jl")
using .Float128FMA          # on Windows this installs Base.fma

include("faa128.jl")
using .Float128FAA          

# Include order: formats → projspec → defaults → decode_encode → project →
# ops_scalar → oracle → tables → kernels → blocks → packed → approx.
# (One deliberate delta from the architecture §11 listing: the evaluation-protocol
# structs BigExactF/EncloseF live in ops_scalar.jl per §6, so ops_scalar precedes
# oracle; oracle's references to them are function-body-late-bound either way, but
# this is the order the harnesses verified.)
include("formats.jl")
include("defaults.jl")
include("projspec.jl")
include("decode_encode.jl")
include("project.jl")
include("ops_scalar.jl")
include("oracle.jl")
include("tables.jl")
include("kernels.jl")
include("blocks.jl")
include("packed.jl")
include("approx.jl")

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------
# the type and every draft §3.2 named format
export Binary
for n in sort!(collect(keys(_NAMED)))
    @eval export $n
end

export binary64, binary32, binary16

# Group M and format introspection
export bitwidth, issigned, isextended, expbias, expbitwidth, trailingsigbits,
       BitwidthOf, PrecisionOf, SignednessOf, DomainOf, ExponentBiasOf,
       ExponentBitwidthOf, TrailingSignificandBitwidthOf,
       MaxFiniteOf, MinFiniteOf, MinPositiveOf, MaxSubnormalOf, MinNormalOf,
       maxfinite_datum, minfinite_datum, formatname, rawvalue, decode   # codepoint extends Base

# projection specifications
export RoundingMode3109, NearestTiesToEven, NearestTiesToAway, TowardPositive,
       TowardNegative, TowardZero, ToOdd, StochasticA, StochasticB, StochasticC,
       SaturationMode, SatFinite, SatPropagate, SatNone,
       ProjSpec, RoundOf, SatOf, roundingmode, saturationmode,
       isstochastic, nrandbits,
       RNE_SatFinite, RNE_SatPropagate, RNE_SatNone,
       RNA_SatFinite, RNA_SatPropagate, RNA_SatNone,
       RTP_SatFinite, RTP_SatPropagate, RTP_SatNone,
       RTN_SatFinite, RTN_SatPropagate, RTN_SatNone,
       RTZ_SatFinite, RTZ_SatPropagate, RTZ_SatNone,
       RTO_SatFinite, RTO_SatPropagate, RTO_SatNone,
       RSA_SatFinite, RSA_SatPropagate, RSA_SatNone,
       RSB_SatFinite, RSB_SatPropagate, RSB_SatNone,
       RSC_SatFinite, RSC_SatPropagate, RSC_SatNone,
       default_projspec, projmode

# session defaults (defaults.jl)
export DefaultType, DefaultType!,
       DefaultReturnType, DefaultReturnType!,
       DefaultAccumulatorType, DefaultAccumulatorType!,
       DefaultRoundingMode, DefaultRoundingMode!,
       DefaultSaturationMode, DefaultSaturationMode!,
       DefaultProjection, DefaultProjection!,
       DefaultRNG, DefaultRNG!,
       DefaultRbits, DefaultRbits!,
       with_default_type, with_default_returntype,
       with_default_accumulatortype, with_default_projection

# comparison, classification, stepping (Groups D/M)
export TotalOrder, Class, FPClass,
       ClassNaN, ClassNegInf, ClassNegNormal, ClassNegSubnormal, ClassZero,
       ClassPosSubnormal, ClassPosNormal, ClassPosInf,
       NextGreaterThan, NextLessThan

# the spec register: every draft operation (scalar + array methods)
for op in OP_REGISTRY
    @eval export $(op.name)
end
export vmap, vmap!

# table cache introspection
export table_bytes, empty_tables!

# blocks and scaled operations
export Block, BlockVector, blocksize, scaleformat, elemformat, PackedVector
for op in OP_REGISTRY
    op.name === :Convert && continue
    @eval export $(Symbol(:Block, op.name))
    @eval export $(Symbol(:Scaled, op.name))
end
export BlockReduceAdd, BlockReduceMultiply, BlockDotProduct,
       ConvertFromBlock, ConvertToBlock, ConvertToBlockMaxAbsFinite

# conformance and κ-approximation
export conformance, conformance_dict, conformance_report, draft_revision,
       ConformanceDeclaration, ApproxImpl,
       measure_kappa, codedistance, register_approx!, unregister_approx!,
       approx, kappa, kappa_measured, list_approx, ftz_variant

# ---------------------------------------------------------------------------
# Tier-1 precompile workload (design §7.2): the standard profile's hot entries
# compile during precompilation; everything else specializes lazily on first use.
# ---------------------------------------------------------------------------
@setup_workload begin
    @compile_workload begin
        T = Binary8p4se; S = Binary8p3se
        a, b = T(1.5), T(0.25)
        Add(T, RNE_SatNone, a, b); Multiply(T, RNE_SatFinite, a, b)
        Exp(T, RNE_SatNone, a); Convert(S, RNE_SatNone, a)
        a + b; exp(b); fma(a, b, a); min(a, b)
        get_table(:Exp, T, T, RNE_SatNone)
        A = [a, b, a, b]; B = [b, a, b, a]; d = similar(A)
        vmap!(d, Val(:Add), T, RNE_SatNone, A, B)
        vmap!(d, Val(:Exp), T, RNE_SatNone, A)
        ScaledAdd(T, RNE_SatNone, one(S), a, one(S), b)
        bx = Block(one(S), (a, b, a, b)); by = Block(one(S), (b, a, b, a))
        BlockDotProduct(T, RNE_SatNone, bx, by)
        BlockAdd(T, RNE_SatNone, bx, by, one(S))
        empty_tables!()          # tables are cheap to rebuild; don't bloat the image
    end
end

function __init__()
    # Float128 revision plan §5: every Float128 path fronts a complete MPFR path
    # with identical semantics; the switch trades speed only, never results.
    if get(ENV, "ByteFloats_Float128", "") == "disable"
        _USE_FLOAT128[] = false
    end
    nothing
end

end # module
