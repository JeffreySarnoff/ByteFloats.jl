# Cheat Sheet

A one-page reference for the ByteFloats.jl operations you are most likely to
write. For explanations and semantics, see the [User Guide](@ref); for
implementation details, see the [Technical Guide](@ref).

## Start here

```julia
using ByteFloats
using Random: Xoshiro

T = Binary8p4se
ρ = RNE_SatNone

x = T(1.6)
y = T(0.25)
z = Add(T, ρ, x, y)
```

The explicit operation shape is always:

```julia
Op(result_format, projection_spec, operands...; rng, R)
```

For same-format operands under the default projection, ordinary Julia syntax is
available:

```julia
z == x + y
exp(x)
fma(x, y, z)
```

## Format names

```text
Binary K p P (s|u) (e|f)
       │   │  │     └─ extended (Inf) or finite domain
       │   │  └────── signed or unsigned
       │   └───────── significand precision, including the implicit bit
       └───────────── total bitwidth, 3 through 8
```

Examples:

| Name | Meaning |
|---|---|
| `Binary8p4se` | 8-bit, precision 4, signed, extended |
| `Binary8p4sf` | 8-bit, precision 4, signed, finite |
| `Binary6p3ue` | 6-bit, precision 3, unsigned, extended |
| `Binary5p5uf` | 5-bit, precision 5, unsigned, finite |

The parametric spelling is equivalent:

```julia
Binary8p4se === Binary{8, 4, true, true}
```

Useful format queries:

```julia
bitwidth(T)       # K
precision(T)      # P
issigned(T)
isextended(T)
expbias(T)
expbitwidth(T)
trailingsigbits(T)

MaxFiniteOf(T)
MinFiniteOf(T)
MinPositiveOf(T)
MinNormalOf(T)
MaxSubnormalOf(T)
```

Every query above also takes a **value** instead of a type — `bitwidth(x)` ≡
`bitwidth(typeof(x))` — and folds to the same constant:

```julia
x = Binary8p4se(1.6)
bitwidth(x)       # 8
MaxFiniteOf(x)    # Binary8p4se(224.0 ≡ 0x7e)
```

## Values and code points

```julia
x = T(1.6)                    # numeric value: project into T
x = Convert(T, ρ, 1.6)       # explicit projection

c = codepoint(x)              # UInt8 code point
x = T(c)                      # UInt8 means validated CODE POINT
x = rawvalue(T, c)            # unchecked code-point construction

d = decode(x)                 # exact Float64 datum
```

!!! warning "`UInt8` means code point"
    `T(0x02)` constructs code point `0x02`; `T(2)` projects the numeric value
    two. Use `Convert` when intent should be unmistakable.

Classification and stepping:

```julia
isnan(x); isinf(x); isfinite(x); iszero(x)
signbit(x); issubnormal(x)
Class(x)
NextGreaterThan(x); NextLessThan(x)
nextfloat(x); prevfloat(x)
TotalOrder(x, y)
```

## Projection specifications

A projection specification is `(rounding mode, saturation mode)`:

```julia
ρ = ProjSpec(TowardPositive(), SatFinite())
roundingmode(ρ)
saturationmode(ρ)
```

### Deterministic projections

| Rounding | `SatFinite` | `SatPropagate` | `SatNone` |
|---|---|---|---|
| nearest, ties even | `RNE_SatFinite` | `RNE_SatPropagate` | `RNE_SatNone` |
| nearest, ties away | `RNA_SatFinite` | `RNA_SatPropagate` | `RNA_SatNone` |
| toward +∞ | `RTP_SatFinite` | `RTP_SatPropagate` | `RTP_SatNone` |
| toward −∞ | `RTN_SatFinite` | `RTN_SatPropagate` | `RTN_SatNone` |
| toward zero | `RTZ_SatFinite` | `RTZ_SatPropagate` | `RTZ_SatNone` |
| round to odd | `RTO_SatFinite` | `RTO_SatPropagate` | `RTO_SatNone` |

`RNE_SatNone` is the package-wide default. Saturation modes mean:

| Mode | Out-of-range behavior |
|---|---|
| `SatFinite` | clamp everything to the finite range |
| `SatPropagate` | preserve representable infinities; clamp other overflow |
| `SatNone` | apply the draft's domain-, signedness-, and direction-dependent rows |

### Stochastic projections

Constructors are grouped by stochastic variant:

| Variant | `SatFinite` | `SatPropagate` | `SatNone` |
|---|---|---|---|
| A | `RSA_SatFinite(N)` | `RSA_SatPropagate(N)` | `RSA_SatNone(N)` |
| B | `RSB_SatFinite(N)` | `RSB_SatPropagate(N)` | `RSB_SatNone(N)` |
| C | `RSC_SatFinite(N)` | `RSC_SatPropagate(N)` | `RSC_SatNone(N)` |

`N` is the random-bit budget, `1 ≤ N ≤ 60`. Omitting it uses `N = 8`.

```julia
σ = RSA_SatNone(8)

Add(T, σ, x, y; rng=Xoshiro(1))  # reproducible stream
Add(T, σ, x, y; R=17)            # exact draw, ideal for tests

isstochastic(σ)                   # true
nrandbits(σ)                      # 8
```

For an `N`-bit mode, explicit `R` must be in `0:(2^N - 1)`.

### Session defaults

Read with `DefaultX()`, set with `DefaultX!(v)`:

```julia
DefaultType()            # Binary8p2se     DefaultType!(Binary8p4se)
DefaultReturnType()      # Binary8p2se     DefaultReturnType!(Binary8p3se)
DefaultAccumulatorType() # binary32        DefaultAccumulatorType!(binary64)
DefaultRoundingMode()    # NearestTiesToEven()
DefaultSaturationMode()  # SatNone()
DefaultProjection()      # RNE_SatNone
DefaultRNG()             # Xoshiro         DefaultRNG!(Xoshiro(42))
DefaultRbits()           # 8               DefaultRbits!(16)
```

Setting a rounding/saturation component rebuilds `DefaultProjection`; setting
`DefaultProjection!` directly updates both components. Always:
`DefaultProjection() === ProjSpec(DefaultRoundingMode(), DefaultSaturationMode())`.

The convenience methods (`a + b`, `Exp(x)`, `T(2.1)`) **do** consult
`DefaultProjection` — `default_projspec` reads it. Changing it changes their
results globally, and costs specialization (ρ is no longer a compile-time
constant, so those calls allocate). Name ρ explicitly — `Add(T, ρ, x, y)` — in
hot code; explicit forms are unaffected by the session default.

Consume a default via the combinators — never by computing on a bare
`DefaultX()` read. No dispatch while the default is unchanged, one barrier
dispatch after a change; zero-alloc when `f`'s result type doesn't depend on
the default (a default-typed result boxes once at escape):

```julia
with_default_type((T, x) -> T(x), 1.5)          # Binary8p2se(1.5 ≡ 0x41)
with_default_projection((ρ, x, y) -> Add(T, ρ, x, y), x, y)
with_default_returntype(f, args...)              # f(DefaultReturnType(), args...)
with_default_accumulatortype(f, args...)
```

## Scalar operation catalog

```julia
# explicit result format and projection
Add(T, ρ, x, y)
FMA(T, ρ, x, y, z)
Convert(T, ρ, external_value)

# same-format default-projection convenience
Add(x, y)
FMA(x, y, z)
```

| Arity | Operations |
|---|---|
| Unary | `Abs`, `Negate`, `Sqrt`, `RSqrt`, `Recip`, `Exp`, `Exp2`, `ExpMinusOne`, `Log`, `Log2`, `LogOnePlus`, `Softplus`, `Sin`, `Cos`, `Tan`, `ArcSin`, `ArcCos`, `ArcTan`, `Sinh`, `Cosh`, `Tanh`, `ArcSinh`, `ArcCosh`, `ArcTanh`, `SinPi`, `CosPi`, `TanPi`, `ArcSinPi`, `ArcCosPi`, `ArcTanPi` |
| Binary | `CopySign`, `Add`, `Subtract`, `Multiply`, `Divide`, `Hypot`, `ArcTan2`, `ArcTan2Pi`, `Maximum`, `Minimum`, `MaximumNumber`, `MinimumNumber`, `MaximumMagnitude`, `MinimumMagnitude`, `MaximumMagnitudeNumber`, `MinimumMagnitudeNumber`, `MinimumFinite`, `MaximumFinite` |
| Ternary | `FMA`, `FAA`, `Clamp` |
| Conversion | `Convert` |

Common Base spellings under `RNE_SatNone`:

```julia
x + y; x - y; x * y; x / y
-x; abs(x); inv(x); sqrt(x)
exp(x); exp2(x); expm1(x)
log(x); log2(x); log1p(x)
sin(x); cos(x); tan(x)
asin(x); acos(x); atan(x); atan(y, x)
sinh(x); cosh(x); tanh(x)
min(x, y); max(x, y); clamp(x, y, z)
fma(x, y, z); muladd(x, y, z)
```

## Arrays

Every registered operation has elementwise array methods:

```julia
A = T.([1.0, 1.5, 2.0])
B = T.([0.25, 0.5, 0.75])

C = Add(T, ρ, A, B)
E = Exp(T, ρ, A)
F = FMA(T, ρ, A, B, C)

C = vmap(:Add, T, ρ, A, B)
vmap!(C, Val(:Add), T, ρ, A, B)
```

Operands and destination must have matching axes. Deterministic unary/binary
operations use cached result tables; affordable ternary signatures may also use
tables. Stochastic operations compute each element and consume one draw per
projection.

```julia
table_bytes()
empty_tables!()
```

## Blocks and scaled operations

```julia
FS = Binary8p1uf
FE = Binary8p4se

b = Block(FS(4.0), FE(1.5), FE(-0.75), FE(2.0), FE(0.5))

blocksize(b)
scaleformat(b)
elemformat(b)

ConvertFromBlock(T, ρ, b)
BlockReduceAdd(T, ρ, b)
BlockReduceMultiply(T, ρ, b)
BlockDotProduct(T, ρ, b, b)
```

Every scalar operation except `Convert` has generated block and scaled forms:

```julia
BlockAdd(T, ρ, b1, b2, result_scale)
BlockExp(T, ρ, b, result_scale)

ScaledAdd(T, ρ, scale1, x1, scale2, x2)
ScaledExp(T, ρ, scale, x)
```

Quantize against a supplied or automatically selected scale:

```julia
ConvertToBlock(FS, FE, ρ, values_tuple, scale)
ConvertToBlockMaxAbsFinite(FS, FE, scale_ρ, element_ρ, values_tuple)
```

`BlockVector` stores many equal-shape blocks in a structure-of-arrays layout.

## Packed storage

```julia
v = Binary5p2se.([0.5, 1.0, 1.5, 2.0])
pv = PackedVector(v)

pv[2]
pv[2] = Binary5p2se(0.75)
collect(pv)

out = vmap(:Exp, Binary5p2se, ρ, pv)
```

`PackedVector` stores each code point in exactly `bitwidth(F)` bits. Computation
unpacks tiles internally; packed arithmetic is deliberately not in-place.

## Conformance and approximations

```julia
conformance()
conformance_dict()
conformance_report()
draft_revision()

κ, exhaustive = measure_kappa(fn, :Exp, T, (T,), ρ)
register_approx!(:fast_exp, :Exp, T, (T,), ρ, fn; κ)
impl = approx(:fast_exp)
kappa(impl)
kappa_measured(impl)
list_approx()
unregister_approx!(:fast_exp)
```

Approximate implementations are never substituted into the default API.

## Common traps

| Trap | Correct pattern |
|---|---|
| Treating `UInt8` as a number | `T(Int(c))` for a numeric integer; `T(c::UInt8)` for a code point |
| Silently mixing byte-float formats | Convert explicitly: `Convert(T, ρ, x)` |
| Assuming `SatNone` always clamps | Choose `SatFinite` when clamping is required |
| Expecting IEEE division-by-zero | P3109 semantics here define `x / 0 → NaN` |
| Expecting negative zero | Every format has one zero; `T(-0.0) === zero(T)` |
| Passing a `Rational` | Convert explicitly to an exact supported carrier or knowingly to `Float64` |
| Reproducibility with stochastic rounding | Supply a seeded `rng`, or an explicit `R` in tests |
| Using `rawvalue` on unchecked input | Prefer validated `T(c::UInt8)` outside kernels |

## Performance checklist

- Keep format types and projection specs in `const` bindings, function arguments,
  or type parameters.
- Put code using runtime-selected formats behind a function barrier.
- Expect the first deterministic array call for a specialization to build a table;
  benchmark warm calls separately.
- Use `PackedVector` when storage bandwidth matters more than direct byte access.
- Use `table_bytes()` to inspect cache footprint and `empty_tables!()` to reset it.
- Do not replace explicit projection semantics with `@fastmath`.

For worked applications, continue to [User Examples](@ref). For correctness and
benchmark methodology, continue to [Technical Examples](@ref).
