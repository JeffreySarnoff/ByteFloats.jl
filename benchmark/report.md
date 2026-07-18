# ByteFloats.jl benchmark report

Generated: 2026-07-18 04:25 UTC  ·  Julia 1.12.6  ·  sapphirerapids (1 threads)  ·  Float128 paths: enabled  ·  Chairmarks 1.3.1

Reference format for per-operation tables: `Binary8p4se` under `(NearestTiesToEven, SatNone)`; operand pools draw uniformly over all code points (an honest NaN/±Inf mix). Times are per call; medians with minima alongside. Methodology per the recorded benchmark doctrine: type-parameterized barriers, untimed setup, specialization preflight.

## Core primitives

The decode/compare/step/classify layer plus the projection engine.

| operation | median | min | allocs |
|---|---|---|---|
| `decode` | 2.7 ns | 2.7 ns | 0 |
| `order_key` | 3.0 ns | 2.7 ns | 0 |
| `x < y` | 2.7 ns | 2.4 ns | 0 |
| `TotalOrder` | 3.3 ns | 3.0 ns | 0 |
| `Class` | 3.6 ns | 2.4 ns | 0 |
| `NextGreaterThan` | 3.3 ns | 2.7 ns | 0 |
| `project (RNE·SatNone)` | 11.4 ns | 2.8 ns | 0 |
| `project (StochasticA[8], R drawn)` | 12.1 ns | 3.4 ns | 0 |

## Scalar operations — unary (30)

Sorted by median. Random operands: transcendental rows mix special-row fast returns with enclosure-path evaluations, so these are *scalar-path* costs; bulk unary work routes through 256-byte tables (see Array kernels).

| operation | median | min | allocs |
|---|---|---|---|
| `ArcTanh` | 2.8 ns | 2.5 ns | 0 |
| `ArcCosh` | 4.2 ns | 3.4 ns | 0 |
| `Log2` | 9.7 ns | 3.3 ns | 0 |
| `Abs` | 11.5 ns | 3.4 ns | 0 |
| `Negate` | 11.6 ns | 3.4 ns | 0 |
| `Sqrt` | 18.4 ns | 6.3 ns | 0 |
| `RSqrt` | 20.4 ns | 4.6 ns | 0 |
| `ArcSinPi` | 22.0 ns | 3.4 ns | 0 |
| `ArcSin` | 307.7 ns | 3.5 ns | 0 |
| `ArcCosPi` | 518.8 ns | 2.7 ns | 0 |
| `CosPi` | 581.5 ns | 4.6 ns | 0 |
| `ArcCos` | 623.9 ns | 3.3 ns | 0 |
| `SinPi` | 692.6 ns | 4.1 ns | 0 |
| `ArcTanPi` | 712.6 ns | 5.0 ns | 0 |
| `TanPi` | 721.0 ns | 5.2 ns | 0 |
| `ArcTan` | 741.1 ns | 4.9 ns | 0 |
| `Cos` | 748.1 ns | 4.3 ns | 0 |
| `Log` | 759.7 ns | 759.7 ns | 0 |
| `Sin` | 831.6 ns | 4.7 ns | 0 |
| `Tan` | 910.4 ns | 4.8 ns | 0 |
| `Recip` | 995.7 ns | 7.6 ns | 0 |
| `ExpMinusOne` | 1.05 μs | 5.6 ns | 0 |
| `Tanh` | 1.11 μs | 5.0 ns | 0 |
| `Exp` | 1.11 μs | 5.9 ns | 0 |
| `Cosh` | 1.15 μs | 6.4 ns | 0 |
| `Sinh` | 1.15 μs | 5.6 ns | 0 |
| `LogOnePlus` | 1.3 μs | 1.3 μs | 0 |
| `ArcSinh` | 1.44 μs | 6.7 ns | 0 |
| `Exp2` | 1.65 μs | 6.2 ns | 0 |
| `Softplus` | 2.13 μs | 6.8 ns | 0 |

## Scalar operations — binary (18)

Sorted by median.

| operation | median | min | allocs |
|---|---|---|---|
| `CopySign` | 12.2 ns | 3.1 ns | 0 |
| `MinimumMagnitude` | 12.3 ns | 3.1 ns | 0 |
| `MaximumMagnitude` | 12.4 ns | 2.8 ns | 0 |
| `MinimumMagnitudeNumber` | 12.6 ns | 9.1 ns | 0 |
| `MaximumMagnitudeNumber` | 12.7 ns | 7.9 ns | 0 |
| `Multiply` | 12.9 ns | 3.1 ns | 0 |
| `Maximum` | 13.9 ns | 2.8 ns | 0 |
| `Minimum` | 14.0 ns | 2.8 ns | 0 |
| `MaximumFinite` | 14.0 ns | 8.6 ns | 0 |
| `MinimumFinite` | 14.0 ns | 8.6 ns | 0 |
| `MaximumNumber` | 14.5 ns | 8.0 ns | 0 |
| `MinimumNumber` | 14.6 ns | 8.1 ns | 0 |
| `Subtract` | 15.8 ns | 6.7 ns | 0 |
| `Add` | 17.7 ns | 4.0 ns | 0 |
| `Hypot` | 522.3 ns | 5.6 ns | 0 |
| `ArcTan2` | 819.8 ns | 6.3 ns | 0 |
| `ArcTan2Pi` | 840.7 ns | 5.2 ns | 0 |
| `Divide` | 1.02 μs | 12.7 ns | 0 |

## Scalar operations — ternary (3)

Sorted by median.

| operation | median | min | allocs |
|---|---|---|---|
| `FMA` | 16.2 ns | 6.1 ns | 0 |
| `FAA` | 16.6 ns | 6.1 ns | 0 |
| `Clamp` | 17.4 ns | 3.4 ns | 0 |

## Format sensitivity

Same three binary ops across formats; `Binary8p1uf` exercises the wide-exponent-spread escalations, small-K formats the tiny-table regime.

| operation | median | min | allocs |
|---|---|---|---|
| `Add⟨Binary8p4se⟩` | 13.5 ns | 4.0 ns | 0 |
| `Divide⟨Binary8p4se⟩` | 990.3 ns | 8.3 ns | 0 |
| `Multiply⟨Binary8p4se⟩` | 13.0 ns | 3.1 ns | 0 |
| `Add⟨Binary8p3sf⟩` | 13.0 ns | 4.0 ns | 0 |
| `Divide⟨Binary8p3sf⟩` | 23.9 ns | 6.7 ns | 0 |
| `Multiply⟨Binary8p3sf⟩` | 11.9 ns | 3.1 ns | 0 |
| `Add⟨Binary8p1uf⟩` | 133.1 ns | 13.0 ns | 0 |
| `Divide⟨Binary8p1uf⟩` | 16.2 ns | 5.8 ns | 0 |
| `Multiply⟨Binary8p1uf⟩` | 12.1 ns | 3.1 ns | 0 |
| `Add⟨Binary5p2se⟩` | 13.6 ns | 4.0 ns | 0 |
| `Divide⟨Binary5p2se⟩` | 17.6 ns | 6.1 ns | 0 |
| `Multiply⟨Binary5p2se⟩` | 12.7 ns | 3.1 ns | 0 |
| `Add⟨Binary3p1se⟩` | 10.3 ns | 3.8 ns | 0 |
| `Divide⟨Binary3p1se⟩` | 12.5 ns | 5.8 ns | 0 |
| `Multiply⟨Binary3p1se⟩` | 9.7 ns | 2.7 ns | 0 |

## Projection by rounding/saturation mode

`project(Binary8p4se, ρ, x)` over the mode vocabulary (stochastic budgets N = 8).

| operation | median | min | allocs |
|---|---|---|---|
| `NearestTiesToEven` | 11.3 ns | 2.8 ns | 0 |
| `NearestTiesToAway` | 13.1 ns | 3.7 ns | 0 |
| `TowardPositive` | 10.1 ns | 4.0 ns | 0 |
| `TowardNegative` | 10.1 ns | 3.7 ns | 0 |
| `TowardZero` | 8.9 ns | 4.1 ns | 0 |
| `ToOdd` | 11.1 ns | 3.2 ns | 0 |
| `StochasticA[8]` | 12.1 ns | 2.8 ns | 0 |
| `StochasticB[8]` | 12.2 ns | 3.1 ns | 0 |
| `StochasticC[8]` | 12.1 ns | 3.1 ns | 0 |
| `RNE · SatFinite` | 10.7 ns | 3.1 ns | 0 |
| `RNE · SatPropagate` | 10.8 ns | 2.8 ns | 0 |

## Array kernels (vmap)

Warm caches: table specializations prebuilt, so table rows measure the gather; scalar-loop rows measure the full compute pipeline per element.

| operation | median | min | allocs | per element |
|---|---|---|---|---|
| `vmap unary (table gather), n=65536` | 15.59 μs | 14.46 μs | 0 | 0.24 ns/elem — 4.2 Gelem/s |
| `vmap binary (table gather), n=65536` | 30.34 μs | 28.6 μs | 0 | 0.46 ns/elem — 2.16 Gelem/s |
| `vmap ternary (scalar loop), n=65536` | 1.84 ms | 1.76 ms | 0 | 28.13 ns/elem — 0.04 Gelem/s |
| `vmap binary stochastic (scalar loop), n=65536` | 1.2 ms | 1.13 ms | 0 | 18.29 ns/elem — 0.05 Gelem/s |
| `vmap unary through PackedVector, n=65536` | 150.95 μs | 133.66 μs | 519 | 2.3 ns/elem — 0.43 Gelem/s |

## Sorting (64 K values)

Counting sort is installed as the default algorithm for `Binary` vectors.

| operation | median | min | allocs |
|---|---|---|---|
| `sort! (counting sort via defalg), n=65536` | 283.68 μs | 266.47 μs | 2 |
| `sort! (stock comparison sort), n=65536` | 2.41 ms | 2.31 ms | 3 |
| `sort! rev=true (counting sort), n=65536` | 285.38 μs | 269.09 μs | 2 |

## Table builds (oracle + projection, Float128-first)

Cold cache per sample (`empty_tables!` in untimed setup); JIT pre-warmed.

| operation | median | min | allocs |
|---|---|---|---|
| `Exp⟨8p4se⟩ (256 entries)` | 341.41 μs | 317.15 μs | 257 |
| `Tanh⟨8p4se⟩ (256 entries)` | 305.16 μs | 290.89 μs | 257 |
| `Add⟨8p4se×8p4se⟩ (64 K entries)` | 42.97 ms | 39.92 ms | 785668 |
| `Divide⟨8p4se×8p4se⟩ (64 K)` | 94.33 ms | 89.48 ms | 784906 |
| `Add⟨8p1uf×8p1uf⟩ (64 K, wide-spread)` | 78.79 ms | 75.33 ms | 1587030 |

## Block and scaled operations

Elements `Binary8p4se`, scales `Binary8p1uf`, B = 32.

| operation | median | min | allocs | per lane |
|---|---|---|---|---|
| `BlockAdd (B=32)` | 1.58 μs | 1.49 μs | 35 | 49.36 ns/lane |
| `BlockDotProduct → 8p4se (B=32)` | 707.0 ns | 361.0 ns | 3 | 22.09 ns/lane |
| `BlockReduceAdd → 8p4se (B=32)` | 783.2 ns | 32.6 ns | 0 | 24.48 ns/lane |
| `ConvertToBlockMaxAbsFinite (B=32)` | 3.25 μs | 3.06 μs | 75 | 101.71 ns/lane |

## Conversions and packed storage

| operation | median | min | allocs | per element |
|---|---|---|---|---|
| `T(::UInt8) code-point constructor` | 2.1 ns | 2.1 ns | 0 |  |
| `rawvalue (unchecked kernel route)` | 2.1 ns | 2.1 ns | 0 |  |
| `T(::Float64) numeric constructor (projects)` | 11.4 ns | 3.1 ns | 0 |  |
| `Convert 8p4se → 8p3se (scalar)` | 12.5 ns | 3.4 ns | 0 |  |
| `Float64 → 8p4se (project)` | 11.3 ns | 3.0 ns | 0 |  |
| `PackedVector pack, n=65536` | 62.36 μs | 57.73 μs | 3 | 0.95 ns/elem |
| `PackedVector unpack (collect), n=65536` | 51.03 μs | 49.81 μs | 3 | 0.78 ns/elem |

---
*All numbers from this machine/run; absolute values vary by host. Regenerate with `julia --project=benchmark benchmark/benchmarking.jl`.*