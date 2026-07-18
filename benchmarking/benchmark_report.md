# ByteFloats.jl benchmark report

Generated: 2026-07-18 05:49 UTC  ·  Julia 1.12.6  ·  cascadelake (1 threads)  ·  Float128 paths: enabled  ·  Chairmarks 1.3.1

Reference format for per-operation tables: `Binary8p4se` under `(NearestTiesToEven, SatNone)`; operand pools draw uniformly over all code points (an honest NaN/±Inf mix). Times are per call; medians with minima alongside. Methodology per the recorded benchmark doctrine: type-parameterized barriers, untimed setup, specialization preflight.

## Core primitives

The decode/compare/step/classify layer plus the projection engine.

| operation | median | min | allocs |
|---|---|---|---|
| `decode` | 2.9 ns | 2.9 ns | 0 |
| `order_key` | 3.3 ns | 3.0 ns | 0 |
| `x < y` | 3.5 ns | 3.5 ns | 0 |
| `TotalOrder` | 4.4 ns | 3.8 ns | 0 |
| `Class` | 3.8 ns | 2.7 ns | 0 |
| `NextGreaterThan` | 3.8 ns | 3.2 ns | 0 |
| `project (RNE·SatNone)` | 23.7 ns | 3.6 ns | 0 |
| `project (StochasticA[8], R drawn)` | 23.3 ns | 3.6 ns | 0 |

## Scalar operations — unary (30)

Sorted by median. Random operands: transcendental rows mix special-row fast returns with enclosure-path evaluations, so these are *scalar-path* costs; bulk unary work routes through 256-byte tables (see Array kernels).

| operation | median | min | allocs |
|---|---|---|---|
| `ArcSinPi` | 7.1 ns | 3.3 ns | 0 |
| `Log` | 7.4 ns | 3.3 ns | 0 |
| `Log2` | 8.1 ns | 3.6 ns | 0 |
| `ArcCos` | 14.5 ns | 14.5 ns | 0 |
| `ArcCosh` | 14.8 ns | 14.8 ns | 0 |
| `ArcSin` | 14.9 ns | 14.9 ns | 0 |
| `Abs` | 20.7 ns | 3.6 ns | 0 |
| `Negate` | 23.0 ns | 3.9 ns | 0 |
| `ArcTanh` | 30.9 ns | 5.3 ns | 0 |
| `Sqrt` | 32.4 ns | 9.3 ns | 0 |
| `Recip` | 51.9 ns | 3.7 ns | 0 |
| `Exp` | 57.7 ns | 3.7 ns | 0 |
| `Tanh` | 58.2 ns | 4.0 ns | 0 |
| `Exp2` | 58.6 ns | 12.7 ns | 0 |
| `Sin` | 59.9 ns | 3.4 ns | 0 |
| `Cosh` | 60.1 ns | 3.7 ns | 0 |
| `Cos` | 60.7 ns | 3.3 ns | 0 |
| `Sinh` | 64.6 ns | 3.4 ns | 0 |
| `LogOnePlus` | 66.6 ns | 3.6 ns | 0 |
| `ArcTan` | 68.0 ns | 4.3 ns | 0 |
| `Tan` | 70.0 ns | 3.1 ns | 0 |
| `Softplus` | 78.9 ns | 3.7 ns | 0 |
| `ArcSinh` | 79.3 ns | 3.7 ns | 0 |
| `SinPi` | 813.4 ns | 5.2 ns | 0 |
| `CosPi` | 859.9 ns | 3.7 ns | 0 |
| `TanPi` | 1.05 μs | 4.0 ns | 0 |
| `ArcTanPi` | 1.09 μs | 5.7 ns | 0 |
| `ArcCosPi` | 1.27 μs | 1.27 μs | 0 |
| `ExpMinusOne` | 1.42 μs | 7.4 ns | 0 |
| `RSqrt` | 1.57 μs | 1.57 μs | 0 |

## Scalar operations — binary (18)

Sorted by median.

| operation | median | min | allocs |
|---|---|---|---|
| `Maximum` | 20.7 ns | 3.3 ns | 0 |
| `Minimum` | 20.7 ns | 3.3 ns | 0 |
| `CopySign` | 21.3 ns | 3.3 ns | 0 |
| `MinimumFinite` | 21.9 ns | 14.0 ns | 0 |
| `MaximumMagnitudeNumber` | 21.9 ns | 13.7 ns | 0 |
| `MinimumMagnitudeNumber` | 22.7 ns | 15.0 ns | 0 |
| `MaximumFinite` | 22.7 ns | 13.9 ns | 0 |
| `MinimumNumber` | 22.9 ns | 13.0 ns | 0 |
| `Multiply` | 22.9 ns | 3.3 ns | 0 |
| `MaximumMagnitude` | 23.0 ns | 3.6 ns | 0 |
| `MinimumMagnitude` | 24.5 ns | 3.6 ns | 0 |
| `MaximumNumber` | 25.0 ns | 13.0 ns | 0 |
| `Add` | 27.1 ns | 4.8 ns | 0 |
| `Subtract` | 31.3 ns | 8.9 ns | 0 |
| `Divide` | 47.5 ns | 4.8 ns | 0 |
| `Hypot` | 739.0 ns | 18.5 ns | 0 |
| `ArcTan2` | 1.19 μs | 7.1 ns | 0 |
| `ArcTan2Pi` | 1.27 μs | 6.0 ns | 0 |

## Scalar operations — ternary (3)

Sorted by median.

| operation | median | min | allocs |
|---|---|---|---|
| `Clamp` | 23.9 ns | 3.6 ns | 0 |
| `FMA` | 31.3 ns | 8.5 ns | 0 |
| `FAA` | 32.5 ns | 9.0 ns | 0 |

## Format sensitivity

Same three binary ops across formats; `Binary8p1uf` exercises the wide-exponent-spread escalations, small-K formats the tiny-table regime.

| operation | median | min | allocs |
|---|---|---|---|
| `Add⟨Binary8p4se⟩` | 27.1 ns | 4.8 ns | 0 |
| `Divide⟨Binary8p4se⟩` | 47.5 ns | 4.8 ns | 0 |
| `Multiply⟨Binary8p4se⟩` | 22.9 ns | 3.4 ns | 0 |
| `Add⟨Binary8p3sf⟩` | 24.0 ns | 4.8 ns | 0 |
| `Divide⟨Binary8p3sf⟩` | 45.1 ns | 5.1 ns | 0 |
| `Multiply⟨Binary8p3sf⟩` | 21.9 ns | 3.3 ns | 0 |
| `Add⟨Binary8p1uf⟩` | 179.7 ns | 7.7 ns | 0 |
| `Divide⟨Binary8p1uf⟩` | 23.9 ns | 4.8 ns | 0 |
| `Multiply⟨Binary8p1uf⟩` | 21.0 ns | 3.3 ns | 0 |
| `Add⟨Binary5p2se⟩` | 26.8 ns | 4.7 ns | 0 |
| `Divide⟨Binary5p2se⟩` | 26.2 ns | 5.0 ns | 0 |
| `Multiply⟨Binary5p2se⟩` | 22.7 ns | 3.6 ns | 0 |
| `Add⟨Binary3p1se⟩` | 17.1 ns | 4.7 ns | 0 |
| `Divide⟨Binary3p1se⟩` | 15.6 ns | 5.0 ns | 0 |
| `Multiply⟨Binary3p1se⟩` | 15.4 ns | 3.6 ns | 0 |

## Projection by rounding/saturation mode

`project(Binary8p4se, ρ, x)` over the mode vocabulary (stochastic budgets N = 8).

| operation | median | min | allocs |
|---|---|---|---|
| `NearestTiesToEven` | 22.4 ns | 4.4 ns | 0 |
| `NearestTiesToAway` | 19.4 ns | 4.7 ns | 0 |
| `TowardPositive` | 18.6 ns | 4.7 ns | 0 |
| `TowardNegative` | 17.7 ns | 4.7 ns | 0 |
| `TowardZero` | 18.3 ns | 4.7 ns | 0 |
| `ToOdd` | 23.3 ns | 3.9 ns | 0 |
| `StochasticA[8]` | 24.2 ns | 3.6 ns | 0 |
| `StochasticB[8]` | 23.0 ns | 3.6 ns | 0 |
| `StochasticC[8]` | 22.1 ns | 3.9 ns | 0 |
| `RNE · SatFinite` | 20.0 ns | 3.9 ns | 0 |
| `RNE · SatPropagate` | 19.7 ns | 3.6 ns | 0 |

## Array kernels (vmap)

Warm caches: table specializations prebuilt, so table rows measure the gather; scalar-loop rows measure the full compute pipeline per element.

| operation | median | min | allocs | per element |
|---|---|---|---|---|
| `vmap unary (table gather), n=65536` | 36.73 μs | 33.29 μs | 0 | 0.56 ns/elem — 1.78 Gelem/s |
| `vmap binary (table gather), n=65536` | 51.88 μs | 48.78 μs | 0 | 0.79 ns/elem — 1.26 Gelem/s |
| `vmap ternary (scalar loop), n=65536` | 2.32 ms | 2.24 ms | 0 | 35.46 ns/elem — 0.03 Gelem/s |
| `vmap binary stochastic (scalar loop), n=65536` | 1.85 ms | 1.77 ms | 0 | 28.18 ns/elem — 0.04 Gelem/s |
| `vmap unary through PackedVector, n=65536` | 242.67 μs | 229.56 μs | 519 | 3.7 ns/elem — 0.27 Gelem/s |

## Sorting (64 K values)

Counting sort is installed as the default algorithm for `Binary` vectors.

| operation | median | min | allocs |
|---|---|---|---|
| `sort! (counting sort via defalg), n=65536` | 331.77 μs | 321.98 μs | 2 |
| `sort! (stock comparison sort), n=65536` | 3.99 ms | 3.93 ms | 3 |
| `sort! rev=true (counting sort), n=65536` | 418.34 μs | 402.14 μs | 2 |

## Table builds (oracle + projection, Float128-first)

Cold cache per sample (`empty_tables!` in untimed setup); JIT pre-warmed. The warm-hit column is the steady-state cost of `get_table` when the specialization is already cached (median / min).

| operation | median | min | allocs | warm hit |
|---|---|---|---|---|
| `Exp⟨8p4se⟩ (256 entries)` | 87.93 μs | 74.11 μs | 257 | 179.5 ns / 175.7 ns |
| `Tanh⟨8p4se⟩ (256 entries)` | 93.28 μs | 85.72 μs | 257 | 179.2 ns / 175.5 ns |
| `Add⟨8p4se×8p4se⟩ (64 K entries)` | 69.74 ms | 65.44 ms | 785668 | 178.6 ns / 176.3 ns |
| `Divide⟨8p4se×8p4se⟩ (64 K)` | 72.93 ms | 65.64 ms | 784906 | 178.8 ns / 176.5 ns |
| `Add⟨8p1uf×8p1uf⟩ (64 K, wide-spread)` | 115.86 ms | 108.59 ms | 1587030 | 179.7 ns / 176.3 ns |

## Block and scaled operations

Elements `Binary8p4se`, scales `Binary8p1uf`, B = 32.

| operation | median | min | allocs | per lane |
|---|---|---|---|---|
| `BlockAdd (B=32)` | 2.31 μs | 2.15 μs | 35 | 72.19 ns/lane |
| `BlockDotProduct → 8p4se (B=32)` | 799.0 ns | 472.0 ns | 3 | 24.97 ns/lane |
| `BlockReduceAdd → 8p4se (B=32)` | 991.1 ns | 56.5 ns | 0 | 30.97 ns/lane |
| `ConvertToBlockMaxAbsFinite (B=32)` | 4.38 μs | 4.17 μs | 75 | 136.97 ns/lane |

## Conversions and packed storage

| operation | median | min | allocs | per element |
|---|---|---|---|---|
| `T(::UInt8) code-point constructor` | 2.9 ns | 2.9 ns | 0 |  |
| `rawvalue (unchecked kernel route)` | 2.9 ns | 2.9 ns | 0 |  |
| `T(::Float64) numeric constructor (projects)` | 23.9 ns | 3.6 ns | 0 |  |
| `Convert 8p4se → 8p3se (scalar)` | 24.5 ns | 3.3 ns | 0 |  |
| `Float64 → 8p4se (project)` | 19.5 ns | 3.9 ns | 0 |  |
| `PackedVector pack, n=65536` | 79.53 μs | 73.22 μs | 3 | 1.21 ns/elem |
| `PackedVector unpack (collect), n=65536` | 59.39 μs | 58.21 μs | 3 | 0.91 ns/elem |

---
*All numbers from this machine/run; absolute values vary by host. Regenerate with `julia --project=benchmark benchmark/benchmarking.jl`.*
