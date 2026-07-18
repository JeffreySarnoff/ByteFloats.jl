# Float128FAA.jl — `faa(x, y, z)`: fused add-add for `Quadmath.Float128`

`faa(x, y, z)` computes the **correctly rounded** (IEEE 754 round-to-nearest,
ties-to-even) value of `x + y + z` with a **single rounding** — the additive
analogue of `fma`. Pure Julia, allocation-free, no calls into libquadmath, so
it works identically on every platform including Windows (where Float128
`fma` is unavailable, see
[JuliaMath/Quadmath.jl#31](https://github.com/JuliaMath/Quadmath.jl/issues/31)).
Julia ≥ 1.6; developed and validated on Julia 1.12.6. Companion package to
Float128FMA.jl (each is self-contained; they can be used together or alone).

```julia
pkg> dev path/to/Float128FAA

julia> using Quadmath, Float128FAA

julia> faa(Float128(1), Float128(2)^-113, Float128(2)^-226)
1.00000000000000000000000000000000019e+00        # nextfloat(1.0)

julia> (Float128(1) + Float128(2)^-113) + Float128(2)^-226
1.0                                              # sequential addition misses it
```

`faa128` is an alias. There is no `Base.faa`, so nothing is pirated; call the
module's function directly. Single-file use also works:
`include("src/Float128FAA.jl"); using .Float128FAA`.

## Design

Soft-float three-way addition on the raw binary128 bit patterns: decompose
into sign / unbiased exponent / 113-bit significand (subnormals normalized),
sort by magnitude |a| ≥ |b| ≥ |c|, accumulate in 256-bit fixed point anchored
at the largest operand (leading bit at position 254), fold in b and c with
shift-right-jam alignment, then round once to nearest-even at 113 bits with
gradual underflow and overflow to ±Inf.

The three-addend hazards, and why the algorithm is exact:

* **Accumulator overflow.** Three same-sign operands can reach ~3·2^255 in
  accumulator units, past 256 bits — so the bit-255 carry is renormalized
  after folding in b, *before* c is aligned (against the updated anchor).
  A carry-jam bit can only exist while the partial sum keeps its leading bit
  at ≥ 253, leaving it ≥ 140 positions below the rounding position.
* **Jam vs. cancellation.** A jam requires an alignment shift > 142, making
  the jammed operand < 2^112 while every reachable partial sum stays
  ≥ 2^141 − 2^112, so jam bits always end ≥ 28 positions below the guard bit
  after normalization. Conversely, whenever deep cancellation is possible
  (exponent gap ≤ 1 between the cancelling pair) the alignment shift loses
  no bits and the subtraction is exact.
* **Total cancellation of the two largest.** If a + b == 0 exactly, the
  answer is c itself — detected and returned directly (trivially correctly
  rounded), rather than reading c's jammed remains out of the accumulator.

## Semantics for special values

IEEE 754 addition applied to the *fused* sum:

* Only actual **input** infinities produce an infinite result; finite
  operands whose partial sums would overflow do not
  (`faa(floatmax, floatmax, -floatmax) == floatmax`, while sequential
  addition overflows to `Inf`).
* `+Inf` and `−Inf` among the inputs → the x86 "indefinite" NaN
  `0xffff8000…0`.
* An exact zero result is `+0` in round-to-nearest; three zeros give `−0`
  only when all three are `−0`.
* NaN propagation mirrors libquadmath's addition chain `(x + y) + z`
  (glibc soft-fp, x86 `_FP_CHOOSENAN`): among stage-1 NaNs (x, y) the larger
  raw 112-bit fraction wins (ties to x); the survivor is quieted and then
  competes with z's raw fraction (ties to the survivor); the winner is
  returned quieted with sign and payload preserved. For any special-value
  inputs — where no rounding is involved — the result is bit-identical to
  native `(x + y) + z`.
* Floating-point exception flags are not modeled.

## Validation

No native `faaq` exists, so three independent oracles are used
(`test/runtests.jl`, run via `Pkg.test()`):

1. **Exact-value oracle**: x + y + z computed in 40000-bit BigFloat, with
   the round-to-nearest / ties-to-even property checked directly against
   both Float128 neighbours (no BigFloat→Float128 conversion involved).
2. **Sequential-native oracle** for special values (no rounding occurs, so
   fused ≡ sequential): a full special-value grid including NaN payload,
   sign, and quiet-bit probes, plus 200k special-mix fuzz cases — all
   bit-identical to native `(x + y) + z`.
3. **Boldo–Melquiond cross-check**: an independent correctly-rounded sum
   built from native Float128 operations via error-free transforms and
   emulated round-to-odd, mass-compared in a safe exponent range with any
   disagreement adjudicated by the BigFloat oracle.

Batteries: broad random; a directed exponent-gap sweep across every
alignment regime (gaps 0, 1, …, 113, 142, 143, …, 30000 — this sweep caught
the accumulator-overflow bug during development); cancellation and residual
patterns including the return-c path under all operand orders; subnormal and
underflow-boundary mixes; overflow boundary; constructed ties at bit 113 and
double-rounding traps. **≈ 6.6 million cases validated with zero failures**,
including 4 million adjudicated Boldo–Melquiond cross-checks (zero
disagreements) and the all-positive same-exponent carry region specifically.

## Performance

Per call, random full-precision operands, Julia 1.12.6, x86-64:

| implementation                              | ns/op | allocations |
|---------------------------------------------|------:|------------:|
| `faa` (this package, pure Julia)            |  ~55  | 0 |
| Boldo–Melquiond from native Float128 ops    | ~380  | 0 |
| `(x + y) + z` native (NOT correctly rounded)|  ~65  | 0 |

Correctly rounded and still faster than two plain libquadmath additions.
Inferred effects are `+c,+e,+t` (consistent, effect-free, terminating), so
constant arguments fold at compile time.

## Note for ByteFloats.jl

`faa` computes exactly the quantity ByteFloats' `ωeval(Val(:FAA), x, y, z)`
defines (with Float64 inputs converted exactly to Float128): a candidate
replacement for the `_bigsum3` MPFR fallback when the result is
representable, or a fast pre-filter in front of it.

## Layout

```
Float128FAA/
├── Project.toml
├── README.md
├── src/Float128FAA.jl        # self-contained; also usable via include()
└── test/runtests.jl
```