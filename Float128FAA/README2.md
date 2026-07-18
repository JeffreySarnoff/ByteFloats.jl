# Float64faa.jl & Float32faa.jl — fused add-add for Float64 and Float32

`faa(x, y, z)` computes the **correctly rounded** (IEEE 754 round-to-nearest,
ties-to-even) value of `x + y + z` with a **single rounding** — the additive
analogue of `fma` — for `Float64` (binary64) and `Float32` (binary32).

Same algorithm as the companion Float128FAA.jl, instantiated per format.
Pure Julia bit-level integer arithmetic, **zero dependencies** (not even
Quadmath), allocation-free, deterministic and identical on every platform.
Julia ≥ 1.6; developed and validated on Julia 1.12.6.

```julia
include("Float64faa.jl"); using .Float64FAA     # exports faa, faa64
include("Float32faa.jl"); using .Float32FAA     # exports faa, faa32

faa(1.0, 2.0^-53, 2.0^-106)      # 1.0000000000000002 == nextfloat(1.0)
(1.0 + 2.0^-53) + 2.0^-106       # 1.0 — sequential addition misses it

faa(floatmax(Float64), floatmax(Float64), -floatmax(Float64))
# == floatmax(Float64); sequential addition overflows to Inf
```

If both modules are loaded, the shared export `faa` becomes ambiguous — use
the sized aliases `faa64` / `faa32` (or qualify: `Float64FAA.faa`). There is
no `Base.faa`, so nothing is pirated.

## Design (shared with Float128FAA, rescaled)

Decompose into sign / unbiased exponent / p-bit significand (p = 53 or 24;
subnormals normalized); sort by magnitude |a| ≥ |b| ≥ |c|; accumulate in a
single-word fixed-point accumulator — `UInt128` for Float64 (leading bit
anchored at 126, guard bit 73), `UInt64` for Float32 (anchor 62, guard 38) —
folding b and c in with shift-right-jam alignment; round once to
nearest-even with gradual underflow and overflow to ±Inf.

The three-addend hazards, with per-format constants re-derived in each file:
the accumulator carry (bit 127 / 63) is renormalized after folding in b and
*before* c is aligned against the updated anchor (three same-sign operands
would otherwise overflow the word); a jam requires an alignment gap ≥ 75 /
40, which keeps every reachable partial sum large enough that jam bits stay
below bit 54 / 26 — strictly under the guard bit — after any normalization,
while deep cancellation (gap ≤ 1) is always bit-exact; and `a + b == 0`
exactly returns c directly.

## Semantics for special values

IEEE 754 addition applied to the fused sum: only actual **input** infinities
produce an infinite result; `+Inf` with `−Inf` gives the x86 indefinite NaN
(`0xfff8000000000000` / `0xffc00000`); an exact zero result is `+0`, and
three zeros give `−0` only when all three are `−0`. For every non-NaN
special input the result is bit-identical to native `(x + y) + z`.

NaN payload propagation is **deterministic**, using the same two-stage rule
as the Float128 family (larger raw fraction field wins, quiet bit included;
ties to the earlier stage; stage-1 survivor quieted before competing with z;
winner returned quieted). Hardware NaN propagation is platform-specific
(x86 and ARM differ, and compilers may commute operands), so payloads here
are well-defined rather than hardware-matching. Exception flags are not
modeled.

## Validation (`test_faa_64_32.jl`)

Per format: exact-value oracle in high-precision BigFloat (4096 / 512 bits —
both exceed the worst-case exact-sum width) checking nearest/ties-to-even
against both neighbours; directed exponent-gap sweeps across every alignment
regime including the all-positive accumulator-carry region; cancellation and
return-c paths under all operand orders; subnormal, underflow, and overflow
boundaries; constructed ties and double-rounding traps; a special-value grid
bit-compared against native sequential addition plus directed NaN-rule
cases; 10 million Boldo–Melquiond round-to-odd cross-checks per format
(built from native hardware ops, adjudicated by the BigFloat oracle); and
zero-allocation / inference checks.

A cross-format diagnostic compares `faa32(x,y,z)` against
`Float32(faa64(x,y,z))`. These deliberately **deviate** (~0.3% of random
cases): a three-addend exact sum can land exactly on a Float32 tie point
after the Float64 rounding, so rounding through a wider format is *not*
single rounding — the very reason `faa` exists. The suite requires the
BigFloat oracle to rule `faa32` correct in every deviation; across 10
million cases it did (31,985 deviations, all adjudicated in `faa32`'s
favor).

Result on Julia 1.12.6 / x86-64: **all 59 testsets pass; ≈ 21 million cases
per run with zero adjudicated failures.**

## Performance

Per call, random full-precision operands, Julia 1.12.6, x86-64:

| implementation                          | ns/op | allocations |
|------------------------------------------|------:|------------:|
| `faa64` (pure Julia, correctly rounded)  |  ~33  | 0 |
| `faa32` (pure Julia, correctly rounded)  |  ~28  | 0 |
| `(x+y)+z` hardware (two roundings)       | ~1–1.7| 0 |

Guaranteed single rounding costs roughly 20–30× two hardware additions —
tens of millions of fused sums per second, with inferred effects
`+c,+e,+t` (consistent, effect-free, terminating), so constant arguments
fold at compile time.

## Files

* `Float64faa.jl` — module `Float64FAA`, exports `faa`, `faa64`
* `Float32faa.jl` — module `Float32FAA`, exports `faa`, `faa32`
* `test_faa_64_32.jl` — joint validation suite (no dependencies beyond
  stdlib; run: `julia test_faa_64_32.jl`)
  