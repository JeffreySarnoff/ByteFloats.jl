# Float128FMA.jl — correctly-rounded `fma` for `Quadmath.Float128` on Windows

Quadmath.jl deliberately does not bind libquadmath's `fmaq` on Windows because
calling it subtly corrupts floating-point state there
([JuliaMath/Quadmath.jl#31](https://github.com/JuliaMath/Quadmath.jl/issues/31)).
As a result, `fma(::Float128, ::Float128, ::Float128)` on Windows throws
`ErrorException("fma not defined for Float128")` via Base's promotion fallback.

This package fills the gap with a **pure-Julia, allocation-free, IEEE 754
correctly-rounded** (round-to-nearest, ties-to-even) fused multiply-add that
is **bit-for-bit compatible with libquadmath's `fmaq`** — including NaN
payload propagation and invalid-operation NaN generation — so results on
Windows are indistinguishable from Linux/macOS. It never calls into
libquadmath (only `reinterpret` and integer arithmetic), so the FP-state
corruption of issue #31 cannot occur. Platform-independent, Julia ≥ 1.6,
developed and validated on Julia 1.12.6.

## Installation & usage

As a package (recommended):

```julia
pkg> dev path/to/Float128FMA      # or: Pkg.develop(path="path/to/Float128FMA")

julia> using Quadmath, Float128FMA   # on Windows this installs Base.fma

julia> fma(Float128(2)^60 + Float128(2)^-52,
           Float128(2)^60 - Float128(2)^-52,
           -Float128(2)^120)         # == -2^-104, the single-rounding answer
```

Or as a single file: `include("src/Float128FMA.jl"); using .Float128FMA`.

* `fma128(x, y, z)` — always this implementation, on every OS.
* `Base.fma` is extended **only when no Float128 method exists** (i.e. on
  Windows), so on Linux/macOS the native `fmaq` keeps priority and no
  method overwrite occurs. `Float128FMA.install!()` forces installation.
* Load order does not matter for callers that invoke `fma` dynamically:
  defining the method any time before first use is sufficient.

### Wiring into ByteFloats.jl on Windows

```julia
using Quadmath, Float128FMA
using ByteFloats            # Float128 Class-R fast paths now work on Windows
```

Without this, ByteFloats' `_f128()` fast paths in `oracle.jl`
(`Divide`/`Recip`/`Sqrt`/`RSqrt`) and `blocks.jl` (`_bp_element`) throw on
the first inexact quotient or root.

## Design

Classic soft-float `mulAdd` on the raw binary128 bit patterns:

1. Decompose x, y, z into sign / unbiased exponent / 113-bit significand,
   normalizing subnormal inputs.
2. Exact 113×113 → 226-bit product with four 64×64→128 multiplies.
3. Product and addend live in a 256-bit fixed-point accumulator (a pair of
   `UInt128`) with the leading significand bit pinned at bit 254.
4. Alignment uses **shift-right-jam** (any shifted-out bit sets the LSB).
   Jamming is exact for round-to-nearest here: a jam can only occur when the
   exponent gap exceeds ~30 (product shifted) or ~142 (addend shifted), in
   which case catastrophic cancellation is impossible and the jam bit stays
   ≥ 100 positions below the guard bit; conversely, whenever cancellation
   can occur (|gap| ≤ 1) the alignment shift loses no bits and the
   subtraction is exact.
5. One round-to-nearest-even at 113 bits (guard bit 141, sticky bits 0–140),
   with gradual underflow (an extra jam-shift into the subnormal grid before
   rounding; the min-normal rounding carry falls out of the bit assembly for
   free) and overflow to ±Inf.

### NaN semantics (bit-exact with libquadmath)

Determined empirically against `fmaq` and matching glibc soft-fp's x86
`_FP_CHOOSENAN`, applied in two stages:

* **Stage 1 (product):** if both x and y are NaN, the one with the larger
  raw 112-bit fraction field wins (quiet bit participates; tie → x). A
  single NaN wins outright. An invalid product (`0 × ∞`) generates the x86
  "indefinite" NaN `0xffff8000…0` (sign set, payload 0).
* **Quieting between stages:** the stage-1 NaN has its quiet bit set
  *before* stage 2 — this matters when both stages hold signaling NaNs.
* **Stage 2 (sum):** the (now quiet) stage-1 NaN competes with z's **raw**
  fraction by the same larger-fraction rule, tie → stage-1. `∞ − ∞` with no
  NaN inputs yields the indefinite NaN. The winner is returned with its
  quiet bit set, sign and payload preserved.

Floating-point exception *flags* are not modeled (Julia does not expose
them for Float128).

## Validation

`test/runtests.jl` (`Pkg.test()`) — on Linux every result is compared
**bit-for-bit** (NaN sign and payload included) against native `fmaq`;
thousands of cases per battery are additionally verified against an
independent oracle computing the exact value in 40000-bit BigFloat and
checking nearest/ties-to-even directly against both Float128 neighbours.
On Windows the BigFloat oracle and the directed NaN-rule tests run.

Batteries: uniform random bit patterns; mid/broad-range normals;
subnormal-heavy mixes plus a directed sweep of products landing around the
underflow threshold with all guard/sticky flavours; overflow boundary;
massive cancellation (`z = −(x ⊗ y)`) and the division-residual identity
`fl(x/y) exact ⇔ fma(q, y, −x) == 0` (exactly ByteFloats' witness);
constructed halfway ties at bit 113; a 5832-combination special-value grid;
scaled double-rounding traps; directed NaN-rule cases; and NaN/special fuzz.

Across development, **≈ 12 million cases were validated with zero
mismatches** against native `fmaq` on Julia 1.12.6 / x86-64 Linux
(including 7 million dedicated NaN/special-value fuzz cases).

## Performance

Per call, random full-precision operands, Julia 1.12.6, x86-64:

| implementation                     | ns/op | allocations |
|------------------------------------|------:|------------:|
| `fma128` (this package, pure Julia)|  ~33  | 0 |
| `fmaq` via ccall (Linux reference) | ~1060 | 0 |
| `x*y + z` in Float128 (2 roundings)|  ~55  | 0 |

~30× faster than round-tripping through libquadmath, and faster than a
single non-fused multiply-add: Float128 arithmetic is soft-float
everywhere, and this implementation stays in registers with native
64/128-bit integer ops. Inferred effects are `+c,+e,+t` (consistent,
effect-free, terminating), so constant arguments fold at compile time.

## Layout

```
Float128FMA/
├── Project.toml
├── README.md
├── src/Float128FMA.jl        # self-contained; also usable via include()
└── test/runtests.jl
```
