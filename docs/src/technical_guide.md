# Technical Guide

How ByteFloats.jl works inside: the encoding, the projection engine, the oracle and its
correctness protocol, the performance layers, and the verification doctrine that holds
it all together. Read the [User Guide](@ref) first; this page assumes its vocabulary.

## Layer map

Source files load in dependency order, each layer speaking only downward:

| layer | file | provides |
|---|---|---|
| formats | `formats.jl` | `Binary{K,P,SGN,EXT}`, the 120 named aliases, Group M queries |
| specs | `projspec.jl` | rounding/saturation singletons, `ProjSpec{R,S}` |
| codec | `decode_encode.jl` | decode (generated tables + bit-composed compute), encode, order keys, counting sort, `Class`, Next ops |
| engine | `project.jl` | `round_to_precision` (mask-based Float64 core + generic core), `saturate`, `project`, `project_interval` |
| ops | `ops_scalar.jl` | result-kind protocol, `apply_op`, the operation registry, both API registers |
| oracle | `oracle.jl` | ω-semantics for all 52 operations |
| tables | `tables.jl` | the pure-ρ result-table cache |
| kernels | `kernels.jl` | Shape-A gathers, Shape-B scalar loops, `vmap` |
| blocks | `blocks.jl` | `Block`, block/scaled ops, exact reductions |
| packed | `packed.jl` | sub-byte `PackedVector` storage |
| approx | `approx.jl` | κ measurement/registry, conformance declaration |

## Encoding and decoding

A value is its code point (`UInt8`). The bit layout is the draft's: sign (signed
formats), biased exponent, trailing significand; one NaN at the negative-zero slot;
no −0; ±Inf adjacent to the extremes in extended formats.

`decode` is a **`@generated` constant-tuple lookup**: per format, a `2^K`-tuple of
`Float64` datums built once from the computational decode (so table and computation
are correct by construction and asserted equivalent exhaustively). Constant inputs
still fold — `maxfinite_datum(T)` is a compile-time constant — while runtime decode
is a single indexed load (≈ 0.7 ns). The computational decode assembles the Float64
**by bits** (normalize with `leading_zeros`, place exponent and mantissa fields,
`reinterpret`), valid because every datum's exponent sits deep inside Float64's
normal range; `Float64` is thereby the *exact* carrier for all datums, an invariant
the suite checks against an independent big-float transliteration of the draft.

Ordering runs on **integer order keys**: a sign-magnitude fold into `UInt16`,
monotone with the total order, NaN mapped to `typemax`. Same-format `TotalOrder`,
`isless`, and the numeric comparisons are key comparisons (~1 ns); since a format has
at most `2^K + 1` distinct keys, vectors sort with an **O(n) counting sort**
installed via `Base.Sort.defalg` (forward and reverse orderings; anything exotic
falls back to the stock algorithm).

## The projection engine

`project(fr, ρ, X; R, sticky)` is the single write path into a code point:

```
RoundToPrecision  →  Saturate  →  Encode
```

**RoundToPrecision** produces a `Rounded(kind, sign, S, Q)` — an exact scaled
significand `S` and exponent `Q` per the draft's `Q = max(⌊log₂|X|⌋, 1−B) − P + 1`.
Two implementations, proven equivalent exhaustively:

- the **generic core** (`_rtp_core`) works on any carrier (`Float64`, `Float128`,
  `BigFloat`) via exact power-of-two scaling and a fraction ν compared against the
  mode's decision points;
- the **mask-based Float64 core** (`_rtp_f64`) extracts sign/exponent/significand
  fields directly, represents ν as a 128-bit fixed-point integer with an OR-mask
  sticky for bits shifted out, and evaluates every mode — including the stochastic
  `⌊ν·2^N⌋` comparisons and `RNITE` ties — as integer field tests. Specials, zeros,
  and (Convert-only-reachable) subnormal Float64 inputs bail to the generic core.

**Symbolic sticky** is how enclosures talk to the engine: `sticky ∈ {−1, 0, +1}`
declares the true value to be the carried value plus an infinitesimal of that sign.
The engine folds it into every comparison; the delicate case — true value just
*below* a representable dyadic — decrements into the previous binade and sets
ν = 1⁻ (encoded in the mask core as an all-ones fixed-point fraction, which
reproduces every predicate including the RNITE tie behavior). This is what lets
directed modes land exactly on asymptotes: `tanh → 1⁻` under `TowardNegative`
projects to `NextLessThan(1)`, automatically.

**Saturate** classifies against the format's range with two integer comparisons per
side: the extremal magnitude in canonical `(S, Q)` form is a constant-folded function
of the type parameters, the rounded value and the extremum share one lexicographic
`(Q, S)` order (subnormals and the lowest normal binade share the same `Q`), signed
formats use sign–magnitude symmetry, unsigned underflow is just the sign bit, and an
internal `HUGEQ` sentinel represents "finite but astronomically large" (the ν = 1⁻
image of an infinite endpoint) so directed `SatNone` clamps it to MaxFinite
correctly. The draft's twenty saturation rows then map the classification to
`as-is / MaxFinite / MinFinite / ±Inf / NaN`.

**Encode** is pure integer bit assembly, including the significand-carry
renormalization and the subnormal/normal field split.

## The oracle and the result-kind protocol

Every operation's defined result is computed by `ωeval`, which returns one of five
result kinds; `apply_op` fast-splits the common one and finishes the rest:

| kind | meaning | finished by |
|---|---|---|
| `Float64` | exact (specials; representable arithmetic) | direct `project` |
| `Float128` | exact by **width analysis** | direct `project` (Float128 carrier) |
| `BigExactF` | exact at 2200-bit precision (wide-spread tail) | `project` on `BigFloat` |
| `Enclose128F` | correctly-rounded Float128 **bracket** | sticky agreement, MPFR fallback |
| `EncloseF` | MPFR directed enclosure `f(prec)`, optional Float128 pre-filter `fq`, optional eager Float64 estimate `yd` | three-stage: `yd` → `fq` → interval protocol |

Two **rigor classes** govern every non-`Float64` path, and their arguments are never
mixed:

**Class R (unconditional).** (a) Sums of decoded datums are *exactly representable*
in `Float128` whenever operand bits + exponent spread fit 113 bits — checked in
advance by integer exponent arithmetic (`Add` at ΔE ≤ 100, `FMA` ≤ 92, `FAA` span
≤ 98); beyond, the 2200-bit path takes over. (b) IEEE mandates correct rounding for
division and `sqrt` at *every* binary width. `Divide`/`Recip` therefore use the
**Float64** quotient directly: it is within half an ulp of truth, so it serves as
`EncloseF`'s eager `yd` estimate with no Float128 arithmetic on the path at all.
`Sqrt`/`RSqrt` use the **Float128** CR result: an inexact nearest-CR value `q`
brackets the truth in the *open* interval `(prevfloat(q), nextfloat(q))` — no
accuracy assumption in either case. Exactness itself is detected by an `fma`
residual test.

**Class E (envelope-conditional).** Faithful-but-not-CR evaluations stand in for an
enclosure only inside a generous envelope, resolved by the two-sided sticky gate.
Two stages, cheapest first. *Float64 stage:* for operations whose Float64 libm is
faithful (≤ 1 ulp ≈ 2⁻⁵² relative — the exp/log families, the hyperbolics, forward
trig inside the |x| ≤ 10¹⁵ reduction window, `atan`, `asinh`, the softplus
composition), the eager estimate `yd` carries envelope `E = |yd|·2⁻⁴⁵`, ≥ 2⁷
slacker than any faithful-libm error (measured margin on this machine: ≥ 2⁶·⁶).
*Float128 stage:* libquadmath's estimate `y = fq()` carries `E = |y|·2⁻⁹⁰`, at
least 2¹⁸ slacker than any published libquadmath bound. In both stages, if the
sticky-projected endpoints agree, that code point is the answer; if not, the next
stage (ultimately the MPFR ladder) decides — so an envelope failure can only cost
speed, never correctness, unless both endpoints agree *and* the envelope is wrong,
which the differential-build tests rule out empirically by building every standard
table twice (Float128-first and pure-MPFR) and byte-diffing, and which the
exhaustive oracle cross-check re-verifies per operation against the 3072-bit
protocol.

**The interval protocol** (`project_interval`) is the termination backbone: evaluate
the MPFR enclosure at precision `p`; if `lo == hi` the value is exactly
representable — project it; otherwise project `lo` with `sticky = +1` and `hi` with
`sticky = −1`; agreement (projection is monotone at fixed `R`) yields the answer;
disagreement doubles `p` up to 4096 bits. Termination requires that the enclosure
not be chasing a value the grid can actually hit — which is why the π-scaled
operations carry **Niven peels**: for dyadic arguments, `tan(πr)` takes rational
values only at the quarter-integers (±1) and `atan2/π` only on the diagonals
(±¼, ±¾); those cases are answered exactly before any enclosure is built, and
Niven's theorem proves the peel set complete.

## Tables and kernels

For pure specs, unary and binary operations are **finite functions** — 256 or 65 536
entries — so the kernel layer materializes them once per `(op, formats, ρ)` into a
locked cache (`Dict{TableKey, Memory{UInt8}}`, double-checked locking, builds outside
the lock) and serves every later array call as a gather: Shape-A, one load per
element, measured 0.27 ns/elem unary and 0.5 ns/elem binary. Tables are built
*through the scalar path*, so they inherit its bit-exactness; the suite asserts
table ≡ scalar over every entry. Ternary and stochastic calls take Shape-B — the
fully specialized scalar pipeline per element, with the RNG resolved once per array.

## Blocks: exactness without a superaccumulator

`blockdecode` produces each lane's `scale × element` exactly (≤ 17-bit significands
in Float64). Reductions then apply **span filters**: one integer pass over lane
exponents decides whether the whole sum (or dot product, with ≤ 33-bit lane
products) is exactly representable in `Float128`; if so, a plain `Float128`
accumulation *is* the exact answer; if not, an exact big-float accumulation takes
over. Either way there is exactly one projection, at the end. `ωBlockProject`
follows the draft's special rows for the scale (NaN, 0, ±Inf) and divides each
element result by the scale through its own cheapest-first CR-bracket / enclosure
cascade (exact Float64 quotient → CR Float128 → bracket/pre-filter → MPFR
interval), mirroring the scalar quotient group's rigor arguments.

## The κ registry and conformance

κ is the maximum code-point distance (along the total order) between an
implementation's result and the defined result, over inputs with finite defined
results; any mismatch on non-finite defined results makes κ = NaN. Because inputs
are enumerable, `register_approx!` *measures* κ exhaustively at registration and
rejects understatement — a declared bound is a verified property, not a promise.
`conformance()` assembles the declaration live from the operation registry, the
table cache, and the approximation registry.

## Verification doctrine

The value sets are small enough that sampling is never necessary, so the suite
enumerates: formats against an independent draft transliteration (14 679); ordering
over all 2.5 M same-format pairs plus Next-op edge tables (7.6 M); every unary
operation on every input against a 3072-bit protocol run; divide and the ternaries
exhaustively; stochastic R-sweeps with directed-asymptote pins; table ≡ scalar;
blocks against a from-scratch reference composition; Float128 carrier ≡ Float64
carrier and Float128-first ≡ MPFR differential builds; the mask rounding core ≡ the
generic core over datums and reachable sums/products at boundary stochastic budgets
N ∈ {45, 60}; packed round-trips; and κ/conformance behavior — ≈ 8.8 M assertions.
Deterministic **specialization regressions** (concrete inferred return types at the
public entry points, zero warm-path allocation) stand in for timing assertions.

## Benchmark doctrine

Recorded after two measurement post-mortems: a benchmark closure over any non-`const`
global measures Julia's dispatch machinery, not the code under test — and it distorts
*ratios*, not just absolutes, because dispatch cost varies with call shape (a dynamic
keyword call costs ~1 µs; a dynamic positional call far less; six unresolved interior
sites cost six times one). The rules the shipped `benchmark/benchmarking.jl`
enforces structurally: format types enter as type parameters; operands come from
untimed setup; functions retrieved reflectively pass through argument barriers to
specialize; and a preflight aborts the run if warm scalar paths allocate. The
table-build section reports both the cold build (cache evicted per sample) and the
steady-state warm cache hit, since callers amortize the former through the latter.
Vary one binding per variant; verify specialization before believing a number.

## Deliberate limitations

No implicit cross-format arithmetic (promotion is to `Float64`, explicitly). No
in-place packed arithmetic. No hidden threading (kernels are single-threaded in this
version). `Irrational`/`Rational` inputs to `Convert` are rejected rather than
double-rounded silently. The `Float128` machinery never changes results — disabling
it (`ByteFloats_Float128=disable`) is a tested no-op semantically.
