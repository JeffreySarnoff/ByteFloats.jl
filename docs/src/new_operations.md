# Adding Operations

How to add a new operation to ByteFloats.jl so that it inherits everything the
existing catalog has: bit-exact defined results, the full projection-mode
vocabulary, generated scalar/array/block surfaces, table caching, exhaustive
test coverage, benchmark pickup, and a truthful conformance declaration. Read
the [Technical Guide](@ref) first — this page assumes its vocabulary (the
result-kind protocol, the two rigor classes, the interval protocol, the `yd`
envelope stage).

The architecture is registry-driven: for the supported arities, **most of the
surface is generated**. Adding an operation is mostly (1) one registry line and
(2) one rigorous `ωeval` method — plus the *duties* that keep the package's
guarantees true. Every subsection below states the duties explicitly; they are
not optional.

## Universal duties (every new operation)

1. **Special rows are spec, not optimization.** Write out the NaN/±Inf/zero
   rows explicitly in `ωeval`, in draft style, before any general evaluation.
   Mark any behavior the draft text under-determines with `[interp]`.
2. **Result-kind protocol.** `ωeval` must return one of
   `Float64 | Float128 | StickyF | BigExactF | EncloseF | Enclose128F`, with
   `Float64`/`Float128` reserved for *provably exact* results and `StickyF`
   for exact-head-plus-tail-direction results whose soundness bound is
   discharged at the emitting site (see the `FMA`/`FAA` wide-spread
   escalations). Keep the return union narrow (ideally `{Float64, EncloseF}`)
   — a wide union defeats inference and costs allocations (measured: a 4-type
   union turned a 50 ns op into 240 ns with 2 allocations).
3. **Termination (the Niven duty).** `project_interval` terminates only if the
   enclosure is not chasing a value the projection grid can actually hit. Any
   input whose *exact* result is a representable dyadic must be **peeled or
   detected** before an enclosure is built. For transcendental ops this is
   usually a finite set of special rows (prove completeness, as the π-family
   does with Niven's theorem); for algebraic ops it is an exactness test.
4. **Envelope justification.** If you supply an eager `yd` estimate, you owe an
   error bound: the true value must lie within `|yd|·2⁻⁴⁵` (`_F64_RELEXP`).
   Acceptable bases, strongest first: IEEE-CR hardware ops (÷, sqrt: ≤ ½ ulp,
   unconditional); short compositions of CR ops (≤ ~1.5 ulp); faithful libm
   (≤ 1 ulp, an empirical-but-published claim); short faithful compositions
   with condition number ≤ 1 (≤ ~3 ulp). **Measure it**: sample your estimator
   against ≥ 256-bit MPFR truth over the op's domain, including the treacherous
   regions (near zeros of the result, domain boundaries), and demand ≥ 2⁶ slack
   under 2⁻⁴⁵. The same duty applies to a Float128 `fq` estimator against its
   2⁻⁹⁰ envelope. Beware **argument-error amplification**: an estimate computed
   as `f(g(x))` where `g` rounds (e.g. `sin(π·r)` in Float64) can be
   catastrophically worse near zeros of `f` than a native result-domain
   evaluation (`sinpi(r)`) — measure the composition you actually ship.
5. **Ladder validity.** The MPFR closure must return a genuine directed
   enclosure. Whole-expression `setrounding` is valid only when the composition
   is monotone in every intermediate (state the argument in a comment, as
   `Softplus` does); otherwise round each step in the correct direction
   explicitly (as `_mpfr_rsqrt` does for the anti-monotone `1/√x`).
6. **Gates.** Run the full suite. The unary exhaustive ×hp gate and the
   coverage check pick up registry ops automatically; extend the hand-curated
   lists where applicable (monotonicity ops, binary exhaustive sweeps). Then
   regenerate the benchmark report — `bench_scalar_ops` iterates the registry
   lists, so your op appears in all four operand-class tables automatically;
   **add a `_SAFE_DOMAINS` entry** in `benchmarking/benchmarking.jl` if the op is
   argument-restricted, or the safe-args table will use the oracle-derived
   fallback pool.

## Scalar operations

The registry lives in `src/ops_scalar.jl`. Adding a name to `_UNARY_OPS` /
`_BINARY_OPS` / `_TERNARY_OPS` causes, with no further work:

- registration in `OP_REGISTRY` (group `:B`/`:C`/`:A` per the loop defaults —
  override with an explicit `register_op!` if the default group is wrong);
- the generated spec-register method `Op(fr, ρ, operands…; rng, R)` and the
  same-format convenience `Op(x…)`;
- generated array methods routing through `vmap` (Shape-A tables for pure ρ —
  always at arity ≤ 2, bitwidth-gated at arity 3 — Shape-B scalar/threaded loops
  otherwise);
- generated `BlockOp` and `ScaledOp` variants (`src/blocks.jl`);
- the export (`src/ByteFloats.jl` loops the registry);
- table cacheability, conformance listing, and test-harness coverage.

What is *not* generated is the mathematics: `ωeval` in `src/oracle.jl`.

### Adding a unary operation

Worked example: **`Cbrt`** (cube root), chosen because it exercises the
hardest duty — exactness detection for an algebraic op.

**Step 1 — registry** (`src/ops_scalar.jl`): append `:Cbrt` to `_UNARY_OPS`.
The loop registers it as group `:B` with arity 1.

**Step 2 — ω-semantics** (`src/oracle.jl`), with every duty discharged:

```julia
# Exact-cube detection (termination duty): x = ±m·2^(3k) with m an odd perfect
# cube ⇒ cbrt(x) is a Float64-exact dyadic that project_interval must never
# chase. Oracle-path cost is irrelevant; correctness of the peel is everything.
function _exact_cbrt(x::Float64)
    fr, e = frexp(abs(x))               # |x| = fr·2^e, fr ∈ [0.5, 1)
    m = Int64(ldexp(fr, 53)); e -= 53   # |x| = m·2^e exactly
    tz = trailing_zeros(m); m >>= tz; e += tz
    e % 3 == 0 || return nothing
    r = round(Int64, cbrt(Float64(m)))
    for c in max(r - 1, 1):(r + 1)      # cbrt is faithful ⇒ true root ∈ r ± 1
        Int128(c)^3 == Int128(m) && return flipsign(ldexp(Float64(c), e ÷ 3), x)
    end
    nothing
end

function ωeval(::Val{:Cbrt}, x::Float64)
    isnan(x) && return NaN                      # special rows first: they are spec
    isinf(x) && return x                        # cbrt(±∞) = ±∞
    iszero(x) && return 0.0
    c = _exact_cbrt(x)
    c === nothing || return c                   # exact dyadic result — peeled
    # yd: Float64 cbrt is faithful (≤ 1 ulp ≈ 2⁻⁵²; measure it) ⇒ ≥ 2⁷ slack
    # under the 2⁻⁴⁵ envelope. fq: libquadmath cbrtq under the 2⁻⁹⁰ envelope.
    # Ladder: cbrt is monotone ⇒ whole-expression directed rounding encloses.
    _encl1(cbrt, x; fq=() -> cbrt(Float128(x)), yd=cbrt(x))
end
```

**Step 3 — optional surface**: a Base veneer (`Base.cbrt(x::Binary) = Cbrt(x)`
next to the others in `ops_scalar.jl`), a docstring, and — since cbrt is
monotone — its symbol in the test suite's monotonicity list.

**Step 4 — validate**: envelope measurement for `cbrt` (include tiny and huge
`x`; cbrt has no dangerous regions — the result is never near zero except at
zero, which is peeled), then the full suite: the exhaustive ×hp gate now checks
`Cbrt` on every code point of the reference format under four modes against the
3072-bit protocol, automatically.

### Adding a binary operation

Worked example: **`LogAddExp`** — `log(eˣ + eʸ)` — chosen because it exercises
composed estimators and the ladder-monotonicity argument.

**Step 1 — registry**: append `:LogAddExp` to `_BINARY_OPS` (default group
`:C`; the special five arithmetic names get `:A`).

**Step 2 — ω-semantics**. One generic evaluator serves all three precisions —
Float64 (`yd`), Float128 (`fq`), and BigFloat (the ladder):

```julia
# Stable form: h + log1p(exp(l − h)), h = max, l = min. Every step is monotone
# nondecreasing in (x, y), so whole-expression directed rounding in _mpfr2 is a
# valid enclosure (the same argument Softplus records).
_logaddexp(a, b) = (h = max(a, b); l = min(a, b); h + log1p(exp(l - h)))

function ωeval(::Val{:LogAddExp}, x::Float64, y::Float64)
    (isnan(x) | isnan(y)) && return NaN
    (x == Inf || y == Inf) && return Inf        # ∞ dominates
    x == -Inf && return y                       # log(0 + eʸ) = y, exactly
    y == -Inf && return x
    # yd: two faithful libm calls plus additions, condition ≤ 1 throughout
    # (result ≥ max(x, y)) ⇒ ≤ ~3 ulp ≈ 2⁻⁵⁰ — 2⁵ slack; measure it, including
    # x ≈ y (where log1p(exp(0)) dominates) and widely separated operands.
    _encl2(_logaddexp, x, y;
           fq=() -> _logaddexp(Float128(x), Float128(y)),
           yd=_logaddexp(x, y))
end
```

**Termination note**: can `logaddexp(x, y)` be exactly dyadic for dyadic
operands with the special rows peeled? `log(eˣ + eʸ) = d` requires
`eˣ + eʸ = e^d` — by Lindemann–Weierstrass this forces relations that dyadic
arguments cannot satisfy except on the peeled rows, so no further peel is
needed; **record that argument in a comment**. When you cannot prove such a
statement for your op, hunt for exact cases the way `_exact_cbrt` does.

**Step 3 — tests**: binary ops are not in the unary exhaustive gate; add a
divide-style exhaustive sweep (256×256 operand pairs against a ladder-only
reference, a handful of modes) to the suite, and an envelope-sanity loop entry
if you introduced a new libm dependence.

### Adding a ternary operation

Worked example: **`FMS`** — fused multiply-subtract, `x·y − z` — chosen because
the cheapest correct implementation is *delegation*, and delegation is a
first-class technique here.

**Step 1 — registry**: append `:FMS` to `_TERNARY_OPS` (group `:A` by the
loop). The generated array method routes through the same ternary bitwidth
policy every other ternary op uses (`tables.jl`'s `_ternary_table_for`): small
operand formats table (eagerly or adaptively), `K = 8` runs the — optionally
threaded — Shape-B scalar loop, and stochastic ρ always draws per element. No
op-specific work needed; the policy keys on formats and ρ, not on which op it is.

**Step 2 — ω-semantics** by delegation to `FMA`'s already-verified analysis
(width thresholds, `_twosum` exactness, wide-spread fallback):

```julia
# x·y − z ≡ FMA(x, y, −z). Negation of z follows the Subtract convention:
# preserve NaN, and map −0 to +0 so the single-zero datum model is respected.
ωeval(::Val{:FMS}, x::Float64, y::Float64, z::Float64) =
    ωeval(Val(:FMA), x, y, isnan(z) ? z : (iszero(z) ? 0.0 : -z))
```

Delegation inherits FMA's result kinds, its exactness thresholds, *and* its
test pedigree — but verify the sign algebra of your reduction on the special
rows (here: `FMS(∞, 1, ∞)` must be NaN, which the delegation produces because
`FMA(∞, 1, −∞)` is `∞ − ∞`). Add the delegated op to the ternary exhaustive
sweep in the suite.

### Adding a quaternary operation

Arity 4 is **not yet plumbed** — the generation loops, kernels, and block layer
handle arities 1–3. Adding a quaternary op is therefore two tasks: extend the
plumbing (once), then add the op. Worked example: **`FMMA`** — `x·y + z·w`, a
two-term dot product.

**Plumbing (one-time), four sites:**

1. `src/ops_scalar.jl` — register manually after the arity loops
   (`register_op!(:FMMA, 4, :A)`), and add an `op.arity == 4` branch to the
   spec-register generation loop, mirroring the ternary branch with a fourth
   operand:

   ```julia
   else # arity 4
       @eval begin
           @inline function $name(fr::Type{<:Binary}, ρ::ProjSpec,
                                  x::Binary, y::Binary, z::Binary, w::Binary;
                                  rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing)
               apply_op($V(), fr, ρ, _drawR(ρ, rng, R),
                        decode(x), decode(y), decode(z), decode(w))
           end
           @inline $name(x::T, y::T, z::T, w::T; kw...) where {T<:Binary} =
               $name(T, default_projspec(T), x, y, z, w; kw...)
       end
   end
   ```

   (`apply_op` and `ωeval` are already variadic — no change needed there.)
2. `src/kernels.jl` — a four-array `vmap!` method. The plain Shape-B loop (copy
   the *body* of the ternary method's compute branch, add operand `D`) is
   enough to get correct results; the ternary method's table-tiering and
   threading are a bitwidth-driven optimization layered on top (`tables.jl`'s
   `_ternary_table_for` and its LRU cache), not required plumbing — extend to
   arity 4 only if a table win at that arity is worth the 2^(4K) growth.
   Also add an arity-4 branch in the registry-generated array-surface loop and
   the matching rng-threading overload (mirroring the three-array one) so
   stochastic ρ dispatches correctly.
3. `src/blocks.jl` — either an arity-4 branch in the Block/Scaled generation
   loop (four operand blocks, four `blockdecode`s) or an explicit skip
   (`op.arity > 3 && continue`) if a block form is not wanted; be deliberate.
4. `src/approx.jl` — `conformance_report` prints arities with `for a in 1:3`;
   widen to `1:4`.

**The op itself** follows the `FAA`/`FMA` width-analysis playbook exactly:

```julia
const _DE_FMMA = 92     # two exact ≤17-bit products: 18 + ΔE ≤ 113, margin 3

function ωeval(::Val{:FMMA}, x::Float64, y::Float64, z::Float64, w::Float64)
    (isnan(x) | isnan(y) | isnan(z) | isnan(w)) && return NaN
    (((iszero(x) && isinf(y)) || (isinf(x) && iszero(y))) ||
     ((iszero(z) && isinf(w)) || (isinf(z) && iszero(w)))) && return NaN
    p = x * y; q = z * w                          # exact: ≤8-bit significands
    if isinf(p) || isinf(q)
        (isinf(p) && isinf(q) && p != q) && return NaN
        return isinf(p) ? p : q
    end
    s, e = _twosum(p, q)
    e == 0.0 && return iszero(s) ? 0.0 : s        # exact in Float64
    (_f128() && !iszero(p) && !iszero(q) && _expdiff(p, q) <= _DE_FMMA) &&
        return Float128(p) + Float128(q)          # exact by width
    BigExactF(() -> setprecision(() ->
        BigFloat(x) * BigFloat(y) + BigFloat(z) * BigFloat(w), BigFloat, _BIGP))
end
```

**Testing duty is heavier at arity 4**: the full cross-product is 2³² tuples,
so the suite entry must be *sampled* (seeded, ladder-referenced, ≥ 10⁵ tuples
across several modes) plus hand-picked exhaustion of the special-row algebra
and the width-threshold boundary (`ΔE ∈ {91, 92, 93}` constructions). Document
the sampling in the test, since it breaks the suite's "enumerate, never
sample" norm — that exception must be visible, not silent.

## Block operations

The block layer (`src/blocks.jl`) implements the draft's elementwise schema
**once** — `blockdecode → ω-op lanewise → blockproject` — and generates every
registry op's `BlockOp`/`ScaledOp` from it. The reductions are hand-written on
a second shared pattern. Which subsection you need depends on whether your
operation fits one of those two patterns.

### Adding a scaled block operation

If the scalar operation is (or can be) a **registry op, the scaled form is
free**: registering `:Cbrt` above already produced `ScaledCbrt(fr, ρ, s, x)`
and `BlockCbrt(fr, ρ, b, sr)` with the draft's §5.4 semantics — decode
`scale × element` exactly per lane, apply the ω-op, `BlockProject` against the
result scale through the division cascade. **Prefer this route**; it is the
package's non-divergence mechanism.

Write a scaled op by hand only when it is *not* an elementwise image of a
registry op. The primitive to compose with is `_bp_element(fr, ρ, R, res, Sdat)`
— it owns the draft's scale-special rows (NaN/0/±Inf scale) and the entire
division-by-scale rigor cascade. Example, a scaled affine `s₁·x + s₂·y`
projected against a *unit* result scale:

```julia
function ScaledAffine(fr::Type{<:Binary}, ρ::ProjSpec,
                      s1::Binary, x::Binary, s2::Binary, y::Binary;
                      rng::MaybeRNG=nothing)
    a = ωeval(Val(:Multiply), decode(s1), decode(x))::Float64   # exact lanes
    b = ωeval(Val(:Multiply), decode(s2), decode(y))::Float64
    _bp_element(fr, ρ, _drawR(ρ, rng, nothing), ωeval(Val(:Add), a, b), 1.0)
end
```

The `::Float64` assertions are load-bearing: lane decodes are exact by width
(≤ 17-bit significands), and the assertion turns that analysis into a runtime
check. The `Add` result may be any kind in the protocol — `_bp_element`
finishes all of them. Dividing by a scale of `1.0` is the identity fast path.

### Adding an elementwise block operation

Same schema, whole-block granularity. The three obligations: decode operand
blocks with `blockdecode` (exact), evaluate lanewise with `ωeval` (never with
ad-hoc Float64 math), and finish with `blockproject` against the supplied
result scale (never with per-lane `project` — the scale-division rows belong
to `blockproject`). Example — `BlockRelu`, an elementwise `max(x, 0)` that
reuses the verified `MaximumNumber` semantics:

```julia
"""BlockRelu(fr, ρ, b, sr): lanewise max(xᵢ, 0) — NaN lanes propagate NaN per
MaximumNumber's number-preferring rows only when both operands are NaN, i.e.
never here; state your lane semantics this explicitly for any new op."""
function BlockRelu(fr::Type{<:Binary}, ρ::ProjSpec, b::Block{B}, sr::Binary;
                   rng::MaybeRNG=nothing) where {B}
    X = blockdecode(b)
    Z = ntuple(i -> ωeval(Val(:MaximumNumber), X[i], 0.0), Val(B))
    blockproject(fr, ρ, sr, Z; rng)
end
```

Add the name to the block-surface export list and to `conformance()`'s
`blocknames` extension vector so the declaration stays truthful. Test against
a from-scratch reference composition (the suite's existing block tests show
the pattern) — element format × scale format × B, exhaustively where feasible.

### Adding a reductive block operation

Reductions do **not** fit the elementwise schema; they follow the second
pattern, visible in `BlockReduceAdd`/`BlockDotProduct`: (1) resolve the ∞/NaN
**fold algebra** on Float64 classifications first; (2) an integer **span
filter** decides whether the whole reduction is exactly representable in
`Float128` — if so, plain `Float128` accumulation *is* the exact answer;
(3) otherwise a `BigExactF` exact accumulator at provably sufficient
precision; (4) exactly **one** projection, at the very end, via `_finish`.
Never round intermediates.

Worked example — **`BlockSumOfSquares`**, `Σ xᵢ²`:

```julia
"""BlockSumOfSquares(fr, ρ, b): project(Σ decode-laneᵢ²) — one rounding total.
Fold algebra: any NaN lane ⇒ NaN; else any ±Inf lane ⇒ +Inf (squares cannot
cancel — no opposite-infinity NaN case exists, unlike ReduceAdd); else exact."""
function BlockSumOfSquares(fr::Type{<:Binary}, ρ::ProjSpec, b::Block{B};
                           rng::MaybeRNG=nothing, R::Union{Nothing,Int}=nothing) where {B}
    X = blockdecode(b)
    res = if any(isnan, X)
        NaN
    elseif any(isinf, X)
        Inf
    elseif all(iszero, X)
        0.0
    else
        # span filter: lanes carry ≤17-bit significands ⇒ squares ≤34 bits and
        # a square's exponent is 2·(lane exponent); the sum is exact in Float128
        # when 34 + 2·span + ⌈log₂B⌉ + 1 ≤ 113  ⇒  2·span + ⌈log₂B⌉ ≤ 78.
        if _f128() && 2 * _expspan(X) + _log2ceil(B) <= 78
            acc = Float128(0)
            for v in X
                q = Float128(v)
                acc += q * q                      # exact squares, exact sum
            end
            acc
        else
            BigExactF(() -> setprecision(BigFloat, _REDPREC) do
                acc = BigFloat(0)
                for v in X
                    acc += BigFloat(v)^2          # exact at _REDPREC
                end
                acc
            end)
        end
    end
    _finish(fr, ρ, _drawR(ρ, rng, R), res)
end
```

The two lines that demand your care in any new reduction are the **fold
algebra** (derive it from the draft's reduce definition on the extended reals —
here the absence of an opposite-infinities NaN case is a *theorem about
squares*, stated in the docstring) and the **span-filter inequality** (derive
it from operand widths, write the derivation in the comment, and add
boundary-condition tests at the threshold, exactly as the suite does for the
existing `_DE_*` constants). Everything else — the exact accumulator, the
single `_finish`, stochastic `R` plumbing — is pattern.

## After any addition: the closing checklist

- Full test suite green, including your new sweep entries.
- Envelope measurements recorded for any new `yd`/`fq` estimator.
- Benchmark report regenerated; `_SAFE_DOMAINS` entry added if the op is
  argument-restricted; sanity-check the op's row in the safe-args table.
- Docstring on the public name; `[interp]` markers where you interpreted.
- `conformance_report()` shows the op (automatic for registry ops; manual
  block additions extend `blocknames`).
