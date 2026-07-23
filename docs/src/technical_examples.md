# Technical Examples

Internals-level recipes: introspecting the pipeline, verifying custom code against
the oracle, measuring κ, and benchmarking without fooling yourself. Outputs are
captured from real sessions. Several examples use unexported internals
(`ByteFloats.round_to_precision`, `ByteFloats.project`, …); those are stable enough to learn
from but are not covered by semantic-versioning guarantees.

## Basic

### Watching RoundToPrecision work

`round_to_precision` returns the draft's exact `(sign, S, Q)` decomposition before
saturation ever sees it:

```julia
using ByteFloats
using ByteFloats: round_to_precision

r = round_to_precision(4, 8, NearestTiesToEven(), 2.30078125, 0, 0)   # P=4, bias=8
(r.sign, r.S, r.Q, r.S * 2.0^r.Q)
```

```
(1, 9, -2, 2.25)          # S·2^Q = 9/4: the nearest P=4 value
```

### Symbolic sticky: projecting "just below" a number

The engine accepts `sticky ∈ {-1, 0, +1}` meaning "the true value is the carrier
plus an infinitesimal of this sign". This is how enclosure endpoints and asymptotes
(`tanh → 1⁻`) are projected exactly:

```julia
using ByteFloats: project
project(Binary8p4se, ProjSpec(TowardNegative(), SatNone()), 1.0; sticky=-1)
```

```
Binary8p4se(0.9375 ≡ 0x3f)     # == NextLessThan(1.0): the engine crossed the binade
```

### Tables are the scalar path, memoized

```julia
using ByteFloats: get_table
empty_tables!()
tbl = get_table(:Exp, Binary8p4se, Binary8p4se, RNE_SatNone)
(tbl[Int(0x45) + 1],
 codepoint(Exp(Binary8p4se, RNE_SatNone, rawvalue(Binary8p4se, 0x45))),
 table_bytes())
```

```
(0x52, 0x52, 256)              # table entry ≡ scalar result; one 256-byte table cached
```

### Order keys

Ordering is integer arithmetic: a sign-magnitude fold, monotone with the total
order, NaN at the top. This is what comparisons and the O(n) counting sort run on:

```julia
using ByteFloats: order_key
[(v, order_key(v)) for v in (Binary8p4se(-1.0), Binary8p4se(0.0),
                             Binary8p4se(1.0), Binary8p4se(NaN))]
```

```
-1.0 → 64      0.0 → 129      1.0 → 193      NaN → 65535
```

## General AI

### Ranking safety: an exhaustive monotonicity audit of score conversion

Decision and ranking systems consume scores through *order*. If scores are
produced in one format and compared in another, the conversion must be monotone
— and formats are small enough to prove it by enumeration rather than trust it.
Walk every `Binary8p3se` datum in total order and check that its `Binary8p4se`
image never goes backward on the integer order keys:

```julia
using ByteFloats
using ByteFloats: order_key

function monotone_conversion(::Type{From}, ::Type{To}) where {From,To}
    prev = nothing
    for v in sort(From.(0x00:UInt8(2^bitwidth(From) - 1)))
        isnan(decode(v)) && continue
        g = Convert(To, RNE_SatNone, v)
        prev !== nothing && order_key(g) < order_key(prev) && return false
        prev = g
    end
    true
end
monotone_conversion(Binary8p3se, Binary8p4se)
```

```
true
```

Run this for every `(From, To, ρ)` triple your system actually uses; it is a few
hundred integer comparisons per triple. A ranking pipeline whose conversions all
pass this audit cannot invert a preference by changing formats — a guarantee no
amount of spot-testing provides.

### κ-safe decision margins for approximate evaluators

Search under a compute budget often wants a cheap, approximate evaluation
function. The κ registry turns "how approximate?" into a *measured* code-point
bound — and code-point bounds compose into a decision rule: **if two defined
evaluations differ by more than 2κ code points along the total order, the
approximate evaluator cannot invert their comparison.** Verify the rule
exhaustively for a κ = 2 evaluator:

```julia
using ByteFloats: order_key

fast(x) = (r = Exp(Binary8p4se, RNE_SatNone, x);
           isfinite(decode(r)) ? NextGreaterThan(NextGreaterThan(r)) : r)
κ, exhaustive = measure_kappa(fast, :Exp, Binary8p4se, (Binary8p4se,), RNE_SatNone)

function margin_audit(fast, κ)
    codes = [rawvalue(Binary8p4se, UInt8(c)) for c in 0:255]
    safe = violations = 0
    for a in codes, b in codes
        da, db = Exp(Binary8p4se, RNE_SatNone, a), Exp(Binary8p4se, RNE_SatNone, b)
        (isnan(decode(da)) || isnan(decode(db))) && continue
        codedistance(da, db) > 2κ || continue
        safe += 1
        (order_key(fast(a)) < order_key(fast(b))) ==
            (order_key(da) < order_key(db)) || (violations += 1)
    end
    (safe, violations)
end
(κ, exhaustive, margin_audit(fast, κ)...)
```

```
(2.0, true, 49522, 0)
```

49,522 operand pairs clear the 2κ margin, and the approximate evaluator agrees
with the defined ordering on every one of them — zero violations, exhaustively.
Inside the margin, comparisons are genuinely undecidable at this κ; a search can
treat sub-margin comparisons as ties to expand, or escalate those few nodes to
the exact evaluator. Either way the pruning is *provably* sound, with κ measured
at registration rather than promised.

## Machine Learning

### Verifying a quantizer exhaustively against an independent reference

Formats are small enough that "test on a few points" is never necessary. Here every
in-range datum of `Binary8p3se` is projected into `Binary8p4se` and checked against a
brute-force nearest search under 256-bit arithmetic (distance agreement; the draft's
tie rule is the projection engine's job and is pinned elsewhere in the suite):

```julia
using ByteFloats

function ref_nearest_distance(::Type{T}, x) where {T}
    fins = [v for v in (rawvalue(T, UInt8(c)) for c in 0:255) if isfinite(decode(v))]
    minimum(abs(setprecision(() -> BigFloat(decode(v)) - BigFloat(x), BigFloat, 256))
            for v in fins)
end

function verify_quantizer()
    ok = true
    for c in 0x00:0xff
        x = decode(rawvalue(Binary8p3se, c))
        (isfinite(x) && abs(x) <= decode(MaxFiniteOf(Binary8p4se))) || continue
        got = Convert(Binary8p4se, RNE_SatNone, x)
        ok &= abs(decode(got) - x) == Float64(ref_nearest_distance(Binary8p4se, x))
    end
    ok
end
verify_quantizer()
```

```
true
```

This pattern — enumerate the inputs, compare against an independently written
reference — is how the entire package is tested, and it is available to *your*
quantization code at trivial cost.

### The κ workflow, including the part where it says no

Register an approximate kernel and the registry measures it; understate the bound
and it refuses:

```julia
step2(x) = (r = Exp(Binary8p4se, RNE_SatNone, x);
            isfinite(decode(r)) ? NextGreaterThan(NextGreaterThan(r)) : r)

measure_kappa(step2, :Exp, Binary8p4se, (Binary8p4se,), RNE_SatNone)
```

```
(2.0, true)                # κ = 2 code points, verified over all 256 inputs
```

```julia
register_approx!(:cheater, :Exp, Binary8p4se, (Binary8p4se,), RNE_SatNone, step2; κ=1)
```

```
ERROR: ArgumentError: declared κ = 1.0 understates measured κ = 2.0 — registration rejected
```

### Exporting the conformance declaration

```julia
d = conformance_dict()          # plain nested Dict{String,Any}
(d["package"], length(d["formats"]), length(d["operations"]),
 [s["op"] for s in d["cached_specializations"]])
```

Serialize `d` with any JSON/TOML writer to attach a machine-readable conformance
statement to experiment artifacts.

## Deep Learning

### Proving your fused kernel exact: BlockDotProduct vs 512-bit truth

The block layer promises *one* projection with exact lane arithmetic. Trust, then
verify — against big-float truth, over random blocks with mixed scales and an honest
special-value mix:

```julia
using ByteFloats, Random

rng = Xoshiro(5)
mk() = Block(Binary8p1uf(2.0^rand(rng, -3:3)),
             ntuple(_ -> rawvalue(Binary8p4se, UInt8(rand(rng, 0:255))), 32))

function verify_dots(trials)
    for _ in 1:trials
        bx, by = mk(), mk()
        lx = [decode(bx.s) * decode(v) for v in bx.x]
        ly = [decode(by.s) * decode(v) for v in by.x]
        (any(!isfinite, lx) || any(!isfinite, ly)) && continue
        truth = setprecision(() -> sum(BigFloat(lx[i]) * BigFloat(ly[i]) for i in 1:32),
                             BigFloat, 512)
        got = BlockDotProduct(Binary8p4se, RNE_SatNone, bx, by)
        ref = setprecision(() -> ByteFloats.project(Binary8p4se, RNE_SatNone, truth), BigFloat, 512)
        codepoint(got) == codepoint(ref) || return false
    end
    true
end
verify_dots(200)
```

```
true
```

The same shape verifies any custom fused kernel you write: compose the truth in
`BigFloat`, project once, compare code points.

### Stochastic rounding, audited: the full-R sweep

Every stochastic projection is a deterministic function of the draw `R`, so its
*distribution* is checkable exactly. For `x = 2 + 3/64` in `Binary8p4se` the fraction
is ν = 3/16 of an ulp, so `StochasticA{4}` must round up for exactly 3 of the 16
draws:

```julia
σ4 = RSA_SatNone(4)                  # ≡ ProjSpec(StochasticA{4}(), SatNone())
x = 2.0 + 3/64
count(decode(ByteFloats.project(Binary8p4se, σ4, x; R)) == 2.25 for R in 0:15)
```

```
3
```

Sweeping `R` like this turns "is my stochastic pipeline unbiased?" from a statistical
question into an exhaustive one — the pattern the shipped test suite uses.

### An accelerator-style activation, κ-measured: hard-tanh vs Tanh

Hardware activation units often ship piecewise approximations. The κ machinery
makes the substitution honest: measure the approximation's worst code-point
deviation from the defined activation, exhaustively, and register it under that
measured bound. Hard-tanh — `clamp(x, −1, 1)` — as a stand-in for `Tanh`:

```julia
one4 = Binary8p4se(1.0)
hardtanh(x) = Clamp(Binary8p4se, RNE_SatNone, x, Negate(one4), one4)

κ, exhaustive = measure_kappa(hardtanh, :Tanh, Binary8p4se, (Binary8p4se,), RNE_SatNone)
(κ, exhaustive)
```

```
(4.0, true)                # worst deviation: 4 code points, verified on all 256 inputs
```

Where is it worst? At `x = 1.0`: hard-tanh returns 1.0 while the defined
`tanh(1.0)` is 0.75 — four code points along the total order. The registration
round-trip, including the conformance surface:

```julia
impl = register_approx!(:hardtanh_act, :Tanh, Binary8p4se, (Binary8p4se,),
                        RNE_SatNone, hardtanh; κ=4)
(kappa(:hardtanh_act), kappa_measured(impl), :hardtanh_act in list_approx())
```

```
(4.0, 4.0, true)           # declared = measured; conformance_report() now lists it
```

(`unregister_approx!(:hardtanh_act)` removes it.) Declaring `κ=3` would be
*rejected* — understatement is impossible by construction. The decision of
whether a κ = 4 activation is acceptable now belongs to your accuracy budget,
not to hope: combined with the margin rule from the General AI section, any two
pre-activations whose defined `Tanh` images sit more than 8 code points apart
keep their order under hard-tanh, provably.

### Benchmarking without measuring the dispatcher

The package's benchmark doctrine, in one snippet (needs the `benchmark/` environment
for Chairmarks). Format types enter as **type parameters**; operands come from
Chairmarks' *untimed* `setup`; and if you retrieve functions reflectively
(`getfield(ByteFloats, op)`), pass them through an argument barrier so they specialize:

```julia
using Chairmarks, ByteFloats, Random
using Statistics: median

function bench_add(::Type{T}) where {T<:Binary}          # T: type parameter, not a global
    pool = [rawvalue(T, rand(UInt8)) for _ in 1:4096]
    @be (rand(pool), rand(pool)) (t -> Add(T, RNE_SatNone, t[1], t[2]))(_)
end

b = bench_add(Binary8p4se)
(round(median(b).time * 1e9; digits=1), median(b).allocs)
```

```
(16.3, 0.0)          # ns per full scalar Add, zero allocations
```

!!! warning "The 60× trap"
    The same call with `T` read from a non-`const` global measures ~1 µs — Julia's
    dynamic keyword dispatch, not this package. Two project post-mortems trace to
    exactly this mistake; the shipped `benchmark/benchmarking.jl` asserts
    specialization (zero warm-path allocation) before it believes any number, and
    so should yours.
