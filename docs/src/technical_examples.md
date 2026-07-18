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
σ4 = ProjSpec(StochasticA{4}(), SatNone())
x = 2.0 + 3/64
count(decode(ByteFloats.project(Binary8p4se, σ4, x; R)) == 2.25 for R in 0:15)
```

```
3
```

Sweeping `R` like this turns "is my stochastic pipeline unbiased?" from a statistical
question into an exhaustive one — the pattern the shipped test suite uses.

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
