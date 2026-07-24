# Regenerating and reading the benchmark report

<!-- updated by /doc-it -->
The runnable suite lives in `benchmark/` (its own environment, `benchmark/Project.toml`).
This directory holds the generated artifacts — `benchmark_report.md`/`.pdf`,
`domain_ops.csv`, `safe_domain_ops.csv` — and this guide. `benchmarking/benchmarking.jl`
and `benchmarking/simple_benchmarking.jl` are earlier copies kept for reference;
`benchmark/benchmarking.jl` is the current script.

## Regenerating

From your repository root, it's two commands.

One-time environment setup:

    julia --project=benchmark -e "using Pkg; Pkg.develop(path=\".\"); Pkg.instantiate()"

Then generate the report:

    julia --project=benchmark benchmark/benchmarking.jl benchmarking/benchmark_report.md

The trailing argument is the output path — omit it and the report lands at
benchmark_report.md in the current directory; pass benchmarking/benchmark_report.md
to match the default location (create the benchmarking/ folder first if it does
not exist — the script does not create directories).

From the REPL:

    include("benchmark/benchmarking.jl");
    generate_report("benchmarking/benchmark_report.md")

Expect roughly 8–14 minutes and ~280 measurement rows (add two more with
`julia -t N`, N > 1 — the threaded-vs-sequential K=8 ternary comparison only
runs with more than one Julia thread). The table-builds section deliberately
measures cold builds (up to ~90 ms each, many samples) plus warm cache hits,
and the scalar-operation tables appear in FOUR variants per arity
(unary 30, binary 18, ternary 3 — each measured four ways).

## How to read the scalar tables

Each scalar section (unary / binary / ternary) contains four tables, ordered
most-filtered first. **If you want one number per operation, read the first
table ("safe args") — it is the honest, fully unmasked cost.** Each later table
adds back one class of cheap fast-row dilution:

| order | title suffix       | operands sampled                                      | what it tells you |
|-------|--------------------|-------------------------------------------------------|-------------------|
| 1st   | "safe args"        | finite operands inside each op's safe domain          | the true per-op cost — no NaN dilution at all |
| 2nd   | "no NaN, Inf args" | finite operands only (zeros/subnormals kept)          | cost with ±Inf/NaN operands removed; domain-restricted ops still hit NaN rows on out-of-domain finite operands |
| 3rd   | "no NaN args"      | everything except the NaN code point (±Inf sampled)   | isolates the effect of the NaN operand alone |
| 4th   | *(no suffix)*      | ALL code points — NaN and ±Inf sampled                | the historical uniform-sweep view; medians of domain-restricted ops are heavily diluted by instant NaN rows |

**Safe-args operand generation.** Eleven argument-restricted operations draw
their arguments directly from explicit safe domains declared in the
`_SAFE_DOMAINS` map in benchmarking.jl:

    Sqrt        0 ≤ x < ∞          ArcSin    −1 ≤ x ≤ 1
    RSqrt       0 < x < ∞          ArcCos    −1 ≤ x ≤ 1
    Log, Log2   0 < x < ∞          ArcCosh    1 ≤ x < ∞
    LogOnePlus −1 < x < ∞          ArcTanh   −1 < x < 1
    Recip       x ≠ 0              Divide     y ≠ 0 (x unrestricted)

Every operation NOT in that map uses an oracle-derived pool instead: finite
operand tuples whose defined result is not NaN — so unlisted ops can never
drift out of sync with the implementation. To adjust a domain or add an op,
edit `_SAFE_DOMAINS` (one predicate per argument position); the report
structure follows automatically.

## The other sections

The non-scalar sampled tables (core primitives, format sensitivity, projection
modes, array kernels, sorting, blocks, conversions) all use the all-code-points
pool — NaN and ±Inf sampled — and each says so in its note. The table-builds
section enumerates every code point by construction (NaN and ±Inf included) and
reports both the cold build and the steady-state warm cache hit.

## Trusting the numbers

The script runs a specialization preflight first and aborts rather than publish
numbers if the warm scalar paths allocate — if that trips on your machine,
something is measuring dispatch, not arithmetic.

Numbers are machine-specific (my figures came from a single-threaded
sapphire-rapids container), so expect different absolutes but similar ratios.

Reading guidance for the middle tables: excluding NaN/±Inf OPERANDS does not
eliminate NaN RESULTS — domain-restricted ops still take instant NaN fast rows
on out-of-domain finite operands (e.g. a negative operand to Sqrt), so their
medians in the "no NaN, Inf args" and "no NaN args" tables remain partially
diluted; the report's notes state this. That residual dilution is precisely
what the leading "safe args" table removes, so cross-checking an op's median
across the four tables shows you exactly how much of its uniform-sweep median
was NaN fast rows versus real work.