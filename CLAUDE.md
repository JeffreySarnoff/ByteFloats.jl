# CLAUDE.md

Guidance for AI coding agents working in this repository.

## What this package is

ByteFloats.jl implements the **IEEE P3109 draft standard** (*Arithmetic Formats
for Machine Learning*): all 120 binary formats at bitwidths 3–8, every
projection (rounding × saturation) mode including three stochastic families, the
full scalar operation catalog, block-scaled (MX-style) operations, sub-byte
packed storage, and the draft's conformance and κ-approximation machinery.

The implemented draft revision is the string in
[`src/approx.jl`](src/approx.jl) (`DRAFT_REVISION`, exported as
`draft_revision()`).

Requires Julia 1.12. Runtime deps: `PrecompileTools`, `Quadmath`, `Random`.

## Commands

```
julia --project=. -e 'using Pkg; Pkg.test()'      # exhaustive suite (~8.9M assertions)
julia docs/builddocs.jl                            # build docs site locally
julia --project=docs docs/make.jl                  # docs build with env already set up
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark benchmark/benchmarking.jl benchmarking/benchmark_report.md
```

`ENV["ByteFloats_Float128"] = "disable"` (set **before** loading) selects the
pure-MPFR configuration. It must produce **bit-identical results** — the switch
trades speed only. CI runs both configurations
([`.github/workflows/CI.yml`](.github/workflows/CI.yml)).

## Architecture

Twelve source layers, included in a fixed order by
[`src/ByteFloats.jl`](src/ByteFloats.jl):

```
formats → projspec → defaults → decode_encode → project → ops_scalar
        → oracle → tables → kernels → blocks → packed → approx
```

`ops_scalar` deliberately precedes `oracle` (the evaluation-protocol structs
`BigExactF`/`EncloseF` live in `ops_scalar.jl`). **Keep the include-order comment
in `src/ByteFloats.jl` in sync when adding a layer.**

| layer | role |
|---|---|
| `formats.jl` | `Binary{K,P,SGN,EXT}` type, Group M traits, format names, Base API |
| `projspec.jl` | rounding × saturation as zero-size *types* so kernels specialize on ρ |
| `defaults.jl` | settable session defaults (`DefaultType`, `DefaultProjection`, …); hot paths never consult them |
| `decode_encode.jl` | ωDecode / ωEncode, ordering keys, `Next*` ops |
| `project.jl` | the projection engine: ωRoundToPrecision → ωSaturate → ωEncode |
| `ops_scalar.jl` | `OP_REGISTRY`, `apply_op`, spec + Base register veneers |
| `oracle.jl` | ω-semantics catalog: rigorous defined results. Never on the hot path |
| `tables.jl` | table lifecycle + cache for pure-ρ specializations |
| `kernels.jl` | array kernels: Shape A (table gather), Shape B (per-element compute) |
| `blocks.jl` | blocks, scaled operations, reductions (draft §5) |
| `packed.jl` | `PackedVector`: sub-byte packed storage, unpack → compute → repack |
| `approx.jl` | conformance declaration + the κ-approximation registry |

Plus two **vendored** soft-float modules, `src/fma128.jl` and `src/faa128.jl`,
providing correctly-rounded `fma`/`faa` for `Float128` (needed on Windows, where
Quadmath does not bind `fmaq`). Upstream sources live in
[Float128FMA/](Float128FMA/) and [Float128FAA/](Float128FAA/) as standalone
packages. **The vendored copies have diverged from the standalone sources** —
they are not byte-identical; do not assume a fix in one propagates to the other.

## Invariants — do not break these

1. **One write path.** `project.jl`'s engine is the *only* way a code point is
   created. Never encode a result by hand in an operation, kernel, or block path.
2. **Code point vs value.** `UInt8` is the one argument type meaning *code
   point*; every other `Real` means *value*. `T(0x08)` and `T(8.0)` are different
   things by design.
3. **Representation invariant.** The code point occupies the low `K` bits of the
   payload byte; the high `8−K` bits are maintained zero.
4. **Stochastic ρ is never tabulable.** The result is a distribution over R.
   Table builders reject stochastic specs loudly — keep it that way.
5. **Nothing approximate is reachable from the default API.** Approximate
   kernels live only in the κ registry, retrieved by explicit name. κ is
   *measured exhaustively* at `register_approx!` time; understated declarations
   are rejected. A declared bound is a verified property, not a promise.
6. **A table entry IS the defined result.** Builders walk every input through the
   oracle-backed scalar path, so use sites carry no residual correctness burden.
7. **Registry-driven codegen.** `OP_REGISTRY` in `ops_scalar.jl` is the single
   source generating the spec-register functions, same-format convenience
   methods, Base veneers, `Block*`/`Scaled*` variants, table enumeration,
   exhaustive test sets, exports, and `conformance()`. Adding an operation means
   adding a registry row — see [docs/src/new_operations.md](docs/src/new_operations.md),
   which is the authoritative procedure. Never hand-write a per-op variant that
   the registry should generate; that is the non-divergence mechanism.

## Oracle result protocol

`ωeval` returns one of five kinds, finished by `apply_op`/`_finish`:

| kind | meaning |
|---|---|
| `Float64` | exact: draft-tabled specials, exactly representable finites |
| `Float128` | exact by width analysis (significand fits 113 bits) |
| `BigExactF` | `f()` → exact `BigFloat` at 2200 bits (wide-spread tail) |
| `Enclose128F` | correctly-rounded `Float128` bracket, MPFR closure as fallback |
| `EncloseF` | `f(prec)` → MPFR directed `(lo, hi)`; optional `fq` Float128 pre-filter and `yd` eager Float64 estimate |

`sticky ∈ {-1,0,+1}` carries symbolic direction information (true value = carried
value + sticky·ε) so a correctly-rounded endpoint can stand in for an irrational
truth without losing what directed/tie/stochastic rounding needs.

## Performance rules

- **Pass format types through `const` bindings, type parameters, or function
  arguments.** A non-`const` global format type forces dynamic dispatch on every
  call (~1 µs per scalar keyword call, vs ≈ 26 ns for a specialized `Add`). A
  single function barrier `f(::Type{T}, args...) where {T}` restores full speed.
- **Zero warm-path allocation is a tested regression**, not an aspiration. The
  suite pins concrete inferred return types at public entry points and zero
  allocation on warm paths. The benchmark script runs a specialization preflight
  and *aborts* rather than publish numbers if warm scalar paths allocate.
- Table getters are `@noinline` and called **once** per array call, hoisted out
  of the loop. No dict lookups, locks, or global loads per element.
- Benchmark doctrine: a closure over any non-`const` global measures Julia's
  dispatch machinery, not the code under test. See the "Benchmark doctrine"
  section of [docs/src/technical_guide.md](docs/src/technical_guide.md).

## Testing

`test/runtests.jl` is the shipped suite — seven consolidated harnesses, run by
`Pkg.test()`. It **enumerates rather than samples**; the value sets are small
enough that sampling is never necessary. Preserve that property in new tests.

`test/quality.jl` (included from `runtests.jl`) is the hygiene gate:

- `Aqua.test_all` — all checks on, no exclusions. Notably, the
  **unbound-type-parameter** check is why `Block`'s inner constructor takes
  `Tuple{FE,Vararg{FE,Bm1}}` rather than `NTuple{B,FE}`: the empty tuple would
  leave `FE` unbound. Keep new signatures free of that pattern.
- `JET.report_package` — whole-package analysis. One report is filtered by name
  (`_vmap_packed`), because JET widens the correlation between a `::Type{fr}`
  argument and the element type of a container built from it. If you add a
  filter, it must come with a concrete-call gate proving the path is clean.
- `JET.@test_call` on concrete entry points — the analysis that matches the
  specialization doctrine. Add a line here for every new public entry point.

The other files in `test/` are not run by `Pkg.test()` and are included/run
manually:

- `refimpl.jl` — `Rational{BigInt}` reference for ωRoundToPrecision/ωSaturate/
  ωEncode, sharing no code with the engine under test
- `ternary_opt.jl` — gates for the ternary table tiers and the sticky-head
  wide-spread FMA/FAA escalation
- `validate_correctness.jl` — envelope measurement, differential ladder
  comparison, and adversarial edges for the `yd`/fast-path refactors

## Repository layout

```
src/           the twelve layers + the two vendored Float128 modules
test/          runtests.jl (shipped) + three manual harnesses
docs/          Documenter site (make.jl, builddocs.jl, src/*.md) and docs/pdf/
benchmark/     current Chairmarks suite + report generator (own environment)
benchmarking/  generated report, domain CSVs, report-reading guide
Float128FAA/   standalone faa(x,y,z) package
Float128FMA/   standalone Float128 fma package
```

## Known repository gotchas

- **`benchmark/` vs `benchmarking/`.** `benchmark/benchmarking.jl` is the current
  script; `benchmarking/benchmarking.jl` and `simple_benchmarking.jl` are older
  copies kept alongside the generated artifacts. Edit the one in `benchmark/`.
- **Docs deploy from CI only.** [docs/make.jl](docs/make.jl) calls `deploydocs`
  targeting the `gh-pages` branch. Local builds log "could not auto-detect the
  building environment. Skipping deployment" — that warning is expected and
  correct, not a failure.
- **Dangling design-document references.** Source headers cite `design §N`,
  `architecture §N`, `Float128 revision plan §N`, and `checkpoint.md`. None of
  those documents are in this repository. Treat the citations as historical
  markers; do not invent their contents, and do not "fix" a reference by
  guessing what the section said.
