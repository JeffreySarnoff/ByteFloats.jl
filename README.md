# ByteFloats.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://JeffreySarnoff.github.io/ByteFloats.jl/dev/)
[![CI](https://github.com/JeffreySarnoff/ByteFloats.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JeffreySarnoff/ByteFloats.jl/actions/workflows/CI.yml)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

A conforming, performance-oriented Julia implementation of the **IEEE P3109 draft
standard** — *Arithmetic Formats for Machine Learning* — covering all 120 binary
formats at bitwidths 3–8, every projection (rounding × saturation) mode including
three stochastic families, the full scalar operation catalog, block-scaled
(MX-style) operations, sub-byte packed storage, and the draft's conformance and
κ-approximation machinery.

**Bit-exact defined results on every default path.** One projection engine is the
single write path into a code point; results are established against a rigorous
oracle (exact arithmetic, IEEE-correctly-rounded `Float128` where mandated, MPFR
directed enclosures with precision escalation elsewhere); approximation exists only
behind an explicit registry whose deviation bounds are *measured exhaustively* at
registration. The test suite enumerates rather than samples: ≈ 8.9 million
assertions.

```julia
using ByteFloats

Binary8p4se(1.6) + Binary8p4se(0.25)        # Binary8p4se(1.875 ≡ 0x47)
Binary5p3sf(0x08) == Binary5p3sf(1.0)       # true — UInt8 constructs from a code point

σ = ProjSpec(StochasticA{8}(), SatNone())
Add(Binary8p4se, σ, Binary8p4se(2.0), Binary8p4se(0.03125); R = 255)

Exp(Binary8p4se, RNE_SatNone, Binary8p4se.(randn(1000)))   # table-gather kernel
```

## Documentation

<!-- updated by /doc-it -->
**[Read the docs](https://JeffreySarnoff.github.io/ByteFloats.jl/dev/)** —
Markdown sources in `docs/src/` (Introduction, Cheat Sheet, User Guide, User
Examples, Technical Guide, Technical Examples, Adding Operations, API
Reference). Build the HTML site locally:

```
julia docs/builddocs.jl
```

## Tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Benchmarks

One-time setup, then generate a machine-specific report:

```
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark benchmark/benchmarking.jl benchmark_report.md
```

## Repository layout

<!-- updated by /doc-it -->
```
src/           fourteen source layers (formats → … → rand) plus the vendored
               Float128 fma/faa soft-float modules; see the Technical Guide
test/          the consolidated exhaustive suite (runtests.jl) plus the
               independent validation harness (validate_correctness.jl,
               refimpl.jl, ternary_opt.jl)
docs/          Documenter site (make.jl, builddocs.jl, src/*.md) and docs/pdf/
benchmark/     Chairmarks suite + report generator (own environment) — current
benchmarking/  generated report (benchmark_report.md/.pdf), operation-domain
               CSVs, and the report-reading guide (README.md)
Float128FAA/   standalone `faa(x,y,z)` package, vendored as src/faa128.jl
Float128FMA/   standalone `fma` package for Float128, vendored as src/fma128.jl
```

## Notes

- Requires Julia 1.12. `Quadmath` is a hard dependency used only inside the
  oracle/fallback paths; `ENV["ByteFloats_Float128"] = "disable"` (before loading)
  selects the pure-MPFR configuration with bit-identical results.
- **License: MIT** (see `LICENSE`).
