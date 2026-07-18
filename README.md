# ByteFloats.jl

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
registration. The test suite enumerates rather than samples: ≈ 8.8 million
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

Markdown sources in `docs/src/` (Introduction, User Guide, User Examples,
Technical Guide, Technical Examples, API Reference). Build the HTML site locally:

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

```
src/         eleven source layers (formats → … → approx) plus the vendored
             Float128 fma/faa soft-float modules; see the Technical Guide
test/        the consolidated exhaustive suite (runtests.jl)
docs/        Documenter site (make.jl, builddocs.jl, src/*.md)
benchmark/   Chairmarks suite + report generator (own environment)
```

## Notes

- Requires Julia 1.12. `Quadmath` is a hard dependency used only inside the
  oracle/fallback paths; `ENV["ByteFloats_Float128"] = "disable"` (before loading)
  selects the pure-MPFR configuration with bit-identical results.
- **License: MIT** (see `LICENSE`).
