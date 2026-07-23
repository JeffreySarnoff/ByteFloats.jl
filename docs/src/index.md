# ByteFloats.jl

*A conforming, performance-oriented Julia implementation of the IEEE P3109 draft
standard — arithmetic formats for machine learning at bitwidths 3–8.*

Bit-exact defined results on every default path; one projection engine as the single
write path into a code point; approximation only behind an explicitly named,
exhaustively measured κ registry; ≈ 8.9 million test assertions; table-gather kernels
at fractions of a nanosecond per element.

```julia-repl
julia> using ByteFloats

julia> Binary8p4se(1.6) + Binary8p4se(0.25)
Binary8p4se(1.875 ≡ 0x47)

julia> Binary5p3sf(0x08) == Binary5p3sf(1.0)     # UInt8 constructs from a code point
true
```

## Documentation map

- **[Introduction](@ref)** — what the package is, the design pillars, a
  thirty-second tour, installation.
- **[Cheat Sheet](@ref)** — one-page lookup for format names, conversion,
  projections, operations, arrays, blocks, packed storage, and common traps.
- **[User Guide](@ref)** — the complete public API in usage order: formats, values,
  projection specifications, the two operation registers, arrays and sorting,
  blocks, packed storage, conformance and κ, performance guidance.
- **[User Examples](@ref)** — runnable examples in three tiers: basic, machine
  learning, deep learning.
- **[Technical Guide](@ref)** — internals: the encoding and projection engine, the
  oracle's rigor classes, tables and kernels, the block layer's exactness filters,
  the verification and benchmark doctrines.
- **[Technical Examples](@ref)** — internals-level recipes: pipeline introspection,
  exhaustive verification of custom code, κ measurement, doctrine-compliant
  benchmarking.
- **[Adding Operations](@ref)** — how to add new scalar (unary through
  quaternary) and block (scaled, elementwise, reductive) operations with every
  guarantee intact: registry mechanics, ω-semantics duties, worked examples.
- **[API Reference](@ref)** — docstrings for the documented public surface.
