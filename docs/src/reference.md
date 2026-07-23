# API Reference

Docstrings for the documented surface of the package. The full export list is far
larger (393 names, including the 120 format aliases, the predefined projection-spec
grid, and the generated operation and block registers); names without docstrings are
covered descriptively in the [User Guide](@ref) and [Technical Guide](@ref), and
completing per-name docstring coverage is tracked work.

## Public API

```@autodocs
Modules = [ByteFloats]
Private = false
```

## Documented internals

Unexported machinery documented in-source and referenced by the
[Technical Guide](@ref) and [Technical Examples](@ref).

```@autodocs
Modules = [ByteFloats]
Public = false
```
