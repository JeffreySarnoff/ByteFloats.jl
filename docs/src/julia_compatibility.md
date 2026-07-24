# Julia Compatibility

`Binary` values speak ordinary Julia. The **Base register**
(`src/juliacompat.jl`) maps every draft operation that has a Base counterpart
onto its Base name, so generic Julia code — and your fingers — can use the
spellings they already know. Every method there is exactly one same-format
spec-register call under the session default projection; there is no third
semantics: `x + y` *is* `Add(x, y)`, which *is*
`Add(T, DefaultProjection(), x, y)`.

## Using the overloads

```julia-repl
julia> x, y, z = Binary8p4se(1.6), Binary8p4se(0.25), Binary8p4se(2.0);

julia> x + y                       # Add
Binary8p4se(1.875 ≡ 0x47)

julia> exp(y)                      # Exp — a 256-byte table lookup once warm
Binary8p4se(1.25 ≡ 0x42)

julia> atan(y, x)                  # ArcTan2, in Base's (y, x) argument order
Binary8p4se(0.15625 ≡ 0x2a)

julia> fma(x, y, z)                # FMA — one rounding; muladd is the same op
Binary8p4se(2.5 ≡ 0x4a)

julia> sincos(y)                   # composite: (Sin(y), Cos(y)) componentwise
(Binary8p4se(0.25 ≡ 0x30), Binary8p4se(1.0 ≡ 0x40))
```

Because the veneers follow the session default projection, changing it changes
them (see the session-defaults section of the [User Guide](@ref)):

```julia-repl
julia> Binary8p4se(200.0) * Binary8p4se(2.0)     # RNE_SatNone: overflow → Inf
Binary8p4se(Inf ≡ 0x7f)

julia> DefaultProjection!(RTZ_SatFinite);

julia> Binary8p4se(200.0) * Binary8p4se(2.0)     # now clamps to MaxFinite
Binary8p4se(224.0 ≡ 0x7e)
```

Code that must be insensitive to the session default names its projection:
`Multiply(T, RNE_SatNone, x, y)`.

## The mapping

| Base spelling | draft operation |
|---|---|
| `+` `-` `*` `/`, unary `-` | `Add`, `Subtract`, `Multiply`, `Divide`, `Negate` |
| `abs`, `inv`, `sqrt` | `Abs`, `Recip`, `Sqrt` |
| `exp`, `exp2`, `expm1`, `log`, `log2`, `log1p` | `Exp`, `Exp2`, `ExpMinusOne`, `Log`, `Log2`, `LogOnePlus` |
| `sin`, `cos`, `tan`, `asin`, `acos`, `atan` | `Sin`, `Cos`, `Tan`, `ArcSin`, `ArcCos`, `ArcTan` |
| `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh` | `Sinh`, `Cosh`, `Tanh`, `ArcSinh`, `ArcCosh`, `ArcTanh` |
| `sinpi`, `cospi`, `tanpi` | `SinPi`, `CosPi`, `TanPi` |
| `atan(y, x)` | `ArcTan2` (Base's argument order) |
| `copysign`, `hypot` | `CopySign`, `Hypot` |
| `min`, `max` | `Minimum`, `Maximum` (NaN-propagating, exactly Base's float semantics) |
| `fma`, `muladd` | `FMA` (both: one rounding) |
| `clamp` | `Clamp` |
| `sincos`, `sincospi`, `minmax` | componentwise composites of the draft ops |

The mapping is a declarative partition over the op lists, and the test suite
asserts it is exhaustive: every operation is either mapped above or listed in
`_NO_BASE_COUNTERPART`. Adding an operation to the registry forces an explicit
decision here.

## What is *not* mapped, and why

**Draft operations with no Base spelling** — call them by their draft names:
`RSqrt`, `Softplus`, `ArcSinPi`/`ArcCosPi`/`ArcTanPi`/`ArcTan2Pi` (Base has
only the `sinpi` family), the NaN-ignoring/magnitude/finite extremum families
(`MinimumNumber`, `MaximumMagnitude`, `MinimumFinite`, …; Base has no
NaN-ignoring pair), and `FAA` (no Base fused add-add).

```julia-repl
julia> min(Binary8p4se(NaN), x)          # Base semantics: NaN propagates
Binary8p4se(NaN ≡ 0x80)

julia> MinimumNumber(Binary8p4se(NaN), x)   # the NaN-ignoring draft op
Binary8p4se(1.625 ≡ 0x45)
```

**Mixed `Binary` formats** — deliberately a `promotion ... failed to change`
error, never a silent widening. Mixing formats is an explicit `Convert`:

```julia-repl
julia> x + Convert(Binary8p4se, RNE_SatNone, Binary5p3sf(1.0))
Binary8p4se(2.5 ≡ 0x4a)
```

**`Binary` with ordinary numbers** — promotes to `Float64` (the exact carrier)
through the rules in `formats.jl`, so the result is a `Float64`, not a
re-projected `Binary`:

```julia-repl
julia> x + 2.0
3.625
```

Project it back explicitly when that is what you mean:
`Convert(Binary8p4se, RNE_SatNone, decode(x) + 2.0)`.

## The rest of the Base surface

Beyond arithmetic, `Binary` integrates where the other layers provide it:
predicates and constants (`isnan`, `isfinite`, `signbit`, `zero`, `one`,
`eps`, `typemin`/`typemax`, `floatmin`/`floatmax`), comparisons and `isless`
on integer order keys, `sort` via an O(n) counting sort,
`nextfloat`/`prevfloat` as the draft Next operations, `codepoint`, and
`rand`/`randn` (see the [User Guide](@ref)). All of it goes through the same
projection engine and decode tables as the draft-named API.
