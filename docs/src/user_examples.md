# User Examples

Worked, runnable examples. Every output shown was captured from a real session; seeds
are fixed so you can reproduce them exactly.

## Basic

### One value, every rounding mode

The value 2.30078125 sits between the `Binary8p4se` grid points 2.25 and 2.5
(ulp = 0.25 in [2, 4)):

```julia
using ByteFloats

for μ in (NearestTiesToEven(), NearestTiesToAway(), TowardPositive(),
          TowardNegative(), TowardZero(), ToOdd())
    v = Convert(Binary8p4se, ProjSpec(μ, SatNone()), 2.30078125)
    println(rpad(string(typeof(μ)), 20), " → ", decode(v))
end
```

```
NearestTiesToEven    → 2.25
NearestTiesToAway    → 2.25
TowardPositive       → 2.5
TowardNegative       → 2.25
TowardZero           → 2.25
ToOdd                → 2.25
```

(`ToOdd` keeps 2.25 because its significand is already odd; try 2.05 to see it move.)

### Enumerating a whole format

Small formats are small enough to look at in full — 16 code points for `Binary4p2se`,
sorted by the total order (single NaN last):

```julia
decode.(sort(Binary4p2se.(0x00:0x0f)))      # broadcast the code-point constructor
```

```
[-Inf, -2.0, -1.5, -1.0, -0.75, -0.5, -0.25, 0.0,
  0.25, 0.5, 0.75, 1.0, 1.5, 2.0, Inf, NaN]
```

This is a good way to build intuition for a format before committing to it: you can
*see* the subnormal spacing, the binade structure, and the dynamic range.

### Saturation in one line each

```julia-repl
julia> w, two = Binary8p4se(200.0), Binary8p4se(2.0);   # MaxFinite is 224

julia> Multiply(Binary8p4se, RNE_SatNone, w, two)
Binary8p4se(Inf ≡ 0x7f)

julia> Multiply(Binary8p4se, RNE_SatFinite, w, two)
Binary8p4se(224.0 ≡ 0x7e)
```

## Machine Learning

### Quantizing a weight tensor and measuring the damage

```julia
using ByteFloats, Random, Statistics

rng = Xoshiro(42)
w = randn(rng, 10_000) .* 0.25          # typical trained-weight scale
q = Binary8p4se.(w)                     # project every weight
back = decode.(q)

mean((w .- back).^2), maximum(abs.(w .- back)), count(isinf, back) / length(back)
```

```
(4.41e-5, 0.0312, 0.0)                  # MSE, max |error|, overflow fraction
```

The max error is exactly half an ulp of the largest binade the data reached — the
worst case nearest rounding permits.

### Picking a format: precision vs range

For the same Gaussian tensor, `K = 8` formats trade significand bits against exponent
range. RMSE doubles per lost significand bit; MaxFinite grows explosively:

```julia
rng = Xoshiro(7); x = randn(rng, 50_000)
for F in (Binary8p5se, Binary8p4se, Binary8p3se, Binary8p2se)
    back = [decode(Convert(F, RNE_SatFinite, xi)) for xi in x]
    println(rpad(formatname(F), 12), " rmse = ", round(sqrt(mean((x .- back).^2)); sigdigits=3),
            "   maxfinite = ", decode(MaxFiniteOf(F)))
end
```

```
Binary8p5se  rmse = 0.0133   maxfinite = 15.0
Binary8p4se  rmse = 0.0266   maxfinite = 224.0
Binary8p3se  rmse = 0.0528   maxfinite = 49152.0
Binary8p2se  rmse = 0.103    maxfinite = 2.147483648e9
```

For unit-scale data, spend bits on precision; buy range only when your data needs it
(or get it from a block scale — see below).

### Stochastic rounding is unbiased where nearest is not

Nearest rounding maps 0.30078125 to 0.3125 every single time — a systematic +0.0117
bias. Stochastic rounding is right *on average*:

```julia
target = 0.30078125                       # ν = 3/16 of an ulp above 0.296875
σ = ProjSpec(StochasticA{16}(), SatNone())
rng = Xoshiro(1)
m = mean(decode(Convert(Binary8p4se, σ, target; rng)) for _ in 1:200_000)
(decode(Binary8p4se(target)), m)
```

```
(0.3125, 0.300765)                        # RNE result vs stochastic mean ≈ target
```

This is the property that makes stochastic rounding matter for accumulating many
small contributions (gradients, activations statistics) in a low-precision format.

## Deep Learning

### MX-style block quantization — and the staging pitfall

Block formats pair a shared scale with narrow elements, tracking dynamic range that a
bare element format cannot. Here each row of `W` lives at a different scale, spanning
~2¹⁶ across the matrix:

```julia
using ByteFloats, Random, Statistics

rng = Xoshiro(3)
W = randn(rng, 64, 64) .* 2 .* (2.0 .^ (collect(0:63) ./ 4))   # per-row scales

function mx_rows(W, ::Type{FST}, ::Type{FS}, ::Type{FE}, ::Val{B}) where {FST,FS,FE,B}
    n, m = size(W)
    [begin
        # stage through a WIDE-RANGE format so nothing clamps before scaling
        seg = ntuple(k -> Convert(FST, RNE_SatFinite, W[i, (j-1)*B + k]), Val(B))
        ConvertToBlockMaxAbsFinite(FS, FE, RNE_SatNone, RNE_SatNone, seg)
    end for i in 1:n, j in 1:m ÷ B]
end

blocks = mx_rows(W, Binary8p2se, Binary8p1uf, Binary8p4se, Val(32))
recon = [let b = blocks[i, (j-1) ÷ 32 + 1]
             decode(b.s) * decode(b.x[(j-1) % 32 + 1])
         end for i in 1:64, j in 1:64]
plain = [decode(Convert(Binary8p4se, RNE_SatFinite, w)) for w in W]

relerr(A) = sqrt(mean(((W .- A) ./ max.(abs.(W), 1e-9)).^2))
relerr(plain), relerr(recon)
```

```
(0.616, 0.109)      # plain 8p4se clamps the large rows; MX tracks them
```

!!! warning "Stage wide, then block-quantize"
    `ConvertToBlockMaxAbsFinite` takes already-`Binary` elements. If you stage your
    `Float64` data through the *element* format first, `SatFinite` clamps the large
    values **before** the block scale can absorb them, and blocks buy you nothing —
    we measured exactly that (0.617 vs 0.616) with `Binary8p4se` staging. Stage
    through a wide-range format (`Binary8p2se` above); the residual 0.109 here is the
    staging format's own precision, which bounds what any downstream scheme can keep.

At `B = 32` with an 8-bit scale the storage cost is 8.25 bits/value.

### Quantized dot products with one final rounding

`BlockDotProduct` computes every lane product and the accumulation *exactly*, then
projects once:

```julia
rng = Xoshiro(9)
a64, b64 = randn(rng, 32), randn(rng, 32)
qb(v) = ConvertToBlockMaxAbsFinite(Binary8p1uf, Binary8p4se, RNE_SatNone, RNE_SatNone,
            ntuple(i -> Convert(Binary8p4se, RNE_SatFinite, v[i]), Val(32)))
dq = BlockDotProduct(Binary8p4se, RNE_SatNone, qb(a64), qb(b64))
(a64'b64, decode(dq))
```

```
(5.0634, 4.5)       # difference is input quantization only, never accumulation error
```

### Why training loops like stochastic rounding: the swamping demo

Accumulate 400 gradient steps of 0.011 into a `Binary8p3se` accumulator. Under
nearest rounding, the moment the accumulator dwarfs the increment, every add rounds
back to the accumulator — it **stalls**. Stochastic rounding keeps absorbing the
increments in expectation:

```julia
function accumulate_demo(nsteps, g)
    σ = ProjSpec(StochasticA{16}(), SatNone())
    rng = Xoshiro(11)
    acc_rne = Binary8p3se(0.0); acc_sto = Binary8p3se(0.0)
    for _ in 1:nsteps
        acc_rne = Add(Binary8p3se, RNE_SatNone, acc_rne, Convert(Binary8p3se, RNE_SatNone, g))
        acc_sto = Add(Binary8p3se, σ, acc_sto, Convert(Binary8p3se, σ, g; rng); rng)
    end
    (nsteps * g, decode(acc_rne), decode(acc_sto))
end
accumulate_demo(400, 0.011)
```

```
(4.4, 0.125, 5.0)   # exact sum, RNE accumulator (stalled!), stochastic accumulator
```

The stochastic result is noisy (5.0 vs 4.4 on this seed) but unbiased; the RNE result
is *wrong by 35×* and no amount of steps will fix it. In practice you keep a wider
accumulator when you can — and use stochastic rounding when you can't.

### Packing a quantized model

```julia
n = 1_000_000
model = [rawvalue(Binary5p2se, UInt8(rand(0:31))) for _ in 1:n]
pv = PackedVector(model)
(sizeof(model), sizeof(pv.data))
```

```
(1000000, 625000)   # 8 bits → 5 bits per value; indexing and vmap work directly on pv
```
