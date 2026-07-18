# Independent validation of the yd/fast-path refactors.
# Part A: empirical envelope — libm error vs MPFR truth must be ≪ 2^-45.
# Part B: differential — apply_op (yd→fq→ladder) vs ladder-only reference,
#         exhaustive over codes × formats × modes. Any mismatch = unsoundness.
# Part C: adversarial edges outside the operational pools.
using ByteFloats, Random
using ByteFloats: Float128
using ByteFloats: apply_op, ωeval, project, project_interval, codepoint,
                  rawvalue, decode, bitwidth, EncloseF, BigExactF, Enclose128F

# ---------------- Part A: envelope validation ----------------
softplus64(x) = x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))
function envelope_check(name, f, gen; n=20_000, prec=256)
    rng = Xoshiro(42)
    worst = 0.0; wx = 0.0
    for _ in 1:n
        x = gen(rng)
        y = f(x)
        (isfinite(y) && abs(y) >= 6.7e-290) || continue
        t = setprecision(() -> f(BigFloat(x)), BigFloat, prec)
        iszero(t) && continue
        rel = Float64(abs((BigFloat(y) - t) / t))
        rel > worst && (worst = rel; wx = x)
    end
    lg = worst == 0 ? -Inf : log2(worst)
    ok = worst < 2.0^-45
    println(rpad(name, 12), "max rel err 2^", round(lg, digits=1),
            "  at x=", wx, ok ? "   OK (≪ 2^-45)" : "   *** ENVELOPE VIOLATION ***")
    ok
end
bigsoftplus(b::BigFloat) = b > 0 ? b + log1p(exp(-b)) : log1p(exp(b))
softname(x::BigFloat) = bigsoftplus(x)
println("== Part A: libm-vs-MPFR envelope (claim: error ≪ 2^-45) ==")
okA = true
okA &= envelope_check("exp",   exp,   r -> randn(r) * 40)
okA &= envelope_check("exp2",  exp2,  r -> randn(r) * 60)
okA &= envelope_check("log",   log,   r -> exp(randn(r) * 40) * rand(r))
okA &= envelope_check("log2",  log2,  r -> exp(randn(r) * 40) * rand(r))
okA &= envelope_check("log1p", log1p, r -> rand(r) < 0.5 ? rand(r) * 1e6 : rand(r) * 1.999 - 0.999)
okA &= envelope_check("sinh",  sinh,  r -> randn(r) * 40)
okA &= envelope_check("cosh",  cosh,  r -> randn(r) * 40)
okA &= envelope_check("tanh",  tanh,  r -> randn(r) * 5)
okA &= envelope_check("asinh", asinh, r -> randn(r) * exp(randn(r) * 20))
okA &= envelope_check("atan",  atan,  r -> randn(r) * exp(randn(r) * 20))
okA &= envelope_check("sin",   sin,   r -> randn(r) * exp(rand(r) * 34))   # up to ~1e15 window
okA &= envelope_check("cos",   cos,   r -> randn(r) * exp(rand(r) * 34))
okA &= envelope_check("tan",   tan,   r -> randn(r) * exp(rand(r) * 34))
okA &= (let ok = true
    rng = Xoshiro(7); worst = 0.0; wx = 0.0
    for _ in 1:20_000
        x = randn(rng) * 40
        y = softplus64(x)
        (isfinite(y) && abs(y) >= 6.7e-290) || continue
        t = setprecision(() -> bigsoftplus(BigFloat(x)), BigFloat, 256)
        iszero(t) && continue
        rel = Float64(abs((BigFloat(y) - t) / t))
        rel > worst && (worst = rel; wx = x)
    end
    println(rpad("softplus", 12), "max rel err 2^", round(log2(worst), digits=1),
            "  at x=", wx, worst < 2.0^-45 ? "   OK (≪ 2^-45)" : "   *** ENVELOPE VIOLATION ***")
    worst < 2.0^-45
end)

# ---------------- Part B: differential vs ladder-only reference ----------------
# Reference deliberately bypasses BOTH the yd stage and the fq (Float128) stage.
function refcode(op, ::Type{T}, ρ, R, d) where {T}
    res = ωeval(Val(op), d)
    res isa Float64     && return codepoint(project(T, ρ, res; R))
    res isa Float128    && return codepoint(project(T, ρ, res; R))
    res isa BigExactF   && return codepoint(project(T, ρ, res.f(); R))
    res isa EncloseF    && return codepoint(project_interval(T, ρ, res.f; R))
    res isa Enclose128F && return codepoint(project_interval(T, ρ, res.f; R))
    error("unknown result type $(typeof(res))")
end

function refcode(op, ::Type{T}, ρ, R, x, y) where {T}
    res = ωeval(Val(op), x, y)
    res isa Float64     && return codepoint(project(T, ρ, res; R))
    res isa Float128    && return codepoint(project(T, ρ, res; R))
    res isa BigExactF   && return codepoint(project(T, ρ, res.f(); R))
    res isa EncloseF    && return codepoint(project_interval(T, ρ, res.f; R))
    res isa Enclose128F && return codepoint(project_interval(T, ρ, res.f; R))
    error("unknown result type $(typeof(res))")
end

const UNOPS = (:Exp, :Exp2, :Log, :Log2, :LogOnePlus, :Softplus, :Sinh, :Cosh,
               :Tanh, :ArcSinh, :Sin, :Cos, :Tan, :ArcTan, :Recip, :RSqrt)
modes() = [
    (RNE_SatNone, 0), (ProjSpec(NearestTiesToAway(), SatNone()), 0),
    (ProjSpec(TowardPositive(), SatNone()), 0), (ProjSpec(TowardNegative(), SatNone()), 0),
    (ProjSpec(TowardZero(), SatNone()), 0), (ProjSpec(ToOdd(), SatNone()), 0),
    (RNE_SatFinite, 0), (ProjSpec(NearestTiesToEven(), SatPropagate()), 0),
    (ProjSpec(TowardNegative(), SatFinite()), 0),
    (ProjSpec(StochasticA{8}(), SatNone()), 0),
    (ProjSpec(StochasticA{8}(), SatNone()), 137),
    (ProjSpec(StochasticC{8}(), SatNone()), 255)]
println("\n== Part B: apply_op vs ladder-only reference ==")
function unary_diff()
    mism = 0; total = 0
    for T in (Binary8p4se, Binary8p3se, Binary8p1uf, Binary5p2se, Binary3p1se)
        codes = [rawvalue(T, UInt8(c)) for c in 0:(1 << bitwidth(T)) - 1]
        for op in UNOPS, (ρ, R) in modes(), v in codes
            d = decode(v)
            got = codepoint(apply_op(Val(op), T, ρ, R, d))
            want = refcode(op, T, ρ, R, d)
            total += 1
            if got != want
                mism += 1
                mism <= 5 && println("MISMATCH: $op $(T) ρ=$ρ R=$R x=$d got=$got want=$want")
            end
        end
    end
    total, mism
end
total, mism = unary_diff()
println("unary differential: $total comparisons, $mism mismatches")
function divide_diff()
    dm = 0; dt = 0
    T = Binary8p4se
    codes = [decode(rawvalue(T, UInt8(c))) for c in 0:255]
    for (ρ, R) in [(RNE_SatNone, 0), (ProjSpec(TowardNegative(), SatFinite()), 0),
                   (ProjSpec(ToOdd(), SatNone()), 0), (ProjSpec(StochasticA{8}(), SatNone()), 200)]
        for x in codes, y in codes
            got = codepoint(apply_op(Val(:Divide), T, ρ, R, x, y))
            want = refcode(:Divide, T, ρ, R, x, y)
            dt += 1
            got != want && (dm += 1; dm <= 5 && println("MISMATCH Divide $x/$y ρ=$ρ"))
        end
    end
    U = Binary8p1uf
    rng = Xoshiro(9); ucodes = [decode(rawvalue(U, UInt8(c))) for c in 0:255]
    for _ in 1:4096
        x = rand(rng, ucodes); y = rand(rng, ucodes)
        got = codepoint(apply_op(Val(:Divide), U, RNE_SatNone, 0, x, y))
        want = refcode(:Divide, U, RNE_SatNone, 0, x, y)
        dt += 1
        got != want && (dm += 1)
    end
    dt, dm
end
dt, dm = divide_diff()
println("Divide differential: $dt comparisons, $dm mismatches")

# ---------------- Part C: adversarial edges ----------------
println("\n== Part C: adversarial edges ==")
function edge_checks()
    ok = true
    r1 = ωeval(Val(:Recip), 5.0e-324)            # q = Inf ⇒ yd = NaN ⇒ fq/ladder path
    c1 = codepoint(ByteFloats._finish(Binary8p4se, RNE_SatNone, 0, r1))
    c1r = codepoint(project_interval(Binary8p4se, RNE_SatNone, r1.f))
    ok &= (r1 isa EncloseF && isnan(r1.yd) && c1 == c1r)
    println("Recip(5e-324): yd=NaN routed, finish==ladder: ", c1 == c1r)
    r2 = ωeval(Val(:Divide), floatmax(), 5.0e-324)   # q = Inf
    r3 = ωeval(Val(:Divide), 5.0e-324, floatmax())   # q subnormal
    ok &= (r2 isa EncloseF && isnan(r2.yd) && r3 isa EncloseF)
    println("Divide over/underflow degenerate routing: ", r2 isa EncloseF && r3 isa EncloseF)
    wok = true
    for x in (0.999e15, 1.0e15, 1.000001e15)
        res = ωeval(Val(:Sin), x)
        inw = x <= 1.0e15
        wok &= (res isa EncloseF) && (inw ? !isnan(res.yd) : (isnan(res.yd) && res.fq === nothing))
    end
    ok &= wok
    println("Sin window edge (yd/fq only inside |x| ≤ 1e15): ", wok)
    hok = true
    for x in (3.141592653589793, 6.283185307179586, 1.5707963267948966,
              1.0e15 - 0.7853981633974483, 708.0, -708.0)
        for op in (:Sin, :Cos, :Tan, :Exp)
            got = codepoint(apply_op(Val(op), Binary8p4se, RNE_SatNone, 0, x))
            want = refcode(op, Binary8p4se, RNE_SatNone, 0, x)
            got == want || (hok = false; println("EDGE MISMATCH $op($x)"))
        end
    end
    ok &= hok
    println("hard reduction / overflow-boundary points: ", hok ? "all agree with ladder" : "MISMATCH")
    dok = all(ωeval(Val(:Log2), ldexp(1.0, e)) == Float64(e) for e in -20:20)
    ok &= dok
    println("Log2 dyadic screen intact: ", dok)
    ok
end
okC = edge_checks()

println("\nVERDICT: ", (okA && mism == 0 && dm == 0 && okC) ?
        "all validations passed" : "*** VALIDATION FAILURES PRESENT ***")