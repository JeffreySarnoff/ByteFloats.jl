# ===== benchmark/benchmarking.jl — Chairmarks benchmark suite and report generator
#
# Run:  julia --project=benchmark benchmark/benchmarking.jl [output.md]
# or:   include(...); generate_report("report.md")
#
# Measurement doctrine (checkpoint.md, "Resolution of the two flagged measurements"):
#   • every benchmarked call sits behind a type-specialized function barrier —
#     format types enter as type parameters, never as non-const globals;
#   • operands come from Chairmarks' *untimed* setup, so values are runtime-varying
#     (no constant folding) without paying generation cost in the timing;
#   • specialization is asserted before anything is measured (preflight): if the
#     warm scalar paths allocate, the numbers would measure dispatch, so we abort;
#   • medians are reported with minima alongside; times are per single call unless
#     a table says per-element.

using Chairmarks
using Statistics: median
using Random
using P3109
using P3109: project, order_key, get_table, empty_tables!, _f128, opinfo, OP_REGISTRY,
             _UNARY_OPS, _BINARY_OPS, _TERNARY_OPS, rawvalue, decode

# ---------------------------------------------------------------------------
# formatting helpers
# ---------------------------------------------------------------------------
function fmt_time(s::Float64)
    s < 1e-6 && return string(round(s * 1e9; digits=1), " ns")
    s < 1e-3 && return string(round(s * 1e6; digits=2), " μs")
    s < 1.0  && return string(round(s * 1e3; digits=2), " ms")
    string(round(s; digits=3), " s")
end
fmt_alloc(x) = x == 0 ? "0" : string(round(Int, x))

struct Row
    name::String
    med::Float64      # seconds
    min::Float64
    allocs::Float64
    extra::String     # optional trailing column (throughput etc.)
end
Row(name, b::Chairmarks.Benchmark; extra="") =
    Row(name, median(b).time, minimum(b).time, median(b).allocs, extra)

function write_table(io, title, note, rows; extra_header="")
    println(io, "\n## ", title, "\n")
    isempty(note) || println(io, note, "\n")
    hdr = "| operation | median | min | allocs |"
    sep = "|---|---|---|---|"
    if !isempty(extra_header)
        hdr *= " $extra_header |"
        sep *= "---|"
    end
    println(io, hdr); println(io, sep)
    for r in rows
        line = "| `$(r.name)` | $(fmt_time(r.med)) | $(fmt_time(r.min)) | $(fmt_alloc(r.allocs)) |"
        isempty(extra_header) || (line *= " $(r.extra) |")
        println(io, line)
    end
end

# random finite-heavy operand pools (all codes ⇒ honest NaN/Inf mix where noted)
codes_pool(::Type{T}, n) where {T<:Binary} =
    [rawvalue(T, UInt8(rand(0:(1 << bitwidth(T)) - 1))) for _ in 1:n]

# ---------------------------------------------------------------------------
# preflight: abort rather than publish dispatch measurements
# ---------------------------------------------------------------------------
function preflight(::Type{T}) where {T<:Binary}
    a, b, c = T(1.5), T(0.25), T(2.0)
    add2(x, y) = Add(T, RNE_SatNone, x, y);  add2(a, b)
    prj(d) = project(T, RNE_SatNone, d);     prj(2.3)
    fma3(x, y, z) = FMA(T, RNE_SatNone, x, y, z); fma3(a, b, c)
    ok = @allocated(add2(a, b)) == 0 && @allocated(prj(2.3)) == 0 && @allocated(fma3(a, b, c)) == 0
    ok || error("preflight failed: warm scalar paths allocate — measurements would reflect dispatch, not arithmetic")
    nothing
end

# ---------------------------------------------------------------------------
# benchmark sections (each returns Vector{Row}; T enters as a type parameter)
# ---------------------------------------------------------------------------
function bench_primitives(::Type{T}) where {T<:Binary}
    pool = codes_pool(T, 4096)
    σ = ProjSpec(StochasticA{8}(), SatNone())
    rows = Row[]
    push!(rows, Row("decode",        @be rand(pool) decode(_)))
    push!(rows, Row("order_key",     @be rand(pool) order_key(_)))
    push!(rows, Row("x < y",         @be (rand(pool), rand(pool)) (t -> t[1] < t[2])(_)))
    push!(rows, Row("TotalOrder",    @be (rand(pool), rand(pool)) (t -> TotalOrder(t[1], t[2]))(_)))
    push!(rows, Row("Class",         @be rand(pool) Class(_)))
    push!(rows, Row("NextGreaterThan", @be rand(pool) NextGreaterThan(_)))
    push!(rows, Row("project (RNE·SatNone)", @be decode(rand(pool)) project(T, RNE_SatNone, _)))
    push!(rows, Row("project (StochasticA[8], R drawn)",
                    @be decode(rand(pool)) project(T, σ, _)))
    rows
end

# `getfield(P3109, op)` infers `::Any`; captured directly it would make every
# benchmarked call a dynamic dispatch (exactly the harness failure the doctrine
# exists to prevent — caught in review when the same op measured 10× slower here
# than in the sensitivity table). Passing `f` through an argument specializes on
# its concrete function type.
_bench_op(f::F, ::Type{T}, pool, ::Val{1}) where {F,T} =
    @be rand(pool) f(T, RNE_SatNone, _)
_bench_op(f::F, ::Type{T}, pool, ::Val{2}) where {F,T} =
    @be (rand(pool), rand(pool)) (t -> f(T, RNE_SatNone, t[1], t[2]))(_)
_bench_op(f::F, ::Type{T}, pool, ::Val{3}) where {F,T} =
    @be (rand(pool), rand(pool), rand(pool)) (t -> f(T, RNE_SatNone, t[1], t[2], t[3]))(_)

function bench_scalar_ops(::Type{T}, names, arity) where {T<:Binary}
    pool = codes_pool(T, 4096)
    rows = Row[]
    for op in names
        b = _bench_op(getfield(P3109, op), T, pool, Val(arity))
        push!(rows, Row(string(op), b))
    end
    sort!(rows; by=r -> r.med)
end

function bench_format_sensitivity(ops)
    rows = Row[]
    for F in (Binary8p4se, Binary8p3sf, Binary8p1uf, Binary5p2se, Binary3p1se), op in ops
        pool = codes_pool(F, 4096)
        b = _bench_op(getfield(P3109, op), F, pool, Val(2))
        push!(rows, Row("$(op)⟨$(P3109.formatname(F))⟩", b))
    end
    rows
end

function bench_modes(::Type{T}) where {T<:Binary}
    pool = [decode(v) for v in codes_pool(T, 4096)]
    rows = Row[]
    for (nm, ρ) in [("NearestTiesToEven", RNE_SatNone),
                    ("NearestTiesToAway", ProjSpec(NearestTiesToAway(), SatNone())),
                    ("TowardPositive", ProjSpec(TowardPositive(), SatNone())),
                    ("TowardNegative", ProjSpec(TowardNegative(), SatNone())),
                    ("TowardZero", ProjSpec(TowardZero(), SatNone())),
                    ("ToOdd", ProjSpec(ToOdd(), SatNone())),
                    ("StochasticA[8]", ProjSpec(StochasticA{8}(), SatNone())),
                    ("StochasticB[8]", ProjSpec(StochasticB{8}(), SatNone())),
                    ("StochasticC[8]", ProjSpec(StochasticC{8}(), SatNone())),
                    ("RNE · SatFinite", RNE_SatFinite),
                    ("RNE · SatPropagate", ProjSpec(NearestTiesToEven(), SatPropagate()))]
        push!(rows, Row(nm, @be rand(pool) project(T, ρ, _)))
    end
    rows
end

function bench_kernels(::Type{T}; n=65536) where {T<:Binary}
    A = codes_pool(T, n); B = codes_pool(T, n); C = codes_pool(T, n)
    σ = ProjSpec(StochasticA{8}(), SatNone())
    get_table(:Exp, T, T, RNE_SatNone)               # warm the caches: measure gather, not build
    get_table(:Add, T, T, T, RNE_SatNone)
    perel(b, m) = string(round(median(b).time / m * 1e9; digits=2), " ns/elem — ",
                         round(m / median(b).time / 1e9; digits=2), " Gelem/s")
    rows = Row[]
    b = @be similar(A) vmap!(_, Val(:Exp), T, RNE_SatNone, A) evals=1
    push!(rows, Row("vmap unary (table gather), n=$n", b; extra=perel(b, n)))
    b = @be similar(A) vmap!(_, Val(:Add), T, RNE_SatNone, A, B) evals=1
    push!(rows, Row("vmap binary (table gather), n=$n", b; extra=perel(b, n)))
    b = @be similar(A) vmap!(_, Val(:FMA), T, RNE_SatNone, A, B, C) evals=1
    push!(rows, Row("vmap ternary (scalar loop), n=$n", b; extra=perel(b, n)))
    b = @be (similar(A), Xoshiro(1)) (t -> vmap!(t[1], Val(:Add), T, σ, A, B, t[2]))(_) evals=1
    push!(rows, Row("vmap binary stochastic (scalar loop), n=$n", b; extra=perel(b, n)))
    pv = PackedVector(A)
    b = @be P3109.vmap(:Exp, T, RNE_SatNone, pv) evals=1
    push!(rows, Row("vmap unary through PackedVector, n=$n", b; extra=perel(b, n)))
    rows
end

function bench_sorting(::Type{T}; n=65536) where {T<:Binary}
    A = codes_pool(T, n)
    rows = Row[]
    b = @be copy(A) sort!(_) evals=1
    push!(rows, Row("sort! (counting sort via defalg), n=$n", b))
    b = @be copy(A) sort!(_; alg=Base.Sort.DEFAULT_UNSTABLE) evals=1
    push!(rows, Row("sort! (stock comparison sort), n=$n", b))
    b = @be copy(A) sort!(_; rev=true) evals=1
    push!(rows, Row("sort! rev=true (counting sort), n=$n", b))
    rows
end

function bench_table_builds()
    T = Binary8p4se; U = Binary8p1uf
    rows = Row[]
    specs = [("Exp⟨8p4se⟩ (256 entries)",        () -> get_table(:Exp, T, T, RNE_SatNone)),
             ("Tanh⟨8p4se⟩ (256 entries)",       () -> get_table(:Tanh, T, T, RNE_SatNone)),
             ("Add⟨8p4se×8p4se⟩ (64 K entries)", () -> get_table(:Add, T, T, T, RNE_SatNone)),
             ("Divide⟨8p4se×8p4se⟩ (64 K)",      () -> get_table(:Divide, T, T, T, RNE_SatNone)),
             ("Add⟨8p1uf×8p1uf⟩ (64 K, wide-spread)", () -> get_table(:Add, U, U, U, RNE_SatNone))]
    for (nm, f) in specs
        f()                                           # JIT warm; builds are cache-evicted per sample
        b = @be empty_tables!() (_ -> f())(_) evals=1 seconds=3
        push!(rows, Row(nm, b))
    end
    empty_tables!()
    rows
end

function bench_blocks(::Type{FE}, ::Type{FS}; B=32) where {FE<:Binary,FS<:Binary}
    mk() = Block(one(FS), ntuple(_ -> rawvalue(FE, UInt8(rand(0:(1 << bitwidth(FE)) - 1))), B))
    x = mk(); y = mk()
    perlane(b) = string(round(median(b).time / B * 1e9; digits=2), " ns/lane")
    rows = Row[]
    b = @be (mk(), mk()) (t -> BlockAdd(FE, RNE_SatNone, t[1], t[2], one(FS)))(_)
    push!(rows, Row("BlockAdd (B=$B)", b; extra=perlane(b)))
    b = @be (mk(), mk()) (t -> BlockDotProduct(Binary8p4se, RNE_SatNone, t[1], t[2]))(_)
    push!(rows, Row("BlockDotProduct → 8p4se (B=$B)", b; extra=perlane(b)))
    b = @be mk() BlockReduceAdd(Binary8p4se, RNE_SatNone, _)
    push!(rows, Row("BlockReduceAdd → 8p4se (B=$B)", b; extra=perlane(b)))
    tup = ntuple(_ -> rawvalue(FE, UInt8(rand(0:(1 << bitwidth(FE)) - 1))), B)
    b = @be ConvertToBlockMaxAbsFinite(FS, FE, RNE_SatNone, RNE_SatNone, tup)
    push!(rows, Row("ConvertToBlockMaxAbsFinite (B=$B)", b; extra=perlane(b)))
    rows
end

function bench_conversions()
    T = Binary8p4se; S = Binary8p3se
    pool = codes_pool(T, 4096)
    dpool = [decode(v) for v in pool]
    A = codes_pool(T, 65536)
    rows = Row[]
    cpool = [UInt8(rand(0:255)) for _ in 1:4096]
    push!(rows, Row("T(::UInt8) code-point constructor",
                    @be rand(cpool) Binary8p4se(_)))
    push!(rows, Row("rawvalue (unchecked kernel route)",
                    @be rand(cpool) rawvalue(Binary8p4se, _)))
    push!(rows, Row("T(::Float64) numeric constructor (projects)",
                    @be rand(dpool) Binary8p4se(_)))
    push!(rows, Row("Convert 8p4se → 8p3se (scalar)",
                    @be rand(pool) Convert(S, RNE_SatNone, _)))
    push!(rows, Row("Float64 → 8p4se (project)",
                    @be rand(dpool) project(T, RNE_SatNone, _)))
    b = @be PackedVector(A) evals=1
    push!(rows, Row("PackedVector pack, n=65536", b;
                    extra=string(round(median(b).time / 65536 * 1e9; digits=2), " ns/elem")))
    pv = PackedVector(A)
    b = @be collect(pv) evals=1
    push!(rows, Row("PackedVector unpack (collect), n=65536", b;
                    extra=string(round(median(b).time / 65536 * 1e9; digits=2), " ns/elem")))
    rows
end

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
function generate_report(path::AbstractString="benchmark_report.md"; seed=2026)
    Random.seed!(seed)
    T = Binary8p4se
    preflight(T)
    open(path, "w") do io
        println(io, "# ByteFloats.jl benchmark report")
        println(io, "\nGenerated: ", string(Dates_now()), "  ·  Julia ", VERSION,
                "  ·  ", Sys.CPU_NAME, " (", Sys.CPU_THREADS, " threads)",
                "  ·  Float128 paths: ", _f128() ? "enabled" : "disabled",
                "  ·  Chairmarks ", pkgversion(Chairmarks))
        println(io, "\nReference format for per-operation tables: `Binary8p4se` under ",
                "`(NearestTiesToEven, SatNone)`; operand pools draw uniformly over all ",
                "code points (an honest NaN/±Inf mix). Times are per call; medians with ",
                "minima alongside. Methodology per the recorded benchmark doctrine: ",
                "type-parameterized barriers, untimed setup, specialization preflight.")

        write_table(io, "Core primitives",
            "The decode/compare/step/classify layer plus the projection engine.",
            bench_primitives(T))
        write_table(io, "Scalar operations — unary (30)",
            "Sorted by median. Random operands: transcendental rows mix special-row " *
            "fast returns with enclosure-path evaluations, so these are *scalar-path* " *
            "costs; bulk unary work routes through 256-byte tables (see Array kernels).",
            bench_scalar_ops(T, _UNARY_OPS, 1))
        write_table(io, "Scalar operations — binary (18)", "Sorted by median.",
            bench_scalar_ops(T, _BINARY_OPS, 2))
        write_table(io, "Scalar operations — ternary (3)", "Sorted by median.",
            bench_scalar_ops(T, _TERNARY_OPS, 3))
        write_table(io, "Format sensitivity",
            "Same three binary ops across formats; `Binary8p1uf` exercises the " *
            "wide-exponent-spread escalations, small-K formats the tiny-table regime.",
            bench_format_sensitivity((:Add, :Divide, :Multiply)))
        write_table(io, "Projection by rounding/saturation mode",
            "`project(Binary8p4se, ρ, x)` over the mode vocabulary (stochastic budgets N = 8).",
            bench_modes(T))
        write_table(io, "Array kernels (vmap)",
            "Warm caches: table specializations prebuilt, so table rows measure the " *
            "gather; scalar-loop rows measure the full compute pipeline per element.",
            bench_kernels(T); extra_header="per element")
        write_table(io, "Sorting (64 K values)",
            "Counting sort is installed as the default algorithm for `Binary` vectors.",
            bench_sorting(T))
        write_table(io, "Table builds (oracle + projection, Float128-first)",
            "Cold cache per sample (`empty_tables!` in untimed setup); JIT pre-warmed.",
            bench_table_builds())
        write_table(io, "Block and scaled operations",
            "Elements `Binary8p4se`, scales `Binary8p1uf`, B = 32.",
            bench_blocks(Binary8p4se, Binary8p1uf); extra_header="per lane")
        write_table(io, "Conversions and packed storage", "",
            bench_conversions(); extra_header="per element")

        println(io, "\n---\n*All numbers from this machine/run; absolute values vary ",
                "by host. Regenerate with `julia --project=benchmark benchmark/benchmarking.jl`.*")
    end
    path
end

# Dates without adding a dependency: ISO-ish stamp from libc
Dates_now() = Libc.strftime("%Y-%m-%d %H:%M UTC", time())

if abspath(PROGRAM_FILE) == @__FILE__
    out = isempty(ARGS) ? "benchmark_report.md" : ARGS[1]
    println("generating ", out, " …")
    generate_report(out)
    println("done: ", out)
end
