# ===== docs/make_latex.jl — LaTeX resolution pass for the PDF pipeline
#
# A copy of make.jl per docs/pdf/generate_pdf.md Stage 1: format swapped to the
# LaTeX writer with platform = "none" (emit resolved .tex, don't compile),
# build directory build_latex, remotes = nothing (local build), deploydocs and
# HTML-only options dropped. pages / modules / checkdocs / warnonly kept
# exactly as make.jl defines them.

using Documenter
using ByteFloats

makedocs(;
    sitename = "ByteFloats.jl",
    modules = [ByteFloats],
    authors = "Jeffrey Sarnoff",
    build = "build_latex",
    format = Documenter.LaTeX(platform = "none", version = v"0.1.0"),
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "Cheat Sheet" => "cheat_sheet.md",
        "User Guide" => "user_guide.md",
        "User Examples" => "user_examples.md",
        "Julia Compatibility" => "julia_compatibility.md",
        "Technical Guide" => "technical_guide.md",
        "Technical Examples" => "technical_examples.md",
        "Adding Operations" => "new_operations.md",
        "External Reference" => "external_reference.md",
        "Internal Reference" => "internal_reference.md",
    ],
    checkdocs = :none,
    warnonly = [:cross_references, :missing_docs],
)
