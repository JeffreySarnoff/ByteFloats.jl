# ===== docs/make.jl — Documenter build definition
#
# Invoked by docs/builddocs.jl for local builds, or directly on CI with the docs
# environment already instantiated:  julia --project=docs docs/make.jl

using Documenter
using ByteFloats

makedocs(;
    sitename = "ByteFloats.jl",
    modules = [ByteFloats],
    authors = "Jeffrey Sarnoff",
    build = get(ENV, "DOCS_BUILD_DIR", "build"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://JeffreySarnoff.github.io/ByteFloats.jl",
        edit_link = "main",
        assets = String[],
    ),
    repo = Documenter.Remotes.GitHub("JeffreySarnoff", "ByteFloats.jl"),
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
    # The guides intentionally reference some names that carry no docstrings yet
    # (docstring coverage is tracked work); keep those as warnings, not failures.
    checkdocs = :none,
    warnonly = [:cross_references, :missing_docs],
)

# Publishes to the gh-pages branch from CI (no-ops on local builds, which lack
# the deploy credentials). Requires GitHub Pages to serve from gh-pages.
deploydocs(;
    repo = "github.com/JeffreySarnoff/ByteFloats.jl",
    devbranch = "main",
    push_preview = true,
)
