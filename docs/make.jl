# ===== docs/make.jl — Documenter build definition
#
# Invoked by docs/builddocs.jl for local builds, or directly on CI with the docs
# environment already instantiated:  julia --project=docs docs/make.jl

using Documenter
using ByteFloats

makedocs(;
    sitename = "ByteFloats.jl",
    modules = [ByteFloats],
    authors = "ByteFloats.jl contributors",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://example.github.io/ByteFloats.jl",   # set to the real Pages URL
        assets = String[],
    ),
    # Local/non-git builds: disable repository detection (Documenter 1.x errors
    # without a git remote). On CI with a real repository, delete this line and
    # add `edit_link = "main"` to the HTML options to restore edit links.
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Introduction" => "introduction.md",
        "User Guide" => "user_guide.md",
        "User Examples" => "user_examples.md",
        "Technical Guide" => "technical_guide.md",
        "Technical Examples" => "technical_examples.md",
        "Adding Operations" => "new_operations.md",
        "API Reference" => "reference.md",
    ],
    # The guides intentionally reference some names that carry no docstrings yet
    # (docstring coverage is tracked work); keep those as warnings, not failures.
    checkdocs = :none,
    warnonly = [:cross_references, :missing_docs],
)

# Uncomment and point at the real repository to publish from CI:
# deploydocs(;
#     repo = "github.com/<org>/ByteFloats.jl",
#     devbranch = "main",
#     push_preview = true,
# )
