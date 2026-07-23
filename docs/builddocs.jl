#!/usr/bin/env julia
# ===== docs/builddocs.jl — one-command local documentation build
#
#   julia docs/builddocs.jl          (from the repository root)
#   julia builddocs.jl               (from docs/)
#
# Activates the docs environment, dev's the package from the parent directory,
# instantiates, runs make.jl, and reports where the site landed. CI should skip
# this file and run `julia --project=docs docs/make.jl` after its own Pkg setup.

import Pkg

const DOCS = @__DIR__
const PKG  = dirname(DOCS)

Pkg.activate(DOCS)
Pkg.develop(Pkg.PackageSpec(path = PKG))
Pkg.instantiate()

include(joinpath(DOCS, "make.jl"))

const BUILD = get(ENV, "DOCS_BUILD_DIR", "build")
const INDEX = joinpath(DOCS, BUILD, "index.html")
if isfile(INDEX)
    @info "Documentation built successfully" site = INDEX
else
    @warn "make.jl completed but no index.html found — check the build log" dir = joinpath(DOCS, BUILD)
end
