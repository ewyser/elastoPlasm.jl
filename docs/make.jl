# docs/make.jl
using Documenter,elastoPlasm

format = Documenter.HTML()
pages  = [
    "Home" => "index.md",
]
@info "Making documentation..."
makedocs(
    sitename = "elastoPlasm.jl Documentation",
    authors  = "madMax",
    modules  = [elastoPlasm],
    format   = format,
    pages    = pages,
)
@info "Deploying documentation..."
deploydocs(
    repo = "github.com/ewyser/elastoPlasm.jl.git",
    branch = "gh-pages",
    devbranch = "main",
    target = "docs",
)