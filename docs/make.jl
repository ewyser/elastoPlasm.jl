# docs/make.jl
using Documenter,elastoPlasm

format = Documenter.HTML()
pages  = [
    "Home" => "index.md",
    "Getting Started" => "getting_started.md",
    "Functions" => "function.md",
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
    devbranch = "$(getBranch())",
    target = "docs",
    forcepush = true,
)