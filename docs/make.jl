push!(LOAD_PATH,"../src/")
using Documenter
using elastoPlasm

DocMeta.setdocmeta!(elastoPlasm, :DocTestSetup, :(using elastoPlasm); recursive=true)

format = Documenter.HTML()
manual = [
    "Home" => "index.md",
    "Getting Started" => "getting_started.md",
    #"Functions" => mdGenerate(),
]
@info "Making documentation..."
makedocs(;
    modules=[elastoPlasm],
    authors="madmax",
    sitename="elastoPlasm.jl 👻",
    format=Documenter.HTML(;
        canonical="https://mewyser.github.io/elastoPlasm.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs= :none,
)
@info "Deploying documentation..."
deploydocs(
    repo="github.com/ewyser/elastoPlasm.jl",
    devbranch="main",
)
