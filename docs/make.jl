using Documenter,ElastoPlasm

DocMeta.setdocmeta!(ElastoPlasm, :DocTestSetup, :(using ElastoPlasm); recursive=true)

format = Documenter.HTML()
manual = [
    "Home"            => "index.md",
    "Getting Started" => "getting_started.md",
    #"Functions" => mdGenerate(),
]
@info "Making documentation..."
makedocs(;
    modules=[ElastoPlasm],
    authors="madmax",
    sitename="ϵlastσPlasm.jl 👻",
    format=Documenter.HTML(;
        repolink="github.com/ewyser/ElastoPlasm.jl",
        canonical="https://ewyser.github.io/ElastoPlasm.jl/",
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
    repo="github.com/ewyser/ElastoPlasm.jl",
    devbranch="main",
)
