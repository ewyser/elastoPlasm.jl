# docs/make.jl
push!(LOAD_PATH,"../src/")
include("./src/fun/utils.jl")
using Documenter,elastoPlasm



function mdGenerate(file="function.md",path="docs/src")
    # concatenate in a single stirng the list all function declaration in api
    libs = Dict()
    for (k,lib) ∈ enumerate(keys(elastoPlasm.sys.lib))
        lib_list = undef
        for (l,val) ∈ enumerate(elastoPlasm.sys.lib[lib])
            if l == 1
                lib_list = "\"$(val)\","
            else
                lib_list = lib_list*"\"$(val)\","
            end
        end
        push!(libs,lib=>lib_list)
    end
    # check if path already exists
    if isfile(joinpath(path,file))
        rm(joinpath(path,file))
    end
    # open a .md file and writes corresponding doc
    open(joinpath(path,file),"w") do f
        content = 
        """
        ```@meta
        CollapsedDocStrings = true
        ```

        # Functions

        ## Public API
        ```@autodocs
            Modules = [elastoPlasm]
            Pages   = [$(libs["api"])]
        ```
        ## Internals
        ```@autodocs
            Modules = [elastoPlasm]
            Pages   = [$(libs["program"])]
        ```
        """
        write(f,content)
    end
    return file
end

format = Documenter.HTML()
manual = [
    "Home" => "index.md",
    "Getting Started" => "getting_started.md",
    #"Functions" => mdGenerate(),
]
@info "Making documentation..."
makedocs(
    sitename = "ϵlastσPlasm.jl Documentation",
    authors  = "madMax",
    modules  = [elastoPlasm],
    format   = format,
    pages    = manual,
    pagesonly= true,
    checkdocs= :none,
)
@info "Deploying documentation..."
deploydocs(
    repo = "github.com/ewyser/elastoPlasm.jl.git",
    branch = "gh-pages",
    devbranch = "main",
    target = "docs",
)

#=
# docs/make.jl
using TOML,Documenter,cORIUm
function cORIUmTOML()
    return TOML.parsefile("Project.toml")
end
function getBranch()
    return readchomp(`git rev-parse --abbrev-ref HEAD`)
end
function mdGenerate(file="function.md",path="docs/src")
    # concatenate in a single stirng the list all function declaration in api
    libs = Dict()
    for (k,lib) ∈ enumerate(keys(cORIUm.info.sys.lib))
        lib_list = undef
        for (l,val) ∈ enumerate(cORIUm.info.sys.lib[lib])
            if l == 1
                lib_list = "\"$(val)\","
            else
                lib_list = lib_list*"\"$(val)\","
            end
        end
        push!(libs,lib=>lib_list)
    end
    # check if path already exists
    if isfile(joinpath(path,file))
        rm(joinpath(path,file))
    end
    # open a .md file and writes corresponding doc
    open(joinpath(path,file),"w") do f
        content = 
        """
        ```@meta
        CollapsedDocStrings = true
        ```

        # Functions

        ## Public API
        ```@autodocs
            Modules = [cORIUm]
            Pages   = [$(libs["api"])]
        ```
        ## Internals
        ```@autodocs
            Modules = [cORIUm]
            Pages   = [$(libs["program"])]
        ```
        """
        write(f,content)
    end
    return file
end



#format = Documenter.HTML()
format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"
)
manual = [
    "Home" => "index.md",
    "Getting Started" => "getting_started.md",
    "Functions" => mdGenerate(),
]
@info "Making documentation..."
makedocs(
    sitename = "cORIUm.jl Documentation",
    authors  = "madMax",
    modules  = [cORIUm],
    format   = format,
    pages    = manual,
    pagesonly= true,
    checkdocs= :none,
)
@info "Deploying documentation..."
deploydocs(
    repo = "github.com/ewyser/cORIUm.git",
    branch = "gh-pages",
    devbranch = "$(getBranch())",
    target = "docs",
    versions = ["$(cORIUmTOML()["version"])"],
    forcepush = true,
)
=#