function superInc(DIR)
	msg,included = ["method(s) sucessfully included:"],[]
	for (root, dirs, files) in walkdir(DIR)
		for file in files
			f = joinpath(root, file) # path to files
            include(f)
            push!(included,f)
            push!(msg   ,"\n\t(✓) "*file)
		end
	end
	return included
end

# activate project
using Pkg
if splitpath(Base.active_project())[end-1]!="elastoPlasm"
    Pkg.activate(".")
end
# include dependencies & function call(s)
try 
    using LinearAlgebra, KernelAbstractions, Plots, LaTeXStrings, Random, Base.Threads,ProgressMeter
    import KernelAbstractions.@atomic as @atom
    import KernelAbstractions.synchronize as sync
catch
    @error "$(Base.active_project()) needs instantiation"
    @warn "By-default instantiation launched...may take a while"
    Pkg.instantiate()
    using LinearAlgebra, KernelAbstractions, Plots, LaTeXStrings, Random, Base.Threads,ProgressMeter
    import KernelAbstractions.@atomic as @atom
end
println(pwd())
# arithmetic precision & relative path for figs & data
const typeD     = Float64  
const path_plot = "./misc/out/"
if isdir(path_plot)==false mkdir(path_plot) end
const path_test = "./misc/test/"
if isdir(path_test)==false mkdir(path_test) end

# include doc for: help?> ϵp2De()
include("./misc/doc.jl")

# include init
include("./misc/types.jl")
include("./misc/utilities.jl")
include("./misc/setup.jl")
include("./misc/physics.jl")
include("./misc/plot.jl")

# include functions
if @isdefined(perf) 
    if perf
        @warn "ϵp23De() init performance mode on"
        include("./misc/perf/setup.jl")
        include("./misc/perf/shpfun.jl")
        include("./misc/perf/mapsto.jl")
        include("./misc/perf/solve.jl")
        include("./fun_fs/solve.jl")
        #include("./fun_fs/elastoplast.jl")
        include("./misc/perf/elastoplast.jl")
            include("./fun_fs/RetMap/J2RetMap.jl")
            include("./fun_fs/RetMap/MCRetMap.jl")
            include("./fun_fs/RetMap/DPRetMap.jl")
            #include("./fun_fs/RetMap/camC/camCmodRetMap.jl")
            #include("./fun_fs/RetMap/camC/camCcohRetMap.jl")
            include("./fun_fs/RetMap/camC/camCgenRetMap.jl")
    else
        @warn "ϵp23De() init performance mode off"
        include("./fun_fs/shpfun.jl")
        include("./fun_fs/mapsto.jl")
        include("./fun_fs/solve.jl")
        include("./fun_fs/elastoplast.jl")
            include("./fun_fs/RetMap/J2RetMap.jl")
            include("./fun_fs/RetMap/MCRetMap.jl")
            include("./fun_fs/RetMap/DPRetMap.jl")
            #include("./fun_fs/RetMap/camC/camCmodRetMap.jl")
            #include("./fun_fs/RetMap/camC/camCcohRetMap.jl")
            include("./fun_fs/RetMap/camC/camCgenRetMap.jl")
    end
else
    @info "ϵp23De() init by-default mode"
    include("./fun_fs/shpfun.jl")
    include("./fun_fs/mapsto.jl")
    include("./fun_fs/solve.jl")
    include("./fun_fs/elastoplast.jl")
        include("./fun_fs/RetMap/J2RetMap.jl")
        include("./fun_fs/RetMap/MCRetMap.jl")
        include("./fun_fs/RetMap/DPRetMap.jl")
        #include("./fun_fs/RetMap/camC/camCmodRetMap.jl")
        #include("./fun_fs/RetMap/camC/camCcohRetMap.jl")
        include("./fun_fs/RetMap/camC/camCgenRetMap.jl")
end
