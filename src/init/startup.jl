####################################################################################################################################
## start-up function definition
####################################################################################################################################
function superInc(DIR::String)
	msg,included = ["method(s) sucessfully included:"],[]
	for (root, dirs, files) in walkdir(DIR)
		for file in files
			f = joinpath(root, file) # path to files
			if !occursin("/mpi/",f) && last(splitext(f)) == ".jl" 
				include(f)
				push!(included,file)
				push!(msg   ,"\n\t(✓) "*file)
			end
		end
	end
	return included
end
function treeLike(sucess, prefix="\n\t", level=0, max_level=1)
    if level > max_level
        return nothing
    end
    n,printout = length(sucess),[]
    for (i, name) in enumerate(sucess)
        connector = i == n ? "└── " : "├── "
		push!(printout,prefix*connector*name)
    end
	return printout
end
####################################################################################################################################
## conditional list of source code include and external packages deps
####################################################################################################################################
# activate project
using Pkg
if splitpath(Base.active_project())[end-1]!="elastoPlasm"
    Pkg.activate(".")
end
# include dependencies & function call(s)
using Revise
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


lists = ["init/scripts","init/misc","init/fun_fs","init/api"]



