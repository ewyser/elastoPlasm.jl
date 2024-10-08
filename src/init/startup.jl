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
# include dependencies & function call(s)
using Revise,Pkg
using LinearAlgebra,SparseArrays, KernelAbstractions, Plots, LaTeXStrings, Random, Base.Threads,ProgressMeter,REPL.TerminalMenus
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.synchronize as sync
# arithmetic precision & relative path for figs & data
const typeD     = Float64  
const path_DAT = "./misc/DAT/"
if isdir(path_DAT)==false mkdir(path_DAT) end
const path_plot = "./misc/out/"
if isdir(path_plot)==false mkdir(path_plot) end
const path_test = "./misc/test/"
if isdir(path_test)==false mkdir(path_test) end
lists = ["init/api","init/program","init/scripts"]