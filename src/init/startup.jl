####################################################################################################################################
## start-up function definition
####################################################################################################################################
function superInc(DIR::String)
	msg,inc = ["method(s) sucessfully included:"],[]
	for (root, dirs, files) ∈ walkdir(DIR)
		for file ∈ files
			f = joinpath(root, file)
			if !occursin("/mpi/",f) && last(splitext(f)) == ".jl" 
				include(f)
				push!(inc,file           )
				push!(msg,"\n\t(✓) "*file)
			end
		end
	end
	return inc
end
function tree(sucess, prefix="\n\t", level=0, max_level=1)
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
## start-up info struct definition and instanciation
####################################################################################################################################
Base.@kwdef mutable struct moduleCore
	cpu::NamedTuple = (name=nothing,lab=nothing,mtp=nothing)
	gpu::NamedTuple = (name=nothing,lab=nothing,mtp=nothing)
	root::String    = ROOT
	init::String    = dirname(@__FILE__)
	lib::Vector     = ["api","program","scripts"]
	method::Vector  = []
	out::String     = joinpath(dirname(ROOT),"out")
end
####################################################################################################################################
## conditional list of source code include and external packages deps
####################################################################################################################################
# include dependencies
using Revise,Pkg,Test
using Plots,LaTeXStrings,ProgressMeter,REPL.TerminalMenus
using LinearAlgebra,SparseArrays,Random
using KernelAbstractions,Adapt,Base.Threads
import KernelAbstractions.@atomic as @atom
import KernelAbstractions.synchronize as sync
import Adapt.adapt as user_adapt
import Adapt.@adapt_structure as @user_struct
# instantiate sys
sys = moduleCore()
# arithmetic precision & relative path for figs & data
const typeD     = Float64  



