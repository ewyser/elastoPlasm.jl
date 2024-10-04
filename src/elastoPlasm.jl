module elastoPlasm
# define module location as const
const ROOT = dirname(@__FILE__)
# include startup file
include(joinpath(ROOT,"init/startup.jl"))
# include .jl files
sucess = ["elastoPlasm: sucessful superInclude()"]
for (k,child) ∈ enumerate(lists)
	list = superInc(joinpath(ROOT,child))
	push!(sucess,"\n\t✓ "*child)
	if haskey(ENV,"TREE_SUPERINC") && ENV["TREE_SUPERINC"]=="true"
		push!(sucess,join(treeLike(list)))
	end
end
@info join(sucess)
end # module elastoPlasm
