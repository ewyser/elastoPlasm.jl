#using REPL.TerminalMenus
#options   = ["standard","performance"]
#select    = request("select mode:",MultiSelectMenu(options))
#println(select)

module elastoPlasm
# define module location as const
const ROOT = dirname(@__FILE__)
# include startup file
include(joinpath(ROOT,"init/startup.jl"))


lists = ["init/scripts","init/misc","init/fun","init/api"]
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
