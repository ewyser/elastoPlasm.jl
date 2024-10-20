module ElastoPlasm
# define module location as const
const ROOT = dirname(@__FILE__)
# include startup file
include(joinpath(ROOT,"init/startup.jl"))
# include .jl files
sucess = ["welcome to ϵlastσPlasm 👻 \nsucessful superInclude()"]
for (k,child) ∈ enumerate(sys.lib)
	list = superInc(joinpath(sys.init,child))
	if isempty(list)
		push!(sucess,"\n✗ "*child)
	else	
		push!(sys.method,("$(child)"=>list))
		push!(sucess    ,"\n✓ "*child      )
	end
end
@info join(sucess)
@info """new comer ?
- copy-paste the followings:
  L,nel = [64.1584,12.80],40
  instr = slump(L,nel);
- wait for the simulation to end
"""
end # module elastoPlasm

