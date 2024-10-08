@kernel inbounds = true function p2e2D!(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.p2e[p] = (floor(Int64,(mpD.x[p,2]-meD.minC[2])*1.0/meD.h[2])+1)+(meD.nel[2])*floor(Int64,(mpD.x[p,1]-meD.minC[1])*1.0/meD.h[1])
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end
@kernel inbounds = true function p2e3D!(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.p2e[p  ] = (floor(Int64,(mpD.x[p,3]-meD.minC[3])*1.0/meD.h[3])+1)+(meD.nel[3])*floor(Int64,(mpD.x[p,1]-meD.minC[1])*1.0/meD.h[1])+(meD.nel[3]*meD.nel[1])*floor(Int64,(mpD.x[p,2]-meD.minC[2])*1.0/meD.h[2])
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end
function twoDtplgy!(mpD,meD)
    @isdefined(tplgy2D!) ? nothing : tplgy2D! = p2e2D!(CPU())
    tplgy2D!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end
function threeDtplgy!(mpD,meD)
    @isdefined(tplgy3D!) ? nothing : tplgy3D! = p2e3D!(CPU())
    tplgy3D!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end
function shpfun!(mpD,meD,instr)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    if meD.nD == 2
        twoDtplgy!(mpD,meD)
    elseif meD.nD == 3
        threeDtplgy!(mpD,meD)
    end
    # calculate shape functions
    if instr[:shpfun] == :bsmpm 
        ϕ∂ϕbsmpm!(mpD,meD) 
    elseif instr[:shpfun] == :gimpm 
        ϕ∂ϕgimpm!(mpD,meD) 
    elseif instr[:shpfun] == :smpm
        ϕ∂ϕsmpm!(mpD,meD)
    end
    return nothing
end