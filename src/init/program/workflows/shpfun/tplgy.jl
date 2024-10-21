@kernel inbounds = true function p2e2n1D!(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.p2e[p] = cld(mpD.x[p,1]-meD.x₀[1],meD.h[1])
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end
@kernel inbounds = true function p2e2n2D!(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        elx        = fld(mpD.x[p,1]-meD.x₀[1],meD.h[1])#floor(Int64,(mpD.x[p,1]-meD.minC[1])*1.0/meD.h[1])
        ely        = cld(mpD.x[p,2]-meD.x₀[2],meD.h[2])#floor(Int64,(mpD.x[p,2]-meD.minC[2])*1.0/meD.h[2])
        mpD.p2e[p] = ely+(meD.nel[2])*elx
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end
@kernel inbounds = true function p2e2n3D!(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.p2e[p] = (floor(Int64,(mpD.x[p,3]-meD.x₀[3])*1.0/meD.h[3])+1)+(meD.nel[3])*floor(Int64,(mpD.x[p,1]-meD.x₀[1])*1.0/meD.h[1])+(meD.nel[3]*meD.nel[1])*floor(Int64,(mpD.x[p,2]-meD.x₀[2])*1.0/meD.h[2])
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end