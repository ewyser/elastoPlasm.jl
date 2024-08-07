@views function solve!(meD,Δt)
    # viscous damping
    η   = 0.1
    # initialize
    meD.f .= 0.0
    meD.a .= 0.0
    meD.v .= 0.0
    # solve momentum equation on the mesh
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.m[n]>0.0 
            m          = 1.0/meD.m[n]
            fx         = meD.fext[n,1]-meD.fint[n,1]              
            fy         = meD.fext[n,2]-meD.fint[n,2]              
            Dx         = η*sqrt(fx^2+fy^2)*sign(meD.p[n,1]*m)  
            Dy         = η*sqrt(fx^2+fy^2)*sign(meD.p[n,2]*m)  
            fx         = fx-Dx                                          
            fy         = fy-Dx                                     
            meD.a[n,1] = fx*m*meD.bc[n,1]                          
            meD.a[n,2] = fy*m*meD.bc[n,2]                          
            meD.v[n,1] = (meD.p[n,1]+Δt*fx)*m*meD.bc[n,1]          
            meD.v[n,2] = (meD.p[n,2]+Δt*fy)*m*meD.bc[n,2]          
        end
    end
    return nothing
end