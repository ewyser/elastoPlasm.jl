function materialGeom(meD,lz,wl,coh0,cohr,ni)
    if meD.nD == 2
        xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
        zL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:lz-0.5*meD.h[2]/ni
        npx,npz     = length(xL),length(zL)
        xp,zp       = ((xL'.*ones(npz,1  )      )),((     ones(npx,1  )'.*zL )) 
        c           = GRFS_gauss(xp,coh0,cohr,ni,meD.h[1])
        xp,zp,c     = vec(xp),vec(zp),vec(c)
        x           = LinRange(minimum(xp),maximum(xp),200)
        a           = -1.25
        x,z         = x.+0.5.*meD.L[1],a.*x
        xlt = Float64[]
        zlt = Float64[]
        clt = Float64[]
        pos = Float64 
        for mp ∈ eachindex(xp)
            for p ∈ eachindex(z)
                Δx = xp[mp]-x[p]
                Δz = zp[mp]-z[p]
                nx = a
                nz = -1.0
                s  = Δx*nx+Δz*nz        
                if s>0
                    pos = 1
                else
                    pos = 0
                end
                if zp[mp]<wl 
                    pos = 1
                end
            end
            if pos==1
                push!(xlt, xp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
                push!(zlt, zp[mp]) # push!(inArray, What), incremental construction of an array of arbitrary size
                push!(clt, c[mp])
            end
        end
    elseif meD.nD == 3
        xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
        yL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:meD.xB[4]
        zL          = meD.xB[5]+(0.5*meD.h[3]/ni):meD.h[3]/ni:lz-0.5*meD.h[3]/ni
        npx,npy,npz = length(xL),length(yL),length(zL)
        xp          = (xL'.*ones(npz,1  )      ).*ones(1,1,npy)
        yp          = (     ones(npz,npx)      ).*reshape(yL,1,1,npy)
        zp          = (     ones(npx,1  )'.*zL ).*ones(1,1,npy)
        c           = GRFS_gauss(xp,coh0,cohr,ni,meD.h[1])   
        xp,yp,zp,c  = vec(xp),vec(yp),vec(zp),vec(c)
        x           = LinRange(minimum(xp),maximum(xp),200)
        a           = -1.25
        x,z         = x.+0.5.*meD.L[1],a.*x
        xlt,ylt,zlt = Float64[],Float64[],Float64[]
        clt         = Float64[]
        pos         = Float64 
        for mp ∈ eachindex(xp)
            for p ∈ eachindex(z)
                Δx = xp[mp]-x[p]
                Δz = zp[mp]-z[p]
                nx = a
                nz = -1.0
                s  = Δx*nx+Δz*nz        
                if s>0.0
                    pos = 1
                else
                    pos = 0
                end
                if zp[mp]<wl 
                    pos = 1
                end
            end
            if pos==1
                push!(xlt, xp[mp]) 
                push!(ylt, yp[mp]) 
                push!(zlt, zp[mp]) 
                push!(clt, c[mp])
            end
        end
    end
    xp = if meD.nD == 2 hcat(xlt,zlt) elseif meD.nD == 3 hcat(xlt,ylt,zlt) end
    id = shuffle(collect(1:size(xp,1)))
    return xp,clt
end
function meshSetup(nel,L,typeD)
    # geometry                                               
    L,h,nD       = meshGeom(L,nel)
    # mesh 
    x,nn,nel,nno = meshCoord(nD,L,h)
    # boundary conditions
    bc,xB        = meshBCs(x,h,nno,nD)
    # constructor
    meD = (
        nD   = nD,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = L,
        h    = h,
        minC = minimum(x,dims=2),
        # nodal quantities
        xn   = x,
        mn   = zeros(typeD,nno[nD+1]             ), # lumped mass vector
        Mn   = zeros(typeD,nno[nD+1],nno[nD+1]   ), # consistent mass matrix
        oobf = zeros(typeD,nno[nD+1],nD          ),
        Dn   = zeros(typeD,nno[nD+1],nD          ),
        fn   = zeros(typeD,nno[nD+1],nD          ),
        an   = zeros(typeD,nno[nD+1],nD          ),
        pn   = zeros(typeD,nno[nD+1],nD          ),
        vn   = zeros(typeD,nno[nD+1],nD          ),
        Δun  = zeros(typeD,nno[nD+1],nD          ),
        ΔJn  = zeros(typeD,nno[nD+1],nD          ),
        bn   = zeros(typeD,nD       ,nD,nno[nD+1]),
        # mesh-to-node topology
        e2n  = e2N(nD,nno,nel,nn),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )
    return meD
end
function pointSetup(meD,L,cmParam,isGRF,typeD)
    coh0,cohr,phi0,phir,rho0 = cmParam[:c0],cmParam[:cr],cmParam[:ϕ0],cmParam[:ϕr],cmParam[:ρ0]
    # non-dimensional constant                                                   
    if meD.nD == 2 ni,nstr = 2,3 elseif meD.nD == 3 ni,nstr = 2,6 end # number of material point along 1d, number of stress components                                                          
    # material geometry
    lz     = L[end]
    wl     = 0.15*lz
    xp,clt = materialGeom(meD,lz,wl,coh0,cohr,ni)
    # scalars & vectors
    nmp    = size(xp,1)
    if meD.nD == 2
        l0,l   = ones(typeD,nmp,meD.nD).*0.5.*(meD.h[1]./ni)  ,ones(typeD,nmp,meD.nD).*0.5.*(meD.h[1]./ni)
        v0,v   = ones(typeD,nmp).*(2.0.*l0[:,1].*2.0.*l0[:,2]),ones(typeD,nmp  ).*(2.0.*l[:,1].*2.0.*l[:,2])
    elseif meD.nD == 3
        l0,l   = ones(typeD,nmp,meD.nD).*0.5.*(meD.h[1]./ni)                ,ones(typeD,nmp,meD.nD).*0.5.*(meD.h[1]./ni)
        v0,v   = ones(typeD,nmp).*(2.0.*l0[:,1].*2.0.*l0[:,2].*2.0.*l0[:,3]),ones(typeD,nmp  ).*(2.0.*l[:,1].*2.0.*l[:,2].*2.0.*l[:,2])
    end
    m      = rho0.*v0
    if isGRF coh = clt else coh = ones(typeD,nmp).*coh0 end
    #coh,phi  = RFS(xp[:,1],xp[:,2],coh0,cohr,phi0,phir)
    cohr   = ones(typeD,nmp).*cohr
    phi    = ones(typeD,nmp).*phi0
    phi[xp[:,end].<=2*wl] .= phir
    # constructor
    mpD = (
        nD   = meD.nD,
        nmp  = nmp,
        x    = xp,
        u    = zeros(typeD,nmp,meD.nD), 
        v    = zeros(typeD,nmp,meD.nD),
        p    = zeros(typeD,nmp,meD.nD),
        l0   = l0,
        l    = l,
        V0   = v0,
        V    = v,
        m    = m,
        c0   = coh,
        cr   = cohr,
        ϕ    = phi,            
        ϵpII = zeros(typeD,nmp),
        ϵpV  = zeros(typeD,nmp), 
        ΔJ   = ones(typeD,nmp),
        J    = ones(typeD,nmp),
        # tensor in matrix notation
        I    = Matrix(1.0I,meD.nD,meD.nD    ),
        ∇u   = zeros(typeD,meD.nD,meD.nD,nmp),
        ΔF   = zeros(typeD,meD.nD,meD.nD,nmp),
        F    = repeat(Matrix(1.0I,meD.nD,meD.nD),1,1,nmp),
        ∇v   = zeros(typeD,meD.nD,meD.nD,nmp),
        ϵ    = zeros(typeD,meD.nD,meD.nD,nmp),
        ω    = zeros(typeD,meD.nD,meD.nD,nmp),
        σJ   = zeros(typeD,meD.nD,meD.nD,nmp),
        b    = repeat(Matrix(1.0I,meD.nD,meD.nD),1,1,nmp),
        # tensor in voigt notation
        σ    = zeros(typeD,nstr,nmp),
        τ    = zeros(typeD,nstr,nmp),
        # additional quantities
        ϕ∂ϕ  = zeros(typeD,meD.nn,nmp ,meD.nD+1   ),
        δnp  = zeros(typeD,meD.nn,meD.nD,nmp      ),
        B    = zeros(typeD,meD.nn.*meD.nD,nstr,nmp),
        # connectivity
        p2e  = zeros(Int64,nmp),
        p2n  = zeros(Int64,meD.nn,nmp),
    )
    # plot initial cohesion field
    plot_coh(mpD.x,mpD.c0,mpD.ϕ,coh0,phi0)
    return mpD 
end
