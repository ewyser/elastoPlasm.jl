function meshGeom(L,nel)
    nD = length(L)
    if nD == 2
        L   = [L[1],ceil(L[2])]
        h   = [L[1]/nel,L[1]/nel]
    elseif nD == 3
        L   = [L[1],L[2],ceil(L[3])]
        h   = [L[1]/nel[1],L[1]/nel[1],L[1]/nel[1]]
    else 
        err_msg = "nD = $(nD), L= $(L): unsupported mesh geometry"
        throw(error(err_msg))
    end
    return L,h,nD
end
function meshCoord(nD,L,h)
    if nD == 2
        xn  = collect((0.0-2*h[1]):h[1]:(L[1]+2.0*h[1])) 
        zn  = reverse(collect((0.0-2*h[2]):h[2]:(L[2]+2.0*h[2])))
        nno = [length(xn),length(zn),length(xn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,(nno[1]-1)*(nno[2]-1)]
        nn  = 16
        xn  = (xn'.*ones(typeD,nno[2],1     ))     
        zn  = (     ones(typeD,nno[1],1     )'.*zn)
        x   = hcat(vec(xn),vec(zn))
    elseif nD == 3
        xn  = collect((0.0-2*h[1]):h[1]:(L[1]+2.0*h[1])) 
        yn  = collect((0.0-2*h[2]):h[2]:(L[2]+2.0*h[2])) 
        zn  = reverse(collect((0.0-2*h[3]):h[3]:(L[3]+2.0*h[3])))        
        nno = [length(xn),length(yn),length(zn),length(xn)*length(yn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]
        nn  = 64
        xn  = (xn'.*ones(typeD,nno[3],1     ))     .*ones(typeD,1,1,nno[2])
        zn  = (     ones(typeD,nno[1],1     )'.*zn).*ones(typeD,1,1,nno[2])
        yn  = (     ones(typeD,nno[3],nno[1]))     .*reshape(yn,1,1,nno[2])
        x   = hcat(vec(xn),vec(yn),vec(zn))
    end
    return x,nn,nel,nno
end
function meshBCs(xn,h,nno,nD)
    if nD == 2
        xB  = [minimum(xn[:,1])+2*h[1],maximum(xn[:,1])-2*h[1],
               0.0                    ,Inf                    ]                                    
        bcx = vcat(findall(x->x<=xB[1], xn[:,1]),findall(x->x>=xB[2], xn[:,1]))
        bcz = findall(x->x<=xB[3], xn[:,2])
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0.0
        #bcX[bcz] .= 0.0
        bc   = hcat(bcX,bcZ)
    elseif nD == 3
        xB  = [minimum(xn[:,1])+2*h[1],maximum(xn[:,1])-2*h[1],
               minimum(xn[:,2])+2*h[2],maximum(xn[:,2])-2*h[2],
               0.0                    ,Inf                    ]       
        bcx = vcat(findall(x->x<=xB[1], xn[:,1]),findall(x->x>=xB[2], xn[:,1]))
        bcy = vcat(findall(x->x<=xB[3], xn[:,2]),findall(x->x>=xB[4], xn[:,2]))
        bcz = findall(x->x<=xB[5], xn[:,3])
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bcY = ones(nno[nD+1],1)
        bcY[bcy] .= 0.0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0.0
        bc   = hcat(bcX,bcY,bcZ)
    end
    return bc,xB
end
function e2N(nD,nno,nel,nn)
	e2n  = zeros(Int64,nel[nD+1],nn)
    if nD == 2
        gnum = reverse(reshape(1:(nno[3]),nno[2],nno[1]),dims=1)
        iel  = 1
        for i ∈ 1:nel[1]#nelx
            for j ∈ 1:nel[2]#nelz
                if i>1 && i<nel[1] && j>1 && j<nel[2]
                    e2n[iel,1 ] = gnum[j-1,i-1]
                    e2n[iel,2 ] = gnum[j-0,i-1]
                    e2n[iel,3 ] = gnum[j+1,i-1]
                    e2n[iel,4 ] = gnum[j+2,i-1]

                    e2n[iel,5 ] = gnum[j-1,i  ]
                    e2n[iel,6 ] = gnum[j-0,i  ]
                    e2n[iel,7 ] = gnum[j+1,i  ]
                    e2n[iel,8 ] = gnum[j+2,i  ]

                    e2n[iel,9 ] = gnum[j-1,i+1]
                    e2n[iel,10] = gnum[j-0,i+1]
                    e2n[iel,11] = gnum[j+1,i+1]
                    e2n[iel,12] = gnum[j+2,i+1]

                    e2n[iel,13] = gnum[j-1,i+2]
                    e2n[iel,14] = gnum[j-0,i+2]
                    e2n[iel,15] = gnum[j+1,i+2]
                    e2n[iel,16] = gnum[j+2,i+2]
                end
                iel = iel+1;
            end
        end
    elseif nD == 3
        gnum = reverse(reshape(1:(nno[4]),nno[3],nno[1],nno[2]),dims=1)
        iel  = 1
        for k ∈ 1:nel[2]#nely
            for i ∈ 1:nel[1]#nelx
                for j ∈ 1:nel[3]#nelz
                    if i>1 && i<nel[1] && j>1 && j<nel[3] && k>1 && k<nel[2]
                        e2n[iel,1 ] = gnum[j-1,i-1,k-1]
                        e2n[iel,2 ] = gnum[j-0,i-1,k-1]
                        e2n[iel,3 ] = gnum[j+1,i-1,k-1]
                        e2n[iel,4 ] = gnum[j+2,i-1,k-1]
                        e2n[iel,5 ] = gnum[j-1,i  ,k-1]
                        e2n[iel,6 ] = gnum[j-0,i  ,k-1]
                        e2n[iel,7 ] = gnum[j+1,i  ,k-1]
                        e2n[iel,8 ] = gnum[j+2,i  ,k-1]
                        e2n[iel,9 ] = gnum[j-1,i+1,k-1]
                        e2n[iel,10] = gnum[j-0,i+1,k-1]
                        e2n[iel,11] = gnum[j+1,i+1,k-1]
                        e2n[iel,12] = gnum[j+2,i+1,k-1]
                        e2n[iel,13] = gnum[j-1,i+2,k-1]
                        e2n[iel,14] = gnum[j-0,i+2,k-1]
                        e2n[iel,15] = gnum[j+1,i+2,k-1]
                        e2n[iel,16] = gnum[j+2,i+2,k-1]
                        
                        e2n[iel,17] = gnum[j-1,i-1,k  ]
                        e2n[iel,18] = gnum[j-0,i-1,k  ]
                        e2n[iel,19] = gnum[j+1,i-1,k  ]
                        e2n[iel,20] = gnum[j+2,i-1,k  ]
                        e2n[iel,21] = gnum[j-1,i  ,k  ]
                        e2n[iel,22] = gnum[j-0,i  ,k  ]
                        e2n[iel,23] = gnum[j+1,i  ,k  ]
                        e2n[iel,24] = gnum[j+2,i  ,k  ]
                        e2n[iel,25] = gnum[j-1,i+1,k  ]
                        e2n[iel,26] = gnum[j-0,i+1,k  ]
                        e2n[iel,27] = gnum[j+1,i+1,k  ]
                        e2n[iel,28] = gnum[j+2,i+1,k  ]
                        e2n[iel,29] = gnum[j-1,i+2,k  ]
                        e2n[iel,30] = gnum[j-0,i+2,k  ]
                        e2n[iel,31] = gnum[j+1,i+2,k  ]
                        e2n[iel,32] = gnum[j+2,i+2,k  ]
                        
                        e2n[iel,33] = gnum[j-1,i-1,k+1]
                        e2n[iel,34] = gnum[j-0,i-1,k+1]
                        e2n[iel,35] = gnum[j+1,i-1,k+1]
                        e2n[iel,36] = gnum[j+2,i-1,k+1]
                        e2n[iel,37] = gnum[j-1,i  ,k+1]
                        e2n[iel,38] = gnum[j-0,i  ,k+1]
                        e2n[iel,39] = gnum[j+1,i  ,k+1]
                        e2n[iel,40] = gnum[j+2,i  ,k+1]
                        e2n[iel,41] = gnum[j-1,i+1,k+1]
                        e2n[iel,42] = gnum[j-0,i+1,k+1]
                        e2n[iel,43] = gnum[j+1,i+1,k+1]
                        e2n[iel,44] = gnum[j+2,i+1,k+1]
                        e2n[iel,45] = gnum[j-1,i+2,k+1]
                        e2n[iel,46] = gnum[j-0,i+2,k+1]
                        e2n[iel,47] = gnum[j+1,i+2,k+1]
                        e2n[iel,48] = gnum[j+2,i+2,k+1]
                            
                        e2n[iel,49] = gnum[j-1,i-1,k+2]
                        e2n[iel,50] = gnum[j-0,i-1,k+2]
                        e2n[iel,51] = gnum[j+1,i-1,k+2]
                        e2n[iel,52] = gnum[j+2,i-1,k+2]
                        e2n[iel,53] = gnum[j-1,i  ,k+2]
                        e2n[iel,54] = gnum[j-0,i  ,k+2]
                        e2n[iel,55] = gnum[j+1,i  ,k+2]
                        e2n[iel,56] = gnum[j+2,i  ,k+2]
                        e2n[iel,57] = gnum[j-1,i+1,k+2]
                        e2n[iel,58] = gnum[j-0,i+1,k+2]
                        e2n[iel,59] = gnum[j+1,i+1,k+2]
                        e2n[iel,60] = gnum[j+2,i+1,k+2]
                        e2n[iel,61] = gnum[j-1,i+2,k+2]
                        e2n[iel,62] = gnum[j-0,i+2,k+2]
                        e2n[iel,63] = gnum[j+1,i+2,k+2]
                        e2n[iel,64] = gnum[j+2,i+2,k+2]
                    end
                    iel = iel+1;
                end
            end
        end
    end
    e2n = permutedims(e2n,(2,1))
	return e2n
end
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
function pointSetup(meD,L,coh0,cohr,phi0,phir,rho0,isGRF,typeD)
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
        # connectivity
        p2e  = zeros(Int64,nmp),
        p2n  = zeros(Int64,meD.nn,nmp),
    )
    # plot initial cohesion field
    plot_coh(mpD.x,mpD.c0,mpD.ϕ,coh0,phi0)
    return mpD 
end
