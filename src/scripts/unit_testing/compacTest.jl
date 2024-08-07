# include("./scripts/unit_testing/compacTest.jl")
# include dependencies
include("../../src/superInclude.jl")
using Test
# main program
function meshGeom(L,nel)
    nD = length(L)
    if nD == 2
        h   = [L[1],L[2]/nel]
    elseif nD == 3
        h   = [L[1],L[2],L[3]/nel]
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
        xn  =         Array(range(0.0-2*h[1],L[1]+2.0*h[1],step=h[1]))
        yn  =         Array(range(0.0-2*h[2],L[2]+2.0*h[2],step=h[2]))
        zn  = reverse(Array(range(0.0-2*h[3],L[3]+2.0*h[3],step=h[3])))        
        nno = [length(xn),length(yn),length(zn),length(xn)*length(yn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]
        nn  = 64
        xn  = (xn'.*ones(typeD,nno[3],1     ))     .*ones(typeD,1,1,nno[2])
        zn  = (     ones(typeD,nno[1],1     )'.*zn).*ones(typeD,1,1,nno[2])
        yn  = (     ones(typeD,nno[3],nno[1]))     .*reshape(yn,1,1,nno[2])
        xn  = vec(xn)
        yn  = vec(yn)
        zn  = vec(zn)
        x   = hcat(xn,yn,zn)
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
        fext = zeros(typeD,nno[nD+1],nD          ), 
        fint = zeros(typeD,nno[nD+1],nD          ),
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
function materialGeomCompact(meD,lz,wl,coh0,cohr,ni)
    if meD.nD == 2
        xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
        zL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:lz-0.5*meD.h[2]/ni
        npx,npz     = length(xL),length(zL)
        xp,zp       = ((xL'.*ones(npz,1  )      )),((     ones(npx,1  )'.*zL )) 
        xp,zp       = vec(xp),vec(zp)
    elseif meD.nD == 3
        xL          = meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2]
        yL          = meD.xB[3]+(0.5*meD.h[2]/ni):meD.h[2]/ni:meD.xB[4]
        zL          = meD.xB[5]+(0.5*meD.h[3]/ni):meD.h[3]/ni:lz-0.5*meD.h[3]/ni
        npx,npy,npz = length(xL),length(yL),length(zL)
        xp          = (xL'.*ones(npz,1  )      ).*ones(1,1,npy)
        yp          = (     ones(npz,npx)      ).*reshape(yL,1,1,npy)
        zp          = (     ones(npx,1  )'.*zL ).*ones(1,1,npy)
        xp,yp,zp    = vec(xp),vec(yp),vec(zp)
    end
    id = shuffle(collect(1:size(xp,1)))
    return if meD.nD == 2 hcat(xp[id],zp[id]) elseif meD.nD == 3 hcat(xp[id],yp[id],zp[id]) end
end
function pointSetup(meD,L,coh0,cohr,phi0,phir,rho0,typeD)
    # non-dimensional constant                                                   
    if meD.nD == 2 ni,nstr = 2,3 elseif meD.nD == 3 ni,nstr = 2,6 end # number of material point along 1d, number of stress components                                                          
    # material geometry
    lz     = L[end]
    wl     = 0.15*lz
    xp     = materialGeomCompact(meD,lz,wl,coh0,cohr,ni)
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
    coh    = ones(typeD,nmp ).*coh0
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
        coh  = coh,
        cohr = cohr,
        phi  = phi,
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
        b    = repeat(Matrix(1.0I,meD.nD,meD.nD),1,1,nmp),
        # tensor in voigt notation
        ω    = zeros(typeD,nmp),
        σR   = zeros(typeD,nstr,nmp),
        σ    = zeros(typeD,nstr,nmp),
        τ    = zeros(typeD,nstr,nmp),
        dev  = zeros(typeD,nstr,nmp),
        ep   = zeros(typeD,nstr,nmp),
        # additional quantities
        ϕ∂ϕ  = zeros(typeD,meD.nn,nmp ,meD.nD+1   ),
        δnp  = zeros(typeD,meD.nn,meD.nD,nmp      ),
        B    = zeros(typeD,meD.nn.*meD.nD,nstr,nmp),
        # connectivity
        p2e  = zeros(Int64,nmp),
        p2n  = zeros(Int64,meD.nn,nmp),
    )
    return mpD 
end
@views function plotStuff(mpD,t,type,ctr,title)
    xlab,ylab = L"$x-$direction",L"$z-$direction"
    gr(size=(2*250,2*125),legend=true,markersize=2.5,markershape=:circle,markerstrokewidth=0.0,markerstrokecolor=:match,)
    temp = title
    if type == "P"
        p = -mpD.σ[end,:]/1e3
        scatter(mpD.x[:,1],mpD.x[:,2],zcolor=p,
            xlabel = xlab,
            ylabel = ylab,
            label=L"$\sigma_{zz}$",
            aspect_ratio=1,
            c=:viridis,
            ylim=(0.0,50.0),
            title=temp,
            show=true,
            )  
    end
    return ctr+=1
end
@views function solve!(meD,Δt)
    # viscous damping
    η      = 0.05
    # initialize
    meD.fn.= 0.0
    meD.an.= 0.0
    meD.vn.= 0.0
    # solve momentum equation on the mesh
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.mn[n]>0.0 
            m            = (1.0/meD.mn[n]).*meD.bc[n,:]                   #(2,)
            meD.Dn[n,:] .= η.*norm(meD.oobf[n,:]).*sign.(meD.pn[n,:].*m)  #(2,)
            meD.fn[n,:] .= meD.oobf[n,:].-meD.Dn[n,:]                     #(2,)
            meD.an[n,:] .= meD.fn[n,:].*m                                 #(2,)
            meD.vn[n,:] .= (meD.pn[n,:].+Δt.*meD.fn[n,:]).*m              #(2,)
        end
    end
    return nothing
end
@views function getVals(it,cmpl,symb)
    # completion [%]
    cmpl = round(100.0*cmpl,digits=1)
    # save vals
    vals = [("iteration(s)",it),
            (symb*" t/T",cmpl)]
    return vals
end

@views function compactTest(dim,nel,varPlot,ν,E,ρ0,l0; kwargs...)
    cmType = "MC"
    ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar,isGRF = getKwargs(kwargs)
    # mesh & mp setup
    if dim == 2
        L       = [l0/nel,l0       ]                                           # domain geometry
    elseif dim == 3
        L       = [l0/nel,l0/nel,l0]                                            # domain geometry
    end
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    K,G,Del = D(E,ν,meD.nD)                                                     # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    tg      = ceil((1.0/yd)*(2.0*l0)*40.0)
    t,te    = 1.25*tg,1.25*tg
    mpD     = pointSetup(meD,L,c0,cr,ϕ0,ϕr,ρ0,typeD)                            # material point geometry setup 
    z0      = copy(mpD.x[:,end])
    Hp      = -60.0e3*meD.h[1]                                                  # softening modulus
    # constitutive model param.
    cmParam = (Kc = K, Gc = G, Del = Del, Hp = Hp,)
    if meD.nD == 2
        @info "** ϵp$(length(L))De v$(getVersion()): compaction of a two-dimensional column under self weight **"
        @info "mesh & mp feature(s):" nel=Int64(meD.nel[2]-4)
    elseif meD.nD == 3
        @info "** ϵp$(length(L))De v$(getVersion()): compaction of a three-dimensional column under self weight **"
        @info "mesh & mp feature(s):" nel=Int64(meD.nel[3]-4)
    end
    # plot & time stepping parameters
    tw,tC,it,ctr,ηmax,ηtot = 0.0,1.0,0,0,0,0    
    # action
    if meD.nD == 2 g = [0.0,0.0] elseif meD.nD == 3 g = [0.0,0.0,0.0] end
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    while tw<=t
        # set clock on/off
        tic = time_ns()
        # adaptative Δt & linear increase in gravity
        Δt,g  = get_Δt(mpD.v,meD.h,yd),get_g(tw,tg,meD.nD)
        # bsmpm cycle
        shpfun!(mpD,meD,ϕ∂ϕType)
        mapsto!(mpD,meD,g,Δt,trsfrAp,"p->n")                  
        solve!(meD,Δt)
        mapsto!(mpD,meD,g,Δt,trsfrAp,"p<-n")
        ηmax = elastoplast!(mpD,meD,cmParam,cmType,Δt,ϕ∂ϕType,isΔFbar,fwrkDeform,tw>te)
        # update sim time
        tw,it,toc,ηtot = tw+Δt,it+1,((time_ns()-tic)),max(ηmax,ηtot)
        next!(prog;showvalues = getVals(it,tw/t,"(✗)"))
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = getVals(it,1.0,"(✓)"))
    # analytics
    if meD.nD==2
        xN,yN = abs.(mpD.σ[2,:]),z0
    elseif meD.nD==3
        xN,yN = abs.(mpD.σ[3,:]),z0
    end
    xA,yA = abs.(ρ0.*g[end].*(l0.-z0)),z0
    err   = sum(sqrt.((xN.-xA).^2).*mpD.V0)/(abs(g[end])*ρ0*l0*sum(mpD.V0))
    return (xN,yN,xA,yA),meD.h,err
end
@views function compacTest(dim,trsfrAp)
    ϕ∂ϕType    = :gimpm
    fwrkDeform = :finite
    store,H,error = [],[],[]
    try
        @testset "convergence using $(ϕ∂ϕType), $(fwrkDeform) deformation, $(trsfrAp) mapping" begin
            # geometry
            n         = [0,1,2,3,4,5,6]
            nel       = 2.0.^n
            # initial parameters 
            l0,ν,E,ρ0 = 50.0,0.0,1.0e4,80.0
            # init error
            ϵ         = 1.0
            for (it,nel) in enumerate(nel)
                #action
                DAT,h,err = compactTest(dim,nel,"P",ν,E,ρ0,l0;shpfun=ϕ∂ϕType,fwrk=fwrkDeform,trsf=trsfrAp,vollock=true)
                push!(store,DAT )
                push!(H ,h[end])
                push!(error,err)
                # test
                @test (err < ϵ )
                ϵ = err
            end 
        end
    catch

    end
    xN,yN,xA,yA = store[end]
    return H,error
end
@views function runCompacTest(DIM,TSF)
    for k ∈ 1:length(DIM)
        H,error = compacTest(DIM[k],TSF[k])
        p2 = gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
        if DIM[k] == 2 shape = :circle elseif DIM[k] == 3 shape = :star end
        if k == 1
            p2 = plot( 1.0./H,error,seriestype=:line,markerize=5.0,markershape=shape,label="$(DIM[1])D, $(TSF[1])") 
        elseif k<length(DIM)
            p2 = plot!(1.0./H,error,seriestype=:line,markerize=5.0,markershape=shape,label="$(DIM[k])D, $(TSF[k])") 
        elseif k == length(DIM)
            p2 = plot!(1.0./H,error,seriestype=:line,markerize=5.0,markershape=shape,label="$(DIM[k])D, $(TSF[k])",xlabel=L"$1/h$ [m$^{-1}$]",ylabel="error",xaxis=:log10,yaxis=:log10) 
            display(plot(p2; layout=(1,1), size=(450,250)))
            savefig(path_test*"23D_convergence_pass_compacTest.png")
        end
    end
    return "all tests passed...exit"
end

if @isdefined(perf)
    if perf
        runCompacTest([2,3],[:mUSL,:mUSL])
    else
        runCompacTest([2,3],[:mUSL,:mUSL])
    end
else
    runCompacTest([2,3,2,3],[:mUSL,:mUSL,:tpicUSL,:tpicUSL])
end
#runCompacTest([2,3,2,3,2,3],[:mUSL,:mUSL,:picflipUSL,:picflipUSL,:tpicUSL,:tpicUSL])
#runCompacTest([2,2,2],[:picflipUSL,:mUSL,:tpicUSL])







































#=
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    p1 = plot(xN.*1e-3,yN,seriestype=:scatter, label="$(dim)D $(ϕ∂ϕType), $(trsfrAp) mapping")
    p1 = plot!(xA.*1e-3,yA,label=L"\sum_{p}\dfrac{||\sigma_{yy}^p-\sigma_{yy}^a(x_p)||V_0^p}{(g\rho_0l_0)V_0}",xlabel=L"$\sigma_{yy}$ [kPa]",ylabel=L"$y-$position [m]") 
    display(plot(p1; layout=(1,1), size=(450,250)))
    savefig(path_test*"$(dim)D_numericVsAnalytic_compacTest_$(ϕ∂ϕType)_$(fwrkDeform)_$(trsfrAp).png")
=#