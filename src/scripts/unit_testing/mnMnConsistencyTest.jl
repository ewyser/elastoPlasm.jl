# include("./scripts/unit_testing/mnMnConsistencyTest.jl")
# include dependencies
using Test, Random, Plots, Pkg, Base.Threads, LaTeXStrings, LinearAlgebra
include("../../src/misc/types.jl")
include("../../src/misc/utilities.jl")
include("../../src/misc/setup.jl")
include("../../src/misc/physics.jl")
include("../../src/misc/plot.jl")
include("../../src/fun_fs/shpfun.jl")
const path_test = "./docs/test/"
const path_plot = "./docs/out/"
const typeD     = Float64 
# main program
@views function mnMn(nel::Int64,varPlot::String,cmType::String; kwargs...)
    ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar,isGRF = getKwargs(kwargs)
    @info "** ϵp2De v$(getVersion()): lumped consistent mass matrix **"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    K,G,Del = D(1.0e6,0.3,2)                                                    # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    ρ0      = 2700.0                                                            # density [kg/m^3]
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    t,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    L       = [64.1584,12.80]                                                   # domain geometry
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,c0,cr,ϕ0,ϕr,ρ0,typeD)                            # material point geometry setup
    Hp      = -60.0e3*meD.h[1]                                                  # softening modulus
    # constitutive model param.
    cmParam = (Kc = K, Gc = G, Del = Del, Hp = Hp,)
    @info "mesh & mp feature(s):" dim=meD.nD nel=Int64(meD.nel[end]) nno=meD.nno[end] nmp=mpD.nmp
    
    # calculate shape functions
    shpfun!(mpD,meD,ϕ∂ϕType)
    # initialize nodal quantities
    meD.Mn  .= 0.0
    meD.mn  .= 0.0
    # mapping back to mesh
    @threads for dim ∈ 1:meD.nD
        @simd for p ∈ 1:mpD.nmp
            # accumulation
            if dim == 1 
                # lumped mass matrix
                meD.mn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p]
                # consistent mass matrix
                meD.Mn[mpD.p2n[:,p],mpD.p2n[:,p]].+= (mpD.ϕ∂ϕ[:,p,1].*mpD.ϕ∂ϕ[:,p,1]').*mpD.m[p] 
            end
        end
    end
    # lumped mass matrix and lumped mass vector
    mn_lumped = sum(meD.Mn,dims=2)
    mn        = meD.mn


    x = reshape(meD.xn[:,1],(meD.nno[2],(meD.nno[1])))
    z = reshape(meD.xn[:,2],(meD.nno[2],(meD.nno[1])))
    x = x[1,:]
    z = sort(z[:,1])

    mn_L_plot = reshape(mn_lumped,(meD.nno[2],(meD.nno[1])))
    mn_plot   = reshape(mn       ,(meD.nno[2],(meD.nno[1])))
    
    
    gr()
    p1 = heatmap(x,z,mn_plot  ,yflip=true, title=L"$m_i = \sum m_{p} \phi_{i}(x_{p}) $")
    p2 = heatmap(x,z,mn_L_plot,yflip=true, title=L"$m_{i,\mathrm{lump}}=\sum m_{ij} = \sum m_p \phi_{i}(x_p)\phi_{j}(x_p)$")
    p3 = heatmap(x,z,log10.(abs.(mn_L_plot-mn_plot)),yflip=true,title=L"\mathrm{log}_{10}(|m_n-m_{i,\mathrm{lump}}|)")
    display(plot(p1,p2,p3; layout=(3,1), size=(450,550)))
    savefig(path_test*"mnVsMn.png")

    return abs.(mn_lumped.-mn)
end

@testset "lumped consistent mass matrix vs. lumped mass vector:" begin
    diff = maximum(mnMn(40,"P","MC"))
    for (it,limit) in enumerate([1e-1,1e-2,1e-4,1e-8,1e-16])
        @test diff < limit 
    end
end
