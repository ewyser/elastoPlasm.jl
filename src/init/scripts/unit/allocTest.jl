# include("./src/scripts/unit_testing/allocTest.jl")
@warn "unit testing"
using BenchmarkTools
@views function allocCheck(L::Vector{Float64},nel::Int64,varPlot::String,cmType::String; kwargs...)
    ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar,isGRF = getKwargs(kwargs)
    @info "** ϵp$(length(L))De v$(getVersion()): $(fwrkDeform) strain formulation **"
    @info "init..."
    # mesh setup
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    E,ν     = 1.0e6,0.3                                                         # Young's mod. [Pa], Poisson's ratio [-]
    K,G,Del = D(E,ν,meD.nD)                                                     # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    ρ0      = 2700.0                                                            # density [kg/m^3]
    yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    t,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mp setup
    mpD     = pointSetup(meD,L,c0,cr,ϕ0,ϕr,ρ0,isGRF,typeD)                      # material point geometry setup
    Hp      = -60.0e3*meD.h[1]                                                  # softening modulus
    # constitutive model param.
    cmParam = (Kc = K, Gc = G, Del = Del, Hp = Hp,)
    @info "mesh & mp feature(s):" ϕ∂ϕType fwrkDeform trsfrAp isΔFbar nel nthreads()
    # plot & time stepping parameters
    tw,tC,it,ctr,toc,flag,ηmax,ηtot,Δt = 0.0,1.0/1.0,0,0,0.0,0,0,0,1.0e-4    
    # action
    @info "Evaluate core functions:"
    println("launch ϕ∂ϕ!()")
    @btime shpfun!($mpD,$meD,$ϕ∂ϕType)
    println("launch mapsto!(p->n)")
    @btime mapsto!($mpD,$meD,vec([0.0,0.0,9.81]),$Δt,$trsfrAp,"p->n")   
    println("launch solve!()")
    @btime solve!($meD,$Δt)
    println("launch mapsto!(p<-n)")
    @btime mapsto!($mpD,$meD,vec([0.0,0.0,9.81]),$Δt,$trsfrAp,"p<-n")
    println("launch elastoplast!()")
    @btime ηmax = elastoplast!($mpD,$meD,$cmParam,$cmType,$Δt,$ϕ∂ϕType,$isΔFbar,$fwrkDeform,true)
    @warn "Digging deeper in elastoplast!(), "
    println("-> launch deform!()")
    @btime deform!($mpD,$meD,$Δt,$ϕ∂ϕType,$isΔFbar)
    println("-> launch ΔFbar!()")
    @btime ΔFbar!($mpD,$meD)
    println("-> launch elast!()")
    @btime elast!($mpD,$cmParam.Del,$fwrkDeform)
    return msg("(✓) Done! exiting...")
end
export allocCheck
#allocCheck([64.1584,64.1584/4,12.80],80,"P","DP";shpfun=:bsmpm,fwrk=:finite,trsf=:mUSL,vollock=true)