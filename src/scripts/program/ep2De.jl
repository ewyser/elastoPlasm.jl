# include dependencies
include("../../src/superInclude.jl")
# main program
@views function ϵp2De(nel::Int64,varPlot::String,cmType::String; kwargs...)
    ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar = getKwargs(kwargs)
    @info "** ϵp2De v$(getVersion()): $(fwrkDeform) strain formulation **"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    E,ν     = 1.0e6,0.3                                                         # Young's mod. [Pa], Poisson's ratio [-]
    K,G,Del = D(E,ν,meD.nD)                                                     # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
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
    cmParam = (E = E, ν = ν, Kc = K, Gc = G, Del = Del, Hp = Hp,)
    @info "mesh & mp feature(s):" nel=Tuple(meD.nel) nno=Tuple(meD.nno) nmp=mpD.nmp
    # plot & time stepping parameters
    tw,tC,it,ctr,ηmax,ηtot = 0.0,1.0,0,0,0,0    
    # action
    @info "launch $(ϕ∂ϕType) calculation cycle..."
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    while tw<=t
        # plot/save
        if tw >= ctr*tC ctr = plotStuff(mpD,tw,varPlot,ctr) end
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
        # update progress bas
        next!(prog;showvalues = getVals(meD,mpD,it,ηmax,ηtot,tw/t,"(✗)"))
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = getVals(meD,mpD,it,ηmax,ηtot,1.0,"(✓)"))
    ctr     = plotStuff(mpD,tw,varPlot,ctr)
    sleep(2.5)
    savefig(path_plot*"$(varPlot)_$(ϕ∂ϕType)_$(fwrkDeform)_$(trsfrAp)_$(isΔFbar)_$(cmType).png")
    @info "Figs saved in" path_plot
    return msg("(✓) Done! exiting...")
end
# include("./scripts/program/ep2De.jl")
# ϵp2De(40,"P","MC";shpfun=:bsmpm,fwrk=:finite,trsf=:mUSL,vollock=true)