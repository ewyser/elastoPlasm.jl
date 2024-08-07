# include dependencies
include("../../src/superInclude.jl")
#include("../../src/fun_fs/RetMap/camC/2camCRetMap.jl")
# independant physical constant
g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
E,ν     = 1.0e6,0.3                                                         # Young's mod. [Pa], Poisson's ratio [-]
K,G,Del = D(E,ν,3)                                                     # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
ρ0      = 2700.0                                                            # density [kg/m^3]
yd      = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
cmParam = (E = E, ν = ν, Kc = K, Gc = G, Del = Del,)

pc0   = -cmParam.Kc/3.0
pt    = 0.25*pc0

include("../../src/fun_fs/RetMap/camC/camCmodRetMap.jl")
# camC param
β   = 1.0/1.0
a0  = (pt+pc0)/(β+1)
ϕcs = 20.0*pi/180.0
M   = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
p1  = camCplotYieldFun(pc0,pt,a0,β)
display(plot(p1;layout=(1,1),size=(500,250)))
sleep(2.5)
savefig(path_plot*"pqSpace_camCYieldFunMod.png")

include("../../src/fun_fs/RetMap/camC/camCcohRetMap.jl")
# camC param
β   = 0.25
ϕcs = 20.0*pi/180.0
M   = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
p1  = camCplotYieldFun(pc0,M,β)
display(plot(p1;layout=(1,1),size=(500,250)))
sleep(2.5)
savefig(path_plot*"pqSpace_camCYieldFunCoh.png")

include("../../src/fun_fs/RetMap/camC/camCgenRetMap.jl")
# camC param
ϕcs   = 20.0*π/180.0
M     = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
ζ,γ   = 0.0,0.0
α,β   = 0.0,0.0
p1    = camCplotYieldFun(pc0,pt,γ,M,α,β)
display(plot(p1;layout=(1,1),size=(500,250)))
sleep(2.5)
savefig(path_plot*"pqSpace_camCYieldFunGen.png")

# include("./scripts/unit_testing/camCRetMapTest.jl")
