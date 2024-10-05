function D(E,ν,nD)
    Gc = E/(2.0*(1.0+ν))                                                   # shear modulus               [Pa]
    Kc = E/(3.0*(1.0-2.0*ν))                                               # bulk modulus                [Pa]
    if nD == 2
        D  = [ 
            Kc+4/3*Gc Kc-2/3*Gc 0.0 ;
            Kc-2/3*Gc Kc+4/3*Gc 0.0 ;
            0.0       0.0       Gc  ]
    elseif nD == 3
        D  = [ 
            Kc+4/3*Gc Kc-2/3*Gc Kc-2/3*Gc 0.0 0.0 0.0;
            Kc-2/3*Gc Kc+4/3*Gc Kc-2/3*Gc 0.0 0.0 0.0;
            Kc-2/3*Gc Kc-2/3*Gc Kc+4/3*Gc 0.0 0.0 0.0;
            0.0       0.0       0.0       Gc  0.0 0.0;
            0.0       0.0       0.0       0.0 Gc  0.0;
            0.0       0.0       0.0       0.0 0.0 Gc ;]
    end
    return Kc,Gc,D
end
function cm(dim,cmType)
    # independant physical constant          
    E,ν     = 1.0e6,0.3                                                         # Young's mod. [Pa], Poisson's ratio [-]
    K,G,Del = D(E,ν,dim)                                                  # elastic matrix D(E,ν) Young's mod. [Pa] + Poisson's ratio [-]    
    ρ0      = 2700.0                                                            # density [kg/m^3]
    c       = sqrt((K+4.0/3.0*G)/ρ0)                                            # elastic wave speed [m/s]
    c0,cr   = 20.0e3,4.0e3                                                      # cohesion [Pa]
    ϕ0,ϕr,ψ0= 20.0*π/180,7.5*π/180,0.0                                          # friction angle [Rad], dilation angle [Rad]                                                              
    Hp      = -60.0e3                                                           # softening modulus
    # constitutive model param.
    cmParam = (;
        cmType   = cmType, 
        nonlocal = (cond=0,ls=2.5,),
        E   = E, 
        ν   = ν, 
        Kc  = K, 
        Gc  = G, 
        Del = Del, 
        Hp  = Hp, 
        c0  = c0,
        cr  = cr,
        ϕ0  = ϕ0,
        ϕr  = ϕr,
        ρ0  = ρ0,
        c   = c,
    )
    return cmParam
end