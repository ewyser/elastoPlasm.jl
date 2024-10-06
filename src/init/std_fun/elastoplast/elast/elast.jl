@views function mutate(ϵ,Χ,type)
    if type == :tensor # α = 1/2 when ϵ := strain, α = 1.0 when ϵ := stress
        if size(ϵ) == (3,)
            ϵmut = [  ϵ[1] Χ*ϵ[3];
                    Χ*ϵ[3]   ϵ[2]]
        elseif size(ϵ) == (6,)
            ϵmut = [  ϵ[1] Χ*ϵ[6] Χ*ϵ[5];
                    Χ*ϵ[6]   ϵ[2] Χ*ϵ[4];
                    Χ*ϵ[5] Χ*ϵ[4]   ϵ[3]]
        end
    elseif type == :voigt # α = 2.0 when ϵ := strain, α = 1.0 when ϵ := stress
        if size(ϵ) == (2,2)
            ϵmut = vcat(ϵ[1,1],ϵ[2,2],Χ*ϵ[1,2]) #xx,yy,zz,xy
        elseif size(ϵ) == (3,3)
            ϵmut = vcat(ϵ[1,1],ϵ[2,2],ϵ[3,3],Χ*ϵ[2,3],Χ*ϵ[1,3],Χ*ϵ[1,2]) #xx,yy,zz,yz,xz,xy
        end
    end
    return ϵmut
end
@views @kernel inbounds = true function ELAST!(mpD,Del,instr)
    p = @index(Global)
    # deformation framework dispatcher
    if instr[:fwrk] == :finite
        if p ≤ mpD.nmp 
            # update left cauchy-green tensor
            mpD.b[:,:,p].= mpD.ΔF[:,:,p]*mpD.b[:,:,p]*mpD.ΔF[:,:,p]'
            # compute logarithmic strain tensor
            λ,n          = eigen(mpD.b[:,:,p],sortby=nothing)
            mpD.ϵ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mpD.τ[:,p]   = Del*mutate(mpD.ϵ[:,:,p],2.0,:voigt)
        end
    elseif instr[:fwrk] == :infinitesimal
        if p ≤ mpD.nmp 
            # calculate elastic strains & spins
            for i ∈ 1:mpD.nD,j ∈ 1:mpD.nD
                mpD.ϵ[i,j,p] = 0.5*(mpD.ΔF[i,j,p]+mpD.ΔF[j,i,p])-mpD.I[i,j]
                mpD.ω[i,j,p] = 0.5*(mpD.ΔF[i,j,p]-mpD.ΔF[j,i,p])
            end
            Δϵ = mutate(mpD.ϵ[:,:,p],2.0,:voigt)
            Δω = mpD.ω[2,1,p]
            σ0 = [Δω*(2.0*mpD.σ[3,p]),Δω*(-2.0*mpD.σ[3,p]),Δω*(mpD.σ[2,p]-mpD.σ[1,p])]
            for k ∈ 1:length(σ0)
                mpD.σ[k,p]+= (Del[k,1]*Δϵ[1]+Del[k,2]*Δϵ[2]+Del[k,3]*Δϵ[3])+σ0[k]
            end
        end   
    end
end
@views @kernel inbounds = true function elast!(mpD,Del,instr)
    p = @index(Global)
    # deformation framework dispatcher
    if instr[:fwrk] == :finite
        if p ≤ mpD.nmp 
            # update left cauchy-green tensor
            mpD.b[:,:,p].= mpD.ΔF[:,:,p]*mpD.b[:,:,p]*mpD.ΔF[:,:,p]'
            # compute logarithmic strain tensor
            λ,n          = eigen(mpD.b[:,:,p],sortby=nothing)
            mpD.ϵ[:,:,p].= 0.5.*(n*diagm(log.(λ))*n')
            # krichhoff stress tensor
            mpD.τ[:,p]   = Del*mutate(mpD.ϵ[:,:,p],2.0,:voigt)
        end
    elseif instr[:fwrk] == :infinitesimal
        if p ≤ mpD.nmp 
            # calculate elastic strains & spins
            mpD.ϵ[:,:,p] .= 0.5.*(mpD.ΔF[:,:,p]+mpD.ΔF[:,:,p]').-mpD.I
            mpD.ω[:,:,p] .= 0.5.*(mpD.ΔF[:,:,p]-mpD.ΔF[:,:,p]')
            # update cauchy stress tensor
            mpD.σJ[:,:,p].= mutate(mpD.σ[:,p],1.0,:tensor)
            mpD.σJ[:,:,p].= mpD.σJ[:,:,p]*mpD.ω[:,:,p]'+mpD.σJ[:,:,p]'*mpD.ω[:,:,p]
            mpD.σ[:,p]  .+= Del*mutate(mpD.ϵ[:,:,p],2.0,:voigt).+mutate(mpD.σJ[:,:,p],1.0,:voigt)
        end   
    end
end
@views @kernel inbounds = true function transform!(mpD)
    p = @index(Global)
    # deformation framework dispatcher
    if p ≤ mpD.nmp 
        mpD.σ[:,p] .= mpD.τ[:,p]./mpD.J[p]
    end   
end
function stress!(mpD,cmParam,instr,type)
    if type == :update
        if @isdefined(stresses!) 
            nothing 
        else
            if instr[:perf]
                stresses! = ELAST!(CPU())
            else
                stresses! = elast!(CPU())
            end
        end
        stresses!(mpD,cmParam.Del,instr; ndrange=mpD.nmp)
        sync(CPU())
    elseif type == :cauchy
        @isdefined(cauchy!) ? nothing : cauchy! = transform!(CPU())
        cauchy!(mpD; ndrange=mpD.nmp)
        sync(CPU())
    end
    return nothing
end