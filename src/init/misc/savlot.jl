@views function savlot(mpD,t,instr)
    if first(instr[:plot])
        type = last(instr[:plot])
        temp = L"$t = $"*string(round(t,digits=1))*" [s]"
        if type == "P"
            if size(mpD.σ,1) == 3
                d   = -(mpD.σ[1,:]+mpD.σ[2,:])/2/1e3
                lab = L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}\right)/2$"
            elseif size(mpD.σ,1) == 6
                d   = -(mpD.σ[1,:]+mpD.σ[2,:]+mpD.σ[3,:])/3/1e3
                lab = L"$p=-\left(\sigma_{xx,p}+\sigma_{yy,p}+\sigma_{zz,p}\right)/3$"
            end            
            tit   = "pressure, "*temp
            cb    = :viridis
            if minimum(d) == maximum(d)
                cblim = (-1.0,1.0)
            else
                cblim = (minimum(d),maximum(d))
            end
        elseif type == "epII"
            d     = mpD.ϵpII
            lab   = L"$\epsilon_{\mathrm{II}}^{\mathrm{acc}}$"
            tit   = "plastic strain, "*temp
            cb    = :viridis
            if minimum(d) == maximum(d)
                cblim = (-1.0,1.0)
            else
                cblim = (0.0,maximum(d))
            end
        elseif type == "epV"
            d     = mpD.ϵpV
            lab   = L"$\epsilon_{p}^{\mathrm{vol}}$"
            tit   = "volumetric plastic strain, "*temp
            cb    = :seismic
            if minimum(d) == maximum(d)
                cblim = (-1.0,1.0)
            else
                cblim = (-maximum(abs.(d)),maximum(abs.(d)))
            end
        elseif type == "du"
            d     = sqrt.(mpD.u[:,1].^2+mpD.u[:,2].^2)
            lab   = L"$\Delta u$"
            tit   = "displacement, "*temp
            cb    = :viridis
            if minimum(d) == maximum(d)
                cblim = (-1.0,1.0)
            else
                cblim = (0.0,maximum(d))
            end
        else
            err_msg = "$(type): plot option undefined"
            throw(error(err_msg))
        end
        # plot
        gr(legend=true,markersize=2.5,markershape=:circle,markerstrokewidth=0.75,)#markerstrokecolor=:match,)
        p1 = scatter(mpD.x[:,1],mpD.x[:,end],zcolor=d,
        xlabel = L"$x-$direction [m]",
        ylabel = L"$z-$direction [m]",
        label  = lab,
        c      = cb,
        clim   = cblim,
        ylim   = (-10.0,20.0),
        title  = tit,
        aspect_ratio=1,
        )
        display(plot(p1;layout=(1,1),size=(500,250))) 
    else
        nothing
    end
    return nothing
end