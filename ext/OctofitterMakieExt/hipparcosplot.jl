
function hipparcosplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains;
    # If less than 500 samples, just show all of them
    # N=min(size(results, 1) * size(results, 3), 500),
    ts,
    # If showing all samples, include each sample once.
    # Otherwise, sample randomly with replacement
    ii=(
        N == size(results, 1) * size(results, 3) ?
        (1:size(results, 1)*size(results, 3)) :
        rand(1:size(results, 1)*size(results, 3), N)
    ),
    axis=(;),
    colormap=:plasma,
    colorbar=true,
    top_time_axis=true,
    bottom_time_axis=true,
    kwargs...
)
    gs = gridspec_or_fig

    # date_pos, date_strs, xminorticks = _date_ticks(ts)
    
    ax_resids = Axis(
        gs[2, 1:3];
        xlabel="MJD",
        ylabel="along-scan residual [mas]",
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=bottom_time_axis,
        xticklabelsvisible=bottom_time_axis,
        xlabelvisible=bottom_time_axis,
        axis...
    )
    xlims!(ax_resids, extrema(ts))

    hip_like = nothing
    for like_obs in model.system.observations
        if like_obs isa HipparcosIADLikelihood
            if !isnothing(hip_like)
                error("more than one HipparcosIADLikelihood present")
            end
            hip_like = like_obs
        end
    end
    if isnothing(hip_like)
        error("No HipparcosIADLikelihood present")
    end

   
    ax_main =  Axis(gs[1, 1:3],
        xlabel="α* [mas]",
        ylabel="δ [mas]",
        autolimitaspect=1
    )


    scatterlines!(
        ax_main,
        hip_like.table.Δα✱, hip_like.table.Δδ,
        label="Hipparcos model",
        linewidth=4,
        alpha=0.5,
        markersize=4,
        color=:grey
    )

    errorbars!(ax_resids, hip_like.table.epoch,  zeros(size(hip_like.table.epoch)), hip_like.table.sres_renorm, color=:black)

    planet_keys = keys(model.system.planets)
    orbits = map(planet_keys) do pk
        Octofitter.construct_elements(results, pk, :)
    end
    nts = Octofitter.mcmcchain2result(model, results)

    I_best = argmax(vec(results[:logpost]))

    for i in [I_best]
        sim = Octofitter.simulate(hip_like, nts[i,], getindex.(orbits,i))
        
        # Model
        scatterlines!(ax_main,
            sim.α✱_model_with_perturbation[:,1],
            sim.δ_model_with_perturbation[:,1],
            label="Our Model",
            color=:red,
            linewidth=1,
            markersize=4,
            alpha=0.85,
        )

        # Data
        for i in axes(hip_like.table.α✱ₘ,1)
            # @show point1 point2
            lines!(ax_main, hip_like.table.α✱ₘ[i][1:2], hip_like.table.δₘ[i][1:2],color=:black,alpha=0.1)
        end

        resid = map(eachindex(hip_like.table.epoch)) do i
            # point = α✱ₘ
            point =  [
                sim.α✱_model_with_perturbation[i],
                sim.δ_model_with_perturbation[i]
            ]
            line_point_1 =  [hip_like.table.α✱ₘ[i][1], hip_like.table.δₘ[i][1]]
            line_point_2 =  [hip_like.table.α✱ₘ[i][2], hip_like.table.δₘ[i][2]]
            Octofitter.distance_point_to_line(point, line_point_1, line_point_2)
        end
        scatter!(ax_resids, hip_like.table.epoch, resid, markersize=4, color=:red, alpha=1)

        for i in 1:length(resid)
            # Need to draw along scan line and residual on the plot
            # start point at data point x and y
            # move at 90 degree angle to the scan direction, for a length of `resid`
            rise = hip_like.table.δₘ[i][2] - hip_like.table.δₘ[i][1]
            run = hip_like.table.α✱ₘ[i][2] - hip_like.table.α✱ₘ[i][1]            
            x0 = sim.α✱_model_with_perturbation[i]
            y0 = sim.δ_model_with_perturbation[i]

            # Two possible directions, check which is closer
            line_point_1 =  [hip_like.table.α✱ₘ[i][1], hip_like.table.δₘ[i][1]]
            line_point_2 =  [hip_like.table.α✱ₘ[i][2], hip_like.table.δₘ[i][2]]
            angle_1 = atan(rise,run) - π/2
            x1_1 = x0 + resid[i]*cos(angle_1)
            y1_1 = y0 + resid[i]*sin(angle_1)
            d1 = Octofitter.distance_point_to_line([x1_1,y1_1], line_point_1, line_point_2)
            angle_2 = atan(rise,run) + π/2
            x1_2 = x0 + resid[i]*cos(angle_2)
            y1_2 = y0 + resid[i]*sin(angle_2)
            d2 = Octofitter.distance_point_to_line([x1_2,y1_2], line_point_1, line_point_2)
            if d1 < d2
                angle = angle_1
                x1 = x1_1 
                y1 = y1_1
            else
                angle = angle_2
                x1 = x1_2
                y1 = y1_2
            end

            # Now we plot the uncertainties along this direction, centred around the intersection point
            unc_x1 = x1 + hip_like.table.sres_renorm[i]*cos(angle)
            unc_y1 = y1 + hip_like.table.sres_renorm[i]*sin(angle)
            unc_x2 = x1 - hip_like.table.sres_renorm[i]*cos(angle)
            unc_y2 = y1 - hip_like.table.sres_renorm[i]*sin(angle)
            lines!(ax_main, [unc_x1,unc_x2], [unc_y1,unc_y2], color=:blue, alpha=0.25, linewidth=4)


            # Plot residual line (from model point to intersection point [x1,y1])
            lines!(ax_main, [x0,x1], [y0,y1], color=:blue, linewidth=1)

        end
    end
    axislegend(ax_main)
    ylims!(ax_resids, low=0)
    return [ax_resids]
end