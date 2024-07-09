
##################################################
# Astrometry vs time plot
function astromtimeplot(
    model,
    results,
    fname="$(model.system.name)-astromtime.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700,600),
        figure...
    )
    astromtimeplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function astromtimeplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains;
    # If less than 500 samples, just show all of them
    N  = min(size(results, 1)*size(results, 3), 500),
    ts,
    # If showing all samples, include each sample once.
    # Otherwise, sample randomly with replacement
    ii = (
        N == size(results, 1)*size(results, 3) ? 
        (1:size(results, 1)*size(results, 3)) :
        rand(1:size(results, 1)*size(results, 3),N)
    ),
    axis=(;),
    colormap=:plasma,
    colorbar=true,
    top_time_axis=true,
    bottom_time_axis=true,
    mark_epochs_mjd=Float64[],
    kwargs...
)
    gs = gridspec_or_fig

    
    date_pos, date_strs, xminorticks = _date_ticks(ts)
    ax_sep = Axis(
        gs[1, 1];
        ylabel="sep [mas]",
        xaxisposition=:top,
        xticks=(date_pos, date_strs),
        xgridvisible=false,
        ygridvisible=false,
        xminorticks,
        xminorticksvisible=top_time_axis,
        xticksvisible=top_time_axis,
        xticklabelsvisible=top_time_axis,
        axis...
    )
    ax_pa = Axis(
        gs[2, 1];
        xlabel="MJD",
        ylabel="PA [°]",
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=bottom_time_axis,
        xticklabelsvisible=bottom_time_axis,
        xlabelvisible=bottom_time_axis,
        axis...
    )
    linkxaxes!(ax_sep, ax_pa)
    xlims!(ax_sep, extrema(ts))
    xlims!(ax_pa, extrema(ts))


    if length(model.system.planets) == 1
        colormaps = Dict(
            first(keys(model.system.planets)) => colormap
        )
    else
        colormaps = Dict(
            begin
                c = Makie.wong_colors()[i]
                planet_key => Makie.cgrad([c, "#FAFAFA"])
            end
            for (i,planet_key) in enumerate(keys(model.system.planets))
        )
    end

    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        # Now time-series
        sols = orbitsolve.(orbs, ts')
        
        sep_model_t = projectedseparation.(sols)
        pa_model_t = rem2pi.(posangle.(sols), RoundDown)
        color_model_t = rem2pi.(
            eccanom.(sols), RoundDown)
        
        lines!(ax_sep,
            concat_with_nan(ts' .+ 0 .* sep_model_t),
            concat_with_nan(sep_model_t);
            alpha=min.(1, 100 / length(ii)),
            color=concat_with_nan(color_model_t .+ 0 .* ii),
            colorrange=(0,2pi),
            colormap=colormaps[planet_key]
        )
        
        # Go through each PA time series and, whenever we detect PA has wrapped past 0 or 360,
        # do some work to make the line not jump back vertically across the plot.
        # Instead, insert a NaN line break and project the line above/below the top/bottom plot limits.
        # Also, set the plot limits to (0,360)
        allx = Float64[]
        ally = Float64[]
        allcolors = Float64[]
        for yi in axes(pa_model_t, 1)
            x = ts
            y = pa_model_t[yi,:]
            colors = color_model_t[yi,:]
            xi = 0
            xi0 = 0
            # TODO: this does fail if the wrapping happened between first and second indices
            direction = mean(sign.(pa_model_t[yi,2:end] - pa_model_t[yi,1:end-1])) > 0 ? +1.0 : -1.0
            while xi<size(x,1)-1
                xi += 1
                xi0 += 1
                direction_next = sign(pa_model_t[yi,xi0+1] - pa_model_t[yi,xi0])
                if direction != direction_next && direction > 0
                    Δ =  pa_model_t[yi,xi0+1]-(pa_model_t[yi,xi0]-2pi)
                    x_midpoint = (x[xi] + x[xi+1])/2
                    newpt_before = y[xi] + Δ/2
                    newpt_after =  y[xi+1] - Δ/2
                    x = [
                        x[1:xi];
                        x_midpoint;
                        NaN;
                        x_midpoint;
                        x[xi+1:end]
                    ]
                    y = [
                        y[1:xi];
                        newpt_before;
                        NaN;
                        newpt_after;
                        y[xi+1:end]
                    ]
                    colors = [colors[[1:xi;xi]]; NaN; colors[[xi+1;xi+1:end]]]
                    ylims!(ax_pa, 0, 360)
                    xi += 3
                elseif direction != direction_next && direction < 0
                    Δ =  (pa_model_t[yi,xi0+1]-2pi)-pa_model_t[yi,xi0]
                    x_midpoint = (x[xi] + x[xi+1])/2
                    newpt_before = y[xi] + Δ/2
                    newpt_after =  y[xi+1] - Δ/2
                    x = [
                        x[1:xi];
                        x_midpoint;
                        NaN;
                        x_midpoint;
                        x[xi+1:end]
                    ]
                    y = [
                        y[1:xi];
                        newpt_before;
                        NaN;
                        newpt_after;
                        y[xi+1:end]
                    ]
                    colors = [colors[[1:xi;xi]]; NaN; colors[[xi+1;xi+1:end]]]
                    ylims!(ax_pa, 0, 360)
                    xi += 3
                end
            end
            append!(allx, x)
            append!(ally, y)
            append!(allcolors, colors)
            push!(allx, NaN)
            push!(ally, NaN)
            push!(allcolors, 0.0)
        end
        lines!(ax_pa, allx, rad2deg.(ally);
            alpha=min.(1, 100 / length(ii)),
            color=allcolors,
            colorrange=(0,2pi),
            colormap=colormaps[planet_key]
        )
    end


    # Now overplot the data points, if any.
    # This is actually a bit of work since we might need
    # to convert RA and DEC +- uncertainties into SEP and PA.
    like_objs = []
    for planet in model.system.planets
        append!(like_objs, planet.observations)
    end
    append!(like_objs, model.system.observations)
    for like_obj in like_objs
        if nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood
            if hasproperty(like_obj.table, :sep)
                epoch = like_obj.table.epoch
                sep = like_obj.table.sep
                pa = like_obj.table.pa
                σ_sep = like_obj.table.σ_sep
                σ_pa = like_obj.table.σ_pa
            elseif hasproperty(like_obj.table, :ra)
                epoch = like_obj.table.epoch

                
                # We have to do a lot more work in this case
                # We need to transform the 2D error ellipse
                # and then find it's greatest extent along both 
                # axes.
                x = like_obj.table.ra
                y = like_obj.table.dec
                if hasproperty(like_obj.table, :cor)
                    cor = like_obj.table.cor
                else
                    cor = zeros(length(x))
                end
                σ₁ = like_obj.table.σ_ra
                σ₂ = like_obj.table.σ_dec

                # Go through each row one at a time
                sep = Float64[]
                pa = Float64[]
                σ_sep = Float64[]
                σ_pa = Float64[]
                for (x,y,σ₁,cor,σ₂) in zip(x,y,σ₁,cor,σ₂)
                    Σ = [
                        σ₁^2        cor*σ₁*σ₂
                        cor*σ₁*σ₂   σ₂^2
                    ]
                    vals, vecs = eigen(Σ) # should be real and sorted by real eigenvalue
                    length_major = sqrt(vals[2])
                    length_minor = sqrt(vals[1])
                    λ = vecs[:,2]
                    α = atan(λ[2],λ[1])
        
                    xvals = [
                        # Major axis
                        x - length_major*cos(α),
                        x + length_major*cos(α),
                        # Minor axis
                        x - length_minor*cos(α+π/2),
                        x + length_minor*cos(α+π/2),
                    ]
                    yvals = [
                        # Major axis
                        y - length_major*sin(α),
                        y + length_major*sin(α),
                        # Minor axis
                        y - length_minor*sin(α+π/2),
                        y + length_minor*sin(α+π/2),
                    ]

                    push!(sep, sqrt(x^2 + y^2))
                    push!(pa, rem2pi(-atan(y,x)+pi/2,RoundDown)) # TODO: plus offset and sign

                    seps = sqrt.(xvals.^2 .+ yvals.^2)
                    pas = rem2pi.(-atan.(yvals, xvals).+pi/2,RoundDown) # TODO: plus offset and sign
                    x1,x2 = extrema(seps)
                    y1,y2 = extrema(pas)
                    push!(σ_sep, (x2-x1)/2)
                    push!(σ_pa, (y2-y1)/2)
                end
            else
                error("invalid astrometry format")
            end
            Makie.errorbars!(
                ax_sep, epoch, sep, σ_sep;
                color=:black,
                linewidth=3,
            )
            Makie.scatter!(
                ax_sep, epoch, sep;
                color=:white,
                strokewidth=2,
                strokecolor=:black,
                markersize=8,
            )
            Makie.errorbars!(
                ax_pa, epoch, rad2deg.(pa), rad2deg.(σ_pa);
                color=:black,
                linewidth=3,
            )
            Makie.scatter!(
                ax_pa, epoch, rad2deg.(pa);
                color=:white,
                strokewidth=2,
                strokecolor=:black,
                markersize=8,
            )
        elseif nameof(typeof(like_obj)) in (:ImageLikelihood, :LogLikelihoodMap, :InterferometryLikelihood, :GRAVITYWideCPLikelihood)
            # In this case, put scatter points from the posterior
            
            for planet_key in keys(model.system.planets)
                orbs = Octofitter.construct_elements(results, planet_key, ii)
                # Now time-series
                sols = orbitsolve.(orbs, like_obj.table.epoch')
                Makie.scatter!(
                    ax_sep,
                    repeat(like_obj.table.epoch, inner=length(orbs)),
                    vec(projectedseparation.(sols));
                    color=:black,
                    markersize=4,
                )
                Makie.scatter!(
                    ax_pa,
                    repeat(like_obj.table.epoch, inner=length(orbs)),
                    vec(rad2deg.(rem2pi.(posangle.(sols),RoundDown)));
                    color=:black,
                    markersize=4,
                )
            end
        end
    end

     # The user can ask us to plot the position at a particular date
     if !isempty(mark_epochs_mjd)
        i = 0
        for epoch_mjd in mark_epochs_mjd
            i += 1
            for planet_key in keys(model.system.planets)
                orbs = Octofitter.construct_elements(results, planet_key, ii)
                color = Makie.wong_colors()[mod1(i,end)]
                sols = orbitsolve.(orbs, epoch_mjd)
                Makie.scatter!(
                    ax_sep,
                    fill(epoch_mjd, length(orbs)),
                    vec(projectedseparation.(sols));
                    color,
                    markersize=6,
                    strokewidth=1,
                    strokecolor=:black,
                    label = replace(string(mjd2date(epoch_mjd)), "T"=>" ") # TODO: better to use a format string
                )
                Makie.scatter!(
                    ax_pa,
                    fill(epoch_mjd, length(orbs)),
                    vec(rad2deg.(rem2pi.(posangle.(sols),RoundDown)));
                    color,
                    markersize=6,
                    strokewidth=1,
                    strokecolor=:black,
                    label = replace(string(mjd2date(epoch_mjd)), "T"=>" ") # TODO: better to use a format string
                )
            end
        end
        if colorbar
            axislegend(ax_sep)
        end
    end

    if colorbar
        if length(colormaps) == 1
            Colorbar(
                gs[1,2];
                colormap,
                label="mean anomaly →",
                colorrange=(0,2pi),
                ticks=(
                    [0,pi/2,pi,3pi/2,2pi],
                    ["0", "π/2", "π", "3π/2", "2π"]
                )
            )
        else
            col = 1
            for planet_key in keys(colormaps)
                col += 1
                Colorbar(
                    gs[1,col];
                    colormap=colormaps[planet_key],
                    label="mean anomaly ($planet_key) →",
                    colorrange=(0,2pi),
                    ticks=(
                        [0,pi/2,pi,3pi/2,2pi],
                        ["0", "π/2", "π", "3π/2", "2π"]
                    ),
                    ticksvisible=col - 1 == length(colormaps),
                    ticklabelsvisible=col - 1 == length(colormaps),
                )
            end
        end
    end
    return [ax_sep, ax_pa]

end
