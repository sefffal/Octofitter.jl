
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
    colormap=Makie.cgrad([Makie.wong_colors()[1], "#DDDDDD"]),
    colormap_instruments=Makie.cgrad(:Egypt,categorical=true),
    colormap_epochs=Makie.cgrad(:Lakota,categorical=true),
    colorbar=true,
    top_time_axis=true,
    bottom_time_axis=true,
    mark_epochs_mjd=Float64[],
    residuals,
    alpha,
    kwargs...
)
    gs = gridspec_or_fig

    # Detect if should use arcseconds instead of mas for plotting
    use_arcsec = false
    axis_mult = 1.0

    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(model, results, planet_key, ii)
        # Now time-series
        sols = orbitsolve.(orbs, ts')
        # sols_0 = orbitsolve.(orbs, epoch_0)
        try
            raoff(first(sols))
        catch
            continue
        end
        s = maximum(projectedseparation, sols)
        if s > 1500
            axis_mult = 1e-3
            use_arcsec = true
        end
    end
    
    date_pos, date_strs, xminorticks = _date_ticks(ts)
    row = 1
    all_axes = Axis[]
    ax_sep = Axis(
        gs[1, 1];
        ylabel=use_arcsec ? "sep [as] " : "sep [mas]",
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
    push!(all_axes, ax_sep)
    if residuals
        row+=1
        ax_sep_resid = Axis(
            gs[row, 1];
            # ylabel=use_arcsec ? "sep [as] " : "sep [mas]",
            xaxisposition=:top,
            xticks=(date_pos, date_strs),
            xgridvisible=false,
            ygridvisible=false,
            xminorticksvisible=false,
            xticksvisible=false,
            xticklabelsvisible=false,
            axis...
        )
        push!(all_axes, ax_sep_resid)
        Makie.rowgap!(gs, row-1, 0)
    end
    row+=1
    ax_pa = Axis(
        gs[row, 1];
        xlabel="MJD",
        ylabel="PA [°]",
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=bottom_time_axis,
        xticklabelsvisible=bottom_time_axis,
        xlabelvisible=bottom_time_axis,
        axis...
    )
    push!(all_axes, ax_pa)
    Makie.rowgap!(gs, row-1, 10)
    xlims!(ax_sep, extrema(ts))
    xlims!(ax_pa, extrema(ts))
    if residuals
        row+=1
        ax_pa_resid = Axis(
            gs[row, 1];
            xlabel="MJD",
            ylabel="PA [°]",
            xgridvisible=false,
            ygridvisible=false,
            xticksvisible=bottom_time_axis,
            xticklabelsvisible=bottom_time_axis,
            xlabelvisible=bottom_time_axis,
            axis...
        )
        Makie.rowgap!(gs, row-1, 0)
        push!(all_axes, ax_pa)
    end
    linkxaxes!(all_axes...)


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
        orbs = Octofitter.construct_elements(model, results, planet_key, ii)

        sols = orbitsolve.(orbs, ts')

        # Some orbit models provide position angles but not projected separation
        # TODO: we had to disable the ability to plot orbits giving PA but 
        # not sep and vise-versa to add the other-planet perturbations
        # We should find a way to reenable it
        # sep_model_t = try
        #     projectedseparation.(sols)
        # catch
        #     fill(NaN, length(sols))
        # end
        
        # # Other orbit models e.g. RadialVelocityOrbit provide neither
        # pa_model_t = try
        #     rem2pi.(posangle.(sols), RoundDown)
        # catch
        #     fill(NaN, length(sols))
        # end
        color_model_t = rem2pi.(
            eccanom.(sols), RoundDown)

        # Account for planet-star interactions from interior planets
        ra_host_perturbation = zeros(size(sols), )
        dec_host_perturbation = zeros(size(sols), )
        for planet_key′ in keys(model.system.planets)
            if !haskey(results, Symbol("$(planet_key′)_mass"))
                continue
            end

            other_planet_mass = results["$(planet_key′)_mass"][ii]
            orbit_other = Octofitter.construct_elements(model, results, planet_key′, ii)
            try
                raoff.(orbit_other, 0)
            catch
                continue
            end

            # Only account for interior planets
            mask = semimajoraxis.(orbit_other) .< semimajoraxis.(orbs)
            sols′ = orbitsolve.(orbit_other, ts')
            
            ra_host_perturbation .+= mask .* raoff.(sols′, other_planet_mass.*Octofitter.mjup2msol)
            dec_host_perturbation .+= mask .* decoff.(sols′, other_planet_mass.*Octofitter.mjup2msol)
        end

        ra_model = (raoff.(sols) .- ra_host_perturbation)
        dec_model = (decoff.(sols) .- dec_host_perturbation)
        sep_model_t = hypot.(ra_model, dec_model)
        pa_model_t = rem2pi.(atan.(ra_model, dec_model), RoundDown)
        
        if any(isfinite, sep_model_t)
            lines!(ax_sep,
                concat_with_nan(ts' .+ 0 .* sep_model_t),
                concat_with_nan(sep_model_t).*axis_mult;
                alpha,
                color=concat_with_nan(color_model_t .+ 0 .* ii),
                colorrange=(0,2pi),
                colormap=colormaps[planet_key],
                rasterize=4,
            )
        end
        
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
            alpha,
            color=allcolors,
            colorrange=(0,2pi),
            colormap=colormaps[planet_key],
            rasterize=4,
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

    # Colour data based on the instrument name
    rel_astrom_likes = filter(like_objs) do like_obj
        nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood 
    end
    rel_astrom_names = sort(unique(getproperty.(rel_astrom_likes, :name)))
    n_rel_astrom = length(rel_astrom_names)
    
    for like_obj in like_objs
        if nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood
            i_like_obj = findfirst(==(likelihoodname(like_obj)), rel_astrom_names)
            if hasproperty(like_obj.table, :sep)
                epoch = like_obj.table.epoch
                sep = like_obj.table.sep
                pa = rem2pi.(like_obj.table.pa,RoundDown)
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
            if n_rel_astrom == 1
                color = :white
            else
                color = colormap_instruments[mod1(i_like_obj,end)]
            end
            Makie.errorbars!(
                ax_sep, epoch, sep .* axis_mult, σ_sep.*axis_mult;
                color=:black,
                linewidth=3,
            )
            Makie.scatter!(
                ax_sep, epoch, sep .* axis_mult;
                color,
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
                color,
                strokewidth=2,
                strokecolor=:black,
                markersize=8,
            )
        elseif nameof(typeof(like_obj)) in (:ImageLikelihood, :LogLikelihoodMap, :InterferometryLikelihood, :GRAVITYWideCPLikelihood)
            # In this case, put scatter points from the posterior
            
            for planet_key in keys(model.system.planets)
                orbs = Octofitter.construct_elements(model, results, planet_key, ii)
                # Now time-series
                sols = orbitsolve.(orbs, like_obj.table.epoch')
                Makie.scatter!(
                    ax_sep,
                    repeat(like_obj.table.epoch, inner=length(orbs)),
                    vec(projectedseparation.(sols)).*axis_mult;
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
                orbs = Octofitter.construct_elements(model, results, planet_key, ii)
                color = colormap_epochs[mod1(i,end)]
                sols = orbitsolve.(orbs, epoch_mjd)
                Makie.scatter!(
                    ax_sep,
                    fill(epoch_mjd, length(orbs)),
                    vec(projectedseparation.(sols)).*axis_mult;
                    color,
                    markersize=6,
                    strokewidth=1,
                    strokecolor=:black,
                    marker=:rect,
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
                    marker=:rect,
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
