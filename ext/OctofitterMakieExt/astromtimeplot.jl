
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
    kwargs...
)
    gs = gridspec_or_fig


    # Find the minimum time step size
    periods = Float64[]
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, :)
        append!(periods, period.(orbs))
    end
    min_period = minimum(periods)
    med_period = median(periods)

    # Always try and show at least one full period
    # TODO: decide t_start and t_end based on epochs of data if astrometry 
    # available, otherwise other data, otherwise arbitrary
    t_start = years2mjd(1990)
    t_stop = years2mjd(2020)
    t_cur = t_stop-t_start
    if t_cur < med_period
        t_extend = med_period - t_cur
        t_start -= t_extend/2
        t_stop += t_extend/2
    end

    ts = range(t_start, t_stop, step=min_period / 150)
    if length(ts) > 1500
        # Things can get very slow if we have some very short period 
        # solutions
        ts = range(t_start, t_stop, length=1500)
    end
    # Make sure we include the approx. data epochs
    ts = sort([ts; years2mjd(1990.25); years2mjd((1990.25 + 2016) / 2); years2mjd(2016)])
    
    date_pos, date_strs = _date_ticks(ts)
    ax_sep = Axis(
        gs[1, 1];
        ylabel="sep [mas]",
        xaxisposition=:top,
        xticks=(date_pos, date_strs),
        xgridvisible=false,
        ygridvisible=false,
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

    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        # Now time-series
        sols = orbitsolve.(orbs, ts')
        
        # TODO: data overplot...
        sep_model_t = projectedseparation.(sols)
        pa_model_t = rem2pi.(posangle.(sols), RoundDown)
        color_model_t = rem2pi.(
            eccanom.(sols), RoundDown)

        # # Try to unwrap position angle if it extends slightly past 180 deg,
        # # but not if there are 2 or more cycles.
        # pa_model_t_unwrapped = copy(pa_model_t)
        # unwrap!(pa_model_t_unwrapped, 2pi)

        # # Don't use the unwrap if it's cycling a bunch of times
        # a,b = extrema(pa_model_t_unwrapped)
        # if abs(a-b) < 4pi
        #     println("unwrapping")
        #     pa_model_t = pa_model_t_unwrapped
        #     # TO unwrap data, find closest point in time and add the difference
        #     # between wrapped and unwrapped
        # end
        
        lines!(ax_sep,
            concat_with_nan(ts' .+ 0 .* sep_model_t),
            concat_with_nan(sep_model_t);
            alpha=min.(1, 100 / length(ii)),
            color=concat_with_nan(color_model_t .+ 0 .* ii),
            colorrange=(0,2pi),
            colormap
        )

        # TODO: add in direction check, orbits are always increasing or
        # decreasing in PA, but sometimes do a big jump
        delta = 2.0
        allx = Float64[]
        ally = Float64[]
        allcolors = Float64[]
        for yi in axes(pa_model_t, 1)
            x = ts
            y = pa_model_t[yi,:]
            colors = color_model_t[yi,:]
            for xi in 2:(size(pa_model_t, 2)-2)
                if (
                    pa_model_t[yi,xi] > 2pi - delta &&
                    pa_model_t[yi,xi+1] < delta)
                    x = [x[1:xi]; NaN; x[xi+1:end]]
                    y = [y[1:xi-1]; 4pi; NaN; 0; y[xi+2:end]]
                    colors = [colors[1:xi]; NaN; colors[xi+1:end]]
                    ylims!(ax_pa, 0, 360)
                elseif (
                    pa_model_t[yi,xi] < delta &&
                    pa_model_t[yi,xi+1] > 2pi - delta)
                    x = [x[1:xi]; NaN; x[xi+1:end]]
                    y = [y[1:xi-1]; 0; NaN; 4pi; y[xi+2:end]]
                    colors = [colors[1:xi]; NaN; colors[xi+1:end]]
                    ylims!(ax_pa, 0, 360)
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
            colormap
        )
    end


    # Now overplot the data points, if any.
    # This is actually a bit of work since we might need
    # to convert RA and DEC +- uncertainties into SEP and PA.
     for planet in model.system.planets
        for like_obj in planet.observations
            if nameof(typeof(like_obj)) != :PlanetRelAstromLikelihood
                continue
            end
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
        end
    end

    if colorbar
        Colorbar(
            gs[1:2,2];
            colormap,
            label="orbit fraction past periastron",
            colorrange=(0,1)
        )
    end


end
