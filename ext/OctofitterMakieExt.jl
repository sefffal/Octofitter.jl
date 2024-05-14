module OctofitterMakieExt
using Octofitter
using Makie
using LinearAlgebra
using Dates
using Printf
using Statistics


function Octofitter.octoplot(
    ::Val{2},
    model::Octofitter.LogDensityModel,
    results::Chains;
    fname="$(model.system.name)-plot-grid.png",
    show_astrom=nothing,
    show_astrom_time=nothing,
    show_hgca=nothing,
    show_mass=false,
    show_rv=false,
    figure=(;),
    # If less than 500 samples, just show all of them
    N  = min(size(results, 1)*size(results, 3), 500),
    # If showing all samples, include each sample once.
    # Otherwise, sample randomly with replacement
    ii = (
        N == size(results, 1)*size(results, 3) ? 
        (1:size(results, 1)*size(results, 3)) :
        rand(1:size(results, 1)*size(results, 3),N)
    )
    # The user can of course just override the above directly.
)

    # Auto-detect if we should include a given plot
    if isnothing(show_astrom)
        show_astrom = false
        for planet in model.system.planets
            show_astrom |= 
                Octofitter.orbittype(planet) <: Visual{KepOrbit} || 
                Octofitter.orbittype(planet) <: AbsoluteVisual{KepOrbit} || 
                Octofitter.orbittype(planet) <: ThieleInnesOrbit
        end
    end

    if isnothing(show_astrom_time)
        show_astrom_time = false
        for planet in model.system.planets
            show_astrom_time |= 
                Octofitter.orbittype(planet) <: Visual{KepOrbit} || 
                Octofitter.orbittype(planet) <: AbsoluteVisual{KepOrbit} || 
                Octofitter.orbittype(planet) <: ThieleInnesOrbit
        end
    end

    if isnothing(show_hgca)
        show_hgca = false
        for like_obj in model.system.observations
            if like_obj isa HGCALikelihood
                show_hgca = true
            end
        end
    end

    if isnothing(show_mass)
        show_mass = false
        for planet_key in keys(model.system.planets)
            show_mass |= haskey(results, Symbol("$(planet_key)_mass"))
        end
    end

    fig = Figure(;
        figure...
    )
    # Show a colorbar for only the first sub-plot, and don't repeat it.
    colorbar = true
    top_time_axis = true
    item = 0
    cols = 1
    if show_astrom
        item += 1
        col = mod1(item, cols)
        row = cld(item, cols)
        gl = GridLayout(
            fig[row,col],
            width=500,
            height=400,
        )
        Octofitter.astromplot!(gl, model, results; ii, colorbar)
        colorbar = false
    end


    if show_astrom_time
        item += 1
        col = mod1(item, cols)
        row = cld(item, cols)
        gl = GridLayout(
            fig[row,col],
            width=500,
            height=300,
        )
        bottom_time_axis = !(show_hgca || show_rv)
        astromtimeplot!(gl, model, results; ii, colorbar, top_time_axis, bottom_time_axis)
        colorbar = false
        top_time_axis = false
    end


    if show_rv
        item += 1
        col = mod1(item, cols)
        row = cld(item, cols)
        gl = GridLayout(
            fig[row,col],
            width=500,
            height=135,
        )
        bottom_time_axis = !show_hgca
        rvtimeplot!(gl, model, results; ii, colorbar, top_time_axis, bottom_time_axis)
        colorbar = false
        top_time_axis = false
    end

    if show_hgca
        item += 1
        col = mod1(item, cols)
        row = cld(item, cols)
        gl = GridLayout(
            fig[row,col],
            width=500,
            height=480,
        )
        Octofitter.hgcaplot!(gl, model, results; ii, colorbar, top_time_axis)
        colorbar = false
        top_time_axis = false
    end


    if show_mass
        item += 1
        col = mod1(item, cols)
        row = cld(item, cols)
        gl = GridLayout(
            fig[row,col],
            width=500,
            height=400,
        )
        Octofitter.masspostplot!(gl, model, results;)
    end

    # hgcaplot
    Makie.resize_to_layout!(fig)

    save(fname, fig)

    return fig
end


#################################################
# Astrometry plot

function Octofitter.astromplot(
    model,
    results,
    fname="$(model.system.name)-astromplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(; size=(700,600), figure...)
    Octofitter.astromplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function Octofitter.astromplot!(
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
    kwargs...
)
    gs = gridspec_or_fig

    ax = Axis(
        gs[1, 1];
        autolimitaspect=1,
        xreversed=true,
        xlabel="Δra [mas]",
        ylabel="Δdec [mas]",
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )

    # Start by plotting the orbits
    EAs = range(0, 2pi, length=150)

    # Find earliest epoch
    epoch_0 = mjd("2020")
    for planet in model.system.planets
        for like_obj in planet.observations
            x = Float64[]
            y = Float64[]
            xs = Float64[]
            ys = Float64[]
            if nameof(typeof(like_obj)) == :PlanetRelAstromLikelihood
                epoch_0 = min(epoch_0, minimum(like_obj.table.epoch))
            end
        end
    end

    if colorbar 
        Colorbar(
            gs[1,2];
            colormap,
            label="orbit fraction past periastron",
            colorrange=(0,1)
        )
    end

    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)

        # Draws from the posterior
        sols = orbitsolve_eccanom.(orbs, EAs')
        # sols_0 = orbitsolve.(orbs, epoch_0)

        lines!(ax,
            concat_with_nan(raoff.(sols)),
            concat_with_nan(decoff.(sols));
            color=concat_with_nan(
                # rem2pi.(EAs' .- eccanom.(sols_0), RoundDown) .+ 0 .* ii
                rem2pi.(meananom.(sols), RoundDown) .+ 0 .* ii
            ),
            alpha=min.(1, 100 / length(ii)),
            colormap
        )
    end

    # Now over plot any astrometry
    for planet in model.system.planets
        for like_obj in planet.observations
            if nameof(typeof(like_obj)) != :PlanetRelAstromLikelihood
                continue
            end

            x = Float64[]
            y = Float64[]
            xs = Float64[]
            ys = Float64[]

            if hasproperty(like_obj.table, :sep)
                # Plot astrometry point
                x = @. like_obj.table.sep * sin(like_obj.table.pa)
                y = @. like_obj.table.sep * cos(like_obj.table.pa)
                # And overplot with uncertainty ellipse
                if hasproperty(like_obj.table, :cor)
                    cor = like_obj.table.cor
                else
                    cor = 0
                end
                x = like_obj.table.sep
                y = pi / 2 .- like_obj.table.pa
                σ₁ = like_obj.table.σ_sep
                σ₂ = like_obj.table.σ_pa
                error_ellipses = broadcast(x, y, σ₁, cor, σ₂) do x, y, σ₁, cor, σ₂
                    Σ = [
                        σ₁^2 cor*σ₁*σ₂
                        cor*σ₁*σ₂ σ₂^2
                    ]
                    vals, vecs = eigen(Σ) # should be real and sorted by real eigenvalue
                    length_major = sqrt(vals[2])
                    length_minor = sqrt(vals[1])
                    λ = vecs[:, 2]
                    α = atan(λ[2], λ[1])

                    xvals = [
                        # Major axis
                        x - length_major * cos(α),
                        x,
                        x + length_major * cos(α),
                        NaN,
                        # Minor axis
                        x - length_minor * cos(α + π / 2),
                        x,
                        x + length_minor * cos(α + π / 2),
                        NaN,
                    ]
                    yvals = [
                        # Major axis
                        y - length_major * sin(α),
                        y,
                        y + length_major * sin(α),
                        NaN,
                        # Minor axis
                        y - length_minor * sin(α + π / 2),
                        y,
                        y + length_minor * sin(α + π / 2),
                        NaN,
                    ]
                    xvals, yvals
                end
                rs = Base.vcat(getindex.(error_ellipses, 1)...)
                thetas = Base.vcat(getindex.(error_ellipses, 2)...)
                xs = rs .* cos.(thetas)
                ys = rs .* sin.(thetas)
                x, y = x .* cos.(y), x .* sin.(y)
            elseif hasproperty(like_obj.table, :ra)
                x = like_obj.table.ra
                y = like_obj.table.dec

                if hasproperty(like_obj.table, :cor)
                    cor = like_obj.table.cor
                else
                    cor = 0
                end

                σ₁ = like_obj.table.σ_ra
                σ₂ = like_obj.table.σ_dec

                error_ellipses = broadcast(x, y, σ₁, cor, σ₂) do x, y, σ₁, cor, σ₂
                    Σ = [
                        σ₁^2 cor*σ₁*σ₂
                        cor*σ₁*σ₂ σ₂^2
                    ]
                    vals, vecs = eigen(Σ) # should be real and sorted by real eigenvalue
                    length_major = sqrt(vals[2])
                    length_minor = sqrt(vals[1])
                    λ = vecs[:, 2]
                    α = atan(λ[2], λ[1])

                    xvals = [
                        # Major axis
                        x - length_major * cos(α),
                        x + length_major * cos(α),
                        NaN,
                        # Minor axis
                        x - length_minor * cos(α + π / 2),
                        x + length_minor * cos(α + π / 2),
                        NaN,
                    ]
                    yvals = [
                        # Major axis
                        y - length_major * sin(α),
                        y + length_major * sin(α),
                        NaN,
                        # Minor axis
                        y - length_minor * sin(α + π / 2),
                        y + length_minor * sin(α + π / 2),
                        NaN,
                    ]
                    xvals, yvals
                end
                xs = Base.vcat(getindex.(error_ellipses, 1)...)
                ys = Base.vcat(getindex.(error_ellipses, 2)...)
            end
            Makie.lines!(ax, xs, ys, color=:white, linewidth=3)
            Makie.lines!(ax, xs, ys, color=:black, linewidth=2)
            Makie.scatter!(
                ax,
                vec(x),
                vec(y),
                color=:white,
                strokewidth=2,
                strokecolor=:black,
                markersize=8,
            )
        end
    end
    scatter!(ax, [0],[0],marker='⭐', markersize=30, color=:black)
end


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



##################################################
# RV vs time plot
# Note that the "full" plot requires to use
# rvpostplot(). It shows the jitter, residuals,
# and phase folded curves properly for a single sample.
# This just gives a global view.


function rvtimeplot(
    model,
    results,
    fname="$(model.system.name)-rvplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700,600),
        figure...
    )
    rvtimeplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function rvtimeplot!(
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
    ax= Axis(
        gs[1, 1];
        ylabel="rv [m/s]",
        xaxisposition=:top,
        xticks=(date_pos, date_strs),
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=top_time_axis,
        xticklabelsvisible=top_time_axis,
        axis...
    )
    ax_secondary = Axis(
        gs[1, 1];
        xlabel="MJD",
        xgridvisible=false,
        ygridvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false,
        ylabelvisible=false,
        xticksvisible=bottom_time_axis,
        xticklabelsvisible=bottom_time_axis,
        xlabelvisible=bottom_time_axis,
        axis...
    )
    linkxaxes!(ax, ax_secondary)
    xlims!(ax, extrema(ts))
    xlims!(ax_secondary, extrema(ts))



    rv_model_t = zeros(length(ii), length(ts))
    color_model_t = zeros(length(ii), length(ts))
    for (i_planet, planet_key) in enumerate(keys(model.system.planets))
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        # Draws from the posterior
        mass = results["$(planet_key)_mass"][ii] .* Octofitter.mjup2msol

        # Now time-series
        sols = orbitsolve.(orbs, ts')
  
        rv_model_t .+= radvel.(sols, mass) 

              
        # We subtract off the "average" instrument RV offset here.
        # This varies by instrument, and isn't centred around a net hierarchical value.
        # For now just take the average of the different zero points per row.
        if i_planet == 1 && haskey(results, :rv0_1)
            rv0s = []
            for i in 1:10
                k = Symbol("rv0_$i")
                if haskey(results, k)
                    push!(rv0s, results[k][ii])
                end
            end
            ave_rv0 = mean(stack(rv0s),dims=2)
            rv_model_t .-= ave_rv0
        end
            

        color_model_t .+= rem2pi.(
            meananom.(sols), RoundDown) .+ 0 .* ii
    end

    lines!(ax,
        concat_with_nan(ts' .+ 0 .* rv_model_t),
        concat_with_nan(rv_model_t);
        alpha=min.(1, 100 / length(ii)),
        color=concat_with_nan(color_model_t),
        colorrange=(0,2pi),
        colormap
    )



    # # Now overplot the data points, if any.
    # for like_obj in model.system.observations
    #     if nameof(typeof(like_obj)) != :StarAbsoluteRVLikelihood
    #         continue
    #     end
    #     epoch = vec(like_obj.table.epoch)
    #     rv = vec(like_obj.table.rv)
    #     σ_rv = vec(like_obj.table.σ_rv)
    #     Makie.errorbars!(
    #         ax, epoch, rv, σ_rv;
    #         color=:black,
    #         linewidth=3,
    #     )
    #     Makie.scatter!(
    #         ax, epoch, rv;
    #         color=:white,
    #         strokewidth=0.5,
    #         strokecolor=:black,
    #         markersize=2,
    #     )
    # end

    if colorbar
        Colorbar(
            gs[1,2];
            colormap,
            label="orbit fraction past periastron",
            colorrange=(0,1)
        )
    end
end


##################################################
# HGCA Plot
const pmra_label = rich("μ", subscript("α*"), " [mas/yr]")
const pmdec_label = rich("μ", subscript("δ"), " [mas/yr]")
function Octofitter.hgcaplot(
    model,
    results,
    fname="$(model.system.name)-hgcaplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700,600),
        figure...
    )
    Octofitter.hgcaplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end
function Octofitter.hgcaplot!(
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

    # Start by plotting the orbits
    # EAs = range(0, 2pi, length=150)
    # ts = range(years2mjd(1990), years2mjd(2020), step=min_period / 150)
    # Always try and show at least one full period
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
    # Find earliest epoch
    # epoch_0 = years2mjd(1991.25)

    date_pos, date_strs = _date_ticks(ts)
    ax_velra = Axis(
        gs[1, 1:3];
        ylabel=pmra_label,
        xaxisposition=:top,
        xticks=(date_pos, date_strs),
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=top_time_axis,
        xticklabelsvisible=top_time_axis,
        axis...
    )
    ax_veldec = Axis(
        gs[2, 1:3];
        xlabel="MJD",
        ylabel=pmdec_label,
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=bottom_time_axis,
        xticklabelsvisible=bottom_time_axis,
        xlabelvisible=bottom_time_axis,
        axis...
    )
    linkxaxes!(ax_velra, ax_veldec)

    ax_dat1 = Axis(
        gs[3,1],
        xlabel=pmra_label,
        ylabel=pmdec_label,
        autolimitaspect=1.0,
        title="H",
        xticklabelrotation=pi/4,
        xgridvisible=false,
        ygridvisible=false,
    )
    ax_dat2 = Axis(
        gs[3,2],
        xlabel=pmra_label,
        ylabel=pmdec_label,
        autolimitaspect=1.0,
        title="G-H",
        xticklabelrotation=pi/4,
        xgridvisible=false,
        ygridvisible=false,
        ylabelvisible=false,
    )
    ax_dat3 = Axis(
        gs[3,3],
        xlabel=pmra_label,
        ylabel=pmdec_label,
        autolimitaspect=1.0,
        title="G",
        xticklabelrotation=pi/4,
        xgridvisible=false,
        ygridvisible=false,
        ylabelvisible=false,
    )
    # linkxaxes!(ax_dat1,ax_dat2,ax_dat3)
    # linkyaxes!(ax_dat1,ax_dat2,ax_dat3)


    xlims!(ax_velra, extrema(ts))
    xlims!(ax_veldec, extrema(ts))

    pmra_model_t = zeros(length(ii), length(ts))
    pmdec_model_t = zeros(length(ii), length(ts))
    color_model_t = zeros(length(ii), length(ts))
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, ii)
        # Draws from the posterior
        mass = results["$(planet_key)_mass"][ii] .* Octofitter.mjup2msol

        # Now time-series
        sols = orbitsolve.(orbs, ts')
        # Can we use the existing simulator for this please?
        pmra_model_t .+= pmra.(sols, mass) .+ results[:pmra][ii]
        pmdec_model_t .+= pmdec.(sols, mass) .+ results[:pmdec][ii]
        color_model_t .+= rem2pi.(
            meananom.(sols), RoundDown) .+ 0 .* ii
    end
    if colorbar
        Colorbar(
            gs[1:2,4];
            colormap,
            label="orbit fraction past periastron",
            colorrange=(0,1)
        )
    end
    lines!(ax_velra,
        concat_with_nan(ts' .+ 0 .* pmra_model_t),
        concat_with_nan(pmra_model_t);
        alpha=min.(1, 100 / length(ii)),
        color=concat_with_nan(color_model_t),
        colorrange=(0,2pi),
        colormap
    )
    lines!(ax_veldec,
        concat_with_nan(ts' .+ 0 .* pmdec_model_t),
        concat_with_nan(pmdec_model_t);
        alpha=min.(1, 100 / length(ii)),
        color=concat_with_nan(color_model_t),
        colorrange=(0,2pi),
        colormap
    )

    # Now over plot any astrometry
    like_objs = filter(model.system.observations) do like_obj
        nameof(typeof(like_obj)) == :HGCALikelihood
    end
    if isempty(like_objs)
        return
    end
    hgca_like = only(like_objs)

    tx = years2mjd.([
        hgca_like.table.epoch_ra_hip
        (hgca_like.table.epoch_ra_hip + hgca_like.table.epoch_ra_gaia) / 2
        hgca_like.table.epoch_ra_gaia
    ])
    ty = years2mjd.([
        hgca_like.table.epoch_dec_hip
        (hgca_like.table.epoch_dec_hip + hgca_like.table.epoch_dec_gaia) / 2
        hgca_like.table.epoch_dec_gaia
    ])
    x = [
        hgca_like.table.pmra_hip
        hgca_like.table.pmra_hg
        hgca_like.table.pmra_gaia
    ]
    y = [
        hgca_like.table.pmdec_hip
        hgca_like.table.pmdec_hg
        hgca_like.table.pmdec_gaia
    ]

    cor = [
        hgca_like.table.pmra_pmdec_hip
        hgca_like.table.pmra_pmdec_hg
        hgca_like.table.pmra_pmdec_gaia
    ]

    σ₁ = [
        hgca_like.table.pmra_hip_error
        hgca_like.table.pmra_hg_error
        hgca_like.table.pmra_gaia_error
    ]
    σ₂ = [
        hgca_like.table.pmdec_hip_error
        hgca_like.table.pmdec_hg_error
        hgca_like.table.pmdec_gaia_error
    ]

    # # One more model-plot: add scatter points to existing lines at the 
    # # data epochs, colored correctly
    # scatx = vec(stack(map(tx) do tx_i
    #     pmra_model_t[:, argmin(abs.(ts .- tx_i))]
    # end))
    # scaty = vec(stack(map(ty) do ty_i
    #     pmdec_model_t[:, argmin(abs.(ts .- ty_i))]
    # end))
    # scatter!(
    #     ax_vel2d,
    #     scatx,
    #     scaty,
    #     color=repeat(Makie.wong_colors()[1:length(x)], outer=length(ii)),
    #     markersize=5,
    # )

    error_ellipses = broadcast(x, y, σ₁, cor, σ₂) do x, y, σ₁, cor, σ₂
        Σ = [
            σ₁^2 cor*σ₁*σ₂
            cor*σ₁*σ₂ σ₂^2
        ]
        vals, vecs = eigen(Σ) # should be real and sorted by real eigenvalue
        length_major = sqrt(vals[2])
        length_minor = sqrt(vals[1])
        λ = vecs[:, 2]
        α = atan(λ[2], λ[1])

        xvals = [
            # Major axis
            x - length_major * cos(α),
            x + length_major * cos(α),
            NaN,
            # Minor axis
            x - length_minor * cos(α + π / 2),
            x + length_minor * cos(α + π / 2),
            NaN,
        ]
        yvals = [
            # Major axis
            y - length_major * sin(α),
            y + length_major * sin(α),
            NaN,
            # Minor axis
            y - length_minor * sin(α + π / 2),
            y + length_minor * sin(α + π / 2),
            NaN,
        ]
        xvals, yvals
    end

    # 1D plots: stroke twice for contrast
    Makie.errorbars!(
        ax_velra,
        tx,
        x,
        σ₁,
        color=:black,
    )
    Makie.errorbars!(
        ax_veldec,
        ty,
        y,
        σ₂,
        color=:black,
    )
    Makie.scatter!(
        ax_velra,
        tx,
        x,
        σ₁,
        color=Makie.wong_colors()[1:length(x)],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )
    Makie.scatter!(
        ax_veldec,
        ty,
        y,
        σ₂,
        color=Makie.wong_colors()[1:length(x)],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )

    colsize!(gs, 1, Auto(1//3))
    colsize!(gs, 2, Auto(1//3))
    colsize!(gs, 3, Auto(1//3))
    rowsize!(gs, 3, Aspect(3,1.0))

    ## Model

    θ_systems_from_chain = Octofitter.mcmcchain2result(model, results[ii])
    for (θ_system, i) in zip(θ_systems_from_chain,ii)
        orbits = map(keys(model.system.planets)) do planet_key
            Octofitter.construct_elements(results, planet_key, i)
        end
        hgca_like_sim = Octofitter.generate_from_params(hgca_like, θ_system, orbits)
        # HIP Epoch
        Makie.scatter!(
            ax_dat1,
            [only(hgca_like_sim.table.pmra_hip)],
            [only(hgca_like_sim.table.pmdec_hip)],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [years2mjd(only(hgca_like_sim.table.epoch_ra_hip))],
            [only(hgca_like_sim.table.pmra_hip)],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [years2mjd(only(hgca_like_sim.table.epoch_dec_hip))],
            [only(hgca_like_sim.table.pmdec_hip)],
            color=:black,
            markersize=3,
        )
        # HG Epoch
        Makie.scatter!(
            ax_dat2,
            [only(hgca_like_sim.table.pmra_hg)],
            [only(hgca_like_sim.table.pmdec_hg)],
            color=:black,
            markersize=2,
        )
        hg_ra_epoch = (only(hgca_like_sim.table.epoch_ra_hip) +
        only(hgca_like_sim.table.epoch_ra_gaia))/2
        hg_dec_epoch = (only(hgca_like_sim.table.epoch_dec_hip) +
        only(hgca_like_sim.table.epoch_dec_gaia))/2
        Makie.scatter!(
            ax_velra,
            [years2mjd(hg_ra_epoch)],
            [only(hgca_like_sim.table.pmra_hg)],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [years2mjd(hg_dec_epoch)],
            [only(hgca_like_sim.table.pmdec_hg)],
            color=:black,
            markersize=3,
        )
        # GAIA Epoch
        Makie.scatter!(
            ax_dat3,
            [only(hgca_like_sim.table.pmra_gaia)],
            [only(hgca_like_sim.table.pmdec_gaia)],
            color=:black,
            markersize=2,
        )
        Makie.scatter!(
            ax_velra,
            [years2mjd(only(hgca_like_sim.table.epoch_ra_gaia))],
            [only(hgca_like_sim.table.pmra_gaia)],
            color=:black,
            markersize=3,
        )
        Makie.scatter!(
            ax_veldec,
            [years2mjd(only(hgca_like_sim.table.epoch_dec_gaia))],
            [only(hgca_like_sim.table.pmdec_gaia)],
            color=:black,
            markersize=3,
        )
    end

    ## Data 
    Makie.lines!(
        ax_dat1,
        error_ellipses[1][1],
        error_ellipses[1][2],
        color=Makie.wong_colors()[1],
        linewidth=3.5,
    )
    Makie.scatter!(
        ax_dat1, x[1], y[1],
        color=Makie.wong_colors()[1],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )
    Makie.lines!(
        ax_dat2,
        error_ellipses[2][1],
        error_ellipses[2][2],
        color=Makie.wong_colors()[2],
        linewidth=3.5,
    )
    Makie.scatter!(
        ax_dat2, x[2], y[2],
        color=Makie.wong_colors()[2],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )
    Makie.lines!(
        ax_dat3,
        error_ellipses[3][1],
        error_ellipses[3][2],
        color=Makie.wong_colors()[3],
        linewidth=3.5,
    )
    Makie.scatter!(
        ax_dat3, x[3], y[3],
        color=Makie.wong_colors()[3],
        markersize=10,
        strokewidth=1.5,
        strokecolor=:black
    )
end



##################################################
# Mass vs. semi-major axis plot
function Octofitter.masspostplot(
    model,
    results,
    fname="$(model.system.name)-masspostplot!.png",
    args...;
    figure=(;),
    # massloglog=nothing,
    kwargs...
)

    # # auto-determine if we should use a log-log plot
    # if isnothing(massloglog)

    # end

    fig = Figure(;
        size=(500,300),
        figure...,
    )
    Octofitter.masspostplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end

const mass_mjup_label = rich("mass [M", subscript("jup"), "]")
function Octofitter.masspostplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains;
    axis=(;),
    kwargs...
)
    gs = gridspec_or_fig

    ax_hist = Axis(gs[1:2,1];
        # ylabel=rich("mass [M", subscript("jup"), "]"),
        xlabel=mass_mjup_label,
        xgridvisible=false,
        ygridvisible=false,
    )
    ylims!(ax_hist, low=0)
    ax_scat_sma = Axis(gs[1,2];
        ylabel=mass_mjup_label,
        xlabel="sma [AU]",
        # xscale=log10,
        # xticks=2 .^ (0:12),
        # yticks=2 .^ (0:12),
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )
    ax_scat_ecc = Axis(gs[2,2];
        ylabel=mass_mjup_label,
        xlabel="eccentricity",
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )
    xlims!(ax_scat_ecc, 0, 1)
    cred_intervals = []
    for planet_key in keys(model.system.planets)
        mk = Symbol("$(planet_key)_mass")
        if !haskey(results, mk)
            continue
        end
        sma = vec(results["$(planet_key)_a"])
        mass = vec(results[mk])
        ecc = vec(results["$(planet_key)_e"])
        stephist!(
            ax_hist,
            mass,
            linewidth=3
        )
        scatter!(ax_scat_sma, sma, mass;
            markersize=2,
        )
        scatter!(ax_scat_ecc, ecc, mass;
            markersize=2,
        )
        low,mid,high = quantile(mass, (0.16, 0.5, 0.84))
        label = margin_confidence_default_formatter(mid-low,mid,high-mid)
        push!(
            cred_intervals,
            Makie.latexstring("$(planet_key)_mass = "*label)
        )
    end
    if !isempty(cred_intervals)
        Label(
            gs[0,1],
            reduce(*, cred_intervals),
            tellwidth=false
        )
    end


end



# Calculate a list of nicely spaced date ticks given a vector
# of MJD.
# Returns the locations of the ticks in MJD and the formatted strings.
function _date_ticks(ts)

    # For secondary date axis on top
    date_start = mjd2date(ts[begin])
    # Remove any dates before 0 BC, the Julia datetime parser
    # fails on negative years...
    date_start = max.(Date("0000-01-01"), date_start)

    date_end = mjd2date(ts[end])
    date_start = Date(Dates.year(date_start), month(date_start))

    date_end = Date(Dates.year(date_end), month(date_end))
    dates = range(date_start, date_end, step=Year(1))
    dates_str = string.(year.(dates))
    if length(dates) == 1
        dates = range(date_start, date_end, step=Month(1))
        dates_str = map(d->string(Dates.year(d),"-",lpad(month(d),2,'0')),dates)
    else
        year_step = 1
        while length(dates) > 8
            year_step += 1
            dates = range(date_start, date_end, step=Year(year_step))
        end
        dates_str = string.(year.(dates))
    end
    return (mjd.(string.(dates)), dates_str)
end


# From PairPlots.jl
function margin_confidence_default_formatter(low,mid,high)
    largest_error = max(abs(high), abs(low))
    # Fallback for series with no variance
    if largest_error == 0
        if mid == 0
            digits_after_dot = 0
        else
            digits_after_dot = max(0, 1 - round(Int, log10(abs(mid))))
        end
        @static if VERSION >= v"1.10"
            title = @sprintf(
                "\$%.*f",
                digits_after_dot, mid,
            )
        else
            title = @eval @sprintf(
                $("\$%.$(digits_after_dot)f\$"),
                $mid,
            )
        end
        return title
    end

    digits_after_dot = max(0, 1 - round(Int, log10(largest_error)))
    use_scientific = digits_after_dot > 4

    if use_scientific
        if round(low, digits=digits_after_dot) == round(high, digits=digits_after_dot)
            title = @sprintf(
                "\$(%.1f \\pm %.1f)\\times 10^{-%d}\$",
                mid*10^(digits_after_dot-1),
                high*10^(digits_after_dot-1),
                (digits_after_dot-1)
            )
        else
            title = @sprintf(
                "\$(%.1f^{+%.1f}_{-%.1f})\\times 10^{-%d}\$",
                mid*10^(digits_after_dot-1),
                high*10^(digits_after_dot-1),
                low*10^(digits_after_dot-1),
                (digits_after_dot-1)
            )
        end
    else
        # '*' format specifier only supported in Julia 1.10+
        @static if VERSION >= v"1.10"
            if round(low, digits=digits_after_dot) == round(high, digits=digits_after_dot)
                title = @sprintf(
                    "\$%.*f \\pm %.*f\$",
                    digits_after_dot, mid,
                    digits_after_dot, high,
                )
            else
                title = @sprintf(
                    "\$%.*f^{+%.*f}_{-%.*f}\$",
                    digits_after_dot, mid,
                    digits_after_dot, high,
                    digits_after_dot, low
                )
            end
        else
            if round(low, digits=digits_after_dot) == round(high, digits=digits_after_dot)
                title = @eval @sprintf(
                    $("\$%.$(digits_after_dot)f \\pm %.$(digits_after_dot)f\$"),
                    $mid,
                    $high,
                )
            else
                title = @eval @sprintf(
                    $("\$%.$(digits_after_dot)f^{+%.$(digits_after_dot)f}_{-%.$(digits_after_dot)f}\$"),
                    $mid,
                    $high,
                    $low
                )
            end
        end
    end

    return title 
end

concat_with_nan(mat) =
    reduce((row_A, row_B) -> [row_A; NaN; row_B], eachrow(mat), init=Float64[])

# https://discourse.julialang.org/t/equivalent-of-matlabs-unwrap/44882/4?
function unwrap!(x, period = 2π)
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
	end
end

end
