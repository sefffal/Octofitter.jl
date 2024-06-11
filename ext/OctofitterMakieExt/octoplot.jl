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
    planet_rv=nothing, # plot planet(s) absolute RV in addition to star
    show_relative_rv=nothing,
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

    if isnothing(show_relative_rv)
        show_relative_rv = false
        for planet in model.system.planets
            for like_obj in planet.observations
                show_relative_rv |=  nameof(typeof(like_obj)) == :PlanetRelativeRVLikelihood
            end
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

    # determine a consistent time axis/scale for all kinds of data
    # Find the minimum time step size
    periods = Float64[]
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(results, planet_key, :)
        append!(periods, period.(orbs))

        for like_obj in model.system.planets[planet_key].observations
            if isnothing(planet_rv) && nameof(typeof(like_obj)) == :PlanetAbsoluteRVLikelihood
                planet_rv = true
            end
        end
    end
    if isnothing(planet_rv)
        planet_rv = false
    end
    
    min_period = minimum(periods)
    med_period = quantile(periods, 0.35)

    # Always try and show at least one full period
    # decide t_start and t_end based on epochs of data if astrometry 
    # available, otherwise other data, otherwise arbitrary
    t_start = Inf
    t_stop = -Inf
    for like_obj in model.system.observations
        if hasproperty(like_obj, :table) && hasproperty(like_obj.table, :epoch)
            t_start′, t_stop′ = extrema(like_obj.table.epoch,)
            d = t_stop′ - t_start′
            t_start′ -= 0.015d
            t_stop′ += 0.015d
            if t_start′ < t_start
                t_start = t_start′
            end
            if t_stop′ > t_stop
                t_stop = t_stop′
            end
        end
    end
    for planet in model.system.planets
        for like_obj in planet.observations
            if hasproperty(like_obj, :table) && hasproperty(like_obj.table, :epoch)
                t_start′, t_stop′ = extrema(like_obj.table.epoch)
                d = t_stop′ - t_start′
                t_start′ -= 0.015d
                t_stop′ += 0.015d
                if t_start′ < t_start
                    t_start = t_start′
                end
                if t_stop′ > t_stop
                    t_stop = t_stop′
                end
            end
        end
    end
    # Show the full range of the HGCA if plotted
    if show_hgca
        t_start′ = years2mjd(1990)
        t_stop′ = years2mjd(2020)
        d = t_stop′ - t_start′
        t_start′ -= 0.015d
        t_stop′ += 0.015d
        if t_start′ < t_start
            t_start = t_start′
        end
        if t_stop′ > t_stop
            t_stop = t_stop′
        end
    end
    # If no data use current date as midpoint
    if !isfinite(t_start) || !isfinite(t_stop)
        t_start = t_stop = mjd()
    end
    # Now make sure we are showing enough time for a typical period
    t_cur = t_stop-t_start
    if t_cur < med_period
        t_extend = med_period - t_cur
        t_start -= t_extend/2
        t_stop += t_extend/2
    end
    ts = range(t_start, t_stop, step=min_period / 150)

    if show_hgca
        # Make sure we include the approx. data epochs
        ts = sort([ts; years2mjd(1990.25); years2mjd((1990.25 + 2016) / 2); years2mjd(2016)])
        # Find earliest epoch
        # epoch_0 = years2mjd(1991.25)
    end

    if length(ts) > 1500
        # Things can get very slow if we have some very short period 
        # solutions
        ts = range(t_start, t_stop, length=1500)
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
        Octofitter.astromplot!(gl, model, results; ii, ts, colorbar,)
        colorbar = false
    end

    axes_to_link = []

    if show_astrom_time
        item += 1
        col = mod1(item, cols)
        row = cld(item, cols)
        gl = GridLayout(
            fig[row,col],
            width=500,
            height=300,
        )
        bottom_time_axis = !(show_hgca || show_rv || show_relative_rv)
        ax = astromtimeplot!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis)
        top_time_axis = false
        append!(axes_to_link,ax)
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
        bottom_time_axis = !(show_hgca || show_relative_rv)
        ax = rvtimeplot!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis, planet_rv)
        top_time_axis = false
        append!(axes_to_link,ax)
    end


    if show_relative_rv
        item += 1
        col = mod1(item, cols)
        row = cld(item, cols)
        gl = GridLayout(
            fig[row,col],
            width=500,
            height=135,
        )
        bottom_time_axis = !show_hgca
        ax = rvtimeplot_relative!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis)
        top_time_axis = false
        append!(axes_to_link,ax)
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
        ax = Octofitter.hgcaplot!(gl, model, results; ii, ts, colorbar, top_time_axis)
        top_time_axis = false
        append!(axes_to_link,ax)
    end

    if !isempty(axes_to_link)
        Makie.linkxaxes!(axes_to_link...)
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
