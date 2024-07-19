"""

    octoplot(model, chain; fname, ...)

Generate a plot of orbits drawn from the posterior `chain`. There are many options described below.

**Output Panels**
`show_astrom`       : Plot relative RA and DEC of planets vs the primary. Only supported for `Visual` orbit types.
`show_astrom_time`  : Plot relative separation and position angle of planets vs the primary. Only supported for `Visual` orbit types.
`show_rv`           : Plot the absolute radial velocity of the primary. Lines are shifted vertically by RV zero-point to match data. 
`show_relative_rv`  : Plot the relative radial velocity between planets and the primary. 
`show_hgca`         : Plot instantaneous proper motion in R.A. and in Dec., and plot the time averaged HGCA quantities. Only supported for models including an `HGCALikelihood`.
`show_mass`         : Plot a mini corner plot of mass 

**Display Options**
`N`                 : Controls how many orbits are sampled with replacement from the posterior. Default is 250. If there are less than 250 samples, just plot each one once (not sampled).
`ii`                : A vector of posterior sample indices to plot (overrides N and the sampling described above)
`mark_epochs_mjd`   : A vector dates in MJD format to mark on the astrometry plots.
`colormap`          : A symbol from ColorSchemes.jl, or a ColorScheme object. Eg. try  `:viridis`.
`planet_rv`         : Applicable to `show_rv`. Show planets' absolute RV instead of star absolute RV.
`figure`            : Plot into an existing Makie Figure object. Useful for customizing overall figure appearance.

Note that the chain sample subset indices `ii` are not used for certain plots where it makes more sense to show all samples,
including the HGCA scatter plots (`show_hgca`) and mass plot histograms (`show_mass`). If you want to control which samples
are displayed in those plots, apply your filtering before passing the chain to `octoplot`:
```julia
chain_filtered = chain[vec(chain[:b_e]) .> 0.5,:,:]
octoplot(model, chain_filtered)
```
"""
function Octofitter.octoplot(
    model::Octofitter.LogDensityModel,
    results::Chains;
    fname="$(model.system.name)-plot-grid.png",
    show_astrom=nothing,
    show_astrom_time=nothing,
    show_hgca=nothing,
    show_mass=false,
    show_rv=nothing,
    planet_rv=nothing, # plot planet(s) absolute RV in addition to star
    show_relative_rv=nothing,
    mark_epochs_mjd=Float64[],
    figure=(;),
    # If less than 500 samples, just show all of them
    N  = min(size(results, 1)*size(results, 3), 250),
    # If showing all samples, include each sample once.
    # Otherwise, sample randomly with replacement
    ii = (
        N == size(results, 1)*size(results, 3) ? 
        (1:size(results, 1)*size(results, 3)) :
        rand(1:size(results, 1)*size(results, 3),N)
    ),
    # The user can of course just override the above directly.
    colormap=Makie.cgrad(:plasma,rev=true),
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

    if isnothing(show_rv)
        show_rv = false
        for like_obj in model.system.observations
            show_rv |=  nameof(typeof(like_obj)) == :StarAbsoluteRVLikelihood
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
            if isnothing(planet_rv) && nameof(typeof(like_obj)) == :StarAbsoluteRVLikelihood
                planet_rv = true
            end
        end
    end
    if isnothing(planet_rv)
        planet_rv = false
    end
    
    min_period = minimum(periods)
    med_period = quantile(periods, 0.35)
    if !isfinite(med_period)
        med_period = 0.0
    end

    # Always try and show at least one full period
    # decide t_start and t_end based on epochs of data if astrometry 
    # available, otherwise other data, otherwise arbitrary
    t_start = Inf
    t_stop = -Inf
    if !isempty(mark_epochs_mjd)
        t_start, t_stop = extrema(mark_epochs_mjd)
    end
    for like_obj in model.system.observations
        if hasproperty(like_obj, :table) && hasproperty(like_obj.table, :epoch)
            t_start′, t_stop′ = extrema(like_obj.table.epoch,)
            t_start = min(t_start′, t_start)
            t_stop = max(t_stop′, t_stop)
        end
    end
    for planet in model.system.planets
        for like_obj in planet.observations
            if hasproperty(like_obj, :table) && hasproperty(like_obj.table, :epoch)
                t_start′, t_stop′ = extrema(like_obj.table.epoch)
                t_start = min(t_start′, t_start)
                t_stop = max(t_stop′, t_stop)
            end
        end
    end
    # Show the full range of the HGCA if plotted
    if show_hgca
        t_start′ = years2mjd(1990)
        t_stop′ = years2mjd(2020)
        t_start = min(t_start′, t_start)
        t_stop = max(t_stop′, t_stop)
    end
    # If no data use current date as midpoint
    if !isfinite(t_start) || !isfinite(t_stop)
        t_start = t_stop = mjd()
    end

    # Put a little padding on either end
    d = t_stop - t_start
    t_start -= 0.015d
    t_stop += 0.015d
    # Now make sure we are showing enough time for a typical period
    t_cur = t_stop-t_start
    if t_cur < med_period
        t_extend = med_period - t_cur
        t_start -= t_extend/2
        t_stop += t_extend/2
    end

    if isfinite(min_period)
        ts = range(t_start, t_stop, step=min_period / 30)
    else
        ts = [t_start-1, t_start, t_start+1]
    end

    if show_hgca
        # Make sure we include the approx. data epochs of the HGCA
        ts = sort([ts; years2mjd(1990.25); years2mjd((1990.25 + 2016) / 2); years2mjd(2016)])
    end

    # TODO: could have a dense grid of Time to allow different resolutions per orbit.
    # then could step in eccanom.
    if length(ts) > 1000
        # Things can get very slow if we have some very short period 
        # solutions. Limit the maximum number of points per orbit.
        ts = range(t_start, t_stop, length=1000)
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
        Octofitter.astromplot!(gl, model, results; ii, ts, colorbar, colormap, mark_epochs_mjd)
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
        ax = astromtimeplot!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis, colormap, mark_epochs_mjd)
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
        ax = rvtimeplot!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis, planet_rv, colormap)
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
        ax = rvtimeplot_relative!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis, colormap)
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
        ax = Octofitter.hgcaplot!(gl, model, results; ii, ts, colorbar, top_time_axis, colormap)
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
