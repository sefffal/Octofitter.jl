"""
    octoplot(model, chain; kwargs...)

Generate publication-quality figures showing orbit fits and data. Creates a multi-panel
visualization that automatically adapts to show the types of data present in your model.

The file is automatically saved to "\$(model.system.name)-plot-grid.png". You can override 
this by passing `fname="fname.png"`

# Panel Types (enabled with boolean flags)
- `show_astrom=true`: Sky-projected orbits (mas)
- `show_physical_orbit=true`: Physical orbits (AU)
- `show_astrom_time=true`: Separation/PA vs time
- `show_rv=true`: Stellar RV curve
- `show_relative_rv=true`: Planet-star relative RV
- `show_pma=true`: Proper motion anomaly
- `show_mass=true`: Mass posterior

# Optional Arguments
- `N=250`: Number of posterior samples to plot
- `ii=nothing`: Specific sample indices to plot (overrides N)
- `ts=nothing`: Time range for plots (MJD values)
- `colormap=:plasma`: Colormap for orbital phase
- `alpha=auto`: Transparency of orbit lines
- `figscale=1.0`: Overall figure size multiplier
- `mark_epochs_mjd=Float64[]`: Epochs to highlight across panels

# Returns
A Makie figure object that can be further customized and saved in various formats:
```julia
save("orbit_plot.pdf", fig)  # PDF (CairoMakie)
save("orbit_plot.png", fig, px_per_unit=3)  # Higher resolution PNG
"""
function Octofitter.octoplot(
    model::Octofitter.LogDensityModel,
    results::Chains;
    figure=(;),
    fname="$(model.system.name)-plot-grid.png",
    kwargs...
)
    fig = Figure(;
        figure...
    )
    Octofitter.octoplot!(fig.layout,model,results;kwargs...)

    try
        Makie.resize_to_layout!(fig)
    catch
        @warn "Erorr occurred with resize_to_layout!"
    end

    save(fname, fig)
    return fig
end
function Octofitter.octoplot!(
    fig::Makie.GridLayout,
    model::Octofitter.LogDensityModel,
    results::Chains;
    show_astrom=nothing,
    show_absastrom=nothing,
    show_physical_orbit=nothing,
    show_astrom_time=nothing,
    show_pma=nothing,
    show_mass=false,
    show_rv=nothing,
    show_relative_rv=nothing,
    show_hipparcos=nothing,
    show_gaia=nothing,
    residuals=false,
    mark_epochs_mjd=Float64[],
    # If less than 500 samples, just show all of them
    N  = min(size(results, 1)*size(results, 3), 250),
    ts = nothing,
    # If showing all samples, include each sample once.
    # Otherwise, sample randomly with replacement
    ii = (
        N == size(results, 1)*size(results, 3) ? 
        (1:size(results, 1)*size(results, 3)) :
        rand(1:size(results, 1)*size(results, 3),N)
    ),
    # The user can of course just override the above directly.
    colormap=Makie.cgrad([Makie.wong_colors()[1], "#FAFAFA"]),
    alpha=min.(1, 100 / length(ii)),
    figscale=0.6,
    colorbar = true,
    show_post_pred_legend=true,
    show_instrument_names=true
)

    defaults_used = all(isnothing, (
        show_astrom,
        show_physical_orbit,
        show_astrom_time,
        show_pma,
        show_absastrom,
        show_mass,
        show_rv,
        show_relative_rv,
        show_hipparcos,
        show_gaia,
    ))

    if colormap isa Symbol || colormap isa String
        colormap = Makie.cgrad([colormap,"#FAFAFA"])
    end

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

    if isnothing(show_physical_orbit)
        show_physical_orbit = false
        for planet in model.system.planets
            show_physical_orbit |= 
                Octofitter.orbittype(planet) <: KepOrbit || 
                Octofitter.orbittype(planet) <: CartesianOrbit
        end
    end

    if isnothing(show_astrom_time)
        show_astrom_time = false
        for planet in model.system.planets
            for like_obj in planet.observations
                show_astrom_time |=  nameof(typeof(like_obj)) == :PlanetRelAstromObs
            end
        end
    end

    if isnothing(show_relative_rv)
        show_relative_rv = false
        for planet in model.system.planets
            for like_obj in planet.observations
                show_relative_rv |=  nameof(typeof(like_obj)) == :PlanetRelativeRVObs
            end
        end
    end

    if isnothing(show_rv)
        show_rv = false
        for like_obj in model.system.observations
            show_rv |=  nameof(typeof(like_obj)) == :StarAbsoluteRVObs
            show_rv |=  nameof(typeof(like_obj)) == :MarginalizedStarAbsoluteRVObs
        end
    end

    if isnothing(show_pma)
        show_pma = false
        for like_obj in model.system.observations
            show_pma |= like_obj isa Octofitter.HGCAObs
            show_pma |= like_obj isa Octofitter.HGCAInstantaneousObs
            show_pma |= like_obj isa Octofitter.GaiaDifferenceObs
        end
    end


    hgca_detected = false
    for like_obj in model.system.observations
        hgca_detected |= like_obj isa HGCAObs
        hgca_detected |= like_obj isa HGCAInstantaneousObs
    end

    if isnothing(show_absastrom)
        show_absastrom = false
        for like_obj in model.system.observations
            show_absastrom |= like_obj isa Octofitter.GaiaHipparcosUEVAJointObs
        end
    end


    if isnothing(show_hipparcos)
        show_hipparcos = false
        for like_obj in model.system.observations
            if like_obj isa HipparcosIADObs
                show_hipparcos = true
            end
        end
    end

    if isnothing(show_gaia)
        show_gaia = false
        for like_obj in model.system.observations
            if like_obj isa Octofitter.GaiaDR4AstromObs
                show_gaia = true
            end
        end
    end

    if isnothing(show_mass)
        show_mass = false
        for planet_key in keys(model.system.planets)
            show_mass |= haskey(results, Symbol("$(planet_key)_mass"))
        end
    end

    if defaults_used
        @info(
            "You can control the panels included by `octoplot` by passing keyword arguments. Currently selected:",
            show_astrom,
            show_physical_orbit,
            show_astrom_time,
            show_pma,
            show_mass,
            show_rv,
            show_relative_rv,
            show_hipparcos,
            show_gaia,
        )
        @info "pass true or false for one of these arguments to suppress this message."
    end

    if isempty(model.system.planets)
        @warn "No planets in this model"
        return Figure()
    end


    # determine a consistent time axis/scale for all kinds of data
    # Find the minimum time step size
    periods = Float64[]
    for planet_key in keys(model.system.planets)
        orbs = Octofitter.construct_elements(model, results, planet_key, ii)
        append!(periods, period.(orbs))
    end

    if isempty(periods)
        @warn "No samples in this chain"
        return Figure()
    end
    min_period = minimum(periods)
    med_period = quantile(periods, 0.35)
    if !isfinite(med_period)
        med_period = 0.0
    end


    if isnothing(ts)
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
        if hgca_detected
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

        t_start = max(t_start,  mjd("1900"))
        t_stop = min(t_stop,  mjd("2100"))

        if isfinite(min_period)
            ts = range(t_start, t_stop, step=min_period / 30)
        else
            ts = [t_start-1, t_start, t_start+1]
        end

        if hgca_detected
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
        if length(ts) < 200
            # Things can get very slow if we have some very short period 
            # solutions. Limit the maximum number of points per orbit.
            ts = range(t_start, t_stop, length=200)
        end
    end
    
    # Show a colorbar for only the first sub-plot, and don't repeat it.
    top_time_axis = true
    item = 0
    cols = 1


    # Prevent figure updates from computing while we are drawing
    # Speeds up figure creation
    fig = update(fig) do fig
    
        if show_astrom
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            gl = GridLayout(
                fig[row,col],
                width=500figscale,
                height=(400+40*length(mark_epochs_mjd)+ 155*(length(mark_epochs_mjd)>0))*figscale,
            )
            Octofitter.astromplot!(gl, model, results; ii, ts, colorbar, colormap, mark_epochs_mjd, alpha, show_post_pred_legend, show_instrument_names)
            colorbar = false
        end
        if show_physical_orbit
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            gl = GridLayout(
                fig[row,col],
                width=500figscale,
                height=(400+40*length(mark_epochs_mjd)+ 155*(length(mark_epochs_mjd)>0))*figscale,
            )
            physorbplot!(gl, model, results; ii, ts, colorbar, colormap, mark_epochs_mjd, alpha)
            colorbar = false
        end

        axes_to_link = []

        if show_astrom_time
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            gl = GridLayout(
                fig[row,col],
                width=500figscale,
                height=300figscale,
            )
            bottom_time_axis = !(show_pma || show_rv || show_relative_rv || show_absastrom || show_hipparcos)
            ax = astromtimeplot!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis, colormap, mark_epochs_mjd, alpha, residuals)
            top_time_axis = false
            append!(axes_to_link,ax)
        end


        if show_rv
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            matching_like_objs = filter(model.system.observations) do like_obj
                nameof(typeof(like_obj)) ∈ (
                    :MarginalizedStarAbsoluteRVObs,
                    :StarAbsoluteRVObs
                )
            end
            gl = GridLayout(
                fig[row,col],
                width=500figscale,
                height=(135 * max(1, length(matching_like_objs)) *figscale),
            )
            bottom_time_axis = !(show_pma || show_relative_rv || show_absastrom || show_hipparcos)
            ax = rvtimeplot!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis, colormap, alpha)
            top_time_axis = false
            Makie.rowgap!(gl, 10.)
            append!(axes_to_link,ax)
        end


        if show_relative_rv
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            gl = GridLayout(
                fig[row,col],
                width=500figscale,
                height=135figscale,
            )
            bottom_time_axis = !(show_pma || show_absastrom || show_hipparcos)
            ax = rvtimeplot_relative!(gl, model, results; ii, ts, colorbar, top_time_axis, bottom_time_axis, colormap, alpha)
            top_time_axis = false
            Makie.rowgap!(gl, 10.)
            append!(axes_to_link,ax)
        end

        # if show_hgca
        #     item += 1
        #     col = mod1(item, cols)
        #     row = cld(item, cols)
        #     gl = GridLayout(
        #         fig[row,col],
        #         width=500figscale,
        #         height=480figscale,
        #     )
        #     ax = Octofitter.hgcaplot!(gl, model, results; ii, ts, colorbar, top_time_axis, colormap, alpha)
        #     top_time_axis = false
        #     Makie.rowgap!(gl, 10.)
        #     append!(axes_to_link,ax)
        # end

        if show_pma
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            height=300figscale
            for like_obj in model.system.observations
                if like_obj isa HGCAObs || like_obj isa HGCAInstantaneousObs
                    height+=180figscale
                end
            end
            for like_obj in model.system.observations
                if like_obj isa Octofitter.GaiaDifferenceObs
                    height+=180figscale
                end
            end
            gl = GridLayout(
                fig[row,col];
                width=500figscale,
                height,
            )
            ax = pmaplot!(gl, model, results; ii, ts, colorbar, top_time_axis, colormap, alpha)
            top_time_axis = false
            Makie.rowgap!(gl, 10.)
            append!(axes_to_link,ax)
        end

         if show_absastrom
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            height=460figscale
            gl = GridLayout(
                fig[row,col];
                width=500figscale,
                height,
            )
            ax = absastromplot!(gl, model, results; ii, ts, colorbar, top_time_axis, colormap, alpha)
            top_time_axis = false
            Makie.rowgap!(gl, 10.)
            append!(axes_to_link,ax)
        end

        if show_hipparcos
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            gl = GridLayout(
                fig[row,col],
                width=500figscale,
                height=480figscale,
            )
            ax = hipparcosplot!(gl, model, results; ii, ts, colorbar, top_time_axis, colormap, alpha)
            top_time_axis = false
            Makie.rowgap!(gl, 10.)
            append!(axes_to_link,ax)
        end

        if show_gaia
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            gl = GridLayout(
                fig[row,col],
                width=500figscale,
                height=350figscale,
            )
            bottom_time_axis = !show_mass
            ax = Octofitter.gaiatimeplot!(gl, model, results; ii, colorbar, top_time_axis, bottom_time_axis, alpha)
            top_time_axis = false
            Makie.rowgap!(gl, 10.)
            append!(axes_to_link, ax)
        end

        if show_mass
            item += 1
            col = mod1(item, cols)
            row = cld(item, cols)
            gl = GridLayout(
                fig[row,col],
                width=500figscale,
                height=400figscale,
            )
            Octofitter.masspostplot!(gl, model, results;)
            Makie.rowgap!(gl, 10.)
        end
        Makie.rowgap!(fig, 10.)


        # Show warning about pathfinder approximation
        if hasproperty(results.info, :sampler) && results.info.sampler == "pathfinder"
            Label(fig[0,:], Makie.rich("pathfinder approximation",color=:darkred, font=:bold), tellwidth=false)
        end


        if !isempty(axes_to_link)
            Makie.linkxaxes!(axes_to_link...)
            yspace = maximum(Makie.tight_yticklabel_spacing!, axes_to_link)
            for ax in axes_to_link
                ax.yticklabelspace = yspace + 30
            end
        end

        return fig
    end
    return fig
end
 

# Helper from AlgebraOfGraphics -- prevent computing layout updates after adding
# each series
get_layout(gl::Makie.GridLayout) = gl
get_layout(f::Union{Makie.Figure, Makie.GridPosition}) = f.layout
get_layout(l::Union{Makie.Block, Makie.GridSubposition}) = get_layout(l.parent)

# Wrap layout updates in an update block to avoid triggering multiple updates
function update(f, fig)
    layout = get_layout(fig)
    block_updates = layout.block_updates
    layout.block_updates = true
    output = f(fig)
    layout.block_updates = block_updates
    block_updates || Makie.GridLayoutBase.update!(layout)
    return output
end