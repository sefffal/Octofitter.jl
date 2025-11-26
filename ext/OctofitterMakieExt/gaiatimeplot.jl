
##################################################
# Gaia DR4 IAD time plot
# Shows along-scan position vs time with residuals

function Octofitter.gaiatimeplot(
    model,
    results,
    fname="$(model.system.name)-gaiatime.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(700, 500),
        figure...
    )
    Octofitter.gaiatimeplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end

function Octofitter.gaiatimeplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains;
    # If less than 500 samples, just show all of them
    N = min(size(results, 1) * size(results, 3), 500),
    # If showing all samples, include each sample once.
    # Otherwise, sample randomly with replacement
    ii = (
        N == size(results, 1) * size(results, 3) ?
        (1:size(results, 1)*size(results, 3)) :
        rand(1:size(results, 1)*size(results, 3), N)
    ),
    axis=(;),
    colormap=Makie.cgrad([Makie.wong_colors()[1], "#DDDDDD"]),
    colorbar=true,
    top_time_axis=true,
    bottom_time_axis=true,
    alpha=min.(1, 100 / length(ii)),
    kwargs...
)
    gs = gridspec_or_fig

    # Find Gaia DR4 likelihood objects
    gaia_likes = filter(model.system.observations) do like_obj
        like_obj isa Octofitter.GaiaDR4AstromObs
    end
    if isempty(gaia_likes)
        @warn "No GaiaDR4AstromObs found in model"
        return Axis[]
    end
    likeobj = first(gaia_likes)

    # Compute simulations for all samples
    θ_systems_from_chain = Octofitter.mcmcchain2result(model, results)

    # Display all points, unless there are more than 5k
    jj = 1:size(results, 1) * size(results, 3)
    if size(results, 1) * size(results, 3) > 5_000
        jj = ii
    end

    sims = []
    for (θ_system, i) in zip(θ_systems_from_chain, jj)
        orbits = map(keys(model.system.planets)) do planet_key
            Octofitter.construct_elements(model, results, planet_key, i)
        end
        solutions = map(orbits) do orbit
            return orbitsolve.(orbit, likeobj.table.epoch)
        end
        θ_obs = θ_system.observations[Octofitter.normalizename(Octofitter.likelihoodname(likeobj))]
        centroid_pos_al_model_buffer = zeros(size(likeobj.table, 1))
        sim = Octofitter.simulate(
            likeobj,
            θ_system,
            θ_obs,
            orbits,
            solutions,
            0,
            centroid_pos_al_model_buffer
        )
        push!(sims, sim)
    end

    # Set up date ticks
    ts = likeobj.table.epoch
    date_pos, date_strs, xminorticks = _date_ticks(ts)

    # Create axes
    ax = Axis(gs[1, 1];
        ylabel="along scan [mas]",
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
    ax_resid = Axis(gs[2, 1];
        xlabel="time [MJD]",
        ylabel="residuals [mas]",
        xgridvisible=false,
        ygridvisible=false,
        xticksvisible=bottom_time_axis,
        xticklabelsvisible=bottom_time_axis,
        xlabelvisible=bottom_time_axis,
        axis...
    )
    linkxaxes!(ax, ax_resid)
    hlines!(ax_resid, 0, color=:black, linewidth=3)

    # Plot model simulations
    scatter!(ax,
        concat_with_nan(stack([likeobj.table.epoch for _ in sims])),
        concat_with_nan(stack(sims)),
        color=(Makie.wong_colors()[1], alpha),
        markersize=3,
        rasterize=4,
    )

    # Plot data points
    scatter!(ax,
        likeobj.table.epoch,
        likeobj.table.centroid_pos_al,
        color=:black,
        markersize=5
    )
    errorbars!(ax,
        likeobj.table.epoch,
        likeobj.table.centroid_pos_al,
        likeobj.table.centroid_pos_error_al,
        color=:black,
    )

    # Plot residuals as boxplots
    boxplot!(
        ax_resid,
        vec(stack(likeobj.table.epoch for sim in sims)),
        vec((stack(sims) .- likeobj.table.centroid_pos_al)),
        markersize=2,
        color=Makie.wong_colors()[1],
        strokecolor=:grey,
        whiskercolor=:grey,
        mediancolor=:grey,
        width=40
    )
    errorbars!(ax_resid,
        likeobj.table.epoch,
        likeobj.table.centroid_pos_al .* 0,
        likeobj.table.centroid_pos_error_al,
        color=:black,
        linewidth=1,
        whiskerwidth=5,
    )

    Makie.rowgap!(gs, 1, 0)

    return [ax, ax_resid]
end
