
##################################################
# Gaia Star Plot
# Shows the star's orbit in RA/Dec space for a single posterior sample
# Similar to astromplot but shows the star's motion due to companions,
# and like rvpostplot, only plots a single draw at a time since
# the detrending is different per draw.

function Octofitter.gaiastarplot(
    model,
    results,
    args...;
    fname="$(model.system.name)-gaiastarplot.png",
    kwargs...,
)
    fig = Figure(
        size=(600, 600)
    )
    Octofitter.gaiastarplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end

function Octofitter.gaiastarplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains,
    sample_idx=argmax(results["logpost"][:]);
    axis=(;),
)
    if gridspec_or_fig isa Figure
        gs = GridLayout(gridspec_or_fig[1, 1])
    else
        gs = gridspec_or_fig
    end

    mjup2msol = Octofitter.mjup2msol

    # Find Gaia DR4 likelihood objects
    gaia_likes = filter(model.system.observations) do like_obj
        like_obj isa Octofitter.GaiaDR4AstromObs
    end
    if isempty(gaia_likes)
        error("No GaiaDR4AstromObs found in model")
    end
    likeobj = first(gaia_likes)

    # Get the sample
    θ_system = Octofitter.mcmcchain2result(model, results, sample_idx)
    θ_obs = θ_system.observations[Octofitter.normalizename(Octofitter.likelihoodname(likeobj))]

    # Construct orbits for this sample
    orbits = map(keys(model.system.planets)) do planet_key
        Octofitter.construct_elements(model, results, planet_key, sample_idx)
    end

    # Compute simulation at data epochs
    solutions = map(orbits) do orbit
        return orbitsolve.(orbit, likeobj.table.epoch)
    end
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

    # Create axis
    ax = Axis(gs[1, 1],
        xreversed=true,
        autolimitaspect=1,
        xgridvisible=false,
        ygridvisible=false,
        xlabel="Δα* [mas]",
        ylabel="Δδ [mas]",
        axis...
    )
    vlines!(ax, 0, color=:grey, linestyle=:dash)
    hlines!(ax, 0, color=:grey, linestyle=:dash)

    # Plot the Keplerian orbit of the star
    EAs = range(0, 2pi, length=150)
    Δα_kep = zeros(150)
    Δδ_kep = zeros(150)

    for planet_i in eachindex(orbits)
        orbit = orbits[planet_i]
        # Orbit is perfectly periodic, so take equal steps in eccentric anomaly
        sols = orbitsolve_eccanom.(orbit, EAs)
        # Add perturbation from planet
        Δα_kep .+= raoff.(sols, θ_system.planets[planet_i].mass * mjup2msol)
        Δδ_kep .+= decoff.(sols, θ_system.planets[planet_i].mass * mjup2msol)
    end

    lines!(ax, Δα_kep, Δδ_kep, color=Makie.wong_colors()[1], linewidth=2)

    # Calculate residuals and project them back on the orbit
    resids = sim .- likeobj.table.centroid_pos_al
    s = sin.(likeobj.table.scan_pos_angle)
    c = cos.(likeobj.table.scan_pos_angle)
    alpha_res = @. resids * s
    delta_res = @. resids * c

    # Calculate model position at data epochs
    Δα_kep_track = zeros(length(likeobj.table.epoch))
    Δδ_kep_track = zeros(length(likeobj.table.epoch))

    for planet_i in eachindex(orbits)
        sol = orbitsolve.(orbits[planet_i], likeobj.table.epoch)
        # Add perturbation from planet
        Δα_kep_track .+= raoff.(sol, θ_system.planets[planet_i].mass * mjup2msol)
        Δδ_kep_track .+= decoff.(sol, θ_system.planets[planet_i].mass * mjup2msol)
    end

    # Plot error bars along scan direction
    σ = likeobj.table.centroid_pos_error_al
    x_centers = Δα_kep_track .+ alpha_res
    y_centers = Δδ_kep_track .+ delta_res

    n = length(σ)
    x_bars = Vector{Float64}(undef, 3n)
    y_bars = Vector{Float64}(undef, 3n)

    for i in 1:n
        x_bars[3i-2] = x_centers[i] - σ[i] * s[i]
        y_bars[3i-2] = y_centers[i] - σ[i] * c[i]
        x_bars[3i-1] = x_centers[i] + σ[i] * s[i]
        y_bars[3i-1] = y_centers[i] + σ[i] * c[i]
        x_bars[3i] = NaN
        y_bars[3i] = NaN
    end

    lines!(ax, x_bars, y_bars, color=:black)

    # Plot data points
    scatter!(ax,
        Δα_kep_track .+ alpha_res,
        Δδ_kep_track .+ delta_res,
        marker=:rect,
        markersize=5,
        color=:black
    )

    # Plot star at origin
    scatter!(ax, [0], [0], marker='★', markersize=20, color=:white, strokecolor=:black, strokewidth=1.5)

    return ax
end
