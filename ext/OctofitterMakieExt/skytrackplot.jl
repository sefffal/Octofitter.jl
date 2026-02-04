
##################################################
# Sky Track Plot
# Shows the full sky track including parallax and proper motion
# for a single posterior sample. Not part of octoplot.

function Octofitter.skytrackplot(
    model,
    results,
    args...;
    fname="$(model.system.name)-skytrack.png",
    kwargs...,
)
    fig = Figure(
        size=(700, 500)
    )
    # Wrap in update() to prevent StackOverflow from circular Observable updates
    update(fig) do fig
        Octofitter.skytrackplot!(fig.layout, model, results, args...; kwargs...)
    end

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end

function Octofitter.skytrackplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains,
    sample_idx=argmax(results["logpost"][:]);
    ts=nothing,
    axis=(;),
    keplerian_mult=1.0,  # Exaggerate the Keplerian signal by this factor
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

    # Default time range if not provided
    if isnothing(ts)
        t_start, t_stop = extrema(likeobj.table.epoch)
        delta = t_stop - t_start
        ts = range(start=t_start - 0.1 * delta, stop=t_stop + 0.1 * delta, length=200)
    end

    # Create axis
    ax = Axis(gs[1, 1];
        xgridvisible=false,
        ygridvisible=false,
        xlabel="Δα* [mas]",
        ylabel="Δδ [mas]",
        autolimitaspect=1,
        xreversed=true,
        axis...
    )

    # Calculate the parallax factors for each epoch in the grid
    earth_pos_vel = Table(Octofitter.geocentre_position_query.(ts))

    # We use a reference position (could be from Gaia catalog)
    α = likeobj.gaia_sol.ra # Reference RA (arbitrary for relative plot)
    δ = likeobj.gaia_sol.dec # Reference Dec (arbitrary for relative plot)

    # Proper motion contribution
    pmα = θ_obs.pmra * ts / 365.25
    pmδ = θ_obs.pmdec * ts / 365.25
    δ_track = @. δ + pmδ / 60 / 60 / 1000
    α_track = @. α + pmα / 60 / 60 / 1000 * cosd(δ_track)

    # Parallax displacement
    Δα = @. θ_obs.plx * (earth_pos_vel.x * sin(α_track) - earth_pos_vel.y * cos(α_track))
    Δδ = @. θ_obs.plx * (earth_pos_vel.x * cos(α_track) * sin(δ_track) + earth_pos_vel.y * sin(α_track) * sin(δ_track) - earth_pos_vel.z * cos(δ_track))

    # Keplerian perturbation from planets
    Δα_kep = zeros(length(ts))
    Δδ_kep = zeros(length(ts))

    for planet_i in eachindex(orbits)
        sol = orbitsolve.(orbits[planet_i], ts)
        # Add perturbation from planet
        Δα_kep .+= keplerian_mult * raoff.(sol, θ_system.planets[planet_i].mass * mjup2msol)
        Δδ_kep .+= keplerian_mult * decoff.(sol, θ_system.planets[planet_i].mass * mjup2msol)
    end

    # Plot the full sky track
    lines!(ax,
        Δα .+ pmα .+ Δα_kep,
        Δδ .+ pmδ .+ Δδ_kep,
        color=Makie.wong_colors()[1],
        linewidth=2
    )

    # Calculate residuals and project them back on the track
    # sim is a NamedTuple with along_scan_residuals_buffer, ra_offset_buffer, dec_offset_buffer
    resids = sim.along_scan_residuals_buffer .- likeobj.table.centroid_pos_al
    s = sin.(likeobj.table.scan_pos_angle)
    c = cos.(likeobj.table.scan_pos_angle)
    alpha_res = @. resids * s
    delta_res = @. resids * c

    # Calculate model position at data epochs
    pmα_dat = θ_obs.pmra * likeobj.table.epoch / 365.25
    pmδ_dat = θ_obs.pmdec * likeobj.table.epoch / 365.25
    δ_track_dat = @. δ + pmδ_dat / 60 / 60 / 1000
    α_track_dat = @. α + pmα_dat / 60 / 60 / 1000 * cosd(δ_track_dat)

    Δα_kep_track = zeros(length(likeobj.table.epoch))
    Δδ_kep_track = zeros(length(likeobj.table.epoch))

    for planet_i in eachindex(orbits)
        sol = orbitsolve.(orbits[planet_i], likeobj.table.epoch)
        # Add perturbation from planet
        Δα_kep_track .+= keplerian_mult * raoff.(sol, θ_system.planets[planet_i].mass * mjup2msol)
        Δδ_kep_track .+= keplerian_mult * decoff.(sol, θ_system.planets[planet_i].mass * mjup2msol)
    end

    # Parallax displacement at data epochs
    # Note: xyz coordinates are nested in likeobj.table.xyz for GaiaDR4AstromObs
    Δα_dat = @. θ_obs.plx * (likeobj.table.xyz.x * sin(α_track_dat) - likeobj.table.xyz.y * cos(α_track_dat))
    Δδ_dat = @. θ_obs.plx * (likeobj.table.xyz.x * cos(α_track_dat) * sin(δ_track_dat) + likeobj.table.xyz.y * sin(α_track_dat) * sin(δ_track_dat) - likeobj.table.xyz.z * cos(δ_track_dat))

    # Plot data points
    scatter!(ax,
        Δα_dat .+ pmα_dat .+ Δα_kep_track .+ alpha_res,
        Δδ_dat .+ pmδ_dat .+ Δδ_kep_track .+ delta_res,
        marker=:rect,
        markersize=5,
        color=:black
    )

    return ax
end
