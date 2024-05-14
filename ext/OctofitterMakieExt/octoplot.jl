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
