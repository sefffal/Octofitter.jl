using StatsBase

##################################################
# Mass vs. separation/period dot plot with marginal histograms

function Octofitter.dotplot(
    model,
    results,
    fname="$(model.system.name)-dotplot.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(500, 400),
        figure...,
    )
    Octofitter.dotplot!(fig.layout, model, results, args...; kwargs...)

    Makie.save(fname, fig, px_per_unit=3)

    return fig
end

function Octofitter.dotplot!(
    gridspec_or_fig,
    model::Octofitter.LogDensityModel,
    results::Chains;
    mode=:separation,
    epoch=nothing,
    colormap=Makie.cgrad(:deep, rev=false),
    colorbar=true,
    axis=(;),
    kwargs...
)
    gs = gridspec_or_fig

    planet_keys = collect(keys(model.system.planets))
    num_planets = length(planet_keys)

    if num_planets == 0
        error("No planets found in the model")
    end

    # Multi-planet colormaps
    colormaps = if num_planets > 1
        [:blues, :reds, :greens, :oranges, :purples]
    else
        [colormap]
    end
    colors = [:blue, :red, :green, :orange, :purple]

    # Set up tick values for log scales
    yticks = (
        2.0 .^ (-10:11),
        [".001"; ".002"; ".004"; ".008"; ".016"; ".031"; ".063"; ".125"; ".25"; ".5"; string.(2 .^ (0:11))]
    )

    if mode == :separation
        xticks = (
            2.0 .^ (-10:10),
            [".001"; ".002"; ".004"; ".008"; ".016"; ".031"; ".063"; ".125"; ".25"; ".5"; string.(2 .^ (0:10))]
        )
        xlabel = isnothing(epoch) ? "separation [AU]" : "separation [AU]\nat $(epoch)"
    elseif mode == :period
        xticks = (
            10.0 .^ (-2:7),
            string.(10.0 .^ (-2:7))
        )
        xlabel = "period [days]"
    else
        error("mode must be :separation or :period")
    end

    # Main scatter axis
    ax = Axis(
        gs[2, 1];
        yscale=log10,
        xscale=log10,
        xlabel=xlabel,
        ylabel=Makie.rich("mass [M", Makie.subscript("jup"), "]"),
        xticklabelrotation=deg2rad(45),
        yticklabelsize=8,
        xticklabelsize=9,
        xticks=xticks,
        yticks=yticks,
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )

    # Top histogram axis
    ta = Axis(
        gs[1, 1];
        xscale=log10,
        xticks=xticks,
        xgridvisible=false,
        ygridvisible=false,
    )

    # Right histogram axis
    ra = Axis(
        gs[2, 2];
        yscale=log10,
        yticks=yticks,
        yticklabelsize=0.1,
        xticksvisible=false,
        xticklabelsvisible=false,
        xgridvisible=false,
        ygridvisible=false,
    )

    local s = nothing
    for (i, planet_key) in enumerate(planet_keys)
        mass_key = Symbol("$(planet_key)_mass")
        if !haskey(results, mass_key)
            continue
        end

        els = Octofitter.construct_elements(model, results, planet_key, :)
        ecc_key = "$(planet_key)_e"
        ecc = vec(results[ecc_key])
        ii = sortperm(ecc, rev=true)

        if isempty(ii)
            continue
        end

        # Compute x-values based on mode
        if mode == :separation
            if isnothing(epoch)
                # Use semi-major axis as default
                xs = semimajoraxis.(els)
            else
                sols = orbitsolve.(els, mjd(epoch))
                xs = hypot.(posx.(sols), posy.(sols), posz.(sols))
            end
        elseif mode == :period
            xs = period.(els)
        end

        mass = vec(results[mass_key])

        cmap = num_planets > 1 ? colormaps[mod1(i, length(colormaps))] : colormap
        color_vals = num_planets > 1 ? colors[mod1(i, length(colors))] : nothing

        s = scatter!(
            ax,
            xs[ii],
            mass[ii];
            color=ecc[ii],
            colorrange=(0, 1),
            colormap=cmap,
            markersize=clamp(2000 / length(ii), 0.9, 5),
        )

        # Marginal histograms
        hist_color = num_planets == 1 ? :grey : (colors[mod1(i, length(colors))], 0.65)
        # hist!(ta, xs[ii]; gap=0, color=hist_color, bins=xticks[1])
        # hist!(ra, mass[ii]; direction=:x, color=hist_color, bins=yticks[1])



        # Marginal histograms - separation (top) as step histograms
        sep_bins = 2.0 .^ (-10:10)
        h1_sep = fit(Histogram, xs[ii], sep_bins)
        stairs!(ta, h1_sep.edges[1], [h1_sep.weights; 0]; color=hist_color, linewidth=2.5, step=:pre)

        # Marginal histograms - mass (right) as step histograms
        mass_bins = 2.0 .^ (-10:11)
        h1_mass = fit(Histogram, mass[ii], mass_bins)
        stairs!(ra, [h1_mass.weights; 0], h1_mass.edges[1]; color=hist_color, linewidth=2.5, step=:pre)

    end

    linkxaxes!(ax, ta)
    linkyaxes!(ax, ra)
    ylims!(ta, low=0)
    xlims!(ra, low=0)
    hidedecorations!(ra, grid=false, ticklabels=false)
    hidedecorations!(ta, grid=false)

    colsize!(gs, 1, Auto(8))
    rowsize!(gs, 2, Auto(8))
    colgap!(gs, 1, 0)
    rowgap!(gs, 1, 4)

    # Set axis limits
    ylims!(ra, 1e-2, 2^11)
    if mode == :separation
        xlims!(ta, 1e-3, 1024)
    elseif mode == :period
        xlims!(ta, 1e-1, 1e7)
    end

    # Add colorbar
    if colorbar && !isnothing(s)
        Colorbar(gs[2, 3], s; label="eccentricity")
    end

    Label(fig[0,1], "$sysname",tellwidth=false, font=:bold)

    return ax
end
