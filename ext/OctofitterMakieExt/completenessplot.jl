##################################################
# Completeness map visualization
#
# Heatmap of detection completeness on a mass × separation grid,
# with marginal completeness curves on the top and right axes.
# Layout inspired by dotplot.jl.

function Octofitter.completenessplot(
    cmap::Octofitter.CompletenessMap,
    fname="completeness-map.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Figure(;
        size=(550, 450),
        figure...,
    )
    Octofitter.completenessplot!(fig.layout, cmap, args...; kwargs...)
    Makie.save(fname, fig, px_per_unit=3)
    return fig
end

function Octofitter.completenessplot!(
    gs,
    cmap::Octofitter.CompletenessMap;
    xlabel="separation [AU]",
    ylabel=Makie.rich("mass [M", Makie.subscript("jup"), "]"),
    colormap=Makie.cgrad(:viridis),
    colorrange=(0.0, 1.0),
    show_counts::Bool=false,
    title="",
    axis=(;),
)
    masses = cmap.masses
    seps = cmap.separations

    # ── Main heatmap axis ──
    ax = Axis(
        gs[2, 1];
        xscale=log10,
        yscale=log10,
        xlabel=xlabel,
        ylabel=ylabel,
        xgridvisible=false,
        ygridvisible=false,
        axis...
    )

    # Compute cell edges for the heatmap (log-spaced midpoints → edges)
    sep_edges = _log_edges(seps)
    mass_edges = _log_edges(masses)

    hm = heatmap!(ax, sep_edges, mass_edges, cmap.completeness';
        colormap=colormap,
        colorrange=colorrange,
    )

    # Optionally overlay detection counts as text
    if show_counts
        for im in eachindex(masses), is in eachindex(seps)
            n_det = cmap.n_detected[im, is]
            n_tot = cmap.n_total[im, is]
            if n_tot > 0
                text!(ax, seps[is], masses[im];
                    text="$(n_det)/$(n_tot)",
                    align=(:center, :center),
                    fontsize=8,
                    color=cmap.completeness[im, is] > 0.5 ? :black : :white,
                )
            end
        end
    end

    # ── Top marginal: completeness vs separation ──
    ta = Axis(
        gs[1, 1];
        xscale=log10,
        ylabel="compl.",
        xgridvisible=false,
        ygridvisible=false,
    )

    # Marginal completeness: average over mass axis
    marginal_sep = vec(sum(cmap.n_detected, dims=1)) ./ vec(max.(sum(cmap.n_total, dims=1), 1))
    lines!(ta, seps, marginal_sep; color=:black, linewidth=2)
    band!(ta, seps, fill(0.0, length(seps)), marginal_sep; color=(:dodgerblue, 0.3))

    # ── Right marginal: completeness vs mass ──
    ra = Axis(
        gs[2, 2];
        yscale=log10,
        xlabel="compl.",
        xgridvisible=false,
        ygridvisible=false,
    )

    # Marginal completeness: average over separation axis
    marginal_mass = vec(sum(cmap.n_detected, dims=2)) ./ vec(max.(sum(cmap.n_total, dims=2), 1))
    lines!(ra, marginal_mass, masses; color=:black, linewidth=2)
    band!(ra, fill(0.0, length(masses)), marginal_mass, masses; color=(:dodgerblue, 0.3))

    # ── Link axes and style ──
    linkxaxes!(ax, ta)
    linkyaxes!(ax, ra)
    ylims!(ta, 0, 1.05)
    xlims!(ra, 0, 1.05)
    hidedecorations!(ta, grid=false, label=false, ticklabels=false, ticks=false)
    hidexdecorations!(ta)
    hidedecorations!(ra, grid=false, label=false, ticklabels=false, ticks=false)
    hideydecorations!(ra)

    colsize!(gs, 1, Auto(8))
    rowsize!(gs, 2, Auto(8))
    colgap!(gs, 1, 4)
    rowgap!(gs, 1, 4)

    # Colorbar
    Colorbar(gs[2, 3], hm; label="detection fraction")

    if !isempty(title)
        Label(gs[1, 1, Top()], title, tellwidth=false, font=:bold, halign=:left)
    end

    return ax
end


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────

"""
Compute bin edges from log-spaced center values.
Each edge is the geometric mean of adjacent centers; outer edges are extrapolated.
"""
function _log_edges(centers::AbstractVector)
    n = length(centers)
    edges = Vector{Float64}(undef, n + 1)
    for i in 1:n-1
        edges[i+1] = sqrt(centers[i] * centers[i+1])
    end
    # Extrapolate outer edges
    if n >= 2
        edges[1] = centers[1]^2 / centers[2]
        edges[end] = centers[end]^2 / centers[end-1]
    else
        edges[1] = centers[1] / 2
        edges[end] = centers[1] * 2
    end
    return edges
end
