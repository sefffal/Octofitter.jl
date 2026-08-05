# ---------------------------------------------------
# OctofitterMakieExt — the v2 plotting extension.
#
# Presentation only: every model/data/residual quantity comes from the core
# plotting API (`PosteriorSeries`, `plotchannels`, `Octofitter.residuals`,
# `modelcurves`, `rowsignal`), which itself defers to the likelihoods.
# Nothing in this file recomputes likelihood math.
#
# Four generic panels — sky, time-series, phase-folded (each with residual
# strip + marginal histogram), and the octoplot assembly. Observation types
# that implement the plotting protocol get their panels with no code here;
# composite/catalog observations opt out via `defaultpanels` and plug their
# bespoke drawing into the same layout slots and the same PosteriorSeries.
#
# Layout invariant: every panel is a two-column GridLayout — content in
# column 1, side elements (colorbar, marginal histogram) in a Fixed-width
# column 2 — so the right edges of all main axes align down the figure.
# ---------------------------------------------------
module OctofitterMakieExt

using Octofitter
using Octofitter: PosteriorSeries, ObservableQuery, PlotChannel, OctoPlotResult,
    plotchannels, obscontext, modelcurves, default_sky_queries,
    plottable_observations, likelihoodname, refspecs, defaultpanels,
    rowsignal, evalsignal, foldablerows, foldephemeris, foldphase,
    _query, _querystr, _refstr, _solvekw, _isprior
using PlanetOrbits
using PlanetOrbits: plotinfo, plotlabel, _names, _setnames
using Makie
using Random: Xoshiro
using Statistics
using StatsBase: fit, Histogram, ProbabilityWeights
using StatsBase
using MCMCChains: Chains

const WONG = Makie.wong_colors()
_rowcolor(k) = WONG[mod1(k, length(WONG))]
# The draw-ramp: row accent colour fading toward near-white with orbital
# phase; the light endpoint is the one v1 used in octoplot.
_phaseramp(c) = Makie.cgrad([Makie.to_color(c), Makie.to_color("#FAFAFA")])

# Side column (colorbars, marginal histograms): one Fixed width everywhere,
# so column-1 boundaries — and with them every main axis edge — line up
# across panels.
const SIDE_W = 80.0
const COLGAP = 8.0

_layout(gl::Makie.GridLayout) = gl
_layout(gp::Union{Makie.GridPosition,Makie.GridSubposition}) = Makie.GridLayout(gp)

_alpha(n) = min(1.0, 100 / max(n, 1))

# NaN-joined concatenation: one lines! call per row instead of one per draw.
function _nanjoin!(xs, ys, cs, x, y, c)
    append!(xs, x); push!(xs, NaN)
    append!(ys, y); push!(ys, NaN)
    append!(cs, c); push!(cs, last(c))
    return nothing
end

_isangular(sys::Octofitter.System) = !isempty(sys.framevars)

const _PHASE_TICKS = ([0, π / 2, π, 3π / 2, 2π], ["0", "π/2", "π", "3π/2", "2π"])

# The accent colour for a channel group: the hierarchy row whose exterior is
# the query's target (a planet's astrometry gets its sky-track colour);
# star-side queries (RV reflex) fall back to the first Wong colour.
function _querycolor(sys, q)
    q === nothing && return WONG[1]
    tn = Octofitter._leafnames(sys, q.target)
    tn === nothing && return WONG[1]
    for (k, (_, ext, _)) in enumerate(sys.rows)
        ext == tn && return _rowcolor(k)
    end
    return WONG[1]
end

# ---------------------------------------------------
# Display wrapping (position angles)
#
# The v1 astromtimeplot treatment, with robust jump detection: values are
# wrapped into [0, w); wherever consecutive samples jump across the
# boundary, the line exits through the axis edge at the midpoint with its
# incoming slope, breaks with a NaN, and re-enters from the opposite edge —
# no vertical jump line, no gap. Callers set the y limits to exactly
# (0, w) when `wrapped` comes back true.
# ---------------------------------------------------
function _wrap_series(x, y, w)
    xo = Float64[]; yo = Float64[]
    sizehint!(xo, length(x) + 8); sizehint!(yo, length(y) + 8)
    wrapped = false
    yprev = NaN
    for i in eachindex(x)
        yi = mod(y[i], w)
        if !isempty(yo) && isfinite(yi) && isfinite(yprev)
            δ = yi - yprev
            if abs(δ) > w / 2
                wrapped = true
                Δ = δ - sign(δ) * w          # the true (small) increment across the boundary
                xm = (x[i-1] + x[i]) / 2
                push!(xo, xm); push!(yo, yprev + Δ / 2)   # exits past the edge
                push!(xo, NaN); push!(yo, NaN)
                push!(xo, xm); push!(yo, yi - Δ / 2)      # re-enters from the other edge
            end
        end
        push!(xo, x[i]); push!(yo, yi)
        yprev = yi
    end
    return xo, yo, wrapped
end

# ---------------------------------------------------
# Sky panel
# ---------------------------------------------------

"""
    skypanel!(gp, series; tracks=nothing, ntrack=150, colorbar=true, ndraws=nothing)

Sky-plane orbit tracks. By default one track per hierarchy row — the
exterior side relative to the interior side, i.e. exactly the relationship
the row parametrizes — sampled uniformly in eccentric anomaly per draw and
coloured by orbital phase in the row's accent colour. Relative-astrometry
data overlay automatically. Returns `(; sky=ax)`.
"""
function Octofitter.skypanel!(gp, series::PosteriorSeries;
                              tracks=nothing, ntrack=150, colorbar=true,
                              ndraws=nothing)
    sys = series.model.system
    gs = _layout(gp)
    angular = _isangular(sys)
    info_x, info_y = angular ? (plotinfo(raoff), plotinfo(decoff)) :
                     (plotinfo(posx), plotinfo(posy))

    nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
    tstart = isempty(series.data_epochs) ? series.ts[1] : series.data_epochs[1]

    # Collect per-row NaN-joined tracks first: the unit auto-switch (mas ↔
    # arcsec) and the equal-span limits need the full extent before anything
    # draws.
    rowdata = Vector{NTuple{3,Vector{Float64}}}()
    rownames = Symbol[]
    qs = tracks === nothing ? default_sky_queries(sys) :
         [(Octofitter._query(t), Symbol("track$k")) for (k, t) in enumerate(tracks)]
    for (k, (q, owner)) in enumerate(qs)
        xs = Float64[]; ys = Float64[]; cs = Float64[]
        for d in 1:nshow
            posys = series.systems[d]
            ts = tracks === nothing ?
                 PlanetOrbits.orbit_track_epochs(posys, k; n=ntrack, tstart) :
                 series.ts
            traj = tracks === nothing ?
                   orbitsolve(posys, ts; _solvekw(sys)...) : series.trajs[d]
            t = Octofitter.resolveref(posys, q.target)
            r = Octofitter.resolveref(posys, q.ref)
            if angular
                x = [raoff(sol, t, r) for sol in traj]
                y = [decoff(sol, t, r) for sol in traj]
            else
                x = [posx(sol, t, r) for sol in traj]
                y = [posy(sol, t, r) for sol in traj]
            end
            rowk = tracks === nothing ? k : min(k, PlanetOrbits.nrows(posys))
            c = [PlanetOrbits.orbit_phase(posys, rowk, tt) for tt in ts]
            _nanjoin!(xs, ys, cs, x, y, c)
        end
        push!(rowdata, (xs, ys, cs))
        push!(rownames, owner)
    end

    # Data overlay coordinates (gathered before drawing so limits cover them).
    datasets = NTuple{4,Vector{Float64}}[]
    for obs in plottable_observations(sys)
        obs isa Octofitter.RelAstromObs || continue
        tab = obs.table
        if hasproperty(tab, :sep)
            x = tab.sep .* sin.(tab.pa)
            y = tab.sep .* cos.(tab.pa)
            σx = hypot.(tab.σ_sep .* abs.(sin.(tab.pa)), tab.sep .* tab.σ_pa .* abs.(cos.(tab.pa)))
            σy = hypot.(tab.σ_sep .* abs.(cos.(tab.pa)), tab.sep .* tab.σ_pa .* abs.(sin.(tab.pa)))
        else
            x, y = tab.ra, tab.dec
            σx, σy = tab.σ_ra, tab.σ_dec
        end
        push!(datasets, (collect(float.(x)), collect(float.(y)),
                         collect(float.(σx)), collect(float.(σy))))
    end

    # mas → arcsec switch, v1 threshold.
    extent = maximum((isempty(xs) ? 0.0 : maximum(abs, filter(!isnan, xs)) for (xs, _, _) in rowdata); init=0.0)
    unitscale = (angular && extent > 1500) ? 1e-3 : 1.0
    unit = angular ? (unitscale == 1.0 ? "mas" : "arcsec") : "au"

    # Equal-span limits about the joint extent (origin included: the
    # reference point's star marker must stay in frame). With equal spans
    # and a square cell (Aspect row sizing below), DataAspect is satisfied
    # with the axis filling its cell — so the colorbar matches its height
    # and the panel's right edge lines up with the time panels below.
    xmin = xmax = ymin = ymax = 0.0
    for (xs, ys, _) in rowdata, i in eachindex(xs)
        isnan(xs[i]) && continue
        xmin = min(xmin, xs[i]); xmax = max(xmax, xs[i])
        ymin = min(ymin, ys[i]); ymax = max(ymax, ys[i])
    end
    for (x, y, σx, σy) in datasets
        isempty(x) && continue
        xmin = min(xmin, minimum(x .- σx)); xmax = max(xmax, maximum(x .+ σx))
        ymin = min(ymin, minimum(y .- σy)); ymax = max(ymax, maximum(y .+ σy))
    end
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    half = 1.05 * max(xmax - xmin, ymax - ymin, 1e-12) / 2

    # DataAspect keeps the sky square in data units whatever the cell shape;
    # right-aligning pins the axis' right edge to the column boundary, so it
    # lines up with the time panels below even when the cell is not square.
    ax = Makie.Axis(gs[1, 1];
        xlabel="$(info_x.label) [$unit]",
        ylabel="$(info_y.label) [$unit]",
        aspect=Makie.DataAspect(), halign=:right,
        xreversed=info_x.flip,
        xgridvisible=false, ygridvisible=false)
    Makie.xlims!(ax, (cx - half) * unitscale, (cx + half) * unitscale)
    Makie.ylims!(ax, (cy - half) * unitscale, (cy + half) * unitscale)

    α = _alpha(nshow)
    for (k, (xs, ys, cs)) in enumerate(rowdata)
        Makie.lines!(ax, xs .* unitscale, ys .* unitscale;
            color=cs, colormap=_phaseramp(_rowcolor(k)), colorrange=(0, 2π),
            alpha=α, linewidth=0.6)
    end

    # Data overlay: uncertainty crosses; ellipses come later.
    for (x, y, σx, σy) in datasets
        Makie.errorbars!(ax, x .* unitscale, y .* unitscale, σy .* unitscale;
            direction=:y, color=:black, linewidth=1)
        Makie.errorbars!(ax, x .* unitscale, y .* unitscale, σx .* unitscale;
            direction=:x, color=:black, linewidth=1)
        Makie.scatter!(ax, x .* unitscale, y .* unitscale;
            color=:white, strokecolor=:black, strokewidth=2, markersize=8)
    end

    # The reference point of every track sits at the origin.
    Makie.scatter!(ax, [0.0], [0.0];
        marker='★', markersize=20, color=:white, strokecolor=:black, strokewidth=1.5)

    if colorbar
        cb = Makie.Colorbar(gs[1, 2];
            colormap=_phaseramp(_rowcolor(1)), colorrange=(0, 2π),
            width=12, halign=:left,
            ticks=_PHASE_TICKS, label="orbital phase (mean anomaly) →")
        # Match the colorbar's height to the axis itself, not the cell —
        # under DataAspect the axis can be smaller than its cell.
        Makie.on(ax.scene.viewport; update=true) do vp
            cb.height[] = Makie.Fixed(Makie.widths(vp)[2])
        end
    else
        Makie.Box(gs[1, 2]; visible=false)
    end
    Makie.colsize!(gs, 2, Makie.Fixed(SIDE_W))
    Makie.colgap!(gs, 1, COLGAP)
    return (; sky=ax)
end

# ---------------------------------------------------
# Residual strip + marginal histogram (shared by time and phase panels)
# ---------------------------------------------------

# Whitened residuals per posterior draw: (resid / σ_eff) with each draw's
# own parameters (jitters vary by draw), from the likelihood's residuals.
function _zscores(series, obs, ch, nshow)
    r1 = Octofitter.residuals(obs, obscontext(series, obs; draw=1))[ch.name]
    zs = Matrix{Float64}(undef, length(r1.resid), nshow)
    zs[:, 1] = r1.resid ./ r1.σ_eff
    for d in 2:nshow
        r = Octofitter.residuals(obs, obscontext(series, obs; draw=d))[ch.name]
        zs[:, d] = r.resid ./ r.σ_eff
    end
    return zs
end

# Draw the residual strip and marginal histogram for `entries` into
# pre-created axes. `xmap(epochs) -> x` positions points (identity for time
# panels, the fold phase for phase panels). In whitened mode each point is
# the median z-score over draws with a 16–84 % bar, ±1σ guides, and the
# histogram is density-normalized with a unit normal overlaid; in raw mode
# it is the MAP residual with measurement/jitter errorbars and a count
# histogram. Histogram bins are shared across instruments either way.
function _residstrip!(ax_resid, ax_hist, series, entries, xmap, nshow, whiten)
    Makie.hlines!(ax_resid, 0.0; color=:black, linewidth=1)
    whiten && Makie.hlines!(ax_resid, [-1.0, 1.0]; color=(:black, 0.3), linewidth=0.5)

    histdata = Tuple{Any,Vector{Float64}}[]   # (colour, values)
    for (j, (obs, ch)) in enumerate(entries)
        obs === nothing && continue
        r = Octofitter.residuals(obs, obscontext(series, obs))[ch.name]
        u = r.use
        x = xmap(r.epoch)
        color = length(entries) == 1 ? :black : WONG[mod1(j + 1, length(WONG))]
        mcolor = length(entries) == 1 ? :white : color
        if whiten
            zs = _zscores(series, obs, ch, nshow)
            med = [median(view(zs, i, :)) for i in axes(zs, 1)]
            lo = [quantile(view(zs, i, :), 0.16) for i in axes(zs, 1)]
            hi = [quantile(view(zs, i, :), 0.84) for i in axes(zs, 1)]
            Makie.rangebars!(ax_resid, x[u], lo[u], hi[u];
                color=length(entries) == 1 ? "#AAAAAA" : color, linewidth=1)
            Makie.scatter!(ax_resid, x[u], med[u];
                color=mcolor, strokecolor=:black, strokewidth=1.5, markersize=6)
            append!(histdata, ((color isa Symbol ? :black : color, vec(zs[u, :])),))
        else
            Makie.errorbars!(ax_resid, x[u], r.resid[u], r.σ_eff[u]; color="#CCCCCC", linewidth=1)
            Makie.errorbars!(ax_resid, x[u], r.resid[u], r.σ[u]; color=:black, linewidth=2)
            Makie.scatter!(ax_resid, x[u], r.resid[u];
                color=mcolor, strokecolor=:black, strokewidth=1.5, markersize=6)
            append!(histdata, ((color isa Symbol ? :black : color, r.resid[u]),))
        end
    end

    if ax_hist !== nothing && !isempty(histdata)
        pooled = reduce(vcat, (v for (_, v) in histdata))
        isempty(pooled) && return nothing
        h_all = fit(Histogram, pooled)
        edges = h_all.edges[1]
        centers = collect(edges[1:end-1] .+ step(edges) / 2)
        Makie.linkyaxes!(ax_resid, ax_hist)
        for (color, v) in histdata
            h = fit(Histogram, v, edges)
            w = whiten ? h.weights ./ (length(v) * step(edges)) : float.(h.weights)
            Makie.stairs!(ax_hist, w, centers; step=:center, color=color, linewidth=1.0)
        end
        if whiten
            # The unit normal the z-scores should follow.
            g = range(min(edges[1], -3), max(edges[end], 3), length=101)
            Makie.lines!(ax_hist, exp.(-g .^ 2 ./ 2) ./ sqrt(2π), g;
                color=(:black, 0.8), linewidth=1.2, linestyle=:dash)
        end
        Makie.hlines!(ax_hist, 0.0; color=:black, linewidth=1)
        Makie.xlims!(ax_hist; low=0)
    end
    return nothing
end

_stripylabel(whiten) = whiten ? "resid / σ" : "residuals"

# ---------------------------------------------------
# Time-series panel
# ---------------------------------------------------

"""
    timeseriespanel!(gp, series, entries; kwargs...)

Generic data-vs-model panel for one channel group. `entries` is a vector of
`(obs, channel)` pairs sharing one quantity (e.g. several RV instruments);
pass `(nothing, channel)` for a pure model panel with no data.

Layout: main axis (posterior model curves + calibrated data), residual
strip glued below, marginal residual histogram at right with bins shared
across instruments. `whiten=true` replaces the strip and histogram with
whitened residuals — `(data − model)/σ_eff` per posterior draw, summarized
per point as median and 16–84 % interval, histogram against a unit normal;
the default whitens whenever more than one draw is shown. The epoch axis is
calendar-dated via `MJDConversion` and stays correct under zoom; annotate
with `Date`/`DateTime` or MJD alike. The time axes clip exactly to the
model grid (no gap after the last segment); wrapped quantities (position
angle) get the v1 treatment — lines exit one axis edge and re-enter at the
other, with limits pinned to the wrap range.

Returns `(; main, resid, hist)` (the latter two `nothing` when no data).
"""
function Octofitter.timeseriespanel!(gp, series::PosteriorSeries, entries;
                                     top_time_axis=true, bottom_time_axis=true,
                                     show_hist=true, ndraws=nothing, whiten=nothing,
                                     curvecolor=nothing, linewidth=nothing)
    gs = _layout(gp)
    entries = [(obs, ch) for (obs, ch) in entries]
    isempty(entries) && error("timeseriespanel! needs at least one (obs, channel) entry")
    ch1 = entries[1][2]
    ylabel = isempty(ch1.unit) ? String(ch1.label) : "$(ch1.label) [$(ch1.unit)]"
    hasdata = any(obs !== nothing for (obs, _) in entries)
    nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
    whiten = something(whiten, nshow > 1) & hasdata
    curvecolor = something(curvecolor, _querycolor(series.model.system, ch1.query))
    linewidth = something(linewidth, nshow == 1 ? 1.3 : 0.4)

    conv() = PlanetOrbits.MJDConversion()
    ax = Makie.Axis(gs[1, 1];
        ylabel, dim1_conversion=conv(),
        xaxisposition=:top,
        xticklabelsvisible=top_time_axis, xticksvisible=top_time_axis,
        xgridvisible=false, ygridvisible=false)

    ax_resid = ax_hist = nothing
    if hasdata
        ax_resid = Makie.Axis(gs[2, 1];
            ylabel=_stripylabel(whiten), dim1_conversion=conv(),
            xlabel=bottom_time_axis ? "epoch" : "",
            xticklabelsvisible=bottom_time_axis, xticksvisible=bottom_time_axis,
            xgridvisible=false, ygridvisible=false)
        Makie.linkxaxes!(ax, ax_resid)
        Makie.rowgap!(gs, 1, 0)
        Makie.rowsize!(gs, 1, Makie.Auto(2.4))
    end

    # Posterior model curves (one NaN-joined lines! call for all draws).
    q = ch1.query
    if q !== nothing
        curves = modelcurves(series, q)
        xs = Float64[]; ys = Float64[]
        sizehint!(xs, (length(series.ts) + 1) * nshow)
        sizehint!(ys, (length(series.ts) + 1) * nshow)
        anywrapped = false
        for d in 1:nshow
            v = curves[d] .* ch1.scale
            if ch1.wrap !== nothing
                xw, yw, wrapped = _wrap_series(series.ts, v, ch1.wrap)
                anywrapped |= wrapped
                append!(xs, xw); append!(ys, yw)
            else
                append!(xs, series.ts); append!(ys, v)
            end
            push!(xs, NaN); push!(ys, NaN)
        end
        Makie.lines!(ax, xs, ys; color=(curvecolor, _alpha(nshow)), linewidth)
        anywrapped && Makie.ylims!(ax, 0.0, ch1.wrap)
    end

    # Data on the main axis, per instrument, MAP-calibrated.
    for (j, (obs, ch)) in enumerate(entries)
        obs === nothing && continue
        r = Octofitter.residuals(obs, obscontext(series, obs))[ch.name]
        u = r.use
        color = length(entries) == 1 ? :black : WONG[mod1(j + 1, length(WONG))]
        # Averaging windows, when the channel declares them. A catalog proper
        # motion is not measured *at* an epoch, it is averaged over a mission
        # span, and the bar is the only thing that says which span — it is why
        # a Hipparcos point may sit far from a curve the Gaia points hug. (v1
        # `absastromplot`'s one genuinely load-bearing idiom.)
        if haskey(r, :epoch_lo) && haskey(r, :epoch_hi)
            Makie.errorbars!(ax, r.epoch[u], r.data[u],
                (r.epoch .- r.epoch_lo)[u], (r.epoch_hi .- r.epoch)[u];
                direction=:x, color=:black, linewidth=0.7, whiskerwidth=4)
        end
        Makie.errorbars!(ax, r.epoch[u], r.data[u], r.σ_eff[u]; color="#CCCCCC", linewidth=1)
        Makie.errorbars!(ax, r.epoch[u], r.data[u], r.σ[u]; color=:black, linewidth=2)
        Makie.scatter!(ax, r.epoch[u], r.data[u];
            color=length(entries) == 1 ? :white : color,
            strokecolor=:black, strokewidth=1.5, markersize=7,
            label=likelihoodname(obs))
    end
    count(e -> e[1] !== nothing, entries) > 1 &&
        Makie.axislegend(ax; framevisible=false, labelsize=9, padding=(4, 4, 2, 2))

    if ax_resid !== nothing && show_hist
        ax_hist = Makie.Axis(gs[2, 2]; xgridvisible=false, ygridvisible=false)
        Makie.hidedecorations!(ax_hist)
        Makie.colsize!(gs, 2, Makie.Fixed(SIDE_W))
        Makie.colgap!(gs, 1, COLGAP)
    elseif hasdata
        Makie.Box(gs[2, 2]; visible=false)
        Makie.colsize!(gs, 2, Makie.Fixed(SIDE_W))
        Makie.colgap!(gs, 1, COLGAP)
    else
        Makie.Box(gs[1, 2]; visible=false)
        Makie.colsize!(gs, 2, Makie.Fixed(SIDE_W))
        Makie.colgap!(gs, 1, COLGAP)
    end
    hasdata && _residstrip!(ax_resid, ax_hist, series, entries, identity, nshow, whiten)

    # Clip exactly to the model grid: no gap between the last line segment
    # and the axis (the v1 rule).
    Makie.xlims!(ax, first(series.ts), last(series.ts))

    return (; main=ax, resid=ax_resid, hist=ax_hist)
end

# ---------------------------------------------------
# Phase-folded panel
# ---------------------------------------------------

"""
    phasefoldpanel!(gp, series, entries; row=nothing, kwargs...)

`entries` folded on hierarchy row `row`'s orbital phase (a row index, a row
owner name, or `nothing` when only one row can fold this quantity). What is
plotted is the *row's own signal* (`rowsignal`): per draw, the isolated
contribution of that row to the channel's query — e.g. one planet's RV
reflex with the other planets' signals and any perspective acceleration
removed exactly, via the reference grammar rather than plot-side arithmetic.

Fold conventions: each draw's model curve is folded on **its own** period
and phase zero, so every curve is a clean single cycle and the family's
spread shows the fold uncertainty honestly; the data (MAP residuals plus
the row's MAP signal) folds once, on the MAP ephemeris. Phase zero is the
signal's upward zero crossing for radial velocity (the v1 rvpostplot
convention), periastron otherwise. Under the N-body propagator the curve is
the actual model over one osculating period anchored at the data midpoint —
periodicity is then an approximation, which is why `octoplot` only adds
phase panels automatically under `KeplerianApprox`.

Noise-weighted binned means (grey-to-red points) pool all instruments;
the residual strip below shows the same residuals as the time panel against
phase (whitened by default when more than one draw is shown).

Returns `(; main, resid, hist)`.
"""
function Octofitter.phasefoldpanel!(gp, series::PosteriorSeries, entries;
                                    row=nothing, ndraws=nothing, whiten=nothing,
                                    show_binned=nothing, nbins=20, nphase=201,
                                    bottom_axis=true, show_hist=true,
                                    curvecolor=nothing, linewidth=nothing)
    sys = series.model.system
    gs = _layout(gp)
    entries = [(obs, ch) for (obs, ch) in entries]
    isempty(entries) && error("phasefoldpanel! needs at least one (obs, channel) entry")
    ch1 = entries[1][2]
    ch1.query === nothing && error(
        "channel :$(ch1.name) has no observable query, so there is no orbit to fold on")
    q = ch1.query

    rows = foldablerows(series.sys_map, q)
    isempty(rows) && error("no hierarchy row's contribution to $(_querystr(q)) can be isolated")
    k = if row === nothing
        length(rows) == 1 || error(
            "several rows can fold $(_querystr(q)); pass row= one of $(join(("$(j) ($(sys.rows[j][1]))" for j in rows), ", "))")
        rows[1]
    elseif row isa Symbol
        j = findfirst(r -> r[1] === row, sys.rows)
        j === nothing && error("no hierarchy row is owned by :$row")
        j
    else
        Int(row)
    end
    k in rows || error("row $k's contribution to $(_querystr(q)) cannot be isolated")
    sig = rowsignal(series.sys_map, q, k)
    rowname = sys.rows[k][1]

    nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
    hasdata = any(obs !== nothing for (obs, _) in entries)
    whiten = something(whiten, nshow > 1) & hasdata
    show_binned = something(show_binned, ch1.wrap === nothing)
    curvecolor = something(curvecolor, _rowcolor(k))
    linewidth = something(linewidth, nshow == 1 ? 1.3 : 0.4)
    kw = _solvekw(sys)

    tmid = if !isempty(series.data_epochs)
        sum(extrema(series.data_epochs)) / 2
    else
        sum(extrema(series.ts)) / 2
    end
    P0, t00 = foldephemeris(sig, series.sys_map, tmid; solvekw=kw)

    ylabel = isempty(ch1.unit) ? String(ch1.label) : "$(ch1.label) [$(ch1.unit)]"
    ax = Makie.Axis(gs[1, 1];
        ylabel, xticks=-0.5:0.25:0.5,
        xticklabelsvisible=false, xticksvisible=false,
        xgridvisible=false, ygridvisible=false)

    ax_resid = ax_hist = nothing
    if hasdata
        ax_resid = Makie.Axis(gs[2, 1];
            ylabel=_stripylabel(whiten),
            xlabel=bottom_axis ? "orbital phase" : "",
            xticks=-0.5:0.25:0.5,
            xticklabelsvisible=bottom_axis, xticksvisible=bottom_axis,
            xgridvisible=false, ygridvisible=false)
        Makie.linkxaxes!(ax, ax_resid)
        Makie.rowgap!(gs, 1, 0)
        Makie.rowsize!(gs, 1, Makie.Auto(2.4))
    end
    # The row this panel folds on, in its accent colour (the v1 rvpostplot
    # per-planet label).
    Makie.Label(gs[1, 1], String(rowname);
        font=:bold, fontsize=14, color=curvecolor,
        halign=:right, valign=:top, padding=(0, 8, 0, 4),
        tellwidth=false, tellheight=false)

    # Per-draw model curves, each on its own fold ephemeris.
    phases = collect(range(-0.5, 0.5, length=nphase))
    xs = Float64[]; ys = Float64[]
    anywrapped = false
    for d in 1:nshow
        posys = series.systems[d]
        P, t0 = foldephemeris(sig, posys, tmid; solvekw=kw)
        traj = orbitsolve(posys, t0 .+ P .* phases; kw...)
        v = evalsignal(sig, posys, traj) .* ch1.scale
        if ch1.wrap !== nothing
            xw, yw, wrapped = _wrap_series(phases, v, ch1.wrap)
            anywrapped |= wrapped
            append!(xs, xw); append!(ys, yw)
        else
            append!(xs, phases); append!(ys, v)
        end
        push!(xs, NaN); push!(ys, NaN)
    end
    Makie.lines!(ax, xs, ys; color=(curvecolor, _alpha(nshow)), linewidth)
    anywrapped && Makie.ylims!(ax, 0.0, ch1.wrap)

    # Data: MAP residuals + the row's MAP signal, folded on the MAP ephemeris.
    xmap(epochs) = foldphase.(epochs, P0, t00)
    if hasdata
        sig_at_data = evalsignal(sig, series.sys_map, series.data_traj_map)
        pooled_x = Float64[]; pooled_y = Float64[]; pooled_w = Float64[]
        for (j, (obs, ch)) in enumerate(entries)
            obs === nothing && continue
            r = Octofitter.residuals(obs, obscontext(series, obs))[ch.name]
            u = r.use
            s = sig_at_data[series.data_maps[obs]] .* ch.scale
            y = r.resid .+ s
            ch.wrap !== nothing && (y = mod.(y, ch.wrap))
            x = xmap(r.epoch)
            color = length(entries) == 1 ? :black : WONG[mod1(j + 1, length(WONG))]
            Makie.errorbars!(ax, x[u], y[u], r.σ_eff[u]; color="#CCCCCC", linewidth=1)
            Makie.errorbars!(ax, x[u], y[u], r.σ[u]; color=:black, linewidth=2)
            Makie.scatter!(ax, x[u], y[u];
                color=length(entries) == 1 ? :white : color,
                strokecolor=:black, strokewidth=1.5, markersize=5)
            append!(pooled_x, x[u]); append!(pooled_y, y[u])
            append!(pooled_w, 1 ./ r.σ_eff[u] .^ 2)
        end

        # Noise-weighted binned means (the rvpostplot idiom, with the bin
        # mask half-width fixed).
        if show_binned && !isempty(pooled_x)
            edges = range(-0.5, 0.5, length=nbins + 1)
            for i in 1:nbins
                m = (edges[i] .<= pooled_x) .& (pooled_x .< edges[i+1])
                count(m) == 0 && continue
                w = ProbabilityWeights(pooled_w[m])
                μ = mean(pooled_y[m], w)
                σ = count(m) > 1 ? std(pooled_y[m], w) : 0.0
                c = (edges[i] + edges[i+1]) / 2
                Makie.errorbars!(ax, [c], [μ], [σ]; color=:black, linewidth=2)
                Makie.scatter!(ax, [c], [μ];
                    color=:red, markersize=9, strokecolor=:black, strokewidth=1.5)
            end
        end
    end

    if ax_resid !== nothing && show_hist
        ax_hist = Makie.Axis(gs[2, 2]; xgridvisible=false, ygridvisible=false)
        Makie.hidedecorations!(ax_hist)
        Makie.colsize!(gs, 2, Makie.Fixed(SIDE_W))
        Makie.colgap!(gs, 1, COLGAP)
    else
        Makie.Box(gs[hasdata ? 2 : 1, 2]; visible=false)
        Makie.colsize!(gs, 2, Makie.Fixed(SIDE_W))
        Makie.colgap!(gs, 1, COLGAP)
    end
    hasdata && _residstrip!(ax_resid, ax_hist, series, entries, xmap, nshow, whiten)

    Makie.xlims!(ax, -0.5, 0.5)
    ax_resid !== nothing && Makie.xlims!(ax_resid, -0.5, 0.5)

    return (; main=ax, resid=ax_resid, hist=ax_hist)
end

# ---------------------------------------------------
# octoplot assembly
# ---------------------------------------------------

# Group plottable observation channels into panels: channels whose smooth
# curve is the same query share a panel (several RV instruments); channels
# without a query group by observation type + channel name. Observations
# with bespoke `defaultpanels` opt out of the generic grouping.

# Does this channel pass a `channels=` restriction? `nothing` keeps
# everything; a function matches the channel's smooth-curve observable
# (`channels=radvel`); a Symbol matches either the channel's own name or its
# observable's; a collection matches any of its elements.
_channelmatch(::PlotChannel, ::Nothing) = true
_channelmatch(ch::PlotChannel, f::Function) = ch.query !== nothing && ch.query.func === f
_channelmatch(ch::PlotChannel, s::Symbol) =
    ch.name === s || (ch.query !== nothing && nameof(ch.query.func) === s)
_channelmatch(ch::PlotChannel, cs) = any(c -> _channelmatch(ch, c), cs)

function _channelgroups(sys, channels=nothing)
    groups = Vector{Pair{Symbol,Vector{Tuple{Any,PlotChannel}}}}()
    for obs in plottable_observations(sys)
        isempty(defaultpanels(obs)) || continue
        for ch in plotchannels(obs)
            _channelmatch(ch, channels) || continue
            key = ch.query === nothing ?
                  Symbol(nameof(typeof(obs)), :_, ch.name) :
                  Symbol(nameof(ch.query.func), :_, _refstr(ch.query.target), :_, _refstr(ch.query.ref))
            i = findfirst(g -> g.first === key, groups)
            if i === nothing
                push!(groups, key => [(obs, ch)])
            else
                push!(groups[i].second, (obs, ch))
            end
        end
    end
    return groups
end

# Short, unique axis-group names for the result NamedTuple.
function _uniquename!(names, base)
    nm = base
    j = 2
    while nm in names
        nm = Symbol(base, :_, j)
        j += 1
    end
    push!(names, nm)
    return nm
end

_panelnames(groups) =
    (names = Symbol[]; [_uniquename!(names, entries[1][2].name) for (_, entries) in groups])

"""
    octoplot(model, chain; N=250, seed=0, show_sky=nothing, show_phase=nothing,
             whiten=nothing, channels=nothing, fname=nothing, figscale=1,
             ndraws=nothing)

See the core docstring (`?octoplot` without Makie loaded). Panels: a sky
panel when the system has angular observables; one time-series panel per
channel group; one phase-folded panel per (radial-velocity channel group ×
foldable hierarchy row) under `KeplerianApprox` (`show_phase` overrides);
then any bespoke `defaultpanels`. All epoch axes are linked; residual
strips are whitened by default when more than one draw is shown
(`whiten=false` for plain residuals).

`channels=` restricts the figure to some of the data: an observable function
(`channels=radvel` — what [`rvpostplot`](@ref) is), a channel or observable
name as a `Symbol`, or a collection of either. Bespoke panels are dropped
along with everything else not asked for.
"""
function Octofitter.octoplot(model::Octofitter.LogDensityModel, chain::Chains;
                             N=250, seed=0, kwargs...)
    series = PosteriorSeries(model, chain; N, seed)
    return Octofitter.octoplot(series; kwargs...)
end

function Octofitter.octoplot(series::PosteriorSeries;
                             show_sky=nothing, show_phase=nothing, whiten=nothing,
                             channels=nothing, figure=nothing,
                             fname=nothing, figscale=1.0, ndraws=nothing)
    sys = series.model.system
    # Say so when an observation's data cannot be drawn, so that an empty panel
    # reads as a gap in the plotting layer rather than as a fit with no data.
    Octofitter._warn_unplottable(sys)
    groups = _channelgroups(sys, channels)
    pnames = _panelnames(groups)
    do_sky = show_sky === nothing ?
             (_isangular(sys) && !isempty(sys.rows)) : show_sky

    # Phase panels: RV lore — every RV channel group, folded on each row
    # whose signal is isolable; only meaningful under a periodic propagator.
    do_phase = show_phase === nothing ?
               sys.method isa PlanetOrbits.KeplerianApprox : show_phase
    phasepanels = Tuple{Symbol,Int,Any}[]   # (name, row, entries)
    if do_phase
        names = copy(pnames)
        for ((key, entries), nm) in zip(groups, pnames)
            q = entries[1][2].query
            (q !== nothing && q.func === PlanetOrbits.radvel) || continue
            for k in foldablerows(series.sys_map, q)
                push!(phasepanels,
                    (_uniquename!(names, Symbol(nm, :_phase_, sys.rows[k][1])), k, entries))
            end
        end
    end
    # A `channels=` restriction is a request for one kind of panel, so the
    # bespoke ones (which are per-observation, not per-channel) drop out with
    # everything else that was not asked for.
    bespoke = channels === nothing ?
              [obs for obs in sys.observations if !isempty(defaultpanels(obs))] : []

    npanels = length(groups)
    W = round(Int, figscale * 620)
    # The sky cell aims square at the column-1 width (figure width minus
    # side column, gaps, and typical protrusions); DataAspect + edge
    # alignment in the panel absorb the estimate's error.
    skyH = W - 175
    H = round(Int, figscale * 40 + do_sky * figscale * skyH +
              figscale * (npanels * 260 + length(phasepanels) * 230 + 200 * length(bespoke)))
    # `figure=` draws into an existing (emptied) Figure instead of a new one —
    # what an animation needs, since a `VideoStream` records one scene.
    fig = figure === nothing ? Makie.Figure(size=(W, max(H, 300))) : figure

    axpairs = Pair{Symbol,Any}[]
    timeaxes = Makie.Axis[]
    row = 1
    if do_sky
        skyaxes = Octofitter.skypanel!(fig[row, 1], series; ndraws)
        push!(axpairs, :sky => skyaxes)
        Makie.rowsize!(fig.layout, row, Makie.Fixed(figscale * skyH))
        row += 1
    end
    for (i, ((key, entries), nm)) in enumerate(zip(groups, pnames))
        axs = Octofitter.timeseriespanel!(fig[row, 1], series, entries;
            top_time_axis=(i == 1), bottom_time_axis=(i == npanels), ndraws, whiten)
        push!(axpairs, nm => axs)
        push!(timeaxes, axs.main)
        axs.resid !== nothing && push!(timeaxes, axs.resid)
        Makie.rowsize!(fig.layout, row, Makie.Auto(1.0))
        row += 1
    end
    for (i, (nm, k, entries)) in enumerate(phasepanels)
        axs = Octofitter.phasefoldpanel!(fig[row, 1], series, entries;
            row=k, ndraws, whiten, bottom_axis=(i == length(phasepanels)))
        push!(axpairs, nm => axs)
        Makie.rowsize!(fig.layout, row, Makie.Auto(230 / 260))
        row += 1
    end
    for obs in bespoke, (nm, build) in defaultpanels(obs)
        axs = build(fig[row, 1], series)
        if axs isa NamedTuple && haskey(axs, :timeaxes)
            append!(timeaxes, collect(axs.timeaxes))
            axs = Base.structdiff(axs, (; timeaxes=nothing))
        end
        names = first.(axpairs)
        push!(axpairs, _uniquename!(names, nm) => axs)
        Makie.rowsize!(fig.layout, row, Makie.Auto(200 / 260))
        row += 1
    end
    length(timeaxes) > 1 && Makie.linkxaxes!(timeaxes...)

    axes = (; axpairs...)
    fname !== nothing && Makie.save(fname, fig)
    return OctoPlotResult(fig, axes, series)
end

# ---------------------------------------------------
# Bespoke panels for observations that are not epoch series
# ---------------------------------------------------

# A channel's per-draw quantity, as a (npoints × ndraws) matrix.
function _drawmatrix(series, obs, name::Symbol, field::Symbol, nshow)
    r1 = getproperty(Octofitter.residuals(obs, obscontext(series, obs; draw=1)), name)
    v1 = getproperty(r1, field)
    m = Matrix{Float64}(undef, length(v1), nshow)
    m[:, 1] = v1
    for d in 2:nshow
        r = getproperty(Octofitter.residuals(obs, obscontext(series, obs; draw=d)), name)
        m[:, d] = getproperty(r, field)
    end
    return m
end

_medband(m) = ([median(view(m, i, :)) for i in axes(m, 1)],
    [quantile(view(m, i, :), 0.16) for i in axes(m, 1)],
    [quantile(view(m, i, :), 0.84) for i in axes(m, 1)])

"""
    photometrypanel!(gp, series, obs; ndraws=nothing)

Photometry has no epoch axis and no orbital observable behind it — the model
is one body flux variable, constant across the table — so this replaces the
generic time-series panel rather than specialising it. The measurements sit
on a measurement-index axis with their errorbars; the posterior of the
modelled flux is the shaded band across them, with its marginal at the right.

Returns `(; main, hist)`.
"""
function Octofitter.photometrypanel!(gp, series::PosteriorSeries, obs; ndraws=nothing)
    gs = _layout(gp)
    ch = first(plotchannels(obs))
    r = Octofitter.residuals(obs, obscontext(series, obs))[ch.name]
    n = length(r.data)
    nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
    fluxes = vec(_drawmatrix(series, obs, ch.name, :model, nshow)[1, :])

    ylabel = isempty(ch.unit) ? String(ch.label) : "$(ch.label) [$(ch.unit)]"
    ax = Makie.Axis(gs[1, 1]; ylabel,
        xlabel="measurement", xticks=1:max(n, 1),
        xgridvisible=false, ygridvisible=false)
    color = _querycolor(series.model.system, nothing)

    lo, mid, hi = quantile(fluxes, 0.16), median(fluxes), quantile(fluxes, 0.84)
    Makie.hspan!(ax, lo, hi; color=(color, 0.25))
    Makie.hlines!(ax, mid; color=color, linewidth=1.5)
    Makie.errorbars!(ax, 1.0:n, r.data, r.σ; color=:black, linewidth=2)
    Makie.scatter!(ax, 1.0:n, r.data;
        color=:white, strokecolor=:black, strokewidth=1.5, markersize=8)
    Makie.xlims!(ax, 0.5, n + 0.5)

    ax_hist = Makie.Axis(gs[1, 2]; xgridvisible=false, ygridvisible=false)
    Makie.hidedecorations!(ax_hist)
    Makie.linkyaxes!(ax, ax_hist)
    h = fit(Histogram, fluxes)
    edges = h.edges[1]
    centers = collect(edges[1:end-1] .+ step(edges) / 2)
    Makie.stairs!(ax_hist, h.weights ./ (length(fluxes) * step(edges)), centers;
        step=:center, color=color, linewidth=1.0)
    Makie.xlims!(ax_hist; low=0)
    Makie.colsize!(gs, 2, Makie.Fixed(SIDE_W))
    Makie.colgap!(gs, 1, COLGAP)
    return (; main=ax, hist=ax_hist)
end

"""
    likemappanel!(gp, series, obs; ndraws=nothing)

A precomputed log-likelihood surface has no `(data, model, σ)` triple — it is
a reduction output, not a measurement — so what is drawn is the one thing
that *is* well defined per epoch: how far below that map's own maximum the
modelled position falls, summarized over the posterior draws (median, 16–84 %).
Zero means the orbit passes through that epoch's localization peak.

Returns `(; main, timeaxes)`.
"""
function Octofitter.likemappanel!(gp, series::PosteriorSeries, obs; ndraws=nothing)
    gs = _layout(gp)
    ch = first(plotchannels(obs))
    nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
    r = Octofitter.residuals(obs, obscontext(series, obs))[ch.name]
    m = _drawmatrix(series, obs, ch.name, :resid, nshow)
    med, lo, hi = _medband(m)

    ax = Makie.Axis(gs[1, 1];
        ylabel="Δ log-likelihood\nbelow map peak",
        xlabel="epoch", dim1_conversion=PlanetOrbits.MJDConversion(),
        xgridvisible=false, ygridvisible=false)
    Makie.hlines!(ax, 0.0; color=:black, linewidth=1)
    Makie.rangebars!(ax, r.epoch, lo, hi; color="#AAAAAA", linewidth=1)
    Makie.scatter!(ax, r.epoch, med;
        color=:white, strokecolor=:black, strokewidth=1.5, markersize=7,
        label=likelihoodname(obs))
    Makie.ylims!(ax; low=0)
    Makie.Box(gs[1, 2]; visible=false)
    Makie.colsize!(gs, 2, Makie.Fixed(SIDE_W))
    Makie.colgap!(gs, 1, COLGAP)
    return (; main=ax, timeaxes=(ax,))
end

# ---------------------------------------------------
# rvpostplot — the RV slice of octoplot
# ---------------------------------------------------

"""
    rvpostplot(model, chain; kwargs...)
    rvpostplot(series::PosteriorSeries; kwargs...)

The v1 radial-velocity summary figure, as a restriction of [`octoplot`](@ref)
rather than a separate implementation: `channels=radvel` keeps the RV channel
groups and drops everything else, `show_sky=false` drops the sky panel, and
what is left is exactly v1's anatomy — one time-series panel with a residual
strip and marginal histogram, then one phase-folded panel per planet.

Every keyword `octoplot` takes works here (`N`, `seed`, `ndraws`, `whiten`,
`show_phase`, `figscale`, `fname`), and the return value is the same
[`OctoPlotResult`](@ref), so `res.axes.rv.main` and `res.axes.rv_phase_b.main`
name the panels.
"""
Octofitter.rvpostplot(model::Octofitter.LogDensityModel, chain::Chains;
                      N=250, seed=0, kwargs...) =
    Octofitter.rvpostplot(PosteriorSeries(model, chain; N, seed); kwargs...)

function Octofitter.rvpostplot(series::PosteriorSeries; kwargs...)
    res = Octofitter.octoplot(series; show_sky=false,
        channels=PlanetOrbits.radvel, kwargs...)
    isempty(keys(res.axes)) && error(
        "rvpostplot: this model has no radial-velocity observations. " *
        "(`octoplot` draws whatever it does have.)")
    return res
end

"""
    rvpostplot_animated(model, chain; N=50, seed=0, framerate=4,
                        fname="rv-posterior.mp4", kwargs...)

[`rvpostplot`](@ref) recorded over `N` single-draw slices of the chain — v1's
"sweep through the posterior" animation. Each frame is the whole figure
rebuilt for one draw, so every panel (including the phase folds, which move
with the drawn period) is that draw's own. Returns `fname`.

The extension is anything Makie can write from a `VideoStream`: `.mp4`,
`.gif`, `.mkv`.
"""
function Octofitter.rvpostplot_animated(model::Octofitter.LogDensityModel, chain::Chains;
                                        N::Integer=50, seed::Integer=0, framerate=4,
                                        fname::AbstractString="rv-posterior.mp4",
                                        kwargs...)
    nsamples = size(chain, 1) * size(chain, 3)
    ii = sort!(StatsBase.sample(Xoshiro(seed), 1:nsamples, min(N, nsamples); replace=false))
    # One figure, reused: a `VideoStream` records a single scene, so each frame
    # empties the layout and rebuilds the panels into the same figure rather
    # than making a new one (which the stream would never see). The y limits
    # are pinned to the first frame's so the data do not jitter between draws.
    res = Octofitter.rvpostplot(PosteriorSeries(model, chain; ii=[first(ii)]); kwargs...)
    fig = res.figure
    stream = Makie.VideoStream(fig; framerate)
    Makie.recordframe!(stream)
    for i in ii[2:end]
        empty!(fig)
        Octofitter.rvpostplot(PosteriorSeries(model, chain; ii=[i]);
            figure=fig, kwargs...)
        Makie.recordframe!(stream)
    end
    Makie.save(fname, stream)
    return fname
end

# ---------------------------------------------------
# dotplot — mass vs separation/period posterior summary
# ---------------------------------------------------

"""
    dotplot(model, chain; mode=:separation, epoch=nothing, N=5000, seed=0,
            colormap=Makie.cgrad(:deep), colorbar=true, fname=nothing,
            figure=(;), axis=(;))

See the core docstring. A posterior summary with no data in it: one point per
draw per hierarchy row, mass against separation (or period), coloured by
eccentricity, with step-histogram marginals on both axes.

Masses are **M⊙**: v2 has one mass unit throughout, where v1 plotted each
planet's mass in Mⱼᵤₚ.
"""
function Octofitter.dotplot(model::Octofitter.LogDensityModel, chain::Chains;
                            mode::Symbol=:separation, epoch=nothing,
                            N::Integer=5000, seed::Integer=0,
                            colormap=Makie.cgrad(:deep), colorbar=true,
                            fname=nothing, figure=(;), axis=(;))
    mode in (:separation, :period) ||
        error("dotplot: `mode` must be :separation or :period, got :$mode")
    sys = model.system
    isempty(sys.rows) && error("dotplot: this system has no orbiting bodies.")

    nsamples = size(chain, 1) * size(chain, 3)
    rng = Xoshiro(seed)
    ii = N >= nsamples ? collect(1:nsamples) :
         sort!(StatsBase.sample(rng, 1:nsamples, N; replace=false))
    thetas = Octofitter.mcmcchain2result(model, chain, ii)
    systems = map(θ -> Octofitter.construct_system(model, θ), thetas)
    kw = _solvekw(sys)
    tmark = epoch === nothing ? nothing :
            (epoch isa Real ? Float64(epoch) : PlanetOrbits.mjd(epoch))

    fig = Makie.Figure(; size=(560, 460), figure...)
    gs = fig.layout
    xscale, yscale = log10, log10
    xlabel = if mode === :period
        "period [days]"
    elseif tmark === nothing
        "separation [au]"
    else
        "separation [au]\nat $(epoch)"
    end
    ax = Makie.Axis(gs[2, 1]; xscale, yscale, xlabel,
        ylabel=Makie.rich("mass [M", Makie.subscript("⊙"), "]"),
        xgridvisible=false, ygridvisible=false, axis...)
    ta = Makie.Axis(gs[1, 1]; xscale, ylabel="density",
        xgridvisible=false, ygridvisible=false)
    ra = Makie.Axis(gs[2, 2]; yscale, xlabel="density",
        xgridvisible=false, ygridvisible=false)

    hm = nothing
    # The instantaneous-separation mode needs the row's own exterior-vs-interior
    # query, which is exactly the sky panel's default track.
    skyq = default_sky_queries(sys)
    for (k, (owner, ext, _)) in enumerate(sys.rows)
        # The row's mass is the total mass of the bodies it places — for a
        # star+planets system that is the planet's own mass, and for a
        # hierarchy it is the mass of the whole exterior side, which is what
        # the row's Kepler problem actually contains.
        idx = [findfirst(==(n), sys.bodynames) for n in ext]
        ms = [sum(Float64(s.masses[j]) for j in idx) for s in systems]
        es = [Float64(s.rows[k].e) for s in systems]
        xs = if mode === :period
            [Float64(PlanetOrbits.period(s, k)) for s in systems]
        elseif tmark === nothing
            [abs(Float64(s.rows[k].a)) for s in systems]
        else
            q = skyq[k][1]
            map(systems) do s
                sol = only(orbitsolve(s, [tmark]; kw...))
                t = Octofitter.resolveref(s, q.target)
                rf = Octofitter.resolveref(s, q.ref)
                hypot(posx(sol, t, rf), posy(sol, t, rf), posz(sol, t, rf))
            end
        end
        keep = findall(j -> isfinite(xs[j]) && isfinite(ms[j]) && xs[j] > 0 && ms[j] > 0,
            eachindex(xs))
        isempty(keep) && continue
        # High-eccentricity draws first, so the low-e core stays visible.
        ord = keep[sortperm(es[keep], rev=true)]
        cm = length(sys.rows) == 1 ? colormap :
             Makie.cgrad([Makie.to_color("#FAFAFA"), Makie.to_color(_rowcolor(k))])
        hm = Makie.scatter!(ax, xs[ord], ms[ord];
            color=es[ord], colormap=cm, colorrange=(0, 1),
            markersize=clamp(2000 / length(ord), 0.9, 5), rasterize=4)
        c = length(sys.rows) == 1 ? :grey30 : (_rowcolor(k), 0.85)
        _stephist!(ta, xs[keep], log10; color=c)
        _stephist!(ra, ms[keep], log10; color=c, flip=true)
        length(sys.rows) > 1 && Makie.text!(ax, xs[ord[1]], ms[ord[1]];
            text=String(owner), color=_rowcolor(k), fontsize=11, font=:bold)
    end

    Makie.linkxaxes!(ax, ta)
    Makie.linkyaxes!(ax, ra)
    Makie.hidexdecorations!(ta)
    Makie.hideydecorations!(ra)
    Makie.colsize!(gs, 1, Makie.Auto(8))
    Makie.rowsize!(gs, 2, Makie.Auto(8))
    Makie.colgap!(gs, 1, 4)
    Makie.rowgap!(gs, 1, 4)
    colorbar && hm !== nothing &&
        Makie.Colorbar(gs[2, 3], hm; label="eccentricity")
    Makie.Label(gs[1, 1, Makie.Top()], String(sys.name);
        tellwidth=false, font=:bold, halign=:left)
    fname !== nothing && Makie.save(fname, fig, px_per_unit=3)
    return fig
end

# Log-spaced step histogram, drawn as an outline (v1's `stairs!` idiom — a
# filled `hist!` hides the other bodies' distributions when they overlap).
function _stephist!(ax, v, ::typeof(log10); color, flip::Bool=false)
    lo, hi = extrema(v)
    lo <= 0 && return nothing
    edges = 10 .^ range(log10(lo), log10(hi) + eps(), length=41)
    h = fit(Histogram, v, edges)
    w = h.weights ./ max(sum(h.weights), 1)
    if flip
        Makie.stairs!(ax, [w; 0], collect(edges); step=:pre, color, linewidth=2)
    else
        Makie.stairs!(ax, collect(edges), [w; 0]; step=:pre, color, linewidth=2)
    end
    return nothing
end

# ---------------------------------------------------
# Gaia DR4 figures
# ---------------------------------------------------

# The (single) observation of a given type in a model, with a message that
# says what to do when there is none or more than one.
function _theobs(sys, T, what)
    os = [o for o in sys.observations if o isa T]
    isempty(os) && error("$what needs a $(nameof(T)) in the model; this one has none.")
    length(os) > 1 && error(
        "$what draws one $(nameof(T)) and this model has $(length(os)). " *
        "Build a `PosteriorSeries` and call the panel directly for the one you want.")
    return only(os)
end

# A one-draw series: the MAP sample by default, so the bespoke single-sample
# figures below all agree on which draw they are showing.
function _onedraw(model, chain, sample_idx)
    i = sample_idx === nothing ?
        (haskey(chain, :logpost) ? argmax(vec(chain[:logpost])) :
         haskey(chain, :loglike) ? argmax(vec(chain[:loglike])) : 1) : Int(sample_idx)
    return PosteriorSeries(model, chain; ii=[i]), i
end

"""
    gaiastarplot(model, chain, sample_idx=MAP; keplerian_mult=1, ntrack=200,
                 fname=nothing, figure=(;), axis=(;))

The source's reflex track in the Gaia frame for one draw, with each transit's
along-scan residual re-projected into the sky plane along that transit's own
scan angle. Gaia constrains one direction per transit, so a measurement is a
*line*, not a point: each datum is drawn at the modelled position displaced
by `(data − model)` along `(sin ψ, cos ψ)`, with its 1σ tick along the same
direction and a dotted connector back to the track.

(v1 displaced by `model − data`, mirroring each point across the track. The
sign here is the one that makes a point sit where the abscissa says the source
was.)
"""
function Octofitter.gaiastarplot(model::Octofitter.LogDensityModel, chain::Chains,
                                 sample_idx=nothing;
                                 keplerian_mult=1.0, ntrack::Integer=200,
                                 fname=nothing, figure=(;), axis=(;))
    sys = model.system
    obs = _theobs(sys, Octofitter.GaiaDR4AstromObs, "gaiastarplot")
    series, _ = _onedraw(model, chain, sample_idx)
    posys = series.systems[1]
    kw = _solvekw(sys)

    q = ObservableQuery(PlanetOrbits.raoff, obs.target, obs.ref)
    qd = ObservableQuery(PlanetOrbits.decoff, obs.target, obs.ref)
    # One closed track, EA-spaced, for the row whose period dominates the
    # source's excursion — row 1 by convention, as in v1.
    ts = PlanetOrbits.orbit_track_epochs(posys, 1; n=ntrack,
        tstart=first(obs.table.epoch))
    traj = orbitsolve(posys, ts; kw...)
    tx = keplerian_mult .* Octofitter.evalquery(q, posys, traj)
    ty = keplerian_mult .* Octofitter.evalquery(qd, posys, traj)

    # The model at the data epochs, and the along-scan residual there.
    rows = series.data_maps[obs]
    dtraj = series.data_trajs[1]
    mx = keplerian_mult .* Octofitter.evalquery(q, posys, dtraj)[rows]
    my = keplerian_mult .* Octofitter.evalquery(qd, posys, dtraj)[rows]
    r = Octofitter.residuals(obs, obscontext(series, obs; draw=1)).along_scan
    s = sin.(obs.table.scan_pos_angle)
    c = cos.(obs.table.scan_pos_angle)
    u = r.use
    dx = mx .+ r.resid .* s
    dy = my .+ r.resid .* c

    fig = Makie.Figure(; size=(600, 600), figure...)
    ax = Makie.Axis(fig[1, 1];
        xlabel="Δα* [mas]", ylabel="Δδ [mas]",
        xreversed=true, aspect=Makie.DataAspect(),
        xgridvisible=false, ygridvisible=false, axis...)
    Makie.vlines!(ax, 0.0; color=(:grey, 0.5), linestyle=:dash, linewidth=1)
    Makie.hlines!(ax, 0.0; color=(:grey, 0.5), linestyle=:dash, linewidth=1)

    # 1σ along the scan direction, as NaN-separated segments.
    sx = Float64[]; sy = Float64[]
    cx = Float64[]; cy = Float64[]
    for i in eachindex(dx)
        u[i] || continue
        σ = r.σ_eff[i]
        append!(sx, (dx[i] - σ * s[i], dx[i] + σ * s[i], NaN))
        append!(sy, (dy[i] - σ * c[i], dy[i] + σ * c[i], NaN))
        append!(cx, (dx[i], mx[i], NaN))
        append!(cy, (dy[i], my[i], NaN))
    end
    Makie.lines!(ax, cx, cy; color=(:grey, 0.5), linestyle=:dot, linewidth=1)
    Makie.lines!(ax, tx, ty; color=WONG[1], linewidth=2)
    Makie.lines!(ax, sx, sy; color=:black, linewidth=1)
    Makie.scatter!(ax, dx[u], dy[u];
        color=WONG[2], strokecolor=:black, strokewidth=1.5, markersize=8)
    Makie.scatter!(ax, [0.0], [0.0];
        marker='★', markersize=20, color=:white, strokecolor=:black, strokewidth=1.5)
    fname !== nothing && Makie.save(fname, fig, px_per_unit=3)
    return fig
end

"""
    gaiatimeplot(model, chain; ndraws=nothing, fname=nothing, figure=(;), axis=(;))

Along-scan abscissa against time for a `GaiaDR4AstromObs`: the posterior cloud
of modelled abscissae over the measurements, and below it a per-epoch boxplot
of the residuals with the quoted formal error drawn at zero for comparison.

`octoplot` shows the same channel with a whitened residual strip, which asks
whether the residuals are normal. This asks a different question — at which
epochs is the *posterior spread* larger than the measurement error — and that
is why it is worth keeping as its own figure.
"""
function Octofitter.gaiatimeplot(model::Octofitter.LogDensityModel, chain::Chains;
                                 N=250, seed=0, ndraws=nothing,
                                 fname=nothing, figure=(;), axis=(;))
    sys = model.system
    obs = _theobs(sys, Octofitter.GaiaDR4AstromObs, "gaiatimeplot")
    series = PosteriorSeries(model, chain; N, seed)
    nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
    mm = _drawmatrix(series, obs, :along_scan, :model, nshow)
    r = Octofitter.residuals(obs, obscontext(series, obs)).along_scan
    u = r.use
    ep = r.epoch

    fig = Makie.Figure(; size=(700, 500), figure...)
    conv() = PlanetOrbits.MJDConversion()
    ax = Makie.Axis(fig[1, 1]; ylabel="along scan [mas]", dim1_conversion=conv(),
        xaxisposition=:top, xgridvisible=false, ygridvisible=false, axis...)
    ax_r = Makie.Axis(fig[2, 1]; ylabel="residuals [mas]", xlabel="epoch",
        dim1_conversion=conv(), xgridvisible=false, ygridvisible=false)
    Makie.linkxaxes!(ax, ax_r)
    Makie.rowgap!(fig.layout, 1, 0)
    Makie.rowsize!(fig.layout, 1, Makie.Auto(2.0))

    xs = Float64[]; ys = Float64[]
    for d in 1:nshow, i in eachindex(ep)
        u[i] || continue
        push!(xs, ep[i]); push!(ys, mm[i, d])
    end
    Makie.scatter!(ax, xs, ys; color=(WONG[1], _alpha(nshow)), markersize=3, rasterize=4)
    Makie.errorbars!(ax, ep[u], r.data[u], r.σ[u]; color=:black, linewidth=1)
    Makie.scatter!(ax, ep[u], r.data[u]; color=:black, markersize=5)

    Makie.hlines!(ax_r, 0.0; color=:black, linewidth=1)
    span = length(ep) > 1 ? (last(ep) - first(ep)) / max(length(ep), 1) : 40.0
    bx = Float64[]; by = Float64[]
    for i in eachindex(ep), d in 1:nshow
        u[i] || continue
        push!(bx, ep[i]); push!(by, mm[i, d] - r.data[i])
    end
    Makie.boxplot!(ax_r, bx, by; width=max(span, 1.0), color=(WONG[1], 0.6),
        strokecolor=:grey30, strokewidth=0.5, whiskercolor=:grey30, mediancolor=:grey30,
        show_outliers=false)
    Makie.errorbars!(ax_r, ep[u], zeros(count(u)), r.σ[u];
        color=:black, linewidth=1, whiskerwidth=5)
    fname !== nothing && Makie.save(fname, fig, px_per_unit=3)
    return fig
end

"""
    skytrackplot(model, chain, sample_idx=MAP; ra=nothing, dec=nothing,
                 gaia_id=nothing, ts=nothing, keplerian_mult=1, npoints=400,
                 fname=nothing, figure=(;), axis=(;))

The source's whole path on the sky for one draw: parallactic loops, proper
motion and the orbital wobble superimposed.

The parallax ellipse has to be projected onto a sky direction. It is taken
from the system's `ra`/`dec` frame variables when the model declares an
absolute frame; otherwise pass `ra=`/`dec=` in degrees, or `gaia_id=`, which
reads the published solution through `Octofitter.gaia_dr3_solution`. (v1 read
it from `obs.gaia_sol`, which the v2 observation no longer carries — it models
a sky path and does not need the catalog row.)

Needs the Earth's position, i.e. the DE440 ephemeris data dependency.
"""
function Octofitter.skytrackplot(model::Octofitter.LogDensityModel, chain::Chains,
                                 sample_idx=nothing;
                                 ra=nothing, dec=nothing, gaia_id=nothing,
                                 ts=nothing, keplerian_mult=1.0, npoints::Integer=400,
                                 fname=nothing, figure=(;), axis=(;))
    sys = model.system
    obs = _theobs(sys, Octofitter.GaiaDR4AstromObs, "skytrackplot")
    series, _ = _onedraw(model, chain, sample_idx)
    posys = series.systems[1]
    θ = series.thetas[1]
    key = Octofitter.normalizename(likelihoodname(obs))
    obsns = hasproperty(θ, :observations) ? θ.observations : (;)
    θ_obs = hasproperty(obsns, key) ? getproperty(obsns, key) : (;)

    α, δ = if ra !== nothing && dec !== nothing
        Float64(ra), Float64(dec)
    elseif gaia_id !== nothing
        cat = Octofitter.gaia_dr3_solution(; gaia_id)
        Float64(cat.ra), Float64(cat.dec)
    elseif hasproperty(θ, :ra) && hasproperty(θ, :dec)
        Float64(θ.ra), Float64(θ.dec)
    else
        error("""
        skytrackplot needs the source's sky direction to project the parallax
        ellipse onto, and this model has no absolute frame to take it from.
        Pass `ra=`/`dec=` in degrees, or `gaia_id=` to read them from the
        published Gaia DR3 solution.""")
    end

    epochs = collect(Float64, obs.table.epoch)
    tgrid = if ts !== nothing
        collect(Float64, ts)
    else
        t0, t1 = extrema(epochs)
        d = max(t1 - t0, 1.0)
        collect(range(t0 - 0.05d, t1 + 0.05d, length=npoints))
    end

    plx = Float64(θ.plx)
    ref_epoch = hasproperty(θ_obs, :ref_epoch) ? Float64(θ_obs.ref_epoch) : 0.0
    pmra_ = hasproperty(θ_obs, :pmra) ? Float64(θ_obs.pmra) : 0.0
    pmdec_ = hasproperty(θ_obs, :pmdec) ? Float64(θ_obs.pmdec) : 0.0
    ra0 = hasproperty(θ_obs, :ra_offset_mas) ? Float64(θ_obs.ra_offset_mas) : 0.0
    dec0 = hasproperty(θ_obs, :dec_offset_mas) ? Float64(θ_obs.dec_offset_mas) : 0.0

    # Parallax factors in RA and Dec, from the Earth's barycentric position.
    # `sind`/`cosd` because α and δ are degrees — v1 applied the radian
    # versions here while the likelihood used the degree ones.
    function pfactors(t)
        e = Octofitter.geocentre_position_query(t)
        return (e.x * sind(α) - e.y * cosd(α),
            e.x * cosd(α) * sind(δ) + e.y * sind(α) * sind(δ) - e.z * cosd(δ))
    end
    linear(t) = ((t - ref_epoch) / PlanetOrbits.year2day_julian)
    q = ObservableQuery(PlanetOrbits.raoff, obs.target, obs.ref)
    qd = ObservableQuery(PlanetOrbits.decoff, obs.target, obs.ref)

    traj = orbitsolve(posys, tgrid; _solvekw(sys)...)
    kx = keplerian_mult .* Octofitter.evalquery(q, posys, traj)
    ky = keplerian_mult .* Octofitter.evalquery(qd, posys, traj)
    px = similar(tgrid); py = similar(tgrid)
    for (i, t) in enumerate(tgrid)
        fa, fd = pfactors(t)
        px[i] = ra0 + pmra_ * linear(t) + plx * fa + kx[i]
        py[i] = dec0 + pmdec_ * linear(t) + plx * fd + ky[i]
    end

    rows = series.data_maps[obs]
    dtraj = series.data_trajs[1]
    dkx = keplerian_mult .* Octofitter.evalquery(q, posys, dtraj)[rows]
    dky = keplerian_mult .* Octofitter.evalquery(qd, posys, dtraj)[rows]
    r = Octofitter.residuals(obs, obscontext(series, obs; draw=1)).along_scan
    s = sin.(obs.table.scan_pos_angle)
    c = cos.(obs.table.scan_pos_angle)
    dx = similar(epochs); dy = similar(epochs)
    for (i, t) in enumerate(epochs)
        fa, fd = pfactors(t)
        dx[i] = ra0 + pmra_ * linear(t) + plx * fa + dkx[i] + r.resid[i] * s[i]
        dy[i] = dec0 + pmdec_ * linear(t) + plx * fd + dky[i] + r.resid[i] * c[i]
    end

    fig = Makie.Figure(; size=(700, 550), figure...)
    ax = Makie.Axis(fig[1, 1];
        xlabel="Δα* [mas]", ylabel="Δδ [mas]",
        xreversed=true, aspect=Makie.DataAspect(),
        xgridvisible=false, ygridvisible=false, axis...)
    Makie.lines!(ax, px, py; color=WONG[1], linewidth=1.5)
    u = r.use
    Makie.scatter!(ax, dx[u], dy[u]; color=:black, marker=:rect, markersize=5)
    fname !== nothing && Makie.save(fname, fig, px_per_unit=3)
    return fig
end

# ---------------------------------------------------
# hipparcosplot
# ---------------------------------------------------

"""
    hipparcosplot(model, chain, sample_idx=MAP; fname=nothing, figure=(;), axis=(;))

Hipparcos intermediate astrometry for one draw, in its own geometry.

Top: the signed along-scan residual against time, with the renormalized
per-transit σ. Bottom: the catalog's five-parameter sky path (grey), the
modelled path including the companion's perturbation (red), each transit's
abscissa line, and the residual segment from the model to that line with its
1σ tick.

Every part of the geometry is reconstructed so that the residual segment lands
*on* the abscissa line by construction, from the same `FrameOffset` the
likelihood uses. (v1 recovered the sign of an unsigned residual by trying both
perpendicular directions and keeping the one nearer the line, and then divided
the offsets by 3.6e6 — an arcsec→degree conversion in a figure whose units are
milliarcseconds.)

Accepts a `HipparcosIADObs`, or a `G23HObs` that keeps its `:iad_hip` channel.
"""
function Octofitter.hipparcosplot(model::Octofitter.LogDensityModel, chain::Chains,
                                  sample_idx=nothing;
                                  fname=nothing, figure=(;), axis=(;))
    sys = model.system
    cands = [o for o in sys.observations
             if o isa Octofitter.HipparcosIADObs ||
                (o isa Octofitter.G23HObs && :iad_hip ∈ o.table.kind)]
    isempty(cands) && error(
        "hipparcosplot needs a HipparcosIADObs, or a G23HObs keeping its " *
        ":iad_hip channel; this model has neither.")
    length(cands) > 1 && error("hipparcosplot draws one Hipparcos observation; " *
                               "this model has $(length(cands)).")
    obs = only(cands)
    series, _ = _onedraw(model, chain, sample_idx)
    ctx = obscontext(series, obs; draw=1)
    sim = Octofitter.simulate(obs, ctx)

    if obs isa Octofitter.HipparcosIADObs
        tab = obs.table
        chname = :along_scan
        Δα_src, Δδ_src = sim.Δα_mas, sim.Δδ_mas
        plx_anchor = hasproperty(ctx.θ_system, :plx) ? Float64(ctx.θ_system.plx) :
                     Float64(obs.hip_sol.plx)
    else
        tab = obs.hip_table
        chname = :along_scan_hip
        Δα_src, Δδ_src = sim.Δα_mas_hip, sim.Δδ_mas_hip
        # `G23HObs` anchors the abscissa channel's parallax on the *catalog*
        # value: there the frame is pure nuisance (see `_g23h_simulate!`).
        plx_anchor = Float64(obs.hip_sol.plx)
    end
    r = Octofitter.residuals(obs, ctx)[chname]
    off = Octofitter.frame_offset(ctx.θ_obs, plx_anchor, Float64)

    α₀, δ₀ = Float64(obs.hip_sol.radeg), Float64(obs.hip_sol.dedeg)
    cϕ = collect(Float64, tab.cosϕ); sϕ = collect(Float64, tab.sinϕ)
    Δt = (collect(Float64, tab.epoch) .- Octofitter.hipparcos_catalog_epoch_mjd) ./
         Octofitter.julian_year
    # The parallax factors in RA and Dec whose scan projection is exactly the
    # table's `parallaxFactorAlongScan`, so the model point below projects onto
    # the abscissa the likelihood computes.
    Pα = @. tab.x * sind(α₀) - tab.y * cosd(α₀)
    Pδ = @. tab.x * cosd(α₀) * sind(δ₀) + tab.y * sind(α₀) * sind(δ₀) - tab.z * cosd(δ₀)
    sx = collect(Float64, Δα_src)
    sy = collect(Float64, Δδ_src)
    mx = @. off.Δra + Δt * off.pmra + off.plx * Pα + sx
    my = @. off.Δdec + Δt * off.pmdec + off.plx * Pδ + sy
    # Foot of the residual on the measured abscissa line: along-scan is
    # (cosϕ, sinϕ), so a signed along-scan residual moves the point there.
    fx = mx .+ r.resid .* cϕ
    fy = my .+ r.resid .* sϕ

    fig = Makie.Figure(; size=(760, 620), figure...)
    ax_r = Makie.Axis(fig[1, 1];
        ylabel="along-scan residual [mas]", xlabel="epoch",
        dim1_conversion=PlanetOrbits.MJDConversion(),
        xgridvisible=false, ygridvisible=false)
    u = r.use
    Makie.hlines!(ax_r, 0.0; color=:black, linewidth=1)
    Makie.errorbars!(ax_r, r.epoch[u], r.resid[u], r.σ_eff[u]; color="#CCCCCC", linewidth=1)
    Makie.errorbars!(ax_r, r.epoch[u], r.resid[u], r.σ[u]; color=:black, linewidth=1)
    Makie.scatter!(ax_r, r.epoch[u], r.resid[u]; color=WONG[2], markersize=5)

    ax = Makie.Axis(fig[2, 1];
        xlabel="Δα* [mas]", ylabel="Δδ [mas]",
        aspect=Makie.DataAspect(), xgridvisible=false, ygridvisible=false, axis...)
    # Each transit's abscissa line: through the foot, along (−sinϕ, cosϕ).
    lx = Float64[]; ly = Float64[]
    ex = Float64[]; ey = Float64[]
    rx = Float64[]; ry = Float64[]
    Lh = 3 * maximum(abs, r.σ; init=1.0)
    for i in eachindex(fx)
        u[i] || continue
        append!(lx, (fx[i] - Lh * sϕ[i], fx[i] + Lh * sϕ[i], NaN))
        append!(ly, (fy[i] + Lh * cϕ[i], fy[i] - Lh * cϕ[i], NaN))
        append!(rx, (mx[i], fx[i], NaN))
        append!(ry, (my[i], fy[i], NaN))
        σ = r.σ[i]
        append!(ex, (fx[i] - σ * cϕ[i], fx[i] + σ * cϕ[i], NaN))
        append!(ey, (fy[i] - σ * sϕ[i], fy[i] + σ * sϕ[i], NaN))
    end
    Makie.lines!(ax, lx, ly; color=(:black, 0.12), linewidth=0.8)
    Makie.lines!(ax, ex, ey; color=(WONG[1], 0.35), linewidth=4)
    Makie.lines!(ax, rx, ry; color=WONG[1], linewidth=1)
    ord = sortperm(collect(Float64, tab.epoch))
    Makie.scatterlines!(ax, collect(Float64, tab.Δα✱)[ord], collect(Float64, tab.Δδ)[ord];
        color=(:grey, 0.6), linewidth=4, markersize=3, label="Hipparcos catalog path")
    Makie.scatterlines!(ax, mx[ord], my[ord];
        color=WONG[2], linewidth=1.5, markersize=3, label="model")
    Makie.axislegend(ax; framevisible=false, labelsize=10)
    # A high-proper-motion star's Hipparcos path is arcseconds long and
    # milliarcseconds wide, and `DataAspect` is what keeps that honest — so the
    # sky row is given most of the height and the strip the rest.
    Makie.rowsize!(fig.layout, 1, Makie.Auto(1))
    Makie.rowsize!(fig.layout, 2, Makie.Auto(1.6))
    fname !== nothing && Makie.save(fname, fig, px_per_unit=3)
    return fig
end

# ---------------------------------------------------
# Completeness maps
#
# A pure function of a `CompletenessMap` — no model, no chain, no likelihood —
# so it stands apart from the panel machinery above. Ported from v1 unchanged
# except for the mass unit: v2 has one mass unit (M⊙) throughout, where v1's
# per-planet masses were Jupiter masses.
# ---------------------------------------------------

function Octofitter.completenessplot(
    cmap::Octofitter.CompletenessMap,
    fname="completeness-map.png",
    args...;
    figure=(;),
    kwargs...
)
    fig = Makie.Figure(; size=(550, 450), figure...)
    Octofitter.completenessplot!(fig.layout, cmap, args...; kwargs...)
    fname !== nothing && Makie.save(fname, fig, px_per_unit=3)
    return fig
end

function Octofitter.completenessplot!(
    gs,
    cmap::Octofitter.CompletenessMap;
    xlabel="separation [AU]",
    ylabel=Makie.rich("mass [M", Makie.subscript("⊙"), "]"),
    colormap=Makie.cgrad(:viridis),
    colorrange=(0.0, 1.0),
    show_counts::Bool=false,
    title="",
    axis=(;),
)
    masses = cmap.masses
    seps = cmap.separations

    ax = Makie.Axis(gs[2, 1];
        xscale=log10, yscale=log10, xlabel, ylabel,
        xgridvisible=false, ygridvisible=false, axis...)

    hm = Makie.heatmap!(ax, _log_edges(seps), _log_edges(masses), cmap.completeness';
        colormap, colorrange)

    if show_counts
        for im in eachindex(masses), is in eachindex(seps)
            cmap.n_total[im, is] > 0 || continue
            Makie.text!(ax, seps[is], masses[im];
                text="$(cmap.n_detected[im, is])/$(cmap.n_total[im, is])",
                align=(:center, :center), fontsize=8,
                color=cmap.completeness[im, is] > 0.5 ? :black : :white)
        end
    end

    # Marginals. Summed counts rather than averaged completeness, so a cell
    # with no trials contributes nothing instead of contributing a zero — see
    # the note on `assemble_completeness` about `n_total == 0`.
    ta = Makie.Axis(gs[1, 1]; xscale=log10, ylabel="compl.",
        xgridvisible=false, ygridvisible=false)
    marginal_sep = vec(sum(cmap.n_detected, dims=1)) ./ vec(max.(sum(cmap.n_total, dims=1), 1))
    Makie.lines!(ta, seps, marginal_sep; color=:black, linewidth=2)
    Makie.band!(ta, seps, fill(0.0, length(seps)), marginal_sep; color=(:dodgerblue, 0.3))

    ra = Makie.Axis(gs[2, 2]; yscale=log10, xlabel="compl.",
        xgridvisible=false, ygridvisible=false)
    marginal_mass = vec(sum(cmap.n_detected, dims=2)) ./ vec(max.(sum(cmap.n_total, dims=2), 1))
    Makie.lines!(ra, marginal_mass, masses; color=:black, linewidth=2)
    Makie.band!(ra, fill(0.0, length(masses)), marginal_mass, masses; color=(:dodgerblue, 0.3))

    Makie.linkxaxes!(ax, ta)
    Makie.linkyaxes!(ax, ra)
    Makie.ylims!(ta, 0, 1.05)
    Makie.xlims!(ra, 0, 1.05)
    Makie.hidexdecorations!(ta)
    Makie.hideydecorations!(ra)

    Makie.colsize!(gs, 1, Makie.Auto(8))
    Makie.rowsize!(gs, 2, Makie.Auto(8))
    Makie.colgap!(gs, 1, 4)
    Makie.rowgap!(gs, 1, 4)
    Makie.Colorbar(gs[2, 3], hm; label="detection fraction")

    isempty(title) ||
        Makie.Label(gs[1, 1, Makie.Top()], title, tellwidth=false, font=:bold, halign=:left)

    return ax
end

# Bin edges from log-spaced centres: each interior edge is the geometric mean
# of its neighbours, the outer two extrapolated by the same ratio.
function _log_edges(centers::AbstractVector)
    n = length(centers)
    n == 1 && return [centers[1] / 2, centers[1] * 2]
    edges = Vector{Float64}(undef, n + 1)
    for i in 1:n-1
        edges[i+1] = sqrt(centers[i] * centers[i+1])
    end
    edges[1] = centers[1]^2 / centers[2]
    edges[end] = centers[end]^2 / centers[end-1]
    return edges
end

end # module
