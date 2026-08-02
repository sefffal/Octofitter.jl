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
function _channelgroups(sys)
    groups = Vector{Pair{Symbol,Vector{Tuple{Any,PlotChannel}}}}()
    for obs in plottable_observations(sys)
        isempty(defaultpanels(obs)) || continue
        for ch in plotchannels(obs)
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
             whiten=nothing, fname=nothing, figscale=1, ndraws=nothing)

See the core docstring (`?octoplot` without Makie loaded). Panels: a sky
panel when the system has angular observables; one time-series panel per
channel group; one phase-folded panel per (radial-velocity channel group ×
foldable hierarchy row) under `KeplerianApprox` (`show_phase` overrides);
then any bespoke `defaultpanels`. All epoch axes are linked; residual
strips are whitened by default when more than one draw is shown
(`whiten=false` for plain residuals).
"""
function Octofitter.octoplot(model::Octofitter.LogDensityModel, chain::Chains;
                             N=250, seed=0, kwargs...)
    series = PosteriorSeries(model, chain; N, seed)
    return Octofitter.octoplot(series; kwargs...)
end

function Octofitter.octoplot(series::PosteriorSeries;
                             show_sky=nothing, show_phase=nothing, whiten=nothing,
                             fname=nothing, figscale=1.0, ndraws=nothing)
    sys = series.model.system
    groups = _channelgroups(sys)
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
    bespoke = [obs for obs in sys.observations if !isempty(defaultpanels(obs))]

    npanels = length(groups)
    W = round(Int, figscale * 620)
    # The sky cell aims square at the column-1 width (figure width minus
    # side column, gaps, and typical protrusions); DataAspect + edge
    # alignment in the panel absorb the estimate's error.
    skyH = W - 175
    H = round(Int, figscale * 40 + do_sky * figscale * skyH +
              figscale * (npanels * 260 + length(phasepanels) * 230 + 200 * length(bespoke)))
    fig = Makie.Figure(size=(W, max(H, 300)))

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

end # module
