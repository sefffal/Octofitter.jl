# ---------------------------------------------------
# OctofitterMakieExt — the v2 plotting extension.
#
# Presentation only: every model/data/residual quantity comes from the core
# plotting API (`PosteriorSeries`, `plotchannels`, `Octofitter.residuals`,
# `modelcurves`), which itself defers to the likelihoods. Nothing in this
# file recomputes likelihood math.
#
# Three generic panels — sky, time-series (with residual strip + marginal
# histogram), and the octoplot assembly. Observation types that implement
# the plotting protocol get their panels with no code here; the long tail of
# bespoke panels (proper-motion vectors, scan geometry, …) plugs into the
# same layout slots later.
# ---------------------------------------------------
module OctofitterMakieExt

using Octofitter
using Octofitter: PosteriorSeries, ObservableQuery, PlotChannel, OctoPlotResult,
    plotchannels, obscontext, modelcurves, default_sky_queries,
    plottable_observations, likelihoodname, refspecs,
    _query, _querystr, _refstr, _solvekw, _isprior
using PlanetOrbits
using PlanetOrbits: plotinfo, plotlabel, _names, _setnames
using Makie
using Statistics
using StatsBase: fit, Histogram
using MCMCChains: Chains

const WONG = Makie.wong_colors()
_rowcolor(k) = WONG[mod1(k, length(WONG))]
# The draw-ramp: row accent colour fading toward near-white with orbital
# phase; the light endpoint is the one v1 used in octoplot.
_phaseramp(c) = Makie.cgrad([Makie.to_color(c), Makie.to_color("#FAFAFA")])

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
    # arcsec) needs the full extent before anything draws.
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

    # mas → arcsec switch, v1 threshold.
    extent = maximum((isempty(xs) ? 0.0 : maximum(abs, filter(!isnan, xs)) for (xs, _, _) in rowdata); init=0.0)
    unitscale = (angular && extent > 1500) ? 1e-3 : 1.0
    unit = angular ? (unitscale == 1.0 ? "mas" : "arcsec") : "au"

    ax = Makie.Axis(gs[1, 1];
        xlabel="$(info_x.label) [$unit]",
        ylabel="$(info_y.label) [$unit]",
        aspect=Makie.DataAspect(),
        xreversed=info_x.flip,
        xgridvisible=false, ygridvisible=false)

    α = _alpha(nshow)
    for (k, (xs, ys, cs)) in enumerate(rowdata)
        Makie.lines!(ax, xs .* unitscale, ys .* unitscale;
            color=cs, colormap=_phaseramp(_rowcolor(k)), colorrange=(0, 2π),
            alpha=α, linewidth=0.6)
    end

    # Data overlay: every relative-astrometry observation, converted to the
    # panel's coordinates. (Uncertainty crosses; ellipses come later.)
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
        Makie.Colorbar(gs[1, 2];
            colormap=_phaseramp(_rowcolor(1)), colorrange=(0, 2π),
            ticks=_PHASE_TICKS, label="orbital phase (mean anomaly) →")
    end
    return (; sky=ax)
end

# ---------------------------------------------------
# Time-series panel
# ---------------------------------------------------

# Apply a display wrap: values into [0, wrap), NaN-ing the point after each
# wraparound so no vertical jump line is drawn.
function _wrapped(v, wrap)
    wrap === nothing && return v
    w = mod.(v, wrap)
    for i in 2:length(w)
        if !isnan(w[i]) && !isnan(w[i - 1]) && abs(w[i] - w[i - 1]) > wrap / 2
            w[i] = NaN
        end
    end
    return w
end

"""
    timeseriespanel!(gp, series, entries; kwargs...)

Generic data-vs-model panel for one channel group. `entries` is a vector of
`(obs, channel)` pairs sharing one quantity (e.g. several RV instruments);
pass `(nothing, channel)` for a pure model panel with no data.

Layout: main axis (posterior model curves + calibrated data), residual
strip glued below, marginal residual histogram at right with bins shared
across instruments. The epoch axis is calendar-dated via `MJDConversion`
and stays correct under zoom; annotate with `Date`/`DateTime` or MJD alike.

Returns `(; main, resid, hist)` (the latter two `nothing` when no data).
"""
function Octofitter.timeseriespanel!(gp, series::PosteriorSeries, entries;
                                     top_time_axis=true, bottom_time_axis=true,
                                     show_hist=true, ndraws=nothing,
                                     curvecolor=WONG[1], linewidth=0.4)
    gs = _layout(gp)
    entries = [(obs, ch) for (obs, ch) in entries]
    isempty(entries) && error("timeseriespanel! needs at least one (obs, channel) entry")
    ch1 = entries[1][2]
    ylabel = isempty(ch1.unit) ? String(ch1.label) : "$(ch1.label) [$(ch1.unit)]"
    hasdata = any(obs !== nothing for (obs, _) in entries)

    conv() = PlanetOrbits.MJDConversion()
    ax = Makie.Axis(gs[1, 1];
        ylabel, dim1_conversion=conv(),
        xaxisposition=:top,
        xticklabelsvisible=top_time_axis, xticksvisible=top_time_axis,
        xgridvisible=false, ygridvisible=false)

    ax_resid = ax_hist = nothing
    if hasdata
        ax_resid = Makie.Axis(gs[2, 1];
            ylabel="residuals", dim1_conversion=conv(),
            xlabel=bottom_time_axis ? "epoch" : "",
            xticklabelsvisible=bottom_time_axis, xticksvisible=bottom_time_axis,
            xgridvisible=false, ygridvisible=false)
        Makie.linkxaxes!(ax, ax_resid)
        Makie.rowgap!(gs, 1, 0)
        Makie.rowsize!(gs, 1, Makie.Auto(2.4))
        Makie.hlines!(ax_resid, 0.0; color=:black, linewidth=1)
    end

    # Posterior model curves (one NaN-joined lines! for all draws).
    q = ch1.query
    if q !== nothing
        nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
        curves = modelcurves(series, q)
        tsn = vcat(series.ts, NaN)
        xs = repeat(tsn, nshow)
        ys = Float64[]
        sizehint!(ys, length(tsn) * nshow)
        for d in 1:nshow
            v = _wrapped(curves[d] .* ch1.scale, ch1.wrap)
            append!(ys, v)
            push!(ys, NaN)
        end
        Makie.lines!(ax, xs, ys;
            color=(curvecolor, _alpha(nshow)), linewidth)
    end

    # Data + residuals, per instrument, MAP-referenced.
    allresids = Float64[]
    perinstr = Tuple{Any,Any,Any}[]  # (label, colour, channel residual NT)
    for (j, (obs, ch)) in enumerate(entries)
        obs === nothing && continue
        r = Octofitter.residuals(obs, obscontext(series, obs))[ch.name]
        color = length(entries) == 1 ? :black : WONG[mod1(j + 1, length(WONG))]
        u = r.use
        # jitter-inflated errorbars behind measurement errorbars
        Makie.errorbars!(ax, r.epoch[u], r.data[u], r.σ_eff[u]; color="#CCCCCC", linewidth=1)
        Makie.errorbars!(ax, r.epoch[u], r.data[u], r.σ[u]; color=:black, linewidth=2)
        Makie.scatter!(ax, r.epoch[u], r.data[u];
            color=length(entries) == 1 ? :white : color,
            strokecolor=:black, strokewidth=1.5, markersize=7,
            label=likelihoodname(obs))
        if ax_resid !== nothing
            Makie.errorbars!(ax_resid, r.epoch[u], r.resid[u], r.σ_eff[u]; color="#CCCCCC", linewidth=1)
            Makie.errorbars!(ax_resid, r.epoch[u], r.resid[u], r.σ[u]; color=:black, linewidth=2)
            Makie.scatter!(ax_resid, r.epoch[u], r.resid[u];
                color=length(entries) == 1 ? :white : color,
                strokecolor=:black, strokewidth=1.5, markersize=6)
            append!(allresids, r.resid[u])
        end
        push!(perinstr, (likelihoodname(obs), color, r))
    end

    # Marginal histogram: bins shared across instruments.
    if ax_resid !== nothing && show_hist && !isempty(allresids)
        h_all = fit(Histogram, allresids)
        edges = h_all.edges[1]
        centers = collect(edges[1:end-1] .+ step(edges) / 2)
        ax_hist = Makie.Axis(gs[2, 2]; xgridvisible=false, ygridvisible=false)
        Makie.hidedecorations!(ax_hist)
        Makie.linkyaxes!(ax_resid, ax_hist)
        for (_, color, r) in perinstr
            h = fit(Histogram, r.resid[r.use], edges)
            Makie.stairs!(ax_hist, h.weights, centers;
                step=:center, color=color isa Symbol ? :black : color, linewidth=1.0)
        end
        Makie.hlines!(ax_hist, 0.0; color=:black, linewidth=1)
        Makie.xlims!(ax_hist; low=0)
        Makie.colsize!(gs, 2, Makie.Relative(0.15))
        Makie.colgap!(gs, 1, 6)
    end

    return (; main=ax, resid=ax_resid, hist=ax_hist)
end

# ---------------------------------------------------
# octoplot assembly
# ---------------------------------------------------

# Group plottable observation channels into panels: channels whose smooth
# curve is the same query share a panel (several RV instruments); channels
# without a query group by observation type + channel name.
function _channelgroups(sys)
    groups = Vector{Pair{Symbol,Vector{Tuple{Any,PlotChannel}}}}()
    for obs in plottable_observations(sys)
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
function _panelnames(groups)
    names = Symbol[]
    for (key, entries) in groups
        base = entries[1][2].name
        nm = base
        k = 2
        while nm in names
            nm = Symbol(base, :_, k)
            k += 1
        end
        push!(names, nm)
    end
    return names
end

"""
    octoplot(model, chain; N=250, seed=0, show_sky=nothing, fname=nothing,
             figscale=1, ndraws=nothing)

See the core docstring (`?octoplot` without Makie loaded). Panels: a sky
panel when the system has angular observables, then one time-series panel
per channel group, all time axes linked.
"""
function Octofitter.octoplot(model::Octofitter.LogDensityModel, chain::Chains;
                             N=250, seed=0, kwargs...)
    series = PosteriorSeries(model, chain; N, seed)
    return Octofitter.octoplot(series; kwargs...)
end

function Octofitter.octoplot(series::PosteriorSeries;
                             show_sky=nothing, fname=nothing, figscale=1.0,
                             ndraws=nothing)
    sys = series.model.system
    groups = _channelgroups(sys)
    pnames = _panelnames(groups)
    do_sky = show_sky === nothing ?
             (_isangular(sys) && !isempty(sys.rows)) : show_sky

    npanels = length(groups)
    H = round(Int, figscale * (do_sky * 500 + npanels * 260 + 40))
    W = round(Int, figscale * 620)
    fig = Makie.Figure(size=(W, max(H, 300)))

    axpairs = Pair{Symbol,Any}[]
    timeaxes = Makie.Axis[]
    row = 1
    if do_sky
        skyaxes = Octofitter.skypanel!(fig[row, 1], series; ndraws)
        push!(axpairs, :sky => skyaxes)
        Makie.rowsize!(fig.layout, row, Makie.Auto(500 / 260))
        row += 1
    end
    for (i, ((key, entries), nm)) in enumerate(zip(groups, pnames))
        axs = Octofitter.timeseriespanel!(fig[row, 1], series, entries;
            top_time_axis=(i == 1), bottom_time_axis=(i == npanels), ndraws)
        push!(axpairs, nm => axs)
        push!(timeaxes, axs.main)
        axs.resid !== nothing && push!(timeaxes, axs.resid)
        row += 1
    end
    length(timeaxes) > 1 && Makie.linkxaxes!(timeaxes...)

    axes = (; axpairs...)
    fname !== nothing && Makie.save(fname, fig)
    return OctoPlotResult(fig, axes, series)
end

end # module
