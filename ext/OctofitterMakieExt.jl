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
    plotchannels, obscontext, modelcurves, default_sky_queries, default_queries,
    plottable_observations, likelihoodname, refspecs, defaultpanels,
    rowsignal, evalsignal, foldablerows, foldephemeris, foldphase, phasebinmeans,
    sharepanel, datacalibration, noisemodel, predictedchannels,
    _query, _querystr, _refstr, _solvekw, _isprior, _calibration,
    _requestedobservables
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
# Instrument accent colours start at Wong 2: Wong 1 is the model curve, and a
# dataset drawn in the model's own colour reads as part of the model.
_instcolor(j) = WONG[mod1(j + 1, length(WONG))]

# `curvecolor=` on an assembled figure (`octoplot`), resolved per panel. Three
# spellings, in increasing specificity:
#
#   nothing            — every panel keeps its own accent colour (the default)
#   a colour           — that colour everywhere
#   a Dict/NamedTuple  — keyed by *body name*: `(; b=:firebrick)`, or
#                        `Dict(:b => :firebrick, :A => :grey20)`
#
# The key is a body because that is what every panel's colour already means:
# the sky track and the astrometry panels of planet `b` share row `b`'s accent
# colour, and a stellar reflex RV panel is the host's — so `:b` recolours the
# planet's orbit wherever it appears, and a body the mapping does not name
# keeps its default. Anything that is not a mapping is treated as a colour, so
# every Makie colour spelling — `(:red, 0.5)`, an `RGBf`, `"#AA3377"` — passes
# through untouched.
_panelcolor(::Nothing, name) = nothing
_panelcolor(d::AbstractDict, name) =
    name === nothing ? nothing : get(d, name, nothing)
_panelcolor(nt::NamedTuple, name) =
    name === nothing ? nothing : get(nt, name, nothing)
_panelcolor(c, name) = c

# The body a panel's curve belongs to: the query's target when that is a
# single body (planet astrometry → `:b`, a stellar reflex RV → `:A`), and
# `nothing` for a composite reference point, which no body name can select.
function _curvebody(sys, q)
    q === nothing && return nothing
    n = Octofitter._leafnames(sys, q.target)
    (n === nothing || length(n) != 1) && return nothing
    return n[1]
end

# Per-instrument marker shapes, so that overlapping datasets stay separable in
# print and for readers who cannot rely on the colour alone (v1's rvpostplot
# idiom, which the first v2 port dropped).
const MARKERS = (:circle, :rect, :diamond, :utriangle, :dtriangle, :pentagon,
    :cross, :xcross, :hexagon, :star5, :ltriangle, :rtriangle, :star4, :star8)
_instmarker(j) = MARKERS[mod1(j, length(MARKERS))]

# How much ink a data point costs, shared by the time-series and phase-folded
# panels. The defaults are octoplot's: one instrument per panel, a handful of
# epochs, so a heavy black measurement bar reads as precision.
#
# `rvplot` overrides them with v8 `rvpostplot`'s conventions, because it packs
# every instrument and every epoch onto one axis — at a few hundred points the
# same bars become a wall of black, and the instrument colours (which are the
# only thing saying *whose* point it is) are what gets buried. `σ_color =
# :instrument` draws the measurement bar in the point's own colour instead,
# which is also v8's answer to "which of these two nested bars is which".
const DATASTYLE = (; markersize=7, resid_markersize=6, strokewidth=1.5,
    σ_linewidth=2.0, σ_color=:black,
    σeff_linewidth=1.0, σeff_color="#CCCCCC")
function _datastyle(s, extra=(;))
    base = merge(DATASTYLE, extra)
    s === nothing && return base
    # A misspelled key would otherwise merge in as a field nothing ever reads,
    # so the override would silently do nothing.
    bad = setdiff(keys(s), keys(DATASTYLE))
    isempty(bad) || error("`datastyle`: unknown $(length(bad) == 1 ? "key" : "keys") " *
                          "$(join(bad, ", ")). Known: $(join(keys(DATASTYLE), ", ")).")
    return merge(base, s)
end
# `:instrument` resolves per entry; anything else is a literal colour.
_σcolor(st, color) = st.σ_color === :instrument ? color : st.σ_color

# The pale jitter-inflated bar goes *behind* the model curve and the residual
# strip's zero line. Drawn in order it would land on top of both, and a light
# grey vertical stroke across a dark line does not read as an error bar — it
# reads as a gap in the line, which is what made v8's curves look dashed
# wherever the data were dense. Still in front of the correlated-noise band
# (that sits at -10).
_behind!(p) = (Makie.translate!(p, 0, 0, -5); p)

# The draw-ramp: row accent colour fading toward near-white with orbital
# phase; the light endpoint is the one v1 used in octoplot.
_phaseramp(c) = Makie.cgrad([Makie.to_color(c), Makie.to_color("#FAFAFA")])

# Side column (colorbars, marginal histograms): one Fixed width everywhere,
# so column-1 boundaries — and with them every main axis edge — line up
# across panels.
const SIDE_W = 80.0
const COLGAP = 8.0
# `rvplot`'s legend column, down the right of every panel (v8's layout).
const LEGEND_W = 150.0

_layout(gl::Makie.GridLayout) = gl
_layout(gp::Union{Makie.GridPosition,Makie.GridSubposition}) = Makie.GridLayout(gp)

# Where a single-axis panel draws. A grid cell (`fig[i, j]`) is the usual
# case; a bare figure or layout means "the only axis in it"; and an `Axis`
# means the caller has already made and styled one, so the panel's own axis
# attributes are theirs to have set — it draws into it and does not restyle.
_panelaxis(ax::Makie.Axis; kwargs...) = ax
_panelaxis(gp::Union{Makie.GridPosition,Makie.GridSubposition}; kwargs...) = Makie.Axis(gp; kwargs...)
_panelaxis(gl::Makie.GridLayout; kwargs...) = Makie.Axis(gl[1, 1]; kwargs...)
_panelaxis(fig::Makie.Figure; kwargs...) = Makie.Axis(fig[1, 1]; kwargs...)

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
    skypanel!(gp, series; tracks=nothing, ntrack=150, colorbar=true,
              legend=true, ndraws=nothing, curvecolor=nothing)

Sky-plane orbit tracks. By default one track per hierarchy row — the
exterior side relative to the interior side, i.e. exactly the relationship
the row parametrizes — sampled uniformly in eccentric anomaly per draw and
coloured by orbital phase in the row's accent colour. Relative-astrometry
data overlay automatically, one colour and marker shape per observation, with
a legend naming them (`legend=false` to drop it).

`curvecolor=` replaces the colour each phase ramp is built *from* — the ramp
still runs from it to near-white with orbital phase, so the panel keeps saying
where in its orbit each piece of track is. Give one colour for every row, or a
`Dict`/`NamedTuple` keyed by body name (`(; b=:firebrick)`) to recolour some
of them; the colorbar follows the first row's.

Right ascension increases to the **left**, as on the sky. Returns `(; sky=ax)`.
"""
function Octofitter.skypanel!(gp, series::PosteriorSeries;
                              tracks=nothing, ntrack=150, colorbar=true,
                              legend=true, ndraws=nothing, curvecolor=nothing)
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
    # Both parametrizations reach the sky panel through the same two channels:
    # `residuals` publishes (Δα*, Δδ) whatever the table holds, with the
    # platescale and north-angle calibration already applied — a sep/pa
    # instrument and a ra/dec one therefore land in one another's frame with no
    # plot-side trigonometry, and no second implementation of it to drift.
    datasets = NTuple{4,Vector{Float64}}[]
    datanames = String[]
    for obs in plottable_observations(sys)
        # `plotobs`, not `obs`: an observation wrapped in `ObsPriorONeil2019`
        # is still relative astrometry and its points still belong on the sky.
        # The wrapper itself is what `obscontext` must be given, though — the
        # calibration parameters are registered under *its* name.
        Octofitter.plotobs(obs) isa Octofitter.RelAstromObs || continue
        r = Octofitter.residuals(obs, obscontext(series, obs))
        u = r.raoff.use
        push!(datasets, (r.raoff.data[u], r.decoff.data[u],
            r.raoff.σ[u], r.decoff.σ[u]))
        push!(datanames, likelihoodname(obs))
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
        xgridvisible=false, ygridvisible=false)
    # `xlims!` *sets* `xreversed` from the order of its arguments — passing
    # low-then-high silently clears an `xreversed=true` given to the
    # constructor, which is how the sky panel lost its east-left convention.
    # Handing it the limits already reversed is the spelling that survives.
    lo, hi = (cx - half) * unitscale, (cx + half) * unitscale
    Makie.xlims!(ax, info_x.flip ? (hi, lo) : (lo, hi))
    Makie.ylims!(ax, (cy - half) * unitscale, (cy + half) * unitscale)

    # The colour each row's phase ramp starts from: its accent colour, unless
    # the caller named one for that body.
    rowbase(k) = something(_panelcolor(curvecolor, get(rownames, k, nothing)), _rowcolor(k))

    α = _alpha(nshow)
    for (k, (xs, ys, cs)) in enumerate(rowdata)
        Makie.lines!(ax, xs .* unitscale, ys .* unitscale;
            color=cs, colormap=_phaseramp(rowbase(k)), colorrange=(0, 2π),
            alpha=α, linewidth=0.6)
    end

    # Data overlay: uncertainty crosses; ellipses come later. One instrument
    # is drawn in white — there is nothing to tell it apart from — and several
    # take an accent colour and a marker shape each.
    single = length(datasets) == 1
    for (j, (x, y, σx, σy)) in enumerate(datasets)
        Makie.errorbars!(ax, x .* unitscale, y .* unitscale, σy .* unitscale;
            direction=:y, color=:black, linewidth=1)
        Makie.errorbars!(ax, x .* unitscale, y .* unitscale, σx .* unitscale;
            direction=:x, color=:black, linewidth=1)
        Makie.scatter!(ax, x .* unitscale, y .* unitscale;
            color=single ? :white : _instcolor(j), marker=_instmarker(j),
            strokecolor=:black, strokewidth=2, markersize=8,
            label=datanames[j])
    end

    # The reference point of every track sits at the origin.
    Makie.scatter!(ax, [0.0], [0.0];
        marker='★', markersize=20, color=:white, strokecolor=:black, strokewidth=1.5)

    legend && !isempty(datasets) &&
        Makie.axislegend(ax, "Instrument"; position=:rb, framevisible=false,
            labelsize=9, titlesize=10, padding=(4, 4, 2, 2), patchlabelgap=3)

    if colorbar
        cb = Makie.Colorbar(gs[1, 2];
            colormap=_phaseramp(rowbase(1)), colorrange=(0, 2π),
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

# The residual a draw is actually left with, and the uncertainty it should be
# judged against. Without a correlated-noise model that is `data − model` over
# `σ_eff`; with one, the GP's predicted mean comes off the residual and its
# predictive variance goes into the uncertainty — otherwise the strip shows
# precisely the correlated structure the GP was fitted to absorb, and no fit,
# however good, produces standard-normal z-scores.
function _netresid(r)
    haskey(r, :gp_mean) || return (r.resid, r.σ_eff)
    return (r.resid .- r.gp_mean, sqrt.(r.σ_eff .^ 2 .+ r.gp_var))
end

# Whitened residuals per posterior draw: (resid / σ_eff) with each draw's
# own parameters (jitters vary by draw), from the likelihood's residuals.
function _zscores(series, obs, ch, nshow)
    zs = nothing
    for d in 1:nshow
        r = Octofitter.residuals(obs, obscontext(series, obs; draw=d))[ch.name]
        resid, σ = _netresid(r)
        zs === nothing && (zs = Matrix{Float64}(undef, length(resid), nshow))
        zs[:, d] = resid ./ σ
    end
    return zs
end

# Whitened or raw? Raw residuals cannot be shown for several draws at once:
# the jitter, the offsets and the other rows' subtracted signals all move
# between draws, so the same measurement has a different residual — and a
# different meaningful *size* — in every one of them. Dividing by that draw's
# own σ_eff is what makes them commensurable, which is why the strip
# summarizes the z-scores per point rather than plotting residuals per draw.
function _whitenmode(whiten, nshow::Integer, hasdata::Bool)
    hasdata || return false
    whiten === nothing && return nshow > 1
    whiten && return true
    nshow > 1 && error("""
    `whiten=false` asks for residuals in data units, but $nshow draws are being
    shown, and a raw residual is not comparable between draws: each carries its
    own jitter, its own instrument offsets, and its own subtracted signals from
    the other bodies. The same point is a different number of σ from the model
    in each — so a strip of raw residuals over many draws is not a plot of
    anything.

    Either show one draw, which pins every nuisance parameter and makes a raw
    residual well defined:

        octoplot(model, chain; ndraws=1, whiten=false)
        rvplot(model, chain)                 # the single-draw RV figure

    or keep the whitened strip, which is the many-draw form of the same thing:
    per point, the median and 16-84 % interval of (data - model)/σ_eff over the
    draws, against a unit normal.""")
    return false
end

# How wide a per-point box may be, in x units, on an axis spanning `xspan`.
#
# Two constraints pull opposite ways — boxes should not swallow the epoch
# structure, and they have to be visible — and at real sampling densities the
# second wins outright. A spectroscopic campaign takes exposures half an hour
# apart within a night on a baseline of months: K2-131's HARPS-N series has a
# median gap of 0.025 d on a 66 d axis, which is 4·10⁻⁴ of the width, a fifth
# of a pixel. Nothing can be drawn there. So the sampling sets the width only
# where the data are sparse enough for that to be possible, and the clamp —
# roughly two pixels to a twenty-fifth of the axis — decides the rest.
#
# The bounds are relative to the **axis**, not to the entry's own extent: RV
# panels share one linked time axis, so an instrument that observed for three
# weeks of a two-month campaign would otherwise size its boxes against a span
# a third of the one they are drawn on, and come out invisible while its
# neighbour's read fine.
#
# Where a night's exposures do overlap they overlap: at that x resolution the
# measurements are already on top of each other (they are on the panel above
# too), and a column of overlapping translucent boxes reads as the cluster it
# is. `boxwidth=` overrides the lot. The gap statistic is a low quantile
# rather than the minimum so that one pair of back-to-back exposures does not
# shrink every box on a sparse panel.
function _autoboxwidth(xs, xspan)
    xspan > 0 || return nothing
    lo, hi = xspan / 200, xspan / 25
    v = sort!(unique(x for x in xs if isfinite(x)))
    length(v) > 1 || return lo
    return clamp(0.8 * quantile(diff(v), 0.1), lo, hi)
end

# Draw the residual strip and marginal histogram for `entries` into
# pre-created axes. `xmap(epochs) -> x` positions points (identity for time
# panels, the fold phase for phase panels). In whitened mode each point
# carries the distribution of its z-score over the draws — a boxplot when
# there are several draws to distribute, otherwise the single value — against
# ±1σ guides, and the histogram is density-normalized with a unit normal
# overlaid; in raw mode it is the MAP residual with measurement/jitter
# errorbars and a count histogram. Histogram bins are shared across
# instruments either way.
function _residstrip!(ax_resid, ax_hist, series, entries, xmap, nshow, whiten,
                      st=DATASTYLE; xspan, boxwidth=nothing)
    Makie.hlines!(ax_resid, 0.0; color=:black, linewidth=1)
    whiten && Makie.hlines!(ax_resid, [-1.0, 1.0]; color=(:black, 0.3), linewidth=0.5)

    # `residuals` conditions the noise model, so it is called once per entry
    # here and reused for the boxes, the marks and the histogram.
    data = [(j, obs, ch, Octofitter.residuals(obs, obscontext(series, obs))[ch.name])
            for (j, (obs, ch)) in enumerate(entries) if obs !== nothing]

    # A box needs a distribution behind it; one draw gives one z-score per
    # point, which is the marker form. One width serves every instrument on
    # the panel — differing widths would read as meaning something.
    boxes = whiten && nshow > 1
    bw = nothing
    if boxes
        bw = boxwidth === nothing ?
             _autoboxwidth(Iterators.flatten(xmap(r.epoch)[r.use] for (_, _, _, r) in data),
                 xspan) : boxwidth
        # A degenerate axis leaves no width to draw at: fall back to the marker
        # form rather than draw a box of undefined width.
        boxes = bw !== nothing && bw > 0
    end

    histdata = Tuple{Any,Vector{Float64}}[]   # (colour, values)
    single = length(data) == 1
    for (j, obs, ch, r) in data
        u = r.use
        x = xmap(r.epoch)
        color = single ? :black : _instcolor(j)
        mcolor = single ? :white : color
        marker = single ? :circle : _instmarker(j)
        if whiten
            zs = _zscores(series, obs, ch, nshow)
            if boxes
                # The point's z-score *distribution* over the draws, not a
                # summary of it. Every nuisance parameter that sets a residual's
                # size — the jitter, the instrument offset, the trend, the other
                # rows' subtracted signals — moves from sample to sample, so the
                # measurement has a different z-score in each draw and a single
                # mark per point cannot say how much of the scatter is the fit's
                # own uncertainty. (v8 drew exactly this for `G23HObs`' catalog
                # channels, on a categorical axis; it is the right mark for any
                # whitened strip, and an epoch axis only adds the width problem
                # `_autoboxwidth` solves.)
                bx = Float64[]; by = Float64[]
                n = count(u)
                sizehint!(bx, n * nshow); sizehint!(by, n * nshow)
                for i in axes(zs, 1)
                    u[i] || continue
                    for d in axes(zs, 2)
                        push!(bx, x[i]); push!(by, zs[i, d])
                    end
                end
                # A box only two or three pixels wide is mostly stroke, so the
                # marks are weighted for that end: a near-solid fill that
                # carries the instrument's colour, a hairline outline, thinner
                # whiskers so the box still reads as the wider mark, and the
                # median in black against both. Outliers are off — 250 draws
                # per point would scatter more dots than there are data.
                boxfill = single ? Makie.to_color("#B4B4B4") : Makie.to_color((color, 0.7))
                boxedge = single ? Makie.to_color(:grey25) : Makie.to_color(color)
                Makie.boxplot!(ax_resid, bx, by;
                    width=bw, gap=0, color=boxfill,
                    strokecolor=boxedge, strokewidth=0.5,
                    whiskercolor=boxedge, whiskerlinewidth=0.8,
                    mediancolor=:black, medianlinewidth=1.0,
                    show_outliers=false)
            else
                med = [median(view(zs, i, :)) for i in axes(zs, 1)]
                lo = [quantile(view(zs, i, :), 0.16) for i in axes(zs, 1)]
                hi = [quantile(view(zs, i, :), 0.84) for i in axes(zs, 1)]
                Makie.rangebars!(ax_resid, x[u], lo[u], hi[u];
                    color=single ? "#AAAAAA" : color, linewidth=st.σeff_linewidth)
                Makie.scatter!(ax_resid, x[u], med[u];
                    color=mcolor, marker, strokecolor=:black,
                    strokewidth=st.strokewidth, markersize=st.resid_markersize)
            end
            append!(histdata, ((color isa Symbol ? :black : color, vec(zs[u, :])),))
        else
            resid, σ_eff = _netresid(r)
            _behind!(Makie.errorbars!(ax_resid, x[u], resid[u], σ_eff[u];
                color=st.σeff_color, linewidth=st.σeff_linewidth))
            Makie.errorbars!(ax_resid, x[u], resid[u], r.σ[u];
                color=_σcolor(st, color), linewidth=st.σ_linewidth)
            Makie.scatter!(ax_resid, x[u], resid[u];
                color=mcolor, marker, strokecolor=:black,
                strokewidth=st.strokewidth, markersize=st.resid_markersize)
            append!(histdata, ((color isa Symbol ? :black : color, resid[u]),))
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

# A residual strip is glued to its main axis with no gap, so the two share a
# spine — and whenever the data run close to it, the main axis' bottom-most
# tick label and the strip's top-most one print on the same line. Blank the
# main axis' lowest label: at that boundary the strip's scale is the one being
# read, and the alternative (a visible gap) breaks the glued look the whole
# figure set uses. Formatting for the rest is Makie's own, so the labels stay
# identical to every other axis.
function _blanklowesttick(vals)
    labels = Makie.get_ticklabels(Makie.automatic, vals)
    return Any[i == 1 ? "" : l for (i, l) in enumerate(labels)]
end

# ---------------------------------------------------
# Time-series panel
# ---------------------------------------------------

"""
    timeseriespanel!(gp, series, entries; kwargs...)

Generic data-vs-model panel for one channel group. `entries` is a vector of
`(obs, channel)` pairs sharing one quantity; pass `(nothing, channel)` for a
pure model panel with no data at all — which is how a fit predicts a quantity
it was never given (see [`predictedchannels`](@ref)).

Layout: main axis (posterior model curves + data), residual strip glued
below, marginal residual histogram at right with bins shared across
instruments. Calendar dates label the top axis via `MJDConversion` and stay
correct under zoom (annotate with `Date`/`DateTime` or MJD alike); the bottom
axis carries the raw MJD numbers you would type back into a script. The time
axes clip exactly to the model grid (no gap after the last segment); wrapped
quantities (position angle) get the v8 treatment — lines exit one axis edge
and re-enter at the other, with limits pinned to the wrap range.

## Which side gets calibrated

`calibrate=:model` leaves the measurements exactly as reported and gives
**each draw's model curve that draw's own instrument offset and trend**
([`datacalibration`](@ref)). This is the only arrangement in which many draws
and one dataset are simultaneously honest, so it is what a per-instrument
panel uses.

`calibrate=:data` instead moves the data onto a single pure-observable curve,
using the maximum-posterior draw's calibration. That is defensible when one
draw is being shown (there is then nothing to be inconsistent with), or when
the calibration is so tightly constrained that every draw agrees — which is
exactly the [`sharepanel`](@ref) question. [`rvplot`](@ref) uses it.

`calibrate=:auto` (the default) reads that decision off the observations:
`:data` when they all share panels, `:model` otherwise. The two coincide for
observations that declare no `datacalibration` at all, except that `:data`
draws one curve family rather than an identical one per instrument.

## Residuals

`whiten=true` shows `(data − model)/σ_eff` per posterior draw. With several
draws each point becomes a **boxplot of its own z-score distribution** —
median, quartiles and 1.5 IQR whiskers over the draws — because every term
that sets a residual's size (jitter, offset, trend, the other rows' signals)
moves from sample to sample, so one mark per measurement cannot say how much
of the scatter is the fit's own uncertainty. `boxwidth=` overrides the
automatic width (`_autoboxwidth`), which is chosen from the spacing of the
epochs; with one draw the points are drawn as markers instead. The histogram
pools every draw against a unit normal.

The default whitens whenever more than one draw is shown, and `whiten=false`
with several draws is an error rather than a plot — see the `Octofitter` issue
text in `_whitenmode`.

Where an observation has a correlated-noise model ([`noisemodel`](@ref)), its
prediction is subtracted from the residuals and added into `σ_eff`, and — for
a **single** draw — drawn as a band around the model curve. `gpband=nothing`
follows the draw count; `gpband=true` forces the band on for many draws too,
which puts each draw's plain orbit curve and its GP-added twin on one axis.

## Marks

`datastyle` is a `NamedTuple` merged over the panel's defaults, carrying
`markersize`, `resid_markersize`, `strokewidth`, `σ_linewidth`, `σ_color`,
`σeff_linewidth` and `σeff_color`. `σ_color = :instrument` draws the inner
measurement bar in each point's own colour rather than black, which is what
lets several instruments share one axis legibly; [`rvplot`](@ref) uses it.

Returns `(; main, resid, hist)` (the latter two `nothing` when no data).
"""
function Octofitter.timeseriespanel!(gp, series::PosteriorSeries, entries;
                                     top_time_axis=true, bottom_time_axis=true,
                                     show_hist=true, ndraws=nothing, whiten=nothing,
                                     calibrate::Symbol=:auto, show_legend=nothing,
                                     gpband=nothing, boxwidth=nothing, datastyle=nothing,
                                     side_w=SIDE_W, colgap=COLGAP, hist_aspect=nothing,
                                     curvecolor=nothing, linewidth=nothing)
    gs = _layout(gp)
    st = _datastyle(datastyle)
    entries = [(obs, ch) for (obs, ch) in entries]
    isempty(entries) && error("timeseriespanel! needs at least one (obs, channel) entry")
    calibrate in (:auto, :model, :data) ||
        error("`calibrate` must be :auto, :model or :data, got :$calibrate")
    if calibrate === :auto
        # Sharing a panel and calibrating the data are the same decision: an
        # observation only shares when its calibration is pinned tightly enough
        # that one draw's version of it stands for all of them.
        calibrate = all(e -> e[1] === nothing || sharepanel(e[1]), entries) ? :data : :model
    end
    ch1 = entries[1][2]
    ylabel = isempty(ch1.unit) ? String(ch1.label) : "$(ch1.label) [$(ch1.unit)]"
    ndata = count(e -> e[1] !== nothing, entries)
    hasdata = ndata > 0
    nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
    whiten = _whitenmode(whiten, nshow, hasdata)
    curvecolor = something(curvecolor, _querycolor(series.model.system, ch1.query))
    linewidth = something(linewidth, nshow == 1 ? 1.3 : 0.4)
    show_legend = something(show_legend, true)
    # The correlated-noise band, by default, only where it can be read. One
    # draw gets a band and the curve it belongs to, which is the picture of
    # what the GP absorbed. Many draws put the same panel's plain orbit curves
    # *and* their GP-added twins on one axis — twice the ink for a band whose
    # envelope no longer means anything, because it is a different draw's
    # activity model everywhere you look.
    gpband = something(gpband, nshow == 1)
    single = ndata == 1

    conv() = PlanetOrbits.MJDConversion()
    ax = Makie.Axis(gs[1, 1];
        ylabel, dim1_conversion=conv(),
        xaxisposition=:top,
        xticklabelsvisible=top_time_axis, xticksvisible=top_time_axis,
        xgridvisible=false, ygridvisible=false)

    ax_resid = ax_hist = nothing
    if hasdata
        # No `dim1_conversion` here: the top axis already says *when*, and the
        # bottom of the figure is where the MJD numbers belong — dates on both
        # edges is the same information printed twice.
        ax_resid = Makie.Axis(gs[2, 1];
            ylabel=_stripylabel(whiten),
            xlabel=bottom_time_axis ? "epoch [MJD]" : "",
            xticklabelsvisible=bottom_time_axis, xticksvisible=bottom_time_axis,
            xgridvisible=false, ygridvisible=false)
        Makie.linkxaxes!(ax, ax_resid)
        ax.ytickformat = _blanklowesttick
        Makie.rowgap!(gs, 1, 0)
        Makie.rowsize!(gs, 1, Makie.Auto(2.4))
    elseif bottom_time_axis
        # A data-free (predictive) panel has no residual strip to carry the
        # MJD axis, so it gets a decoration-only companion instead.
        PlanetOrbits.add_mjd_axis!(gs[1, 1], ax; label="epoch [MJD]", ticklabelscale=1.0)
    end

    # Posterior model curves. In `:model` mode the curve belongs to an
    # instrument — it carries that instrument's zero point and trend — so
    # there is one family per entry; in `:data` mode one family serves them
    # all, because the data have been moved onto it instead.
    q = ch1.query
    curves = q === nothing ? nothing : modelcurves(series, q)
    families = if q === nothing
        Tuple{Any,PlotChannel,Any}[]
    elseif calibrate === :data || !hasdata
        Tuple{Any,PlotChannel,Any}[(nothing, ch1, curvecolor)]
    else
        [(obs, ch, single ? curvecolor : _instcolor(j))
         for (j, (obs, ch)) in enumerate(entries)]
    end

    anywrapped = false
    for (obs, ch, color) in families
        xs = Float64[]; ys = Float64[]
        sizehint!(xs, (length(series.ts) + 1) * nshow)
        sizehint!(ys, (length(series.ts) + 1) * nshow)
        for d in 1:nshow
            v = curves[d]
            if obs !== nothing
                cal = datacalibration(obs, ch, obscontext(series, obs; draw=d), series.ts)
                cal === nothing || (v = v .+ cal)
            end
            v = v .* ch1.scale
            if ch1.wrap !== nothing
                xw, yw, wrapped = _wrap_series(series.ts, v, ch1.wrap)
                anywrapped |= wrapped
                append!(xs, xw); append!(ys, yw)
            else
                append!(xs, series.ts); append!(ys, v)
            end
            push!(xs, NaN); push!(ys, NaN)
        end
        Makie.lines!(ax, xs, ys; color=(color, _alpha(nshow)), linewidth)
    end
    anywrapped && Makie.ylims!(ax, 0.0, ch1.wrap)

    # Correlated-noise bands, per instrument and per draw, clipped to that
    # instrument's own baseline: a GP reverts to its prior outside the data,
    # and a band ballooning across the rest of the figure says nothing about
    # the fit.
    if gpband && curves !== nothing
        for (j, (obs, ch)) in enumerate(entries)
            obs === nothing && continue
            ep = Octofitter.epochs(obs)
            isempty(ep) && continue
            t0, t1 = extrema(ep)
            m = findall(t -> t0 <= t <= t1, series.ts)
            isempty(m) && continue
            tsub = series.ts[m]
            color = single ? curvecolor : _instcolor(j)
            α = min(0.35, 3.5 / max(nshow, 1))
            for d in 1:nshow
                ctx = obscontext(series, obs; draw=d)
                nm = noisemodel(obs, ctx, tsub)
                nm === nothing && break     # this observation has no noise model
                base = curves[d][m]
                if calibrate === :model
                    cal = datacalibration(obs, ch, ctx, tsub)
                    cal === nothing || (base = base .+ cal)
                end
                mid = (base .+ nm.mean) .* ch1.scale
                sd = sqrt.(nm.var) .* ch1.scale
                b = Makie.band!(ax, tsub, mid .- sd, mid .+ sd; color=(color, α))
                Makie.translate!(b, 0, 0, -10)     # behind the data
                Makie.lines!(ax, tsub, mid; color=(color, _alpha(nshow)), linewidth)
            end
        end
    end

    # Data on the main axis, per instrument.
    for (j, (obs, ch)) in enumerate(entries)
        obs === nothing && continue
        ctx = obscontext(series, obs)
        r = Octofitter.residuals(obs, ctx)[ch.name]
        u = r.use
        y = r.data
        if calibrate === :model
            cal = datacalibration(obs, ch, ctx, r.epoch)
            cal === nothing || (y = y .+ cal .* ch.scale)
        end
        _, σ_eff = _netresid(r)
        color = single ? :black : _instcolor(j)
        # Averaging windows, when the channel declares them. A catalog proper
        # motion is not measured *at* an epoch, it is averaged over a mission
        # span, and the bar is the only thing that says which span — it is why
        # a Hipparcos point may sit far from a curve the Gaia points hug. (v1
        # `absastromplot`'s one genuinely load-bearing idiom.)
        if haskey(r, :epoch_lo) && haskey(r, :epoch_hi)
            Makie.errorbars!(ax, r.epoch[u], y[u],
                (r.epoch .- r.epoch_lo)[u], (r.epoch_hi .- r.epoch)[u];
                direction=:x, color=:black, linewidth=0.7, whiskerwidth=4)
        end
        _behind!(Makie.errorbars!(ax, r.epoch[u], y[u], σ_eff[u];
            color=st.σeff_color, linewidth=st.σeff_linewidth))
        Makie.errorbars!(ax, r.epoch[u], y[u], r.σ[u];
            color=_σcolor(st, color), linewidth=st.σ_linewidth)
        Makie.scatter!(ax, r.epoch[u], y[u];
            color=single ? :white : color, marker=single ? :circle : _instmarker(j),
            strokecolor=:black, strokewidth=st.strokewidth, markersize=st.markersize,
            label=likelihoodname(obs))
    end
    # Say whose data these are. With one instrument a legend is four times the
    # ink of the answer, so the panel takes the corner label `phasefoldpanel!`
    # already uses for its row name — and a per-instrument panel *needs* it,
    # since nothing else on the panel distinguishes HARPS from HIRES.
    if hasdata && show_legend
        if single
            j = findfirst(e -> e[1] !== nothing, entries)
            Makie.Label(gs[1, 1], likelihoodname(entries[j][1]);
                font=:bold, fontsize=11, color=:grey30,
                halign=:right, valign=:top, padding=(0, 8, 0, 4),
                tellwidth=false, tellheight=false)
        else
            # Translucent backing rather than a frame: the legend sits inside
            # the axis and the model curves run underneath it.
            Makie.axislegend(ax, "Instrument"; framevisible=false, labelsize=9,
                titlesize=10, padding=(4, 4, 2, 2), backgroundcolor=(:white, 0.75))
        end
    end

    if ax_resid !== nothing && show_hist
        # `hist_aspect=1` makes the histogram box square whatever the strip's
        # height works out to, and `halign=:left` keeps it glued to the strip
        # rather than floating in the middle of the side column (v8's look;
        # `rvplot` pairs it with `colgap=0`).
        ax_hist = Makie.Axis(gs[2, 2]; xgridvisible=false, ygridvisible=false,
            aspect=hist_aspect, halign=hist_aspect === nothing ? :center : :left)
        Makie.hidedecorations!(ax_hist)
        Makie.colsize!(gs, 2, Makie.Fixed(side_w))
        Makie.colgap!(gs, 1, colgap)
    elseif hasdata
        Makie.Box(gs[2, 2]; visible=false)
        Makie.colsize!(gs, 2, Makie.Fixed(side_w))
        Makie.colgap!(gs, 1, colgap)
    else
        Makie.Box(gs[1, 2]; visible=false)
        Makie.colsize!(gs, 2, Makie.Fixed(side_w))
        Makie.colgap!(gs, 1, colgap)
    end
    # The strip's x axis is linked to the main one, which clips to the model
    # grid below — and that grid is shared by every panel in the figure, so
    # this is the span the boxes are actually drawn on.
    hasdata && _residstrip!(ax_resid, ax_hist, series, entries, identity, nshow, whiten, st;
        xspan=last(series.ts) - first(series.ts), boxwidth)

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
signal's upward zero crossing for radial velocity (the v8 rvpostplot
convention), periastron otherwise. Under the N-body propagator the curve is
the actual model over one osculating period anchored at the data midpoint —
periodicity is then an approximation, which is why `octoplot` only adds
phase panels automatically under `KeplerianApprox`.

Noise-weighted binned means (grey-to-red points) pool all instruments, over
`nbins` bins of phase; a bin has to hold at least two points to be drawn,
since a "mean" of one is the measurement itself, already on the axis with its
own error bar. A fold with fewer points than bins therefore shows no binned
series at all, and `show_binned=false` turns them off outright.
The residual strip below shows the same residuals as the time panel against
phase (whitened by default when more than one draw is shown, and then as a
per-point boxplot — see [`timeseriespanel!`](@ref)). `show_resid=false`
drops it — folding does not change a residual's value, so where a time panel
is already on the figure the strip is that panel's redrawn once per row; the
main axis then carries the phase ticks itself. `xticks` sets their spacing and
`datastyle` the marks (see [`timeseriespanel!`](@ref)).

Returns `(; main, resid, hist)`.
"""
function Octofitter.phasefoldpanel!(gp, series::PosteriorSeries, entries;
                                    row=nothing, ndraws=nothing, whiten=nothing,
                                    show_binned=nothing, nbins=20, nphase=201,
                                    bottom_axis=true, show_hist=true, show_resid=true,
                                    xticks=-0.5:0.25:0.5, boxwidth=nothing, datastyle=nothing,
                                    side_w=SIDE_W, colgap=COLGAP, hist_aspect=nothing,
                                    labelposition::Symbol=:inside, labelcolor=nothing,
                                    curvecolor=nothing, linewidth=nothing)
    sys = series.model.system
    gs = _layout(gp)
    st = _datastyle(datastyle, (; markersize=5))
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
    ndata = count(e -> e[1] !== nothing, entries)
    hasdata = ndata > 0
    single = ndata == 1
    whiten = _whitenmode(whiten, nshow, hasdata)
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
    # The phase axis normally lives on the residual strip glued underneath.
    # Without one — no data, or `show_resid=false` — the main axis has to carry
    # it or the panel ends up with no x axis at all.
    strip = hasdata && show_resid
    ax = Makie.Axis(gs[1, 1];
        ylabel, xticks,
        xlabel=(!strip && bottom_axis) ? "orbital phase" : "",
        xticklabelsvisible=(!strip && bottom_axis),
        xticksvisible=(!strip && bottom_axis),
        xgridvisible=false, ygridvisible=false)

    ax_resid = ax_hist = nothing
    if strip
        ax_resid = Makie.Axis(gs[2, 1];
            ylabel=_stripylabel(whiten),
            xlabel=bottom_axis ? "orbital phase" : "",
            xticks,
            xticklabelsvisible=bottom_axis, xticksvisible=bottom_axis,
            xgridvisible=false, ygridvisible=false)
        Makie.linkxaxes!(ax, ax_resid)
        ax.ytickformat = _blanklowesttick
        Makie.rowgap!(gs, 1, 0)
        Makie.rowsize!(gs, 1, Makie.Auto(2.4))
    end
    # The row this panel folds on. `:inside` puts it in the axis' top-right
    # corner in the row's accent colour, which is what a stacked octoplot
    # wants — the panels there have no spare column. `:side` is v8's
    # `rvpostplot` placement: hanging off the top right *outside* the axis, in
    # the side column, plain black. Colour defaults to the row's own rather
    # than `curvecolor`, so a caller drawing every curve in one colour
    # (`rvplot` does) still has the labels telling the planets apart.
    labelcolor = something(labelcolor,
        labelposition === :side ? Makie.to_color(:black) : _rowcolor(k))
    labelposition in (:inside, :side) ||
        error("`labelposition` must be :inside or :side, got :$labelposition")
    if labelposition === :side
        Makie.Label(gs[1, 2], String(rowname);
            font=:bold, fontsize=16, color=labelcolor,
            halign=:left, valign=:top, padding=(5, 0, 0, 0),
            tellwidth=false, tellheight=false)
    else
        Makie.Label(gs[1, 1], String(rowname);
            font=:bold, fontsize=14, color=labelcolor,
            halign=:right, valign=:top, padding=(0, 8, 0, 4),
            tellwidth=false, tellheight=false)
    end

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
            # What is folded is the residual *after* every nuisance term the
            # fit models — including the correlated-noise prediction, which is
            # stellar activity, not orbit, and would otherwise scatter the
            # folded points around a curve it has nothing to do with.
            resid, σ_eff = _netresid(r)
            s = sig_at_data[series.data_maps[obs]] .* ch.scale
            y = resid .+ s
            ch.wrap !== nothing && (y = mod.(y, ch.wrap))
            x = xmap(r.epoch)
            color = single ? :black : _instcolor(j)
            _behind!(Makie.errorbars!(ax, x[u], y[u], σ_eff[u];
                color=st.σeff_color, linewidth=st.σeff_linewidth))
            Makie.errorbars!(ax, x[u], y[u], r.σ[u];
                color=_σcolor(st, color), linewidth=st.σ_linewidth)
            Makie.scatter!(ax, x[u], y[u];
                color=single ? :white : color, marker=single ? :circle : _instmarker(j),
                strokecolor=:black, strokewidth=st.strokewidth, markersize=st.markersize)
            append!(pooled_x, x[u]); append!(pooled_y, y[u])
            append!(pooled_w, 1 ./ σ_eff[u] .^ 2)
        end

        # Noise-weighted binned means (the rvpostplot idiom, with the bin
        # mask half-width fixed).
        #
        # A bin holding a single point is not a mean of anything, and drawing
        # it as one says three false things at once: the mark lands on the bin
        # *centre* rather than on the measurement's own phase, it is bigger
        # and redder than the measurement, and its error bar is the scatter of
        # one value — zero — where the measurement's own is several m/s. In
        # the sparse limit every bin is that bin (ten points against twenty
        # bins gave ten zero-width red marks, each shouldered off its own
        # datum), so the binned series buried exactly the data it was made
        # from. Skip them: the points are already drawn, with their real
        # error bars, and a fold with fewer points than bins simply gets no
        # binned series — which is the honest answer. Bins with two or more
        # points are untouched, so a well-sampled figure is unchanged.
        if show_binned && !isempty(pooled_x)
            b = phasebinmeans(pooled_x, pooled_y, pooled_w, nbins)
            for i in eachindex(b.centre)
                c, μ, σ = b.centre[i], b.mean[i], b.sigma[i]
                Makie.errorbars!(ax, [c], [μ], [σ]; color=:black, linewidth=2)
                Makie.scatter!(ax, [c], [μ];
                    color=:red, markersize=9, strokecolor=:black, strokewidth=1.5)
            end
        end

        # The one thing on this panel that belongs to no row.
        #
        # A multi-row fold decomposes the *kinematic* signal, because `radvel`'s
        # Einstein term is quadratic in velocity and 1/r in separation and does
        # not telescope through A⁻¹ (see `_rowfunc`). It is therefore in the
        # plotted points but not in the plotted curve. That is bounded and
        # correct — the term is common-mode, like an instrument offset — but it
        # must not be *silent*, because a reader is entitled to assume the
        # scatter about the curve is noise. So when it rises above a tenth of
        # the residual scatter, the panel says how big it is.
        if sig.scaled && sig.query.func === Octofitter._kinrv && length(pooled_y) > 1
            kinq = Octofitter.ObservableQuery(Octofitter._kinrv, q.target, q.ref)
            ein = (Octofitter.evalquery(q, series.sys_map, series.data_traj_map) .-
                   Octofitter.evalquery(kinq, series.sys_map, series.data_traj_map)) .* ch1.scale
            span = isempty(ein) ? 0.0 : maximum(ein) - minimum(ein)
            scatter = std(pooled_y)
            if isfinite(span) && isfinite(scatter) && span > 0.1 * scatter
                unit = isempty(ch1.unit) ? "" : " " * ch1.unit
                Makie.Label(gs[1, 1],
                    "+ Einstein term, $(round(span; sigdigits=2))$unit peak-to-peak (not folded)";
                    fontsize=10, color=Makie.to_color(:grey35),
                    halign=:left, valign=:bottom, padding=(8, 0, 0, 4),
                    tellwidth=false, tellheight=false)
            end
        end
    end

    if ax_resid !== nothing && show_hist
        ax_hist = Makie.Axis(gs[2, 2]; xgridvisible=false, ygridvisible=false,
            aspect=hist_aspect, halign=hist_aspect === nothing ? :center : :left)
        Makie.hidedecorations!(ax_hist)
        Makie.colsize!(gs, 2, Makie.Fixed(side_w))
        Makie.colgap!(gs, 1, colgap)
    else
        # Reserve the side column even when nothing is drawn in it, so this
        # panel's main axis lines up with every other panel's. A `:side` row
        # label lives in the same cell; it declares no size, so it neither
        # widens the column nor replaces this.
        Makie.Box(gs[strip ? 2 : 1, 2]; visible=false)
        Makie.colsize!(gs, 2, Makie.Fixed(side_w))
        Makie.colgap!(gs, 1, colgap)
    end
    strip && _residstrip!(ax_resid, ax_hist, series, entries, xmap, nshow, whiten, st;
        xspan=1.0, boxwidth)      # a fold axis is always one full cycle

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

# The quantity a channel measures, independent of who measured it: two
# instruments' `sep` channels share this key, and so does the `sep` a ra/dec
# table exposes by projection.
_quantitykey(obs, ch) = ch.query === nothing ?
                        Symbol(nameof(typeof(Octofitter.plotobs(obs))), :_, ch.name) :
                        Symbol(nameof(ch.query.func), :_,
                            _refstr(ch.query.target), :_, _refstr(ch.query.ref))
# (The query-less branch keys by the observation *type*, so it has to be the
# type that measured it — `plotobs` — or a wrapped and an unwrapped Gaia
# observation would put the same `:along_scan` on two different axes. The
# `nothing` observation `_predictivegroups` passes only ever reaches the other
# branch: a predicted channel always carries a query.)

function _channelgroups(sys, channels=nothing)
    obss = [obs for obs in plottable_observations(sys) if isempty(defaultpanels(obs))]
    # Which quantities some observation measures *natively*. A derived channel
    # joins one of those panels rather than opening one of its own: a fit given
    # only ra/dec should not sprout separation and position-angle panels it was
    # never asked for, but the moment one instrument reports sep/pa, every
    # other instrument's points belong on that panel too.
    native = Set{Symbol}()
    for obs in obss, ch in plotchannels(obs)
        ch.derived || push!(native, _quantitykey(obs, ch))
    end

    groups = Vector{Pair{Symbol,Vector{Tuple{Any,PlotChannel}}}}()
    for obs in obss
        for ch in plotchannels(obs)
            _channelmatch(ch, channels) || continue
            key = _quantitykey(obs, ch)
            # An explicit `channels=` naming a derived channel is a request for
            # it; otherwise it only rides along with a native panel.
            (ch.derived && channels === nothing && !(key in native)) && continue
            # Only observation types that say so share a panel — see
            # `Octofitter.sharepanel`. Everything else is keyed by identity and
            # gets its own, with its data left uncalibrated.
            gkey = sharepanel(obs) ? key : Symbol(key, :_, likelihoodname(obs))
            i = findfirst(g -> g.first === gkey, groups)
            if i === nothing
                push!(groups, gkey => [(obs, ch)])
            else
                push!(groups[i].second, (obs, ch))
            end
        end
    end
    return groups
end

# Model-only groups for observables the model has no data for, so that
# `channels=radvel` on an astrometry-only fit predicts the RV curve instead of
# drawing an empty figure. Only an observable *function or its name* triggers
# this: `channels=:sep` names the separation **data**, and asking for data a
# fit does not have is not a request to invent it.
function _predictivegroups(sys, channels, groups)
    reqs = _requestedobservables(channels)
    isempty(reqs) && return Pair{Symbol,Vector{Tuple{Any,PlotChannel}}}[]
    have = Set{Any}()
    for (_, entries) in groups
        q = entries[1][2].query
        q === nothing || push!(have, q.func)
    end
    out = Pair{Symbol,Vector{Tuple{Any,PlotChannel}}}[]
    for f in reqs
        f in have && continue
        for (_, ch) in predictedchannels(sys, f)
            push!(out, _quantitykey(nothing, ch) => Tuple{Any,PlotChannel}[(nothing, ch)])
        end
    end
    return out
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

# `res.axes.rv` for a single RV instrument, `res.axes.rv_HARPS` /
# `res.axes.rv_HIRES` when there are several — the channel name alone stops
# identifying a panel as soon as each instrument has its own.
function _panelnames(groups)
    base = [entries[1][2].name for (_, entries) in groups]
    names = Symbol[]
    out = Symbol[]
    for (i, (_, entries)) in enumerate(groups)
        b = base[i]
        if count(==(b), base) > 1 && entries[1][1] !== nothing
            b = Symbol(b, :_, Octofitter.normalizename(likelihoodname(entries[1][1])))
        end
        push!(out, _uniquename!(names, b))
    end
    return out
end

"""
    octoplot(model, chain; N=250, seed=0, show_sky=nothing, show_phase=nothing,
             whiten=nothing, channels=nothing, gpband=nothing, boxwidth=nothing,
             tmin=nothing, tmax=nothing, ts=nothing,
             curvecolor=nothing, datastyle=nothing,
             fname=nothing, figscale=1, ndraws=nothing)

See the core docstring (`?octoplot` without Makie loaded). Panels: a sky
panel when the system has angular observables; one time-series panel per
channel group; phase-folded panels; then any bespoke `defaultpanels`. All
epoch axes are linked — calendar dates along the top of the figure, MJD along
the bottom.

**The time axis** is the [`PosteriorSeries`](@ref) epoch grid, which every
panel clips to. `tmax="2035-01-01"` (an MJD number, a date string, or a
`Date`) extends it past the data — `xlims!` alone cannot, since there is no
curve out there to show — `tmin=` moves the other end, and `ts=` replaces the
grid with epochs of your own. They belong to the series, so pass them to
`octoplot(model, chain; …)` or to the `PosteriorSeries` constructor, not to
`octoplot(series)`.

**Colours.** `curvecolor=` sets the model curves': one colour for the whole
figure, or a `Dict`/`NamedTuple` keyed by body name (`(; b=:firebrick)`) for
one orbit at a time — the sky panel builds its phase ramp from it, and every
panel of that body's signal follows. `datastyle=` overrides the marks (see
[`timeseriespanel!`](@ref)). Instrument colours and marker shapes stay
automatic: on a shared panel they are the only thing saying whose point is
whose.

**One panel per instrument.** Channels group by the quantity they measure
*and* by which observation reported it, unless the observation type opts into
sharing ([`sharepanel`](@ref)). Relative astrometry does — `platescale` and
`northangle` are pinned tightly enough that every instrument's points land in
the same place, so they merge onto one separation, position-angle, Δα* and Δδ
panel apiece. Radial velocity does not: a zero point is a free parameter of
order the data range, so each spectrograph gets its own panel, its data are
drawn exactly as reported, and each draw's model curve carries that draw's own
offset and trend. [`rvplot`](@ref) is the single-draw figure that puts them
back on one axis.

**Every measurement on every panel it belongs on.** (sep, pa) and (Δα*, Δδ)
are the same measurement in two bases, so a fit mixing the two conventions
shows all of its astrometry on all four panels, projected deterministically.

**Residuals** are whitened by default when more than one draw is shown, and
that is not optional: `whiten=false` needs `ndraws=1`, because a raw residual
is not comparable between draws. Each point's z-score varies from draw to draw
as well, so a whitened strip over many draws draws the *distribution* per
point as a boxplot rather than a single mark; `boxwidth=` sets the box width
in x units when the automatic choice is wrong for a dataset.

**Correlated-noise bands** ([`noisemodel`](@ref)) are drawn only when a single
draw is shown: 250 draws of "orbit" and 250 of "orbit + activity" on one axis
is twice the ink for an envelope that is a different draw's activity model at
every epoch. `gpband=true` overrides that; `rvplot` is the single-draw figure
where the band is the point.

`channels=` restricts the figure: an observable function (`channels=radvel`),
a channel or observable name as a `Symbol`, or a collection of either.
Naming an *observable* the model has no data for draws it as a **prediction** —
model curves alone, over the queries [`default_queries`](@ref) picks — which
is how you ask a fit what RV signal it implies before taking a spectrum.
Bespoke panels are dropped along with everything else not asked for.

**Phase folds are off** unless a single draw is being shown. A fold puts the
data on one ephemeris — the MAP's — and a posterior with a broad or
multi-modal period does not have one; see the `do_phase` comment below.
`show_phase=true` asks for them anyway and then folds *every* channel that can
fold, not just radial velocity. [`rvplot`](@ref) shows one draw and folds by
default, which is the conventional figure.
"""
function Octofitter.octoplot(model::Octofitter.LogDensityModel, chain::Chains;
                             N=250, seed=0, tmin=nothing, tmax=nothing, ts=nothing,
                             kwargs...)
    series = PosteriorSeries(model, chain; N, seed, tmin, tmax, ts)
    return Octofitter.octoplot(series; kwargs...)
end

function Octofitter.octoplot(series::PosteriorSeries;
                             show_sky=nothing, show_phase=nothing, whiten=nothing,
                             channels=nothing, figure=nothing, legend=true,
                             gpband=nothing, boxwidth=nothing,
                             curvecolor=nothing, datastyle=nothing,
                             fname=nothing, figscale=1.0, ndraws=nothing)
    sys = series.model.system
    nshow = ndraws === nothing ? length(series) : min(ndraws, length(series))
    # Say so when an observation's data cannot be drawn, so that an empty panel
    # reads as a gap in the plotting layer rather than as a fit with no data.
    Octofitter._warn_unplottable(sys)
    groups = _channelgroups(sys, channels)
    predictive = _predictivegroups(sys, channels, groups)
    if !isempty(predictive)
        @info("Nothing in this model measures " *
              join(unique(String(e[1][2].name) for (_, e) in predictive), ", ") *
              ": drawing what the fitted orbits predict, with no data to compare against.")
        append!(groups, predictive)
    end
    pnames = _panelnames(groups)
    do_sky = show_sky === nothing ?
             (_isangular(sys) && !isempty(sys.rows)) : show_sky

    # Phase panels, and by default only for a single draw.
    #
    # A fold collapses the epoch axis through *one* ephemeris. That is exact
    # for one draw, and for many draws it is a claim about the posterior that
    # the posterior need not support: the data can only be folded on one
    # period (here the MAP's), so a chain with two period modes — or simply a
    # period whose spread is a fair fraction of the baseline — has its points
    # placed at phases most of the drawn orbits disagree with, and the panel
    # reads as scatter about a curve rather than as the ambiguity it is. Some
    # posteriors are tight enough for it to look fine, which is the trap.
    # `rvplot` is the single-draw figure and folds by default; `show_phase=true`
    # asks for it here regardless, and then folds anything foldable, i.e. any
    # channel whose query one hierarchy row's signal can be isolated from.
    # Only meaningful under a periodic propagator either way.
    do_phase = show_phase === nothing ?
               (nshow == 1 && sys.method isa PlanetOrbits.KeplerianApprox) : show_phase
    foldall = show_phase === true
    phasepanels = Tuple{Symbol,Int,Any}[]   # (name, row, entries)
    if do_phase
        names = copy(pnames)
        for ((key, entries), nm) in zip(groups, pnames)
            q = entries[1][2].query
            q === nothing && continue
            (foldall || q.func === PlanetOrbits.radvel) || continue
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
        skyaxes = Octofitter.skypanel!(fig[row, 1], series; ndraws, legend, curvecolor)
        push!(axpairs, :sky => skyaxes)
        Makie.rowsize!(fig.layout, row, Makie.Fixed(figscale * skyH))
        row += 1
    end
    # A shared panel's instrument legend repeats verbatim on every panel the
    # same instruments appear on — four identical legends over a sep/pa/Δα*/Δδ
    # stack. Draw it on the first, and let the rest inherit it by eye.
    legendseen = Set{Any}()
    for (i, ((key, entries), nm)) in enumerate(zip(groups, pnames))
        insts = Tuple(likelihoodname(o) for (o, _) in entries if o !== nothing)
        showleg = legend && (length(insts) <= 1 || !(insts in legendseen))
        push!(legendseen, insts)
        axs = Octofitter.timeseriespanel!(fig[row, 1], series, entries;
            top_time_axis=(i == 1), bottom_time_axis=(i == npanels), ndraws, whiten,
            gpband, boxwidth, show_legend=showleg, datastyle,
            curvecolor=_panelcolor(curvecolor, _curvebody(sys, entries[1][2].query)))
        push!(axpairs, nm => axs)
        push!(timeaxes, axs.main)
        axs.resid !== nothing && push!(timeaxes, axs.resid)
        Makie.rowsize!(fig.layout, row, Makie.Auto(1.0))
        row += 1
    end
    for (i, (nm, k, entries)) in enumerate(phasepanels)
        axs = Octofitter.phasefoldpanel!(fig[row, 1], series, entries;
            row=k, ndraws, whiten, boxwidth, bottom_axis=(i == length(phasepanels)),
            datastyle, curvecolor=_panelcolor(curvecolor, sys.rows[k][1]))
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
# rvplot — the single-draw radial-velocity figure
#
# The one figure that deliberately breaks octoplot's one-panel-per-instrument
# rule, and the single draw is what pays for it: with every nuisance parameter
# pinned, a calibrated RV series is a well-defined quantity and every
# instrument can be drawn on one axis in the same frame. Assembled from the
# same generic panels as everything else — `timeseriespanel!` in `:data`
# calibration mode, then one `phasefoldpanel!` per hierarchy row that moves
# the star — so the error-bar conventions, the binned means, the GP band and
# the returned axes all match the rest of the figure set.
# ---------------------------------------------------

"""
    rvplot(model, chain, [sample_idx]; kwargs...)
    rvplot(series::PosteriorSeries; kwargs...)

The radial-velocity summary figure for a **single** posterior draw: one time
series panel carrying every instrument at once, with a residual strip and a
marginal histogram, then one phase-folded panel per hierarchy row that moves
the star, each with noise-weighted binned means.

`sample_idx` defaults to the highest-posterior-density sample. Pass an index
to render a different draw, and `rvplot_animated` to sweep through many.

The single draw is the point of this figure, not an economy. A calibrated RV
series is only defined *per draw* — the instrument zero points, the jitters,
the trend and the other bodies' signals that are subtracted before folding all
move from sample to sample — so it is exactly by fixing one draw that several
instruments can honestly share an axis. [`octoplot`](@ref) is the many-draw
view and gives each instrument a panel of its own instead:

    octoplot(model, chain; show_sky=false, channels=radvel)

## Keywords

| Keyword | Meaning |
|---|---|
| `show_phase=true` | phase-folded panels |
| `show_hist=true` | marginal residual histograms |
| `show_legend=true` | the instrument + symbol legend, in a column down the right |
| `show_phase_resid=false` | a residual strip under each phase panel too |
| `show_perspective=false` | the frame's own radial-velocity drift (see below) |
| `gpband=true` | the correlated-noise band, where a `gaussian_process` was fitted |
| `whiten=nothing` | `false` (raw residuals) by default — one draw makes them well defined |
| `curvecolor=:blue` | one colour for every panel's model curve; `nothing` gives each row its own |
| `datastyle=nothing` | marker and error-bar overrides, merged over `_rvdatastyle` |
| `tmin`, `tmax`, `ts` | the epoch grid the curve is drawn over, as `octoplot` (the `model, chain` form only) |
| `figscale`, `fname`, `figure` | as `octoplot` |

The figure follows v8 `rvpostplot`'s layout: landscape, legend in a column at
the right, and the phase panels showing the fold on its own. Folding
rearranges the residuals along x but leaves their values alone, so a strip
under each phase panel is the time panel's residuals redrawn once per planet
— `show_phase_resid=true` if you want them regardless.

`show_perspective=true` overlays the secular drift of the system barycentre's
own radial velocity, `PlanetOrbits.frame_rv`, for a model with an absolute
frame. Whether that drift is part of the *fit* is each dataset's own
declaration (`secular_acceleration=:model`, the default for an absolute
`RadialVelocityObs`, or `:data_corrected`) — but either way it is not in this
figure's model curve, which is the pure `radvel` query, while the plotted
points are calibrated only by `offset` and `trend_function`. The overlay is
off by default and drawn dashed for that reason: it shows the size of a term
the curve leaves out. Residuals and phase folds go through the full forward
model and are unaffected.

Returns an [`OctoPlotResult`](@ref); `res.axes.rv.main` and
`res.axes.rv_phase_b.main` name the panels.
"""
function Octofitter.rvplot(model::Octofitter.LogDensityModel, chain::Chains,
                           sample_idx::Integer=_max_logpost_index(chain);
                           tmin=nothing, tmax=nothing, ts=nothing, kwargs...)
    nsamples = size(chain, 1) * size(chain, 3)
    1 <= sample_idx <= nsamples || throw(ArgumentError(
        "sample_idx=$sample_idx is outside the chain's 1:$nsamples samples."))
    return Octofitter.rvplot(
        PosteriorSeries(model, chain; ii=[sample_idx], tmin, tmax, ts); kwargs...)
end

# v1 selected the draw with `argmax(results["logpost"][:])`. Chains that were
# not produced by a sampler carrying a log posterior (a prior draw assembled by
# hand, say) have no such column, and falling back to the first sample is
# better than failing to plot at all.
function _max_logpost_index(chain::Chains)
    :logpost in keys(chain) || return 1
    return argmax(vec(chain[:logpost]))
end

# Every radial-velocity channel in the model, grouped by *quantity* — a
# stellar reflex `radvel(A, Barycentre)` and a relative `radvel(b, A)` are
# different measurements and never share an axis, however many instruments
# report each.
function _rvgroups(sys)
    groups = Vector{Pair{Symbol,Vector{Tuple{Any,PlotChannel}}}}()
    for obs in plottable_observations(sys)
        isempty(defaultpanels(obs)) || continue
        for ch in plotchannels(obs)
            (ch.query !== nothing && ch.query.func === PlanetOrbits.radvel) || continue
            key = _quantitykey(obs, ch)
            i = findfirst(g -> g.first === key, groups)
            i === nothing ? push!(groups, key => Tuple{Any,PlotChannel}[(obs, ch)]) :
            push!(groups[i].second, (obs, ch))
        end
    end
    return groups
end

# Through `plotobs`: a wrapper carries the wrapped table but not its
# `gaussian_process` field, and this only decides a legend entry for a band
# `noisemodel` (which *is* delegated) would draw anyway.
_hasgp(obs) = let o = Octofitter.plotobs(obs)
    hasproperty(o, :gaussian_process) && o.gaussian_process !== nothing
end

# v8 `rvpostplot`'s marks: small filled points, hairline strokes, and the
# measurement bar in the instrument's own colour inside the grey
# jitter-inflated one. See `DATASTYLE` for why this figure needs its own.
_rvdatastyle(user) = merge((; markersize=5, resid_markersize=5, strokewidth=0.4,
        σ_linewidth=0.7, σ_color=:instrument, σeff_linewidth=0.7),
    something(user, (;)))

function Octofitter.rvplot(series::PosteriorSeries;
                           show_phase=true, show_hist=true, show_legend=true,
                           show_phase_resid=false,
                           show_perspective=false, whiten=nothing, gpband=true,
                           curvecolor=:blue, datastyle=nothing,
                           figure=nothing, fname=nothing, figscale=1.0)
    sys = series.model.system
    length(series) == 1 || error("""
    `rvplot` draws one posterior sample, and this series carries $(length(series)).
    Sharing one axis between instruments is only honest when every offset,
    jitter and trend is pinned, which is what a single draw does.

        rvplot(model, chain)          # the highest-posterior-density sample
        rvplot(model, chain, 42)      # a sample you pick
        rvplot(PosteriorSeries(model, chain; ii=[42]))

    For the posterior spread, `octoplot(model, chain; channels=radvel)` gives
    each instrument its own panel and every draw its own calibrated curve.""")

    groups = _rvgroups(sys)
    isempty(groups) && error(
        "`rvplot`: this model has no radial-velocity observations. " *
        "(`octoplot` draws whatever it does have.)")

    allobs = [obs for (_, entries) in groups for (obs, _) in entries]
    anygp = any(_hasgp, allobs)
    absframe = series.sys_map.frame isa PlanetOrbits.AbsoluteFrame
    show_perspective &= absframe

    nphase = show_phase ?
             sum(length(foldablerows(series.sys_map, e[1][2].query)) for (_, e) in groups) : 0

    # v8 `rvpostplot`'s proportions: a wide figure with the legend in a column
    # down the right, a tall time panel over its glued residual strip (3.4 row
    # units between them), and squatter phase panels below (2 each).
    #
    # The row unit is v8's `305 + 135·n_planets` solved for one row, plus ~11%:
    # v8's y label was "RV [m/s]" and v9's is the resolver's "radial velocity
    # [m/s]", which at v8's exact heights is long enough to run into the
    # residual strip's own label on a one-planet figure.
    #
    # A phase panel showing the fold alone is squat; one that has been given a
    # residual strip is the same two-axis stack as the time panel, and needs
    # the same room for the same reason.
    phaseunits = show_phase_resid ? 3.4 : 2.0
    rowunits = 3.4 * length(groups) + phaseunits * nphase
    W = round(Int, figscale * (480 + (show_legend ? LEGEND_W + 10 : 0)))
    H = round(Int, figscale * (75 + 75.0 * rowunits))
    # The side column holds a square marginal histogram glued to the residual
    # strip, and the phase panels' row labels. One row unit is the strip's
    # nominal height, less what the axis decorations take off it — sized so the
    # square is limited by the strip's height rather than by this width, and
    # so little is left over that the histogram reads as flush.
    sidew = figscale * 64.0
    fig = figure === nothing ? Makie.Figure(size=(W, max(H, 300))) : figure

    axpairs = Pair{Symbol,Any}[]
    timeaxes = Makie.Axis[]
    pnames = _panelnames(groups)
    names = copy(pnames)
    row = 1
    for ((_, entries), nm) in zip(groups, pnames)
        # The epoch axis belongs to the time panel whether or not phase panels
        # follow it: the phase panels have their own x axis, so deferring the
        # MJD labels to "the bottom of the figure" would lose them entirely.
        axs = Octofitter.timeseriespanel!(fig[row, 1], series, entries;
            top_time_axis=true, bottom_time_axis=(row == length(groups)),
            show_hist, whiten, gpband, calibrate=:data, show_legend=false,
            curvecolor, linewidth=1.0, datastyle=_rvdatastyle(datastyle),
            side_w=sidew, colgap=0, hist_aspect=1)
        push!(axpairs, nm => axs)
        push!(timeaxes, axs.main)
        axs.resid === nothing || push!(timeaxes, axs.resid)
        Makie.rowsize!(fig.layout, row, Makie.Auto(3.4))
        row += 1

        # The barycentre's own radial velocity, as a drift about its value at
        # the first plotted epoch (the constant is degenerate with every
        # instrument offset and would only shift the line off the panel).
        if show_perspective
            fr = [PlanetOrbits.frame_rv(sol) for sol in series.trajs[1]]
            Makie.lines!(axs.main, series.ts, fr .- first(fr);
                color=:darkorange, linewidth=1.5, linestyle=:dash)
        end
    end

    if show_phase
        phasepanels = Tuple{Symbol,Int,Any}[]
        for ((_, entries), nm) in zip(groups, pnames)
            for k in foldablerows(series.sys_map, entries[1][2].query)
                push!(phasepanels,
                    (_uniquename!(names, Symbol(nm, :_phase_, sys.rows[k][1])), k, entries))
            end
        end
        for (i, (nm, k, entries)) in enumerate(phasepanels)
            # No marginal histogram down here, and by default no residual strip
            # either: folding rearranges the residuals along x and leaves their
            # distribution — and, point for point, their values — alone, so
            # both would be the time panel's residuals redrawn once per planet.
            # v8 showed the fold on its own for exactly that reason.
            axs = Octofitter.phasefoldpanel!(fig[row, 1], series, entries;
                row=k, whiten, show_hist=false, show_resid=show_phase_resid,
                bottom_axis=(i == length(phasepanels)),
                xticks=-0.5:0.1:0.5, curvecolor, linewidth=4.0,
                datastyle=_rvdatastyle(datastyle),
                side_w=sidew, colgap=0, labelposition=:side)
            push!(axpairs, nm => axs)
            Makie.rowsize!(fig.layout, row, Makie.Auto(phaseunits))
            row += 1
        end
    end
    length(timeaxes) > 1 && Makie.linkxaxes!(timeaxes...)

    # The legend is a column down the right of every panel, not a band under
    # them: it is what lets the figure keep v8's landscape proportions instead
    # of growing a row nobody reads first.
    if show_legend
        # `curvecolor=nothing` hands each panel back its own accent colour, so
        # there is no single "orbit model" colour to swatch; the time panel's
        # is the representative one, as it carries the whole signal.
        lc = something(curvecolor, _querycolor(sys, groups[1][2][1][2].query))
        _rvlegend!(fig[1:row-1, 2], allobs, anygp, show_perspective, lc)
        Makie.colsize!(fig.layout, 2, Makie.Fixed(figscale * LEGEND_W))
    end

    axes = (; axpairs...)
    fname !== nothing && Makie.save(fname, fig, px_per_unit=3)
    return OctoPlotResult(fig, axes, series)
end

# v1's two-group rvpostplot legend: which symbol is which spectrograph, and
# what each of the figure's other marks means. Both halves earn their space —
# the error bars are two nested segments whose meanings differ, and nothing
# else on the figure says which is which.
function _rvlegend!(gp, allobs, anygp, show_perspective, curvecolor)
    n = length(allobs)
    single = n == 1
    instruments = [Makie.MarkerElement(
        color=single ? :white : _instcolor(j), marker=single ? :circle : _instmarker(j),
        strokecolor=:black, strokewidth=1, markersize=11) for j in 1:n]
    # The measurement bar carries the instrument's colour on this figure, so
    # its legend entry has to be all of them at once: v8's composite element,
    # one short stroke per instrument side by side in a single legend slot.
    σelem = single ? Makie.LineElement(color=:black, linewidth=3) :
            [Makie.LineElement(color=_instcolor(j), linewidth=2,
                points=Makie.Point2f[((j - 0.5) / n, 0), ((j - 0.5) / n, 1)]) for j in 1:n]
    marks = Any[
        σelem,
        # Vertical, like the bar it stands for and like the entry above it —
        # a horizontal grey rule reads as a line style, not an error bar.
        Makie.LineElement(color="#CCCCCC", linewidth=2,
            points=Makie.Point2f[(0.5, 0), (0.5, 1)]),
        Makie.LineElement(color=curvecolor, linewidth=4),
        Makie.MarkerElement(color=:red, marker=:circle, markersize=11,
            strokecolor=:black, strokewidth=1.5),
    ]
    labels = Any[
        "measurement σ",
        anygp ? "σ + jitter\nand GP" : "σ + jitter",
        "orbit model",
        "binned",
    ]
    if anygp
        push!(marks, Makie.PolyElement(color=(curvecolor, 0.35)))
        push!(labels, "activity model")
    end
    if show_perspective
        push!(marks, Makie.LineElement(color=:darkorange, linewidth=3, linestyle=:dash))
        push!(labels, "frame RV drift\n(not in the curve)")
    end
    return Makie.Legend(gp, [instruments, marks],
        [[likelihoodname(o) for o in allobs], labels],
        ["instrument", "legend"];
        framevisible=false, labelsize=10, titlesize=11,
        rowgap=4, groupgap=18, patchsize=(20.0f0, 14.0f0),
        tellheight=false, tellwidth=false, valign=:top, halign=:left)
end

"""
    rvplot_animated(model, chain; N=50, seed=0, framerate=4,
                    fname="rv-posterior.mp4", kwargs...)

[`rvplot`](@ref) recorded over `N` single-draw slices of the chain — v8's
"sweep through the posterior" animation. Each frame is the whole figure
rebuilt for one draw, so every panel (including the phase folds, which move
with the drawn period) is that draw's own. Returns `fname`.

The extension is anything Makie can write from a `VideoStream`: `.mp4`,
`.gif`, `.mkv`.
"""
function Octofitter.rvplot_animated(model::Octofitter.LogDensityModel, chain::Chains;
                                    N::Integer=50, seed::Integer=0, framerate=4,
                                    fname::AbstractString="rv-posterior.mp4",
                                    tmin=nothing, tmax=nothing, ts=nothing,
                                    kwargs...)
    nsamples = size(chain, 1) * size(chain, 3)
    ii = sort!(StatsBase.sample(Xoshiro(seed), 1:nsamples, min(N, nsamples); replace=false))
    # One figure, reused: a `VideoStream` records a single scene, so each frame
    # empties the layout and rebuilds the panels into the same figure rather
    # than making a new one (which the stream would never see).
    res = Octofitter.rvplot(
        PosteriorSeries(model, chain; ii=[first(ii)], tmin, tmax, ts); kwargs...)
    fig = res.figure
    stream = Makie.VideoStream(fig; framerate)
    Makie.recordframe!(stream)
    for i in ii[2:end]
        empty!(fig)
        Octofitter.rvplot(PosteriorSeries(model, chain; ii=[i], tmin, tmax, ts);
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

Masses are **M⊙**: v9 has one mass unit throughout, where v8 plotted each
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
#
# Selected by `plotobs`, so a wrapped observation is still found, but the
# *wrapper* is returned: it is what `obscontext` and `series.data_maps` are
# keyed by. Same rule as the sky panel's overlay.
function _theobs(sys, T, what)
    os = [o for o in sys.observations if Octofitter.plotobs(o) isa T]
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

(v8 displaced by `model − data`, mirroring each point across the track. The
sign here is the one that makes a point sit where the abscissa says the source
was.)

`sample_idx` selects the draw — the maximum-posterior one by default, any row
index otherwise. Comparing draws is the point of [`gaiastarplot!`](@ref),
which draws into a cell of a figure you already have.
"""
function Octofitter.gaiastarplot(model::Octofitter.LogDensityModel, chain::Chains,
                                 sample_idx=nothing;
                                 fname=nothing, figure=(;), kwargs...)
    fig = Makie.Figure(; size=(600, 600), figure...)
    Octofitter.gaiastarplot!(fig[1, 1], model, chain, sample_idx; kwargs...)
    fname !== nothing && Makie.save(fname, fig, px_per_unit=3)
    return fig
end

"""
    gaiastarplot!(gp, model, chain, sample_idx=MAP; keplerian_mult=1,
                  ntrack=200, axis=(;))

[`gaiastarplot`](@ref) into `gp`: a grid cell (`fig[i, j]`), a whole figure or
layout — one axis at `[1, 1]` — or an `Axis` the caller has already made and
styled. Returns the `Axis`.

A grid of cells each showing a different `sample_idx` is the honest picture of
an astrometric posterior: several quite different orbits, each passing through
the same transits.
"""
function Octofitter.gaiastarplot!(gp, model::Octofitter.LogDensityModel, chain::Chains,
                                  sample_idx=nothing;
                                  keplerian_mult=1.0, ntrack::Integer=200,
                                  axis=(;))
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
    # `obs.table.scan_pos_angle` is in degrees (the archive's unit); the
    # radian projection is precomputed on the observation.
    s = obs.sinψ
    c = obs.cosψ
    u = r.use
    dx = mx .+ r.resid .* s
    dy = my .+ r.resid .* c

    ax = _panelaxis(gp;
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
    return ax
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
reads the published solution through `Octofitter.gaia_dr3_solution`. (v8 read
it from `obs.gaia_sol`, which the v9 observation no longer carries — it models
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
    # Degrees in the table, radians in the precomputed projection — see
    # `GaiaDR4AstromObs`.
    s = obs.sinψ
    c = obs.cosψ
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
likelihood uses. (v8 recovered the sign of an unsigned residual by trying both
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
             if Octofitter.plotobs(o) isa Octofitter.HipparcosIADObs ||
                (Octofitter.plotobs(o) isa Octofitter.G23HObs && :iad_hip ∈ o.table.kind)]
    isempty(cands) && error(
        "hipparcosplot needs a HipparcosIADObs, or a G23HObs keeping its " *
        ":iad_hip channel; this model has neither.")
    length(cands) > 1 && error("hipparcosplot draws one Hipparcos observation; " *
                               "this model has $(length(cands)).")
    obs = only(cands)
    series, _ = _onedraw(model, chain, sample_idx)
    ctx = obscontext(series, obs; draw=1)
    sim = Octofitter.simulate(obs, ctx)

    if Octofitter.plotobs(obs) isa Octofitter.HipparcosIADObs
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
