# ---------------------------------------------------
# Plotting API (backend-agnostic)
#
# Everything a plot needs that is not drawing lives here, in core, so that
# the Makie extension is presentation-only and other consumers (tables,
# goodness-of-fit summaries) can reuse it:
#
#   - `ObservableQuery` — "observable f of target vs ref", the unit of the
#     plotting grammar. Built from a function (or its name as a Symbol) and
#     two reference specs; never by overloading the observable itself.
#
#   - The observation plotting protocol: `plotchannels(obs)` declares the
#     1-D data channels an observation exposes, and `residuals(obs, ctx)`
#     computes calibrated data, model values, residuals and effective
#     uncertainties — with the *likelihood's* math (jitter, platescale,
#     northangle, offsets, outlier masks), never a plot-side reimplementation.
#
#   - `PosteriorSeries` — one bundle of (draw selection, rebuilt systems,
#     solved trajectories, per-observation contexts) shared by every panel,
#     so all panels show the same draws and the solve configuration always
#     matches what the likelihood saw.
# ---------------------------------------------------

using StatsBase: sample as _sb_sample

export refspecs, epochs

# ---------------------------------------------------
# Observable queries
# ---------------------------------------------------

const _OBSERVABLE_FUNCS = (
    PlanetOrbits.posx, PlanetOrbits.posy, PlanetOrbits.posz,
    PlanetOrbits.velx, PlanetOrbits.vely, PlanetOrbits.velz,
    PlanetOrbits.radvel, PlanetOrbits.raoff, PlanetOrbits.decoff,
    PlanetOrbits.pmra, PlanetOrbits.pmdec,
    PlanetOrbits.projectedseparation, PlanetOrbits.posangle,
)

function _obsfunc(s::Symbol)
    for f in _OBSERVABLE_FUNCS
        nameof(f) === s && return f
    end
    error("unknown observable :$s; expected one of " *
          join((":$(nameof(f))" for f in _OBSERVABLE_FUNCS), ", "))
end
_obsfunc(f::Function) = f

"""
    ObservableQuery(f, target, ref)

The quantity `f(sol, target, ref)` as a value: an observable function (or
its name as a `Symbol`) plus a target and a reference in the usual grammar
(a `Body` node, a `Symbol`, `Barycentre`/`Barycentre(A, b)`, `Photocentre`).

    ObservableQuery(radvel, :A, Barycentre)     # stellar reflex RV
    ObservableQuery(:radvel, :c, :b)            # planet–planet relative RV
    ObservableQuery(raoff, Barycentre(:Aa, :Ab), Barycentre)

Anywhere a query is accepted, a plain tuple `(f, target, ref)` works too.
"""
struct ObservableQuery{F,TT,TR}
    func::F
    target::TT
    ref::TR
end
function ObservableQuery(f::Union{Function,Symbol}, target, ref)
    fn, t, r = _obsfunc(f), refspec(target), refspec(ref)
    # Explicit type parameters: the bare 3-arg call would re-enter this
    # convenience method (a normalized Function is still a Function).
    return ObservableQuery{typeof(fn),typeof(t),typeof(r)}(fn, t, r)
end
export ObservableQuery

_query(q::ObservableQuery) = q
_query(t::Tuple{Any,Any,Any}) = ObservableQuery(t[1], t[2], t[3])

_querystr(q::ObservableQuery) =
    string(nameof(q.func), "(", _refstr(q.target), ", ", _refstr(q.ref), ")")

Base.show(io::IO, q::ObservableQuery) = print(io, "ObservableQuery ", _querystr(q))

"""
    evalquery(q, posys, traj) -> Vector

Evaluate a query on a solved trajectory: `q.func(sol, target, ref)` at every
epoch, with the references resolved against `posys` once.
"""
function evalquery(q::ObservableQuery, posys, traj)
    t = resolveref(posys, q.target)
    r = resolveref(posys, q.ref)
    return [q.func(sol, t, r) for sol in traj]
end

# ---------------------------------------------------
# Observation plotting protocol
# ---------------------------------------------------

"""
    PlotChannel

One 1-D data channel of an observation, as declared by
[`plotchannels`](@ref):

  - `name` — key into the NamedTuple [`residuals`](@ref) returns.
  - `label`, `unit` — axis labelling.
  - `scale` — display scale factor (e.g. `rad2deg(1)` for position angles);
    `residuals` already applies it, and a smooth model curve from `query`
    must be multiplied by it too.
  - `wrap` — display period the quantity wraps at (in *display* units, e.g.
    `360` for position angle in degrees), or `nothing`.
  - `query` — an [`ObservableQuery`](@ref) whose curve is this channel's
    smooth model prediction, or `nothing` when the channel has no meaning
    off the data epochs (e.g. Gaia along-scan abscissae, which depend on
    per-transit scan angles).
"""
struct PlotChannel{TQ}
    name::Symbol
    label::String
    unit::String
    scale::Float64
    wrap::Union{Nothing,Float64}
    query::TQ
end
PlotChannel(name, label, unit; scale=1.0, wrap=nothing, query=nothing) =
    PlotChannel{typeof(query)}(name, label, unit, scale, wrap, query)
export PlotChannel

"""
    plotchannels(obs) -> Tuple{Vararg{PlotChannel}}

The 1-D data channels this observation exposes for plotting; `()` if it is
not plottable as data-vs-model series (prior terms return the default).

Implement this alongside `ln_like` and [`residuals`](@ref) when adding an
observation type; a generic time-series panel then works with no plot-side
code.
"""
plotchannels(::AbstractObs) = ()
export plotchannels

"""
    residuals(obs, ctx::ObsContext) -> NamedTuple

Calibrated data, model predictions, residuals and uncertainties for each
channel of [`plotchannels`](@ref), at the observation's epochs, under the
parameters in `ctx` — using exactly the likelihood's math (jitters in
quadrature, platescale/northangle applied, instrument offsets removed,
outlier masks respected).

Returns `(; <channel name> = (; epoch, data, model, resid, σ, σ_eff, use), …)`
where `data` is calibrated into the model frame (so it overlays a pure-orbit
model curve), `resid` is `data - model` in the same (display-scaled) units,
`σ` is the measurement uncertainty and `σ_eff` includes fitted jitter.
`use` is a Bool vector; `false` marks points the likelihood excluded.

Not exported to avoid clashing with `StatsBase.residuals`; call it as
`Octofitter.residuals`.
"""
function residuals end

# --- RelAstromObs -------------------------------------------------------------

function plotchannels(obs::RelAstromObs)
    if hasproperty(obs.table, :sep)
        return (
            PlotChannel(:sep, "separation", "mas";
                query=ObservableQuery(PlanetOrbits.projectedseparation, obs.target, obs.ref)),
            PlotChannel(:pa, "position angle", "°"; scale=rad2deg(1.0), wrap=360.0,
                query=ObservableQuery(PlanetOrbits.posangle, obs.target, obs.ref)),
        )
    else
        return (
            PlotChannel(:raoff, "Δα*", "mas";
                query=ObservableQuery(PlanetOrbits.raoff, obs.target, obs.ref)),
            PlotChannel(:decoff, "Δδ", "mas";
                query=ObservableQuery(PlanetOrbits.decoff, obs.target, obs.ref)),
        )
    end
end

function residuals(obs::RelAstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    θ_obs = ctx.θ_obs
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)
    platescale = hasproperty(θ_obs, :platescale) ? θ_obs.platescale : one(T)
    northangle = hasproperty(θ_obs, :northangle) ? θ_obs.northangle : zero(T)

    sim = simulate(obs, ctx)
    ep = collect(Float64, obs.table.epoch)
    n = length(ep)
    use = trues(n)

    if hasproperty(obs.table, :sep)
        sep_model = hypot.(sim.ra_model, sim.dec_model)
        pa_model = atan.(sim.ra_model, sim.dec_model)
        sep_data = obs.table.sep .* platescale
        pa_data = obs.table.pa .+ northangle
        # Wrapped PA difference, exactly as ln_like computes it.
        pa_resid = @. mod(pa_data - pa_model + π, 2π) - π
        r2d = rad2deg(1.0)
        return (;
            sep=(; epoch=ep, data=collect(sep_data), model=sep_model,
                resid=collect(sep_data .- sep_model),
                σ=collect(float.(obs.table.σ_sep)),
                σ_eff=collect(hypot.(obs.table.σ_sep, jitter)), use),
            pa=(; epoch=ep, data=collect(rem2pi.(pa_data, RoundDown) .* r2d),
                model=rem2pi.(pa_model, RoundDown) .* r2d,
                resid=collect(pa_resid .* r2d),
                σ=collect(float.(obs.table.σ_pa) .* r2d),
                σ_eff=collect(hypot.(obs.table.σ_pa, jitter ./ max.(sep_model, one(T))) .* r2d), use),
        )
    else
        # Mirrors ln_like's ra/dec branch exactly: its pa_dat is measured from
        # the RA axis (atan(dec, ra), *not* the position-angle convention), so
        # cos pairs with RA and sin with Dec.
        pa_dat = atan.(obs.table.dec, obs.table.ra) .+ northangle
        sep_dat = hypot.(obs.table.dec, obs.table.ra) .* platescale
        ra_data = sep_dat .* cos.(pa_dat)
        dec_data = sep_dat .* sin.(pa_dat)
        return (;
            raoff=(; epoch=ep, data=collect(ra_data), model=collect(sim.ra_model),
                resid=collect(ra_data .- sim.ra_model),
                σ=collect(float.(obs.table.σ_ra)),
                σ_eff=collect(hypot.(obs.table.σ_ra, jitter)), use),
            decoff=(; epoch=ep, data=collect(dec_data), model=collect(sim.dec_model),
                resid=collect(dec_data .- sim.dec_model),
                σ=collect(float.(obs.table.σ_dec)),
                σ_eff=collect(hypot.(obs.table.σ_dec, jitter)), use),
        )
    end
end

# --- RadialVelocityObs ---------------------------------------------------------

plotchannels(obs::RadialVelocityObs) = (
    PlotChannel(:rv, "radial velocity", "m/s";
        query=ObservableQuery(PlanetOrbits.radvel, obs.target, obs.ref)),
)

function residuals(obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    θ_obs = ctx.θ_obs
    offset = hasproperty(θ_obs, :offset) ? θ_obs.offset : zero(T)
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)
    sim = simulate(obs, ctx)               # includes the offset
    ep = collect(Float64, obs.table.epoch)
    model_pure = sim.rv_model .- offset    # matches the radvel query curve
    data_cal = obs.table.rv .- offset
    return (;
        rv=(; epoch=ep, data=collect(data_cal), model=collect(model_pure),
            resid=collect(obs.table.rv .- sim.rv_model),
            σ=collect(float.(obs.table.σ_rv)),
            σ_eff=collect(hypot.(obs.table.σ_rv, jitter)),
            use=trues(length(ep))),
    )
end

# --- GaiaDR4AstromObs ----------------------------------------------------------
# One composite quantity per transit: the along-scan abscissa. It has no
# meaning off the data epochs (the scan angle is per transit), so there is no
# smooth-curve query; the panel shows model points, residuals, and the
# histogram.

plotchannels(::GaiaDR4AstromObs) = (
    PlotChannel(:along_scan, "along-scan abscissa", "mas"),
)

function residuals(obs::GaiaDR4AstromObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :astrometric_jitter) ? ctx.θ_obs.astrometric_jitter : zero(T)
    sim = simulate(obs, ctx)
    tab = obs.table
    ep = collect(Float64, tab.epoch)
    use = hasproperty(tab, :outlier_flag) ? collect(tab.outlier_flag .== 0) : trues(length(ep))
    return (;
        along_scan=(; epoch=ep, data=collect(float.(tab.centroid_pos_al)),
            model=collect(sim.along_scan),
            resid=collect(tab.centroid_pos_al .- sim.along_scan),
            σ=collect(float.(tab.centroid_pos_error_al)),
            σ_eff=collect(hypot.(tab.centroid_pos_error_al, jitter)), use),
    )
end

# ---------------------------------------------------
# PosteriorSeries
# ---------------------------------------------------

"""
    PosteriorSeries(model, chain; N=250, seed=0, ii=nothing,
                    points_per_period=30, max_points=1000)

Everything the plot layer needs from a fit, computed once and shared by
every panel:

  - a draw selection `ii` (default: `N` samples without replacement, seeded,
    so panels agree and reruns reproduce), plus the maximum-a-posteriori
    draw;
  - the rebuilt `PlanetOrbits.System` per draw (`construct_system`);
  - a dense epoch grid `ts` sized per hierarchy row
    (`PlanetOrbits.plot_epochs`), and each draw's trajectory solved over it;
  - each draw's trajectory at the union of the data epochs, and
    [`ObsContext`](@ref)s for every observation — with `method`,
    `observing_geometry` and `barycentric_lighttime` forwarded from the
    model, so curves cannot silently disagree with the likelihood.

Accepts a chain from `octofit` or any `MCMCChains.Chains` with matching
parameter names.
"""
struct PosteriorSeries{TM,TC,TTh,TSys,TTr,TTrd,TThm,TSysm,TTrdm}
    model::TM
    chain::TC
    ii::Vector{Int}
    thetas::TTh              # nested-θ per draw
    systems::TSys            # PlanetOrbits.System per draw
    ts::Vector{Float64}      # dense plotting epochs
    trajs::TTr               # solved at `ts`, per draw
    data_epochs::Vector{Float64}
    data_maps::Dict{Any,Vector{Int}}
    data_trajs::TTrd         # solved at `data_epochs`, per draw
    i_map::Int               # flat chain index of the MAP sample
    θ_map::TThm
    sys_map::TSysm
    data_traj_map::TTrdm
end
export PosteriorSeries

_solvekw(sys::System) = (; method=sys.method,
    observing_geometry=sys.observing_geometry,
    barycentric_lighttime=sys.barycentric_lighttime)

function PosteriorSeries(model::LogDensityModel, chain::MCMCChains.Chains;
                         N::Integer=250, seed::Integer=0, ii=nothing,
                         points_per_period::Integer=30, max_points::Integer=1000)
    sys = model.system
    nsamples = size(chain, 1) * size(chain, 3)
    if ii === nothing
        rng = Random.Xoshiro(seed)
        ii = sort!(_sb_sample(rng, 1:nsamples, min(N, nsamples); replace=false))
    else
        ii = collect(Int, ii)
    end

    thetas = mcmcchain2result(model, chain, ii)
    systems = map(θ -> construct_system(model, θ), thetas)

    i_map = haskey(chain, :logpost) ? argmax(vec(chain[:logpost])) :
            haskey(chain, :loglike) ? argmax(vec(chain[:loglike])) : last(ii)
    θ_map = mcmcchain2result(model, chain, i_map)
    sys_map = construct_system(model, θ_map)

    data_epochs, data_maps = epoch_plan(sys)

    # Dense grid: data span padded 1.5 %, widened to at least the
    # 35th-percentile period over draws and rows, then per-row point
    # allocation. (The quantile matches v1's behaviour; see octoplot.jl's
    # `med_period`.)
    periods = Float64[]
    for s in systems, k in 1:PlanetOrbits.nrows(s)
        p = PlanetOrbits.period(s, k)
        isfinite(p) && push!(periods, p)
    end
    if isempty(data_epochs)
        t0 = isempty(periods) ? 50000.0 : PlanetOrbits.periastron(sys_map, 1)
        span = isempty(periods) ? 365.0 : quantile(periods, 0.35)
        t_start, t_stop = t0, t0 + span
    else
        t_start, t_stop = extrema(data_epochs)
        d = t_stop - t_start
        t_start -= 0.015d
        t_stop += 0.015d
        if !isempty(periods)
            med_period = quantile(periods, 0.35)
            if t_stop - t_start < med_period
                t_extend = med_period - (t_stop - t_start)
                t_start -= t_extend / 2
                t_stop += t_extend / 2
            end
        end
        t_start = max(t_start, mjd("1900"))
        t_stop = min(t_stop, mjd("2100"))
    end
    if t_stop <= t_start
        t_stop = t_start + 1.0
    end
    ts = PlanetOrbits.plot_epochs(sys_map, t_start, t_stop; points_per_period, max_points)
    # Model curves must hit the data epochs exactly (fine structure between
    # grid points otherwise slips through — the rvpostplot lesson).
    ts = sort!(unique!(vcat(ts, data_epochs)))

    kw = _solvekw(sys)
    trajs = map(s -> orbitsolve(s, ts; kw...), systems)
    data_trajs = isempty(data_epochs) ? map(_ -> nothing, systems) :
                 map(s -> orbitsolve(s, data_epochs; kw...), systems)
    data_traj_map = isempty(data_epochs) ? nothing : orbitsolve(sys_map, data_epochs; kw...)

    return PosteriorSeries(model, chain, collect(Int, ii), thetas, systems, ts, trajs,
        data_epochs, data_maps, data_trajs, i_map, θ_map, sys_map, data_traj_map)
end

Base.length(s::PosteriorSeries) = length(s.ii)

function Base.show(io::IO, ::MIME"text/plain", s::PosteriorSeries)
    println(io, "PosteriorSeries: ", length(s.ii), " draws of \"",
        s.model.system.name, "\", ", length(s.ts), " plot epochs")
end

"""
    obscontext(series, obs; draw=nothing) -> ObsContext

The evaluation context for `obs` under posterior draw index `draw` (into
`series.ii`), or under the MAP sample when `draw === nothing`. This is the
same context the likelihood saw, so `simulate`, `residuals` and `ln_like`
all evaluate consistently.
"""
function obscontext(series::PosteriorSeries, obs; draw=nothing)
    θ, traj = draw === nothing ? (series.θ_map, series.data_traj_map) :
              (series.thetas[draw], series.data_trajs[draw])
    posys = draw === nothing ? series.sys_map : series.systems[draw]
    traj === nothing && error("this model has no observation epochs to solve at")
    key = normalizename(likelihoodname(obs))
    obsns = hasproperty(θ, :observations) ? θ.observations : (;)
    θ_obs = hasproperty(obsns, key) ? getproperty(obsns, key) : (;)
    return ObsContext(θ, θ_obs, posys, traj, series.data_maps[obs])
end
export obscontext

"""
    modelcurves(series, query) -> Vector{Vector{Float64}}

The query evaluated over the dense epoch grid `series.ts`, one vector per
posterior draw. Display scaling (e.g. rad → deg for position angles) is the
caller's job, via the channel's `scale`.
"""
modelcurves(series::PosteriorSeries, query) =
    map((s, tr) -> evalquery(_query(query), s, tr), series.systems, series.trajs)
export modelcurves

"""
    mapcurve(series, query) -> Vector{Float64}

The query evaluated over `series.ts` for the MAP draw.
"""
function mapcurve(series::PosteriorSeries, query)
    kw = _solvekw(series.model.system)
    traj = orbitsolve(series.sys_map, series.ts; kw...)
    return evalquery(_query(query), series.sys_map, traj)
end
export mapcurve

# ---------------------------------------------------
# Default panel derivation
# ---------------------------------------------------

"""
    default_sky_queries(sys) -> Vector{(ObservableQuery, Symbol)}

One (query, rowname) pair per hierarchy row: the exterior side relative to
the interior side — exactly the relationship each row parametrizes. For a
star + planets system this is each planet about the star (matching v1's
octoplot); for hierarchies it generalizes with no special cases (a moon
about its planet, an inner pair's barycentre about the outer body, …).
"""
function default_sky_queries(sys::System)
    out = Tuple{ObservableQuery,Symbol}[]
    for (owner, ext, int) in sys.rows
        t = length(ext) == 1 ? refspec(ext[1]) : BarycentreSpec{ext}()
        r = length(int) == 1 ? refspec(int[1]) : BarycentreSpec{int}()
        push!(out, (ObservableQuery(PlanetOrbits.raoff, t, r), owner))
    end
    return out
end

"""
    plottable_observations(sys)

The observations with at least one declared plot channel, in declaration
order (prior-shaped terms and unported types drop out naturally).
"""
plottable_observations(sys::System) =
    [obs for obs in sys.observations if !isempty(plotchannels(obs)) && !_isprior(obs)]

# ---------------------------------------------------
# Result type
# ---------------------------------------------------

"""
    OctoPlotResult

What [`octoplot`](@ref) returns: the `figure`, the named `axes` (a nested
NamedTuple — `res.axes.sky.sky`, `res.axes.rv.main`, `res.axes.rv.resid`, …)
for direct annotation with ordinary Makie calls, and the underlying
[`PosteriorSeries`](@ref) for further panels. Displays as its figure.
"""
struct OctoPlotResult{TF,TA,TS}
    figure::TF
    axes::TA
    series::TS
end
Base.display(res::OctoPlotResult) = display(res.figure)
Base.show(io::IO, mime::MIME"image/png", res::OctoPlotResult) = show(io, mime, res.figure)
Base.show(io::IO, mime::MIME"image/svg+xml", res::OctoPlotResult) = show(io, mime, res.figure)
function Base.show(io::IO, ::MIME"text/plain", res::OctoPlotResult)
    println(io, "OctoPlotResult with axes: ", keys(res.axes))
    show(io, MIME"text/plain"(), res.figure)
end
export OctoPlotResult

# ---------------------------------------------------
# Makie-extension stubs (methods added by OctofitterMakieExt)
# ---------------------------------------------------

@noinline _require_makie(name) = error(
    "$name requires a Makie backend; run `using CairoMakie` (or GLMakie) first.")

"""
    octoplot(model, chain; kwargs...) -> OctoPlotResult

One figure summarizing a fit: a sky panel of every orbit (when angular
observables exist) and one time-series panel per data channel, with
residuals and marginal histograms, all sharing a calendar-date epoch axis.

Requires a Makie backend to be loaded (e.g. `using CairoMakie`).

Returns a result whose fields are the `figure`, the named `axes` (for direct
annotation: `text!(res.axes.sky, ...)`), and the underlying
[`PosteriorSeries`](@ref). Keywords are documented in the extension method;
`fname="..."` saves the figure (nothing is written by default).
"""
function octoplot end
octoplot(args...; kwargs...) = _require_makie("octoplot")
export octoplot

"""
    timeseriespanel!(gridposition, series, channelgroup; kwargs...)

Generic data-vs-model time-series panel: posterior model curves, calibrated
data with errorbars, a residual strip, and a marginal residual histogram.
Requires Makie. See [`octoplot`](@ref) for the assembled default.
"""
function timeseriespanel! end
export timeseriespanel!

"""
    skypanel!(gridposition, series; kwargs...)

Sky-plane panel: phase-coloured orbit tracks for the default row queries
(or explicit ones), overlaid relative-astrometry data, star marker.
Requires Makie. See [`octoplot`](@ref).
"""
function skypanel! end
export skypanel!
