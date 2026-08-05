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
# Row signals — the part of a query attributable to one hierarchy row
# ---------------------------------------------------

# Observables linear in the target−ref separation (or its velocity), so a
# query over them telescopes exactly into per-row contributions. (Nonlinear
# ones — projectedseparation, posangle — are still foldable when a single
# row is the whole signal.)
const _LINEAR_OBSERVABLES = (
    PlanetOrbits.posx, PlanetOrbits.posy, PlanetOrbits.posz,
    PlanetOrbits.velx, PlanetOrbits.vely, PlanetOrbits.velz,
    PlanetOrbits.radvel, PlanetOrbits.raoff, PlanetOrbits.decoff,
    PlanetOrbits.pmra, PlanetOrbits.pmdec,
)

# A leaf-name view of a spec, used only for colour matching in the Makie
# extension (a planet's panels take its sky-track accent colour).
_leafnames(sys::System, ::BodyRefSpec{Name}) where {Name} = (Name,)
_leafnames(sys::System, ::BarycentreSpec{()}) = Tuple(sys.bodynames)
_leafnames(sys::System, ::BarycentreSpec{Names}) where {Names} = Names
_leafnames(sys::System, @nospecialize _) = nothing

# The body-weight vector of a resolved reference point (a body's indicator,
# a barycentre's mass weights, a photocentre's flux weights), as a tuple.
function _pointweights(posys::PlanetOrbits.System{NB}, spec) where {NB}
    p = resolveref(posys, spec)
    p isa PlanetOrbits.WeightedPoint && return Tuple(p.w)
    return ntuple(j -> j == p.idx ? 1.0 : 0.0, NB)
end

# Δ_k for every row: how much of each row's relative coordinate enters
# `target − ref`. Positions are linear in the row coordinates
# (`pos_i = Σ_k Ainv[i,k]·r_k`), so this is exact under every hierarchy
# convention — Jacobi, astrocentric, mixed — with no set bookkeeping. The
# coefficients depend on the draw's masses, so they are recomputed per draw.
function _rowcoeffs(posys::PlanetOrbits.System{NB,NR}, tspec, rspec) where {NB,NR}
    wt = _pointweights(posys, tspec)
    wr = _pointweights(posys, rspec)
    dw = wt .- wr
    return ntuple(k -> sum(dw .* Tuple(posys.Ainv[:, k])), NR)
end

"""
    RowSignal

The component of an [`ObservableQuery`](@ref) attributable to one hierarchy
row: built by [`rowsignal`](@ref), evaluated by [`evalsignal`](@ref).
Because the signal is referenced within the row itself, frame effects
(perspective acceleration, proper motion) drop out — it is purely orbital.
"""
struct RowSignal{Q,TT,TR}
    query::Q      # the observable evaluated per epoch
    k::Int
    scaled::Bool  # multiply by the draw's Δ_k (linear multi-row case)
    target::TT    # the original endpoints, for recomputing Δ_k per draw
    ref::TR
end

# The draw-dependent scale factor of the row's relative-coordinate signal.
signalcoeff(sig::RowSignal, posys) =
    sig.scaled ? _rowcoeffs(posys, sig.target, sig.ref)[sig.k] : 1.0

evalsignal(sig::RowSignal, posys, traj) =
    signalcoeff(sig, posys) .* evalquery(sig.query, posys, traj)

"""
    rowsignal(posys, query, k) -> RowSignal | nothing

The part of `query` contributed by hierarchy row `k` of a
`PlanetOrbits.System` (e.g. a `PosteriorSeries`' `sys_map`), or `nothing`
when it cannot be isolated. When row `k` is the only row moving the query,
the signal is the query itself — exact for any observable. Otherwise, for
observables linear in the separation, the contribution is
`Δ_k · f(Bc(ext_k), Bc(int_k))` with `Δ_k` from the hierarchy's topology
matrix — e.g. the row-1 part of `radvel(:A, Barycentre)` in a two-planet
system is planet b's reflex alone, with planet c's signal and any
perspective acceleration removed exactly.
"""
function rowsignal(posys::PlanetOrbits.System{NB,NR}, q::ObservableQuery, k::Integer) where {NB,NR}
    1 <= k <= NR || return nothing
    Δ = _rowcoeffs(posys, q.target, q.ref)
    scale = maximum(abs, Δ; init=0.0)
    scale > 0 || return nothing
    affecting = findall(j -> abs(Δ[j]) > 1e-12 * scale, 1:NR)
    k in affecting || return nothing
    affecting == [k] && return RowSignal(q, k, false, q.target, q.ref)
    q.func in _LINEAR_OBSERVABLES || return nothing
    names = PlanetOrbits._names(posys)
    ext = PlanetOrbits._setnames(names, posys.specs[k].ext)
    int = PlanetOrbits._setnames(names, posys.specs[k].int)
    rowq = ObservableQuery(q.func,
        length(ext) == 1 ? refspec(ext[1]) : BarycentreSpec{ext}(),
        length(int) == 1 ? refspec(int[1]) : BarycentreSpec{int}())
    return RowSignal(rowq, k, true, q.target, q.ref)
end

"""
    foldablerows(posys, query) -> Vector{Int}

The hierarchy rows a phase-fold panel of `query` can be made for (those
with an isolable [`rowsignal`](@ref)).
"""
foldablerows(posys::PlanetOrbits.System{NB,NR}, q) where {NB,NR} =
    [k for k in 1:NR if rowsignal(posys, _query(q), k) !== nothing]

"""
    foldephemeris(sig, posys, tmid; solvekw=(;)) -> (P, t0)

The fold ephemeris of a [`RowSignal`](@ref) under one posterior draw: the
row period and the epoch of phase zero, chosen in the cycle containing
`tmid`. Radial velocity keeps the v1 rvpostplot convention — phase 0 at the
signal's upward zero crossing; everything else puts periastron at phase 0.
Under the N-body propagator these are the drawn (osculating) elements — see
`phasefoldpanel!` for what a fold means there.
"""
function foldephemeris(sig::RowSignal, posys, tmid::Real; solvekw=(;))
    P = PlanetOrbits.period(posys, sig.k)
    tp = posys.rows[sig.k].tp
    t0 = tp + floor((tmid - tp) / P) * P
    if sig.query.func === PlanetOrbits.radvel
        ts = collect(range(t0, t0 + P, length=201))
        traj = orbitsolve(posys, ts; solvekw...)
        v = evalsignal(sig, posys, traj)
        for i in 2:length(v)
            if v[i-1] <= 0 <= v[i]
                return (P, (ts[i-1] + ts[i]) / 2)
            end
        end
    end
    return (P, t0)
end

"""
    foldphase(t, P, t0) -> Float64

Epoch `t` folded on period `P` about `t0`, in [-0.5, 0.5).
"""
foldphase(t, P, t0) = mod((t - t0) / P + 0.5, 1.0) - 0.5

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

"""
    defaultpanels(obs) -> Tuple

The escape hatch for observations that are not epoch series (an HGCA row, a
catalog solution): bespoke panels instead of the generic time-series ones.
Return `()` (the default) to use the generic panels derived from
[`plotchannels`](@ref); return `(name => build, …)` pairs to opt out —
`octoplot` then calls each `build(gridposition, series)` in the panel stack
and merges its returned NamedTuple of axes under `name`. Include a
`timeaxes` key (a tuple of epoch axes) in that NamedTuple to have them
linked with the figure's shared time axis.

A bespoke panel special-cases the *drawing*, not the plumbing: it still
receives the shared [`PosteriorSeries`](@ref) (same draws as every other
panel) and should source its residuals/whitening from
[`residuals`](@ref)/`ln_like` machinery, never re-derive it.
"""
defaultpanels(::AbstractObs) = ()
export defaultpanels

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
    sim = simulate(obs, ctx)               # includes the offset and the trend
    ep = collect(Float64, obs.table.epoch)
    # Subtract the trend as well as the offset. The `radvel` query the panel
    # draws its curve from knows nothing about either, so a trended model
    # whose points kept the trend would not overlay its own curve. v1's
    # `rvpostplot` subtracted both; `trend_function` now lives on the core
    # type, so both belong here.
    trend = obs.trend_function.(Ref(θ_obs), obs.table.epoch)
    model_pure = sim.rv_model .- offset .- trend
    data_cal = obs.table.rv .- offset .- trend
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

# --- HipparcosIADObs -----------------------------------------------------------
# The same shape as the DR4 channel: one along-scan abscissa per transit, with
# no meaning off the data epochs (the scan angle is per transit), so no query.

plotchannels(::HipparcosIADObs) = (
    PlotChannel(:along_scan, "Hipparcos abscissa residual", "mas"),
)

# The along-scan channel of a Hipparcos IAD table, shared by `HipparcosIADObs`
# and `G23HObs` (they read the same transits).
#
# The datum plotted is `res`, the catalog's own abscissa residual, not the
# absolute abscissa `proj_meas_alongscan` the likelihood compares against.
# They differ by the catalog five-parameter path projected on scan, which for
# a high-proper-motion star is *arcseconds* — it would set the axis range and
# bury the milliarcsecond signal the fit is about. Subtracting it from the
# data and the model alike is the same calibration an RV instrument offset
# gets, and it leaves the residual (and therefore the likelihood) untouched.
function _hip_alongscan(tab, resid, σ_infl, jitter)
    res = collect(Float64, tab.res)
    r = collect(Float64, resid)
    # The per-transit σ `ln_like` uses: the renormalized formal error scaled by
    # the BINARYS first-harmonic inflation factor, which grows where the binary
    # modulation reduces the signal amplitude.
    σ = collect(Float64, tab.sres_renorm .* σ_infl)
    return (; epoch=collect(Float64, tab.epoch), data=res, model=res .- r, resid=r,
        σ, σ_eff=collect(hypot.(σ, jitter)),
        # `reject` is van Leeuwen's own per-transit rejection flag, and
        # `ln_like` skips those rows.
        use=collect(Bool, .!tab.reject))
end

function residuals(obs::HipparcosIADObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :hip_iad_jitter) ? ctx.θ_obs.hip_iad_jitter : zero(T)
    sim = simulate(obs, ctx)                        # `resid` is measured − model
    return (; along_scan=_hip_alongscan(obs.table, sim.resid, sim.σ_inflation, jitter))
end

# --- G23HObs -------------------------------------------------------------------
#
# A composite catalog observation (design §5.1): five *time-averaged* proper
# motions per axis rather than an epoch series, plus a UEVA datum and — when
# `:iad_hip` is kept — the Hipparcos per-transit abscissae.
#
# What it exposes, and why that shape:
#
#   * `:pmra` / `:pmdec`. Each catalog channel is a proper motion of the same
#     source over a different window, so all five belong on one axis against
#     one curve — which is what v1's `hgcaplot`, `pmaplot` and `absastromplot`
#     all drew in their top two panels. The curve is the *reflex* proper
#     motion `pmra(host, ref)`, and the data are calibrated onto it by
#     removing the reference frame's own proper motion, exactly as an RV
#     instrument offset is removed. Each point carries `epoch_lo`/`epoch_hi`:
#     a PM is an average over its mission window, and the generic panel draws
#     that window as a horizontal bar (v1 `absastromplot`'s one genuinely
#     load-bearing idiom).
#   * `:along_scan_hip`, when the Hipparcos abscissae are included — the same
#     channel `HipparcosIADObs` exposes, from the same transits.
#
# Deliberately *not* exposed as channels: `:ueva_dr3` (an astrometric
# excess-noise variance, not a proper motion — it shares no axis with anything
# and a one-point panel says less than the corner plot does) and `:rv_dr3` (a
# radial-velocity *variability* statistic, likewise). Both still enter
# `ln_like`; they simply have no data-vs-model series to draw.

const _G23H_PM_KINDS = (
    pmra=(:ra_hip, :ra_hg, :ra_dr2, :ra_dr32, :ra_dr3),
    pmdec=(:dec_hip, :dec_hg, :dec_dr2, :dec_dr32, :dec_dr3),
)

function plotchannels(obs::G23HObs)
    kinds = obs.table.kind
    chs = Any[]
    if any(k -> k ∈ kinds, _G23H_PM_KINDS.pmra)
        push!(chs, PlotChannel(:pmra, "μα*", "mas/yr";
            query=ObservableQuery(PlanetOrbits.pmra, obs.host, obs.ref)))
    end
    if any(k -> k ∈ kinds, _G23H_PM_KINDS.pmdec)
        push!(chs, PlotChannel(:pmdec, "μδ", "mas/yr";
            query=ObservableQuery(PlanetOrbits.pmdec, obs.host, obs.ref)))
    end
    :iad_hip ∈ kinds &&
        push!(chs, PlotChannel(:along_scan_hip, "Hipparcos abscissa residual", "mas"))
    return Tuple(chs)
end

# One axis' worth of catalog proper motions, in table (i.e. epoch) order.
function _g23h_pmseries(obs::G23HObs, mom, kinds, offset)
    lut = Dict{Symbol,Int}(l => k for (k, l) in enumerate(mom.labels))
    ep = Float64[]; lo = Float64[]; hi = Float64[]
    data = Float64[]; model = Float64[]; σ = Float64[]
    for k in kinds
        j = get(lut, k, 0)
        j == 0 && continue
        row = findfirst(==(k), obs.table.kind)
        row === nothing && continue
        push!(ep, obs.table.epoch[row])
        push!(lo, obs.table.start_epoch[row])
        push!(hi, obs.table.stop_epoch[row])
        push!(data, mom.catalog[j] - offset)
        push!(model, mom.model[j] - offset)
        push!(σ, mom.sigma[j])
    end
    # σ is the *marginal* σ of the coupled block. The likelihood scores the
    # whole 11-vector against its full covariance (the DR2/DR3 proper motions
    # are correlated through `rho_dr2_dr3`), so a per-point residual over σ is
    # a per-channel diagnostic, not the quantity that enters `ln_like`. It is
    # the same convention v1 used, and the per-channel jitter and the DR3
    # deflation are already in it — so there is no second, looser σ to report:
    # `σ_eff == σ`.
    return (; epoch=ep, epoch_lo=lo, epoch_hi=hi, data, model,
        resid=data .- model, σ, σ_eff=copy(σ), use=trues(length(ep)))
end

function residuals(obs::G23HObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    kinds = obs.table.kind
    out = NamedTuple()

    if any(k -> k ∈ kinds, (_G23H_PM_KINDS.pmra..., _G23H_PM_KINDS.pmdec...))
        mom = _g23h_catalog_moments(obs, ctx)
        sim = simulate(obs, ctx)
        pmra_sys, pmdec_sys = _g23h_pm(ctx.θ_system, ctx.θ_obs, T)
        # The likelihood expresses every channel against the *reference point's*
        # proper motion and then shifts the whole frame by −Δpm_dr3 so that the
        # frame refers to the primary rather than the barycentre (see
        # `_g23h_simulate!`'s closing comment). Removing that constant from both
        # the catalog and the model values is what puts them on the pure
        # `pmra(host, ref)` curve the panel draws — the proper-motion analogue
        # of subtracting an RV instrument's zero point.
        off_ra = Float64(pmra_sys - sim.Δpmra_dr3)
        off_dec = Float64(pmdec_sys - sim.Δpmdec_dr3)
        any(k -> k ∈ kinds, _G23H_PM_KINDS.pmra) && (out = merge(out,
            (; pmra=_g23h_pmseries(obs, mom, _G23H_PM_KINDS.pmra, off_ra))))
        any(k -> k ∈ kinds, _G23H_PM_KINDS.pmdec) && (out = merge(out,
            (; pmdec=_g23h_pmseries(obs, mom, _G23H_PM_KINDS.pmdec, off_dec))))
    end

    if :iad_hip ∈ kinds
        sim = simulate(obs, ctx)
        jitter = hasproperty(ctx.θ_obs, :hip_iad_jitter) ? ctx.θ_obs.hip_iad_jitter : zero(T)
        out = merge(out, (; along_scan_hip=_hip_alongscan(
            obs.hip_table, sim.iad_resid_signed, sim.σ_inflation_hip, jitter)))
    end
    return out
end

# --- PhotometryObs -------------------------------------------------------------
#
# Photometry has no epoch axis at all (`epochs(obs)` is empty by design) and no
# `PlanetOrbits` observable behind it — the forward model is one body variable,
# constant across the table. A time-series panel is therefore the wrong shape,
# so the type declares its channel (which is what makes `Octofitter.residuals`,
# goodness-of-fit tables and the "no plot channels" audit work) and opts out of
# the generic panels with a small bespoke one.

plotchannels(obs::PhotometryObs) = (
    PlotChannel(:phot, "flux (" * string(_photvar(obs)) * ")", ""),
)

function residuals(obs::PhotometryObs, ctx::ObsContext)
    flux = Float64(simulate(obs, ctx).phot_model)
    n = length(obs.table.phot)
    data = collect(Float64, obs.table.phot)
    σ = collect(Float64, obs.table.σ_phot)
    return (;
        # There is no epoch: the row index stands in, and the bespoke panel
        # puts it on the x axis. (A table may carry an `:epoch` column for the
        # user's own bookkeeping, but the model ignores it — see `epochs`.)
        phot=(; epoch=collect(1.0:n), data, model=fill(flux, n),
            resid=data .- flux, σ, σ_eff=copy(σ), use=trues(n)),
    )
end

defaultpanels(obs::PhotometryObs) =
    (:phot => (gp, series) -> photometrypanel!(gp, series, obs),)

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
order (prior-shaped terms and types with no channels yet drop out naturally).
"""
plottable_observations(sys::System) =
    [obs for obs in sys.observations if !isempty(plotchannels(obs)) && !_isprior(obs)]

"""
    unplottable_observations(sys) -> Vector{String}

Observations that carry data but declare no [`plotchannels`](@ref), described
for a log message.

Worth announcing rather than dropping silently: a proper-motion-anomaly or
Hipparcos-only figure with no data on it looks like a *result* — a model with
nothing to constrain it — rather than like a missing feature.

Every observation type Octofitter ships now declares channels, so this is
normally empty; it fires for user-defined observation types that have not
implemented the protocol yet.
"""
unplottable_observations(sys::System) =
    String["$(likelihoodname(o)) ($(nameof(typeof(o))))"
           for o in sys.observations if !_isprior(o) && isempty(plotchannels(o))]

function _warn_unplottable(sys::System)
    missing_ = unplottable_observations(sys)
    isempty(missing_) && return nothing
    @info("Some observations declare no plot channels, so their data are not overlaid " *
          "— only the modelled orbits are drawn: $(join(missing_, ", ")). " *
          "See `Octofitter.plotchannels`.")
    return nothing
end

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
`fname="..."` saves the figure (nothing is written by default), and
`channels=` restricts it to some of the data — an observable function, a
channel or observable name, or a collection of either. [`rvpostplot`](@ref)
is `channels=radvel` with the sky panel off.
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

"""
    octocorner(model, chains...; small=false, includecols=[], excludecols=[],
               labels=Dict(), truth=(), fname=nothing, kwargs...)

Corner (pair) plot of the fit parameters. Labels, units, and radian→degree
conversions come from `PlanetOrbits.paraminfo` — the same resolver table
the axis labels use — keyed by the flat `<owner>_<var>` chain naming, so
custom parameters simply show their column name. `small=true` keeps only
each body's `a`, `e`, `i`, `mass`. `UniformCircular` helper pairs, fixed
values, and `tp` duplicated by a sampled `θ`/`M0` are dropped;
`includecols` forces columns in, `excludecols` out. Extra keywords pass
through to `PairPlots.pairplot`. Nothing is written unless `fname` is set.

Requires both a Makie backend and PairPlots to be loaded.
"""
function octocorner end
octocorner(args...; kwargs...) = _require_makie("octocorner (also load PairPlots)")
export octocorner

"""
    dotplot(model, chain; mode=:separation, epoch=nothing, kwargs...)

Mass against separation (or period) for every body in the fit, coloured by
eccentricity, with marginal histograms. A posterior summary — no data, no
observations — so it works for any model whose bodies have a sampled `mass`
and `e`.

`mode=:separation` uses each draw's semi-major axis, or the instantaneous
3-D separation at `epoch=` when one is given; `mode=:period` uses the orbital
period. Masses are M⊙ (v1 plotted Mⱼᵤₚ; v2 has one mass unit throughout).
Requires a Makie backend.
"""
function dotplot end
dotplot(args...; kwargs...) = _require_makie("dotplot")
export dotplot

"""
    gaiastarplot(model, chain, sample_idx=MAP; kwargs...)

The host's reflex orbit in the Gaia frame for one posterior draw, with each
transit's along-scan residual re-projected into the sky plane along its own
scan angle and drawn as a segment through the modelled track. This is the
"is there a wobble, and does the orbit fit it?" picture for a
[`GaiaDR4AstromObs`](@ref); the along-scan-versus-time half of it is the
generic panel [`octoplot`](@ref) already draws. Requires a Makie backend.
"""
function gaiastarplot end
gaiastarplot(args...; kwargs...) = _require_makie("gaiastarplot")
export gaiastarplot

"""
    gaiatimeplot(model, chain; kwargs...)

Along-scan abscissa against time for a [`GaiaDR4AstromObs`](@ref): the
posterior cloud of modelled abscissae over the measurements, with a per-epoch
boxplot of the residuals against the quoted formal errors below.

This is the same data as `octoplot`'s generic `:along_scan` panel, drawn in
v1's per-epoch-boxplot idiom, which answers a different question: not "are
the residuals normal" but "at which epochs is the posterior spread larger
than the measurement error". Requires a Makie backend.
"""
function gaiatimeplot end
gaiatimeplot(args...; kwargs...) = _require_makie("gaiatimeplot")
export gaiatimeplot

"""
    skytrackplot(model, chain, sample_idx=MAP; ra=nothing, dec=nothing,
                 gaia_id=nothing, ts=nothing, keplerian_mult=1, kwargs...)

The star's whole path on the sky for one draw: parallactic loops, proper
motion, and the orbital wobble superimposed — the picture of *why* the wobble
is hard to extract. `keplerian_mult` exaggerates the orbital term.

The parallax ellipse needs a sky direction to project onto. It is taken from
the system's own `ra`/`dec` frame variables when the model declares an
absolute frame; otherwise give `ra=`/`dec=` in degrees, or `gaia_id=`, which
reads them from the published solution via
[`gaia_dr3_solution`](@ref). Requires a Makie backend, and the DE440
ephemeris data dependency for the Earth's position.
"""
function skytrackplot end
skytrackplot(args...; kwargs...) = _require_makie("skytrackplot")
export skytrackplot

"""
    hipparcosplot(model, chain, sample_idx=MAP; kwargs...)

Hipparcos intermediate astrometry for one draw, in its own geometry: the
catalog's five-parameter sky path, the modelled path with the companion's
perturbation, each transit's abscissa line, and the perpendicular residual
and formal error drawn against it — plus a residual-versus-time strip.

Works with a [`HipparcosIADObs`](@ref) or with a [`G23HObs`](@ref) that keeps
its `:iad_hip` channel. Requires a Makie backend.
"""
function hipparcosplot end
hipparcosplot(args...; kwargs...) = _require_makie("hipparcosplot")
export hipparcosplot

"""
    photometrypanel!(gridposition, series, obs; kwargs...)

Bespoke panel for a [`PhotometryObs`](@ref): the posterior of the modelled
flux as a band, and the measurements with their errorbars. Photometry has no
epoch axis, so the x axis is the measurement index. Requires Makie.
"""
function photometrypanel! end
export photometrypanel!

"""
    likemappanel!(gridposition, series, obs; kwargs...)

Bespoke panel for a `LogLikelihoodMapObs`: per epoch, how far below that
epoch's map maximum the modelled position falls, over the posterior draws.
Requires Makie (and OctofitterImages). See [`defaultpanels`](@ref).
"""
function likemappanel! end
export likemappanel!

"""
    rvpostplot(model, chain; kwargs...) -> OctoPlotResult

The radial-velocity summary figure: one time-series panel with a residual
strip and marginal histogram, plus one phase-folded panel per hierarchy row
that moves the star. Sugar over [`octoplot`](@ref) restricted to the model's
radial-velocity channels — the same panels, without the sky panel and without
any non-RV data.

Requires a Makie backend. `rvpostplot_animated` records the same figure over
successive single-draw slices of the chain.
"""
function rvpostplot end
rvpostplot(args...; kwargs...) = _require_makie("rvpostplot")
function rvpostplot_animated end
rvpostplot_animated(args...; kwargs...) = _require_makie("rvpostplot_animated")
export rvpostplot, rvpostplot_animated

"""
    phasefoldpanel!(gridposition, series, entries; row, kwargs...)

Data-vs-model panel folded on hierarchy row `row`'s orbital phase: the
row's isolated signal ([`rowsignal`](@ref)) per posterior draw, calibrated
data with the other rows' signals removed, noise-weighted binned means, and
a phase-folded residual strip. Requires Makie. See [`octoplot`](@ref).
"""
function phasefoldpanel! end
export phasefoldpanel!
