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

# The kinematic half of `radvel`, in m/s.
#
# `radvel` is the *spectroscopic* velocity: it adds the Einstein term, which is
# quadratic in velocity and 1/r in separation, and therefore does **not**
# telescope into per-row contributions the way the kinematic part does. A
# multi-row fold decomposes the kinematic signal; the Einstein term stays
# common-mode and lands in the residual beside the instrument offset, which is
# where a term no single row owns belongs. (For planetary companions its
# *varying* part is sub-cm/s; see PlanetOrbits' "Precision opt-outs" page for
# the cases where it is not, all of which are single-row.)
#
# The single-row case is untouched: there the signal is the query itself, so
# it keeps `radvel` and stays exact.
_kinrv(sol, t, r) = PlanetOrbits.velz(sol, t, r) *
                    PlanetOrbits.au2m / PlanetOrbits.year2sec_julian

_rowfunc(f) = f
_rowfunc(::typeof(PlanetOrbits.radvel)) = _kinrv

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
row: built by `rowsignal`, evaluated by `evalsignal`.
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
    rowq = ObservableQuery(_rowfunc(q.func),
        length(ext) == 1 ? refspec(ext[1]) : BarycentreSpec{ext}(),
        length(int) == 1 ? refspec(int[1]) : BarycentreSpec{int}())
    return RowSignal(rowq, k, true, q.target, q.ref)
end

"""
    foldablerows(posys, query) -> Vector{Int}

The hierarchy rows a phase-fold panel of `query` can be made for (those
with an isolable `rowsignal`).
"""
foldablerows(posys::PlanetOrbits.System{NB,NR}, q) where {NB,NR} =
    [k for k in 1:NR if rowsignal(posys, _query(q), k) !== nothing]

"""
    foldephemeris(sig, posys, tmid; solvekw=(;)) -> (P, t0)

The fold ephemeris of a [`RowSignal`](@ref) under one posterior draw: the
row period and the epoch of phase zero, chosen in the cycle containing
`tmid`. Radial velocity keeps the v8 rvpostplot convention — phase 0 at the
signal's upward zero crossing; everything else puts periastron at phase 0.
Under the N-body propagator these are the drawn (osculating) elements — see
`phasefoldpanel!` for what a fold means there.
"""
function foldephemeris(sig::RowSignal, posys, tmid::Real; solvekw=(;))
    P = PlanetOrbits.period(posys, sig.k)
    tp = posys.rows[sig.k].tp
    t0 = tp + floor((tmid - tp) / P) * P
    if sig.query.func === PlanetOrbits.radvel || sig.query.func === _kinrv
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
  - `derived` — this channel is a deterministic re-expression of the
    observation's *native* measurement rather than a measurement in its own
    right: a sep/pa table's (ra, dec), say. The likelihood scores the native
    channels; a derived one exists so that a mixed dataset can put every
    point on every panel. `octoplot` only draws one when some other
    observation declares the same channel natively, or when `channels=` asks
    for it by name, and a consumer that wants the measurements themselves
    (a goodness-of-fit table) should filter these out.
"""
struct PlotChannel{TQ}
    name::Symbol
    label::String
    unit::String
    scale::Float64
    wrap::Union{Nothing,Float64}
    query::TQ
    derived::Bool
end
PlotChannel(name, label, unit; scale=1.0, wrap=nothing, query=nothing, derived=false) =
    PlotChannel{typeof(query)}(name, label, unit, scale, wrap, query, derived)
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
    plotobs(obs) -> AbstractObs

The observation whose *measurements* a plot draws for `obs` — itself for
every observation that carries its own data, and the wrapped observation for
a wrapper like [`ObsPriorONeil2019`](@ref).

A wrapper delegates the whole plotting protocol ([`plotchannels`](@ref),
[`residuals`](@ref), [`sharepanel`](@ref), [`datacalibration`](@ref),
[`noisemodel`](@ref)) to what it wraps, so the generic panels need no
unwrapping. This exists for the few places that still ask *what kind of
observation is this* — the sky panel's relative-astrometry overlay, the
bespoke `hipparcosplot`/`gaiastarplot` selectors — where a `isa` test on the
wrapper would silently drop the dataset. Test `plotobs(obs) isa …`, never
`obs isa …`.

The wrapper itself, not the result of this, stays the object handed to
[`obscontext`](@ref) and looked up in a series' `data_maps`: the fitted
`jitter`/`platescale`/`northangle` are registered under the *wrapper's*
`likelihoodname`, so unwrapping before building a context would silently
fall back to the defaults.
"""
plotobs(obs::AbstractObs) = obs
export plotobs

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

Two optional keys carry per-channel extras when the observation has them:
`epoch_lo`/`epoch_hi` bound an averaging window (a catalog proper motion is
not measured *at* an epoch), and `gp_mean`/`gp_var` are this draw's
correlated-noise prediction at the data epochs — see [`noisemodel`](@ref).
`resid` is always the plain `data − model`; a consumer that wants the
residual the fit is left with subtracts `gp_mean` itself.

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

"""
    sharepanel(obs) -> Bool

May this observation's data share a panel with another observation measuring
the same quantity?

Only one draw's parameters can calibrate a panel, so several instruments on
one axis means every instrument but that draw's is drawn slightly wrong. The
question is whether "slightly" is small enough to be worth the far more
readable figure, and the answer is a property of the observation type:

  - **`false` (the default).** Radial velocity: an instrument zero point is
    an unconstrained free parameter of order the data range, so the
    calibrated series moves visibly from draw to draw. Each instrument gets
    its own panel, its data are drawn *uncalibrated*, and each draw's model
    curve carries that draw's own offset and trend
    ([`datacalibration`](@ref)) — which is the only way many draws and one
    dataset can appear together without misrepresenting either.
  - **`true`.** Relative astrometry: `platescale` and `northangle` are
    calibration constants pinned to within a fraction of a percent, so every
    instrument's points land in the same place under any draw. Merging them
    onto one sky, separation and position-angle panel — calibrated by the
    maximum-posterior draw — is what makes the figure legible, and it is what
    Octofitter has always done.

[`rvplot`](@ref) is the deliberate exception on the other side: it puts every
RV instrument on one panel *because* it shows a single draw, so there is no
inconsistency to hide.
"""
sharepanel(::AbstractObs) = false
export sharepanel

"""
    datacalibration(obs, ch::PlotChannel, ctx::ObsContext, epochs) -> Vector | nothing

The additive term that carries channel `ch`'s pure-observable model curve
into this observation's *raw measured* frame under the parameters in `ctx` —
an RV instrument's zero point plus its trend, evaluated at `epochs`.
`nothing` (the default) means the channel needs no additive calibration.

This is the inverse of what [`residuals`](@ref) applies to the data, and the
two must agree: `residuals` subtracts it so the points lie on a pure model
curve, and a panel drawing uncalibrated data adds it to the curve instead.
Both spellings come from here, so they cannot drift apart.
"""
datacalibration(::AbstractObs, ::PlotChannel, ::ObsContext, epochs) = nothing
export datacalibration

# The calibration as a plain vector, zeros when the observation declares none.
_calibration(obs, ch, ctx, epochs) =
    something(datacalibration(obs, ch, ctx, epochs), zeros(length(epochs)))

"""
    noisemodel(obs, ctx::ObsContext, epochs) -> (; mean, var) | nothing

The observation's correlated-noise model — a Gaussian process fitted to the
residuals — conditioned on this draw's residuals and evaluated at `epochs`;
`nothing` (the default) when the observation has none.

Plots use it twice: the band it predicts is drawn around the model curve, and
the residual strip subtracts its mean and adds its variance to `σ_eff`. That
second use is what makes a whitened residual meaningful for a GP fit at all —
without it the strip shows exactly the correlated structure the GP was fitted
to explain, and the z-scores are not standard normal even for a perfect fit.
"""
noisemodel(::AbstractObs, ::ObsContext, epochs) = nothing
export noisemodel

# --- RelAstromObs -------------------------------------------------------------

# Relative astrometry is measured in one of two parametrizations and is
# meaningful in both: (sep, pa) and (Δα*, Δδ) are related by a rotation, with
# no free parameter and no model in between. A survey that switched
# conventions mid-campaign is one dataset, so all four channels are declared
# whatever the table holds, and the two the table does not carry are marked
# `derived` — the likelihood still scores only the native pair.
_relastrom_channels(obs, native_seppa::Bool) = (
    PlotChannel(:sep, "separation", "mas"; derived=!native_seppa,
        query=ObservableQuery(PlanetOrbits.projectedseparation, obs.target, obs.ref)),
    PlotChannel(:pa, "position angle", "°"; scale=rad2deg(1.0), wrap=360.0,
        derived=!native_seppa,
        query=ObservableQuery(PlanetOrbits.posangle, obs.target, obs.ref)),
    PlotChannel(:raoff, "Δα*", "mas"; derived=native_seppa,
        query=ObservableQuery(PlanetOrbits.raoff, obs.target, obs.ref)),
    PlotChannel(:decoff, "Δδ", "mas"; derived=native_seppa,
        query=ObservableQuery(PlanetOrbits.decoff, obs.target, obs.ref)),
)

plotchannels(obs::RelAstromObs) = _relastrom_channels(obs, hasproperty(obs.table, :sep))

sharepanel(::RelAstromObs) = true

# (sep, pa) → (Δα*, Δδ) and back, with first-order error propagation. The
# rotation is exact and parameter-free, so a point's *position* is the same
# measurement either way; only the uncertainty needs care, and the
# uncorrelated-in-the-native-basis assumption below is the same one the
# likelihood makes.
function _seppa_to_radec(sep, pa, σ_sep, σ_pa)
    s, c = sin.(pa), cos.(pa)
    return (ra=sep .* s, dec=sep .* c,
        σ_ra=hypot.(σ_sep .* s, sep .* σ_pa .* c),
        σ_dec=hypot.(σ_sep .* c, sep .* σ_pa .* s))
end

function _radec_to_seppa(ra, dec, σ_ra, σ_dec)
    sep = hypot.(ra, dec)
    # A point at the origin has no defined position angle; guard the
    # denominators rather than emit a NaN into an axis limit.
    d = max.(sep, eps(one(eltype(sep))))
    return (sep=sep, pa=atan.(ra, dec),
        σ_sep=hypot.(ra .* σ_ra, dec .* σ_dec) ./ d,
        σ_pa=hypot.(dec .* σ_ra, ra .* σ_dec) ./ d .^ 2)
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

    r2d = rad2deg(1.0)
    sep_model = hypot.(sim.ra_model, sim.dec_model)
    pa_model = atan.(sim.ra_model, sim.dec_model)

    if hasproperty(obs.table, :sep)
        sep_data = collect(float.(obs.table.sep .* platescale))
        pa_data = collect(float.(obs.table.pa .+ northangle))
        σ_sep = collect(float.(obs.table.σ_sep))
        σ_pa = collect(float.(obs.table.σ_pa))
        σe_sep = collect(hypot.(σ_sep, jitter))
        # The jitter is an angular scatter in mas; on the position-angle axis
        # it is that scatter divided by the separation, exactly as `ln_like`
        # spells it.
        σe_pa = collect(hypot.(σ_pa, jitter ./ max.(sep_model, one(T))))
        p = _seppa_to_radec(sep_data, pa_data, σ_sep, σ_pa)
        pe = _seppa_to_radec(sep_data, pa_data, σe_sep, σe_pa)
        ra_data, dec_data = p.ra, p.dec
        σ_ra, σ_dec, σe_ra, σe_dec = p.σ_ra, p.σ_dec, pe.σ_ra, pe.σ_dec
    else
        # Mirrors ln_like's ra/dec branch exactly: its pa_dat is measured from
        # the RA axis (atan(dec, ra), *not* the position-angle convention), so
        # cos pairs with RA and sin with Dec.
        pa_dat = atan.(obs.table.dec, obs.table.ra) .+ northangle
        sep_dat = hypot.(obs.table.dec, obs.table.ra) .* platescale
        ra_data = collect(float.(sep_dat .* cos.(pa_dat)))
        dec_data = collect(float.(sep_dat .* sin.(pa_dat)))
        σ_ra = collect(float.(obs.table.σ_ra))
        σ_dec = collect(float.(obs.table.σ_dec))
        σe_ra = collect(hypot.(σ_ra, jitter))
        σe_dec = collect(hypot.(σ_dec, jitter))
        p = _radec_to_seppa(ra_data, dec_data, σ_ra, σ_dec)
        pe = _radec_to_seppa(ra_data, dec_data, σe_ra, σe_dec)
        sep_data, pa_data = p.sep, p.pa
        σ_sep, σ_pa, σe_sep, σe_pa = p.σ_sep, p.σ_pa, pe.σ_sep, pe.σ_pa
    end

    # Wrapped PA difference, exactly as ln_like computes it.
    pa_resid = @. mod(pa_data - pa_model + π, 2π) - π
    return (;
        sep=(; epoch=ep, data=sep_data, model=collect(sep_model),
            resid=collect(sep_data .- sep_model),
            σ=σ_sep, σ_eff=σe_sep, use),
        pa=(; epoch=ep, data=collect(rem2pi.(pa_data, RoundDown) .* r2d),
            model=collect(rem2pi.(pa_model, RoundDown) .* r2d),
            resid=collect(pa_resid .* r2d),
            σ=σ_pa .* r2d, σ_eff=σe_pa .* r2d, use),
        raoff=(; epoch=ep, data=collect(ra_data), model=collect(sim.ra_model),
            resid=collect(ra_data .- sim.ra_model),
            σ=σ_ra, σ_eff=σe_ra, use),
        decoff=(; epoch=ep, data=collect(dec_data), model=collect(sim.dec_model),
            resid=collect(dec_data .- sim.dec_model),
            σ=σ_dec, σ_eff=σe_dec, use),
    )
end

# --- RadialVelocityObs ---------------------------------------------------------

plotchannels(obs::RadialVelocityObs) = (
    PlotChannel(:rv, "radial velocity", "m/s";
        query=ObservableQuery(PlanetOrbits.radvel, obs.target, obs.ref)),
)

# The zero point and the trend, which is everything standing between the
# `radvel` query and the number the spectrograph reported.
function datacalibration(obs::RadialVelocityObs, ch::PlotChannel, ctx::ObsContext, epochs)
    ch.name === :rv || return nothing
    T = _system_number_type(ctx.θ_system)
    θ_obs = ctx.θ_obs
    offset = hasproperty(θ_obs, :offset) ? θ_obs.offset : zero(T)
    return [Float64(offset + obs.trend_function(θ_obs, t)) for t in epochs]
end

function noisemodel(obs::RadialVelocityObs, ctx::ObsContext, epochs)
    isnothing(obs.gaussian_process) && return nothing
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    sim = simulate(obs, ctx)
    resid = collect(Float64, obs.table.rv .- sim.rv_model)
    σ² = collect(Float64, obs.table.σ_rv .^ 2 .+ jitter^2)
    ep = collect(Float64, obs.table.epoch)
    # The same three hooks `ln_like` goes through, so the band a plot draws is
    # the noise model the fit actually used.
    fx = gp_condition(obs.gaussian_process(ctx.θ_obs), ep, σ²)
    m, v = gp_predict(fx, resid, collect(Float64, epochs))
    return (; mean=collect(Float64, m), var=max.(collect(Float64, v), 0.0))
end

function residuals(obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    θ_obs = ctx.θ_obs
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)
    sim = simulate(obs, ctx)               # includes the offset and the trend
    ep = collect(Float64, obs.table.epoch)
    ch = first(plotchannels(obs))
    # Subtract the trend as well as the offset. The `radvel` query the panel
    # draws its curve from knows nothing about either, so a trended model
    # whose points kept the trend would not overlay its own curve. v1's
    # `rvpostplot` subtracted both; `trend_function` now lives on the core
    # type, so both belong here — and `datacalibration` is the one place that
    # spells the sum, so an uncalibrated panel adds back exactly this.
    cal = _calibration(obs, ch, ctx, ep)
    out = (; epoch=ep, data=collect(obs.table.rv .- cal),
        model=collect(sim.rv_model .- cal),
        resid=collect(obs.table.rv .- sim.rv_model),
        σ=collect(float.(obs.table.σ_rv)),
        σ_eff=collect(hypot.(obs.table.σ_rv, jitter)),
        use=trues(length(ep)))
    gp = noisemodel(obs, ctx, ep)
    gp === nothing || (out = merge(out, (; gp_mean=gp.mean, gp_var=gp.var)))
    return (; rv=out)
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

# --- ObsPriorONeil2019 ---------------------------------------------------------
#
# The wrapper delegates the whole protocol to what it wraps, exactly as it
# already delegates `simulate`, `epochs` and the correction hooks. The
# Jacobian reweights the prior; it does not change what the instrument
# measured, so every plotting question about the wrapper is a question about
# the wrapped observation.
#
# Without these it answered the `AbstractObs` defaults — `plotchannels` is
# `()` — and `plottable_observations` dropped it: `octoplot` on a model whose
# relative astrometry was wrapped drew the orbit tracks with no data points on
# them and no separation/position-angle panels at all. (`_warn_unplottable`
# did say so, but an @info is not what a missing dataset should look like.)
#
# The *wrapper* stays the object the plot layer passes around: it is the one
# in `sys.observations`, so it is what `epoch_plan` keys its row maps by and
# what `θ.observations` registers the fitted jitter/platescale/northangle
# under. Every method here takes a context built for the wrapper and hands it
# to the wrapped observation's math, which is sound because the two share
# `table`, `priors`, `derived` and therefore `epochs`. Substituting the inner
# observation into `obscontext` instead would rediscover the v8.3.0 bug from
# the other side: `θ_obs` silently `(;)`, calibration silently back to its
# defaults. That is what `plotobs` is for — `isa` tests, and nothing else.
plotobs(obs::ObsPriorONeil2019) = plotobs(obs.wrapped_like)
plotchannels(obs::ObsPriorONeil2019) = plotchannels(obs.wrapped_like)
residuals(obs::ObsPriorONeil2019, ctx::ObsContext) = residuals(obs.wrapped_like, ctx)
sharepanel(obs::ObsPriorONeil2019) = sharepanel(obs.wrapped_like)
datacalibration(obs::ObsPriorONeil2019, ch::PlotChannel, ctx::ObsContext, epochs) =
    datacalibration(obs.wrapped_like, ch, ctx, epochs)
noisemodel(obs::ObsPriorONeil2019, ctx::ObsContext, epochs) =
    noisemodel(obs.wrapped_like, ctx, epochs)

# `defaultpanels` is deliberately *not* delegated. Its contract is a tuple of
# `build(gridposition, series)` closures, and every implementation closes over
# the observation it was asked about — `PhotometryObs`' panel over that
# `PhotometryObs`, `LogLikelihoodMapObs`' over that map. Forwarding the call
# would hand back closures holding the *inner* observation, which the panel
# then looks up in `series.data_maps` (keyed by the wrapper): a `KeyError`, or
# a context built under the wrong name. A wrapped observation with a bespoke
# panel therefore falls back to the generic time-series panels its own
# `plotchannels`/`residuals` already support — a legible figure rather than a
# wrong one. Nothing Octofitter ships reaches that: the two types with bespoke
# panels are photometry (no epochs at all, so the Jacobian has nothing to sum
# over) and likelihood maps.

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
    return _obsctx(series.model.system, θ, θ_obs, posys, traj, series.data_maps[obs])
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
star + planets system this is each planet about the star (matching v8's
octoplot); for hierarchies it generalizes with no special cases (a moon
about its planet, an inner pair's barycentre about the outer body, …).
"""
default_sky_queries(sys::System) = default_queries(sys, PlanetOrbits.raoff)

# Observables describing how a body *reflexes* about the system barycentre,
# rather than how two bodies are separated from each other. Asked for one of
# these with no observation to hang it on, the panel a user means is the
# host's own signal — nobody predicting "the RV of this system" wants the
# planet's velocity relative to its star.
const _REFLEX_OBSERVABLES = (
    PlanetOrbits.radvel, PlanetOrbits.velx, PlanetOrbits.vely, PlanetOrbits.velz,
    PlanetOrbits.pmra, PlanetOrbits.pmdec,
)

"""
    default_queries(sys, f) -> Vector{(ObservableQuery, Symbol)}

The natural queries of observable `f` for a system that has no observation
declaring them — what [`octoplot`](@ref) draws when `channels=` names a
quantity the model was not fitted to, so that a fit can predict a
not-yet-observed signal.

Two conventions, by the kind of observable:

  - **separations** (`raoff`, `projectedseparation`, `posx`, …): one query per
    hierarchy row, the exterior side about the interior side — exactly the
    relationship that row parametrizes, and the same set the sky panel draws.
  - **reflex signals** (`radvel`, `pmra`, the velocities): one query per root
    body — the bodies no row places — against the whole-system barycentre.

The second element of each pair names the row (or body) the query belongs to,
which is what the panel labels and colours itself by.
"""
function default_queries(sys::System, f)
    fn = _obsfunc(f)
    if fn in _REFLEX_OBSERVABLES
        return [(ObservableQuery(fn, refspec(h), Barycentre), h) for h in _rootbodies(sys)]
    end
    out = Tuple{ObservableQuery,Symbol}[]
    for (owner, ext, int) in sys.rows
        t = length(ext) == 1 ? refspec(ext[1]) : BarycentreSpec{ext}()
        r = length(int) == 1 ? refspec(int[1]) : BarycentreSpec{int}()
        push!(out, (ObservableQuery(fn, t, r), owner))
    end
    return out
end
export default_queries

# The bodies no hierarchy row places: the host(s) everything else orbits. A
# system whose every body is placed by a row has no root (it is described
# entirely relatively), and then the first body stands in.
function _rootbodies(sys::System)
    placed = Set{Symbol}()
    for (_, ext, _) in sys.rows, n in ext
        push!(placed, n)
    end
    roots = Symbol[n for n in sys.bodynames if !(n in placed)]
    return isempty(roots) ? Symbol[first(sys.bodynames)] : roots
end

"""
    predictedchannels(sys, f) -> Vector{Tuple{Nothing,PlotChannel}}

Model-only channels for observable `f` — one per query
[`default_queries`](@ref) picks, labelled from `PlanetOrbits.plotinfo`, with
no observation behind them.

This is how a fit draws a quantity it has no data for: `channels=radvel` on a
relative-astrometry fit predicts the reflex RV curve the orbit implies, which
is the figure you want when deciding whether a target is worth spectroscopic
time. Every observable is plottable this way whether or not it was observed;
what makes a panel is the model, not the dataset.
"""
function predictedchannels(sys::System, f)
    fn = _obsfunc(f)
    info = PlanetOrbits.plotinfo(fn)
    # `plotinfo` reports angles in radians. Every Octofitter panel shows
    # degrees, which is exactly what a channel's `scale`/`wrap` are for.
    isangle = info.unit == "rad"
    scale = isangle ? rad2deg(1.0) : 1.0
    return Tuple{Nothing,PlotChannel}[
        (nothing, PlotChannel(nameof(fn), info.label, isangle ? "°" : info.unit;
            scale, wrap=(info.wrap === nothing ? nothing : info.wrap * scale), query=q))
        for (q, _) in default_queries(sys, fn)]
end
export predictedchannels

# The observables a `channels=` restriction names outright. Only these can
# raise a predicted panel: a *channel* name (`:sep`) says "the separation data",
# and asking for data a model does not have is not a request to invent it.
_requestedobservables(::Nothing) = Any[]
_requestedobservables(f::Function) = f in _OBSERVABLE_FUNCS ? Any[f] : Any[]
function _requestedobservables(s::Symbol)
    for f in _OBSERVABLE_FUNCS
        nameof(f) === s && return Any[f]
    end
    return Any[]
end
_requestedobservables(cs) =
    reduce(vcat, (_requestedobservables(c) for c in cs); init=Any[])

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
# Defer *whether* an image format is offered to the figure, not just how it is
# written. Makie gates a bare `Figure` on the active backend
# (`showable(mime, ::FigureLike)` → `_backend_showable`), which is how
# `CairoMakie.activate!(type="png")` — the default — suppresses SVG. A wrapper
# type gets Julia's fallback instead, where `showable` is true for any MIME that
# has a `show` method, so the two methods above would advertise SVG
# unconditionally and defeat that gate. Documenter prefers SVG when it is
# offered, which turned every orbit plot in the manual into a ~28 MB vector file
# (main's PNGs are ~100-300 KB) and made the render stage take hours.
#
# Strictly these two MIMEs, matching the `show` methods above. Delegating *every*
# MIME is wrong: CairoMakie's supported set also contains `Makie.WEB_MIMES`, and
# `type="png"` disables only svg and pdf, so a blanket rule would claim
# `text/html` — which this type cannot render, and which Documenter probes first
# (`Documenter/src/utilities/utilities.jl:617`). That is a `MethodError` on the
# first plot in the manual, not a fallback.
Base.showable(mime::Union{MIME"image/png",MIME"image/svg+xml"}, res::OctoPlotResult) =
    showable(mime, res.figure)
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
channel or observable name, or a collection of either. A `channels=` the
model has no data for is drawn as a **prediction**: the model curves alone,
over the queries [`default_queries`](@ref) picks.

This is the **many-draw** figure, and three of its defaults follow from that:
residuals are whitened and each point carries the boxplot of its z-score over
the draws, phase-folded panels are off (a fold needs one ephemeris, and a
posterior need not have one), and correlated-noise bands are off (they would
double every curve). Each turns back on for `ndraws=1`, or explicitly via
`show_phase=`/`gpband=`. [`rvplot`](@ref) is the single-draw figure where all
three are the point.
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
period. Masses are M⊙ (v8 plotted Mⱼᵤₚ; v9 has one mass unit throughout).
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

Draw several draws side by side with [`gaiastarplot!`](@ref), which takes a
grid cell instead of making its own figure.
"""
function gaiastarplot end
gaiastarplot(args...; kwargs...) = _require_makie("gaiastarplot")
export gaiastarplot

"""
    gaiastarplot!(gridposition_or_axis, model, chain, sample_idx=MAP; kwargs...)

[`gaiastarplot`](@ref) into a cell of a figure you already have —
`gaiastarplot!(fig[i, j], model, chain, idx)` — or into an axis you made
yourself. Returns the `Axis`, so a grid of draws can be linked and have its
interior decorations hidden in the usual Makie way. Requires a Makie backend.
"""
function gaiastarplot! end
gaiastarplot!(args...; kwargs...) = _require_makie("gaiastarplot!")
export gaiastarplot!

"""
    gaiatimeplot(model, chain; kwargs...)

Along-scan abscissa against time for a [`GaiaDR4AstromObs`](@ref): the
posterior cloud of modelled abscissae over the measurements, with a per-epoch
boxplot of the residuals against the quoted formal errors below.

This is the same data as `octoplot`'s generic `:along_scan` panel, drawn in
v8's per-epoch-boxplot idiom, which answers a different question: not "are
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
    rvplot(model, chain, [sample_idx]; kwargs...) -> OctoPlotResult

The radial-velocity summary figure for a single posterior draw
(`sample_idx`, by default the highest-posterior-density sample): one
time-series panel carrying **every instrument at once**, with a residual strip
and marginal histogram, plus one phase-folded panel per hierarchy row that
moves the star.

This is the one figure allowed to put several RV instruments on one axis (see
[`sharepanel`](@ref)), and showing a single draw is what buys that. A
calibrated RV series is only defined per draw — the zero points, the jitters,
the trend and the other rows' subtracted signals all move between samples — so
everything here belongs to one sample and nothing is misrepresented.
[`octoplot`](@ref) is the many-draws view and gives each instrument its own
panel instead, with the data left uncalibrated and each draw's own offset and
trend carried by its model curve.

Requires a Makie backend. `rvplot_animated` records the same figure over
successive single-draw slices of the chain.

Called `rvpostplot` before v9. The old name still works and forwards here:
what the figure shows is one draw, not the posterior.
"""
function rvplot end
rvplot(args...; kwargs...) = _require_makie("rvplot")
function rvplot_animated end
rvplot_animated(args...; kwargs...) = _require_makie("rvplot_animated")
export rvplot, rvplot_animated

"""
    phasefoldpanel!(gridposition, series, entries; row, kwargs...)

Data-vs-model panel folded on hierarchy row `row`'s orbital phase: the
row's isolated signal (`rowsignal`) per posterior draw, calibrated
data with the other rows' signals removed, noise-weighted binned means, and
a phase-folded residual strip. Requires Makie. See [`octoplot`](@ref).
"""
function phasefoldpanel! end
export phasefoldpanel!
