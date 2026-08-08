# ---------------------------------------------------
# `ObsPriorONeil2019` — the observable-based prior
#
# Wraps another likelihood (`ObsPriorONeil2019(astrom_like)`) and contributes
# O'Neil et al. 2019's Jacobian over that likelihood's epochs. It wraps data,
# so it is an observation with `_isprior = false`.
#
# The orbit is named explicitly, and epoch lookup goes through
# `ctx.epoch_index` rather than by matching `sol.t == epoch`: that equality
# does not hold once light-travel time is modelled, so a scan would silently
# find the wrong solution. One method covers astrometry and RV alike.
#
# Included after relative-astrometry.jl and radial-velocity.jl: it wraps them.
# See `design/observation-types-migration.md` §3.7.
# ---------------------------------------------------

"""
    ObsPriorONeil2019(wrapped; orbit=nothing, name=…)

The observable-based priors of K. O'Neil et al. 2019, *"Improving Orbit
Estimates for Incomplete Orbits with a New Approach to Priors: with
Applications from Black Holes to Planets"*, applied on top of `wrapped`.

    astrom = RelAstromObs(tab; target=b, ref=A, name="GPI")
    System(…, observations=[ObsPriorONeil2019(astrom)], …)

**Pass only the wrapper**, not the wrapped likelihood as well: this object
evaluates `wrapped`'s own log-likelihood and adds the Jacobian term to it, so
listing both double-counts the data.

The correction is only correct if you put `Uniform` priors on all Campbell
orbital elements **and** a `Uniform` prior on period (not on semi-major axis).
The period prior's range has a large effect on the fit and no recommendation
was published with the paper.

# The orbit

The Jacobian is a property of one orbit evaluated at the wrapped likelihood's
epochs. In v8 that orbit was implicit — the likelihood was attached to a
planet — and in v9 nothing attaches, so it is named:

  - `orbit=nothing` (the default) takes the wrapped observation's own
    `target`, which is right for relative astrometry (`target=b, ref=A` → `b`'s
    orbit) and for relative RVs.
  - `orbit=b` names it explicitly. **Stellar-reflex radial velocities need
    this**: a `RadialVelocityObs(…; target=A, ref=Barycentre)` measures the
    host, whose motion is caused by the companion whose orbit the prior is
    about, and `A` has no orbit of its own.
  - `orbit=(b, c)` sums the term over several orbits, which is what v8's
    system-attached method did over every planet.

The orbit meant by a body is the hierarchy row that *places* it, so any
hierarchy convention may be named.

# What changed inside

v8 recovered each epoch's solution by scanning the planet's solution vector
for `sol.t == epoch`, falling back to re-solving the orbit when no match was
found — O(N) per epoch, and a silent fallback for `AbsoluteVisual` orbits,
whose stored solution times are light-travel corrected and so never compare
equal. v9 reads the row's own elements and the trajectory's emission epoch,
which is exactly what the propagator's Kepler solve used, at O(1).

Note the anomalies come from the row's Keplerian elements rather than from an
integrated trajectory, deliberately: this is a prior over the *sampled
elements*, so it means the same thing under `KeplerianApprox` and under a
full N-body propagator, where the row only sets initial conditions.
"""
struct ObsPriorONeil2019{TL<:AbstractObs,TTable,TO<:Tuple} <: AbstractObs
    wrapped_like::TL
    table::TTable
    priors::Priors
    derived::Derived
    orbit::TO
end

function ObsPriorONeil2019(wrapped::AbstractObs; orbit=nothing)
    specs = _oneil_orbitspecs(wrapped, orbit)
    return ObsPriorONeil2019{typeof(wrapped),typeof(wrapped.table),typeof(specs)}(
        wrapped, wrapped.table, wrapped.priors, wrapped.derived, specs)
end

export ObsPriorONeil2019

# Default: the wrapped observation's own target. It is a real default rather
# than a guess for the astrometric case — "the orbit of the body whose position
# this measures" is what v1's per-planet attachment meant — but it is wrong for
# a stellar-reflex RV, whose target is the host, so an unplaceable default is
# reported at model-build time (through `refspecs`) rather than silently.
function _oneil_orbitspecs(wrapped::AbstractObs, orbit)
    specs = if orbit === nothing
        hasproperty(wrapped, :target) || error(
            "ObsPriorONeil2019 could not infer which orbit it applies to: " *
            "$(typeof(wrapped).name.name) has no `target`. Pass `orbit=b`.")
        (refspec(wrapped.target),)
    elseif orbit isa Tuple || orbit isa AbstractVector
        map(refspec, Tuple(orbit))
    else
        (refspec(orbit),)
    end
    if !all(s -> s isa BodyRefSpec, specs)
        got = join(map(_refstr, specs), ", ")
        error(
            "ObsPriorONeil2019's `orbit=` names the body (or bodies) whose orbits " *
            "the Jacobian is computed for, so each entry must be a `Body` model " *
            "node or a `Symbol` naming one; got $got. A `Barycentre`/`Photocentre` " *
            "has no orbit of its own — name the companion instead, e.g. " *
            "`ObsPriorONeil2019(rvs; orbit=b)`.")
    end
    return specs
end

likelihoodname(obs::ObsPriorONeil2019) = "obspri_" * likelihoodname(obs.wrapped_like)

# It wraps data, and it only reaches the log density through calculations
# involving those data points, so it is *not* one of the prior-shaped terms —
# it counts as a likelihood and participates in cross-validation.
_isprior(::ObsPriorONeil2019) = false

# The wrapper's epochs are the wrapped likelihood's epochs, so `ctx` built for
# the wrapper indexes the trajectory correctly for both.
epochs(obs::ObsPriorONeil2019) = epochs(obs.wrapped_like)
refspecs(obs::ObsPriorONeil2019) = (refspecs(obs.wrapped_like)..., obs.orbit...)
_refdesc(obs::ObsPriorONeil2019) = _refdesc(obs.wrapped_like) *
                                   "  [Jacobian over " *
                                   join(map(_refstr, obs.orbit), ", ") * "]"

likeobj_from_epoch_subset(obs::ObsPriorONeil2019, inds) = ObsPriorONeil2019(
    likeobj_from_epoch_subset(obs.wrapped_like, inds); orbit=obs.orbit)

# Simulation and data generation are the wrapped likelihood's business: the
# Jacobian reweights the prior, it does not change what the instrument saw.
simulate(obs::ObsPriorONeil2019, ctx::ObsContext) = simulate(obs.wrapped_like, ctx)

# ...but the *wrapper* has to survive, exactly as it does in
# `likeobj_from_epoch_subset` above. Returning the bare inner likelihood would
# hand back a different model: the Jacobian term is gone, and
# `likelihoodname` changes from `obspri_<name>` to `<name>`, so a chain fit to
# the simulated system does not line up column-for-column with one fit to the
# original — which is precisely what posterior-predictive checks, SBC and
# completeness compare.
generate_from_params(obs::ObsPriorONeil2019, ctx::ObsContext; add_noise) =
    ObsPriorONeil2019(generate_from_params(obs.wrapped_like, ctx; add_noise);
                      orbit=obs.orbit)

function ln_like(obs::ObsPriorONeil2019, ctx::ObsContext)
    ll = ln_like(obs.wrapped_like, ctx)
    ks = _selected_rows(ctx, obs.orbit)
    posys = ctx.system
    nep = length(ctx.epoch_index)
    for k in ks
        row = @inbounds posys.rows[k]
        e = row.e
        # `P` in years; the Jacobian's cbrt(P) factor is the paper's, and the
        # unit matters because the prior is improper up to a constant only
        # within one convention.
        P = PlanetOrbits._period(row) / PlanetOrbits.year2day_julian
        jac = zero(typeof(row.a))
        for i in 1:nep
            M, E = _row_anomalies(row, solutionat(ctx, i))
            jac += abs(3M * (e + cos(E)) + 2 * (-2 + e^2 + e * cos(E)) * sin(E))
        end
        jac *= cbrt(P) / sqrt(1 - e^2)
        ll += 2log(jac)
    end
    return ll
end

"""
    _row_anomalies(row, sol) -> (M, E)

Mean and eccentric anomaly [rad] of hierarchy `row` at the epoch of `sol`.

`M` is the *wrapped* mean anomaly `E - e sin E`, which is what v8's
`meananom(sol)` returned (it reconstructed it from the solved `E` rather than
from `n·(t − tp)`), so the Jacobian sees the same value it always did.

The epoch is the trajectory's emission epoch, not the observation epoch: that
is the time the propagator's own Kepler solve used, and taking it from there
is what fixes v8's `AbsoluteVisual` hole — v8 matched solutions by
`sol.t == epoch`, which light-travel time makes false, and silently re-solved
at the *uncorrected* epoch instead. For every other frame the two are the same
number.
"""
@inline function _row_anomalies(row, sol)
    t_em = @inbounds sol.traj.t_em[sol.k]
    MA = row.n / PlanetOrbits.year2day_julian * (t_em - row.tp)
    E = PlanetOrbits.kepler_solver(MA, row.e)
    return (E - row.e * sin(E), E)
end

# The O'Neil prior reweights another observation; the corrections reach it
# only through that one.
has_correction_impact(::Type{<:ObsPriorONeil2019{TL}}) where {TL} =
    has_correction_impact(TL)
correction_impact(obs::ObsPriorONeil2019, a::ObsContext, b::ObsContext) =
    correction_impact(obs.wrapped_like, a, b)
correction_advisories(obs::ObsPriorONeil2019, ctx::ObsContext) =
    correction_advisories(obs.wrapped_like, ctx)
