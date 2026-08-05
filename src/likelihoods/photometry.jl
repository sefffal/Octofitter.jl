# ---------------------------------------------------
# Photometry — `PhotometryObs`
#
# Compares a body's own `flux_<band>` variable against photometric
# measurements: `PhotometryObs(tab; target=b, band=:H, name="NIRC2")`. It
# touches no epoch and no orbit, and a `~` line in the body's `@variables`
# block would express the same comparison — but only an *observation* can be
# held out by cross-validation or simulated from a fitted model, which is why
# it stays a type. See `design/observation-types-migration.md` §3.4.
# ---------------------------------------------------

const phot_cols = (:phot, :σ_phot)

"""
    PhotometryObs(data; target, band=:default, name, variables=@variables begin end)

Photometry of one body in one band: compares `target`'s `flux_<band>`
variable against the measurements in `data`.

    phot = PhotometryObs(tab; target=b, band=:H, name="NIRC2")

`data` needs `:phot` and `:σ_phot` columns, in whatever flux units the model's
`flux_H` variables are in (setting the host's flux to 1.0 makes every other
body's a contrast ratio). Any Tables.jl source works. One band per
observation; for a second band, add a second `PhotometryObs`.

`band=:H` reads the variable `flux_H` from `target`'s `@variables` block;
`band=:default` (the default) reads a plain `flux`. This is the same
variable-name → band mapping [`Photocentre`](@ref) weights are built from, so
photometry and photocentres cannot disagree about what a body's brightness is.

# When to reach for this, and when not to

v9's `@variables` blocks take `~` over derived quantities, so the comparison
itself is a one-liner in `target`'s own block:

    flux_H ~ Normal(15.0, 3.0)      # an ad-hoc constraint, not data

That is the right spelling for an *ad-hoc constraint*. Use `PhotometryObs`
when the numbers are **data**, because of the general criterion this refactor
sorts every borderline case by:

> If it carries data you might want to hold out or simulate, it is an
> observation. If it only reshapes the prior, it is a `~` line or an
> `_isprior` term.

A `~` line becomes a `UserLikelihood` with `_isprior = true`: it is
excluded from likelihood counts, has no table to subset, and has no
`generate_from_params`. So it can never be held out by cross-validation and
you can never simulate photometry from a fitted model — both of which people
routinely do with real photometry.

# Variables

None are required. The forward model *is* the body's flux variable, which
lives in the body's block rather than here — that is the flux/band
unification: v8 declared `flux` inside the `PhotometryObs`'s own `variables`
block, so two instruments observing the same body in the same band each
carried their own independent flux parameter.
"""
struct PhotometryObs{TTable<:Table,TT,Var} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    target::TT
    name::String
end

function PhotometryObs(observations;
                       target,
                       band::Symbol=:default,
                       name,
                       variables::Tuple{Priors,Derived}=(Priors(), Derived()))
    (priors, derived) = variables
    # Collect the columns: a multithreaded `CSV.read` returns `ChainedVector`s,
    # which the per-epoch indexing in the likelihood cannot take a row view of.
    table = materialize_cols(Table(observations))
    equal_length_cols(table) ||
        error("The columns in the input data do not all have the same length")
    issubset(phot_cols, Tables.columnnames(table)) ||
        error("Expected columns $phot_cols")

    t = refspec(target)
    t isa BodyRefSpec || error(
        "PhotometryObs measures the flux of a single body, so `target` must be a " *
        "`Body` model node or a `Symbol` naming one; got $(_refstr(t)). " *
        "(A `Barycentre`/`Photocentre` is a point on the sky, not a brightness — " *
        "for the combined flux of a blended pair, give each member its own " *
        "`PhotometryObs`, or derive the sum in a system variable and constrain " *
        "that with a `~` line.)")

    # The band is carried as the *variable name* rather than the band symbol,
    # in a type parameter, so the lookup in `ln_like` is a compile-time
    # `getproperty` on the body's namespace with no runtime `Symbol` building.
    var = _flux_varname(band)
    return PhotometryObs{typeof(table),typeof(t),var}(
        table, priors, derived, t, String(name))
end

export PhotometryObs

# Inverse of `_flux_band` in model/nodes.jl, which maps a body's variable name
# to the band it declares. Keeping the two adjacent in meaning is what lets a
# `PhotometryObs` and a `Photocentre(:H, …)` read the same number.
_flux_varname(band::Symbol) = band === :default ? :flux : Symbol(:flux_, band)

_photband(::PhotometryObs{<:Any,<:Any,Var}) where {Var} = _flux_band(Var)
_photvar(::PhotometryObs{<:Any,<:Any,Var}) where {Var} = Var

refspecs(obs::PhotometryObs) = (obs.target,)

# No epoch, no orbit — so it contributes nothing to the epoch union even if the
# user's photometry table happens to carry an `:epoch` column for bookkeeping.
# (The generic fallback would read that column and force a trajectory solve at
# epochs nothing ever queries.)
epochs(::PhotometryObs) = Float64[]

# Not "target vs reference": photometry differences nothing. Say which
# variable it constrains instead, which is the whole content of the type.
_refdesc(obs::PhotometryObs) = _refstr(obs.target) * "." * string(_photvar(obs))

likeobj_from_epoch_subset(obs::PhotometryObs, inds) = PhotometryObs(
    obs.table[inds, :, 1]; target=obs.target, band=_photband(obs), obs.name,
    variables=(obs.priors, obs.derived))

"""
    simulate(obs::PhotometryObs, ctx) -> (; phot_model)

The modelled flux — `target`'s `flux_<band>` variable for this sample. Scalar,
repeated per row, since photometry has no epoch dependence; returned as a
NamedTuple for consistency with the other `simulate` methods.
"""
function simulate(obs::PhotometryObs, ctx::ObsContext)
    return (; phot_model=_targetflux(obs, ctx))
end

function ln_like(obs::PhotometryObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    flux = _targetflux(obs, ctx)
    # A non-finite flux is reachable from a legal draw (a derived flux out of
    # an atmosphere-model interpolation off the end of its grid), and rejecting
    # the sample is much better than propagating a NaN into the log density.
    isfinite(flux) || return convert(T, -Inf)

    ll = zero(T)
    for i in eachindex(obs.table.phot)
        resid = flux - obs.table.phot[i]
        σ² = obs.table.σ_phot[i]^2
        # v1's exact arithmetic, deliberately: this port is checked against v1
        # bit-for-bit, and `-log(sqrt(2πσ²))` is not the same float as
        # `-log(2πσ²)/2`.
        ll += -(1 / 2) * resid^2 / σ² - log(sqrt(2π * σ²))
    end
    return ll
end

function generate_from_params(obs::PhotometryObs, ctx::ObsContext; add_noise)
    flux = _targetflux(obs, ctx)
    σ_phot = obs.table.σ_phot
    phot = add_noise ? flux .+ randn.() .* σ_phot : fill(flux, length(σ_phot))
    return PhotometryObs(Table(; phot, σ_phot); target=obs.target,
        band=_photband(obs), obs.name, variables=(obs.priors, obs.derived))
end

# The body's flux variable, read from the model namespace rather than from
# `PlanetOrbits.fluxes(ctx.system, band)`. The two are the same number by
# construction (codegen builds the latter from the former), but `fluxes`
# reports a body that declares no flux in the band as **zero** — which is the
# right default for a photocentre weight and a silent disaster here, where it
# would compare 0.0 against the data and leave the sampler with a constant
# offset it cannot move. Reading the variable directly turns that into an
# error naming the variable the model is missing.
@inline function _targetflux(::PhotometryObs{<:Any,BodyRefSpec{Name},Var},
                             ctx::ObsContext) where {Name,Var}
    θ_body = getproperty(ctx.θ_system.bodies, Name)
    hasproperty(θ_body, Var) || _err_nofluxvar(Name, Var)
    return getproperty(θ_body, Var)
end

@noinline _err_nofluxvar(name::Symbol, var::Symbol) = error(
    "PhotometryObs: body :$name declares no `$var` variable, so there is nothing " *
    "for the photometry to constrain. Add it to :$name's `@variables` block " *
    "(`$var ~ Uniform(0, 10)`, or a derived expression from an evolutionary " *
    "model), or point the observation at a different `target=`/`band=`.")
