

const phot_cols = (:phot, :σ_phot)

"""
    data = Table(
        (phot=15.0, σ_phot=3.0),
        (phot=14.8, σ_phot=0.5),
    )
    PhotometryObs(
        data,
        name="INSTRUMENT",
        variables=@variables begin
            flux ~ Uniform(0, 10)
        end
    )

An observation type for comparing measured photometry points in a single
filter band to data (provided here). Requires the `:phot` and `:σ_phot`
columns. Can be provided with any Tables.jl compatible data source.

For multiple bands, create separate PhotometryObs objects.

The flux variable should be defined in the `variables` block rather than
in the planet definition. This can be derived from physical models that
take planet mass and other system parameters as input.

The `name` is used for variable naming in the chain output.
"""
struct PhotometryObs{TTable<:Table} <: AbstractObs
    table::TTable
    name::String
    priors::Priors
    derived::Derived
    function PhotometryObs(
            observations;
            name="PHOTOMETRY",
            variables::Tuple{Priors,Derived}=(@variables begin;end)
        )
        (priors,derived)=variables
        table = Table(observations)
        if !equal_length_cols(table)
            error("The columns in the input data do not all have the same length")
        end
        if !issubset(phot_cols, Tables.columnnames(table))
            error("Expected columns $phot_cols")
        end
        return new{typeof(table)}(table, name, priors, derived)
    end
end

# Backwards compatibility alias
const PhotometryLikelihood = PhotometryObs

export PhotometryObs, PhotometryLikelihood

function likeobj_from_epoch_subset(obs::PhotometryObs, obs_inds)
    return PhotometryObs(obs.table[obs_inds,:,1]; name=obs.name, variables=(obs.priors, obs.derived))
end

# PhotometryObs: attached to a system
function ln_like(photometry::PhotometryObs, ctx::SystemObservationContext)
    (; θ_obs) = ctx
    T = _system_number_type(ctx.θ_system)
    ll = zero(T)

    for i in eachindex(photometry.table.phot)
        flux_param = θ_obs.flux
        phot_meas = photometry.table.phot[i]
        if !isfinite(flux_param)
            return -Inf
        end
        σ_phot = photometry.table.σ_phot[i]

        resid = flux_param - phot_meas
        σ² = σ_phot^2
        χ² = -(1/2)*resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    return ll
end

# PhotometryObs: attached to a planet
function ln_like(photometry::PhotometryObs, ctx::PlanetObservationContext)
    (; θ_obs) = ctx
    T = _system_number_type(ctx.θ_system)
    ll = zero(T)

    for i in eachindex(photometry.table.phot)
        flux_param = θ_obs.flux
        phot_meas = photometry.table.phot[i]
        if !isfinite(flux_param)
            return -Inf
        end
        σ_phot = photometry.table.σ_phot[i]

        resid = flux_param - phot_meas
        σ² = σ_phot^2
        χ² = -(1/2)*resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    return ll
end
