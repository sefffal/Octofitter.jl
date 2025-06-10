

const phot_cols = (:phot, :σ_phot)

"""
    data = Table(
        (phot=15.0, σ_phot=3.0),
        (phot=14.8, σ_phot=0.5),
    )
    PhotometryLikelihood(
        data,
        instrument_name="INSTRUMENT",
        variables=@variables begin
            flux ~ Uniform(0, 10)
        end
    )

A likelihood for comparing measured photometry points in a single
filter band to data (provided here). Requires the `:phot` and `:σ_phot`
columns. Can be provided with any Tables.jl compatible data source.

For multiple bands, create separate PhotometryLikelihood objects.

The flux variable should be defined in the `variables` block rather than 
in the planet definition. This can be derived from physical models that 
take planet mass and other system parameters as input.

The `instrument_name` is used for variable naming in the chain output.
"""
struct PhotometryLikelihood{TTable<:Table} <: AbstractLikelihood
    table::TTable
    instrument_name::String
    priors::Priors
    derived::Derived
    function PhotometryLikelihood(
            observations;
            instrument_name="PHOTOMETRY",
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
        return new{typeof(table)}(table, instrument_name, priors, derived)
    end
end
# PhotometryLikelihood(observations::NamedTuple...; kwargs...) = PhotometryLikelihood(observations; kwargs...)
export PhotometryLikelihood

function likeobj_from_epoch_subset(obs::PhotometryLikelihood, obs_inds)
    return PhotometryLikelihood(obs.table[obs_inds,:,1]; instrument_name=obs.instrument_name, variables=(obs.priors, obs.derived))
end

# PhotometryLikelihood: attached to a system
function ln_like(photometry::PhotometryLikelihood, θ_system, θ_obs, orbits::NTuple{N,<:AbstractOrbit}, args...) where N
    T = _system_number_type(θ_system)
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

# PhotometryLikelihood: attached to a planet
function ln_like(photometry::PhotometryLikelihood, θ_system, θ_planet, θ_obs, args...)
    T = _system_number_type(θ_system)
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
