

const phot_cols = (:band, :phot, :σ_phot)

"""
    PhotometryLikelihood(
        (band = :H, phot=15.0, σ_phot=3.),
        (band = :J, phot=13.5, σ_phot=0.5),
        (band = :K, phot=11.0, σ_phot=1.0);
        variables=@variables begin
            H ~ Uniform(0, 10)
            J ~ Uniform(0, 10)
            K ~ Uniform(0, 10)
        end
    )

A likelihood for comparing measured photometry points in one or more 
filter bands to data (provided here). Requires the `:band`, `:phot`,
and `:σ_phot` columns. Can be provided with any Tables.jl compatible
data source.

Band variables (e.g., H, J, K) should be defined in the `variables` block
rather than in the planet definition. These can be derived from physical
models that take planet mass and other system parameters as input.
"""
struct PhotometryLikelihood{TTable<:Table} <: AbstractLikelihood
    table::TTable
    priors::Priors
    derived::Derived
    function PhotometryLikelihood(
            observations;
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
        ii = sortperm(table.epoch)
        table = table[ii]
        return new{typeof(table)}(table, priors, derived)
    end
end
PhotometryLikelihood(observations::NamedTuple...; kwargs...) = PhotometryLikelihood(observations; kwargs...)
export PhotometryLikelihood

function likeobj_from_epoch_subset(obs::PhotometryLikelihood, obs_inds)
    return PhotometryLikelihood(obs.table[obs_inds,:,1]; variables=(obs.priors, obs.derived))
end

# PhotometryLikelihood: attached to a system
function ln_like(photometry::PhotometryLikelihood, θ_system, θ_obs, orbits::NTuple{N,<:AbstractOrbit}, args...) where N
    T = _system_number_type(θ_system)
    ll = zero(T)

    for i in eachindex(photometry.table.band)
        band = photometry.table.band[i]
        phot_param = getproperty(θ_obs, band)
        phot_meas = photometry.table.phot[i]
        if !isfinite(phot_param)
            return -Inf
        end
        # Experimenting with fitting sigma phot
        σ_phot = photometry.table.σ_phot[i]

        resid = phot_param - phot_meas
        σ² = σ_phot^2
        χ² = -(1/2)*resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    return ll
end

# PhotometryLikelihood: attached to a planet
function ln_like(photometry::PhotometryLikelihood, θ_system, θ_planet, θ_obs, orbits::AbstractOrbit, args...)
    T = _system_number_type(θ_system)
    ll = zero(T)

    for i in eachindex(photometry.table.band)
        band = photometry.table.band[i]
        phot_param = getproperty(θ_obs, band)
        phot_meas = photometry.table.phot[i]
        if !isfinite(phot_param)
            return -Inf
        end
        # Experimenting with fitting sigma phot
        σ_phot = photometry.table.σ_phot[i]

        resid = phot_param - phot_meas
        σ² = σ_phot^2
        χ² = -(1/2)*resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    return ll
end
