

const phot_cols = (:band, :phot, :σ_phot)

"""
    PhotometryLikelihood(
        (band = :Z, phot=15.0, σ_phot=3.),
        (band = :J, phot=13.5, σ_phot=0.5),
        (band = :L, phot=11.0, σ_phot=1.0)
    )

A likelihood for comparing measured photometry points in one or more 
filter bands to data (provided here). Requires the `:band`, `:phot',
and `:σ_phot` columns. Can be provided with any Tables.jl compatible
data source.
"""
struct PhotometryLikelihood{TTable<:Table} <: AbstractLikelihood
    table::TTable
    function PhotometryLikelihood(observations...)
        table = Table(observations...)
        if !issubset(phot_cols, Tables.columnnames(table))
            error("Expected columns $phot_cols")
        end
        return new{typeof(table)}(table)
    end
end
PhotometryLikelihood(observations::NamedTuple...) = PhotometryLikelihood(observations)
export PhotometryLikelihood

function likeobj_from_epoch_subset(obs::PhotometryLikelihood, obs_inds)
    return PhotometryLikelihood(obs.table[obs_inds,:,1]...)
end

# PhotometryLikelihood: attached to a system
function ln_like(photometry::PhotometryLikelihood, θ_system, orbits::NTuple{N,<:AbstractOrbit}, args...) where N
    T = _system_number_type(θ_system)
    ll = zero(T)

    for i in eachindex(photometry.table.band)
        band = photometry.table.band[i]
        phot_param = getproperty(θ_system, band)
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
function ln_like(photometry::PhotometryLikelihood, θ_planet, orbits::AbstractOrbit, args...)
    T = _system_number_type(θ_planet)
    ll = zero(T)

    for i in eachindex(photometry.table.band)
        band = photometry.table.band[i]
        phot_param = getproperty(θ_planet, band)
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
