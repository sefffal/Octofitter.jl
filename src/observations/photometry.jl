

const phot_cols = (:band, :phot, :σ_phot)
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
PhotometryLikelihood(observations::NamedTuple) = PhotometryLikelihood(observations)
export PhotometryLikelihood


# PhotometryLikelihood
function ln_like(photometry::PhotometryLikelihood, θ_planet, _elements=nothing, _interior_planets=nothing)
    ll = 0.0
    for i in eachindex(photometry.table.band)
        band = photometry.table.band[i]
        phot_param = getproperty(θ_planet, band)
        phot_meas = photometry.table.phot[i]
        if !isfinite(phot_param)
            return -Inf
        end
        # Experimenting with fitting sigma phot
        σ_phot = photometry.table.σ_phot[i]
        # σₓ = θ_planet_i.σ

        resid = phot_param - phot_meas
        σ² = σ_phot^2
        χ² = -(1/2)*resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    return ll
end
