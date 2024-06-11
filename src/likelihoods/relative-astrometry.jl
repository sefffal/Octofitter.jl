
# PlanetRelAstromLikelihood Data type
const astrom_cols1 = (:epoch, :ra, :dec, :σ_ra, :σ_dec)
const astrom_cols3 = (:epoch, :pa, :sep, :σ_pa, :σ_sep)

"""
    PlanetRelAstromLikelihood(
        (epoch = 5000, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 5050, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 5100, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
    )

Represents a likelihood function of relative astometry between a host star and a secondary body.
`:epoch` is a required column, in addition to either `:ra`, `:dec`, `:σ_ra`, `:σ_dec` or `:pa`, `:sep`, `:σ_pa`, `:σ_sep`.
All units are in **milliarcseconds** or **radians** as appropriate.

In addition to the example above, any Tables.jl compatible source can be provided.
"""
struct PlanetRelAstromLikelihood{TTable<:Table} <: AbstractLikelihood
    table::TTable
    function PlanetRelAstromLikelihood(observations...)
        table = Table(observations...)
        if !issubset(astrom_cols1, Tables.columnnames(table)) && 
           !issubset(astrom_cols3, Tables.columnnames(table))
            error("Expected columns $astrom_cols1 or $astrom_cols3")
        end
        return new{typeof(table)}(table)
    end
end
PlanetRelAstromLikelihood(observations::NamedTuple...) = PlanetRelAstromLikelihood(observations)
export PlanetRelAstromLikelihood


# Plot recipe for astrometry data
using LinearAlgebra

# PlanetRelAstromLikelihood likelihood function
function ln_like(astrom::PlanetRelAstromLikelihood, θ_planet, orbit,)

    # Note: since astrometry data is stored in a typed-table, the column name
    # checks using `hasproperty` ought to be compiled out completely.

    ll = 0.0
    for i in eachindex(astrom.table.epoch)

        # Covariance between the two dimensions
        cor = 0.0

        local o
        try
            o = orbitsolve(orbit, astrom.table.epoch[i])
        catch err
            @warn "Error during orbit solve (maxlog=10)" maxlog=10 err
            return -Inf
        end
        # PA and Sep specified
        if hasproperty(astrom.table, :pa) && hasproperty(astrom.table, :sep)
            ρ = projectedseparation(o)
            pa = posangle(o)

            pa_diff = ( astrom.table.pa[i] - pa + π) % 2π - π;
            pa_diff = pa_diff < -π ? pa_diff + 2π : pa_diff;
            resid1 = pa_diff
            resid2 = astrom.table.sep[i] - ρ
            σ₁ = astrom.table.σ_pa[i ]
            σ₂ = astrom.table.σ_sep[i]

        # RA and DEC specified
        else
            x = raoff(o)
            y = decoff(o)
            resid1 = astrom.table.ra[i] - x
            resid2 = astrom.table.dec[i] - y
            σ₁ = astrom.table.σ_ra[i ]
            σ₂ = astrom.table.σ_dec[i]

            # Add non-zero correlation if present
            if hasproperty(astrom.table, :cor)
                cor = astrom.table.cor[i]
            end
        end

        # Dispatch to a slightly faster (4%) version when
        # correlation between x and y are zero.
        if cor == 0
            # Manual definition:
            # χ²1 = -(1/2)*resid1^2 / σ₁^2 - log(sqrt(2π * σ₁^2))
            # χ²2 = -(1/2)*resid2^2 / σ₂^2 - log(sqrt(2π * σ₂^2))
            # ll += χ²1 + χ²2

            # Leveraging Distributions.jl to make this clearer:
            χ²1 = logpdf(Normal(0, σ₁), resid1)
            χ²2 = logpdf(Normal(0, σ₂), resid2)
            ll += χ²1 + χ²2

        else
            # Same as above, with support for covariance:
            Σ = @SArray[
                σ₁^2        cor*σ₁*σ₂
                cor*σ₁*σ₂   σ₂^2
            ]
            dist = MvNormal(Σ)
            resids = @SArray[resid1, resid2]
            ll += logpdf(dist, resids)
        end

    end
    return ll
end

# Generate new astrometry observations
function generate_from_params(like::PlanetRelAstromLikelihood, θ_planet, orbit::PlanetOrbits.AbstractOrbit)

    # Get epochs and uncertainties from observations
    epoch = like.table.epoch

    if hasproperty(like.table, :pa) && hasproperty(like.table, :sep)

        σ_sep = like.table.σ_sep 
        σ_pa = like.table.σ_pa

        # Generate now astrometry data
        sep = projectedseparation.(orbit, epoch)
        pa = posangle.(orbit, epoch)
        if hasproperty(like.table, :cov)
            astrometry_table = Table(;epoch, sep, pa, σ_sep, σ_pa, like.table.cov)
        else
            astrometry_table = Table(;epoch, sep, pa, σ_sep, σ_pa)
        end
    else
        σ_ra = like.table.σ_ra 
        σ_dec = like.table.σ_dec

        # Generate now astrometry data
        ra = raoff.(orbit, epoch)
        dec = decoff.(orbit, epoch)
        if hasproperty(like.table, :cov)
            astrometry_table = Table(;epoch, ra, dec, σ_ra, σ_dec, like.table.cov)
        else
            astrometry_table = Table(;epoch, ra, dec, σ_ra, σ_dec)
        end
    end

    return PlanetRelAstromLikelihood(astrometry_table)
end
