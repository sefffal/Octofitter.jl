
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
struct PlanetRelAstromLikelihood{TTable<:Table,TDistTuple} <: AbstractLikelihood
    table::TTable
    precomputed_pointwise_distributions::TDistTuple
    function PlanetRelAstromLikelihood(observations...)
        table = Table(observations...)
        if !equal_length_cols(table)
            error("The columns in the input data do not all have the same length")
        end
        if !issubset(astrom_cols1, Tables.columnnames(table)) && 
           !issubset(astrom_cols3, Tables.columnnames(table))
            error("Expected columns $astrom_cols1 or $astrom_cols3")
        end


        # For data points with a non-zero correlation, we can speed things up slightly
        # by pre-computing the appropriate 2x2 matrix factorizations. 
        # Distributions.jl handles this all for us--we just create N MvNormal distributions
        # PA and Sep specified
        if hasproperty(table, :pa) && hasproperty(table, :sep)
            σ₁ = table.σ_pa
            σ₂ = table.σ_sep
        # RA and DEC specified
        else
            σ₁ = table.σ_ra
            σ₂ = table.σ_dec
        end

        # Add non-zero correlation if present
        if hasproperty(table, :cor)
            cor = table.cor

            precomputed_pointwise_distributions = broadcast(σ₁, σ₂, cor) do σ₁, σ₂, cor
                Σ = @SArray[
                    σ₁^2        cor*σ₁*σ₂
                    cor*σ₁*σ₂   σ₂^2
                ]
                dist = MvNormal(Σ)
                return dist
            end
        else
            precomputed_pointwise_distributions = broadcast(σ₁, σ₂) do σ₁, σ₂
                Σ = Diagonal(@SArray[σ₁^2, σ₂^2])
                dist = MvNormal(Σ)
                return dist
            end
        end

        precomputed_pointwise_distributions_tuple = (precomputed_pointwise_distributions...,)
        return new{typeof(table),typeof(precomputed_pointwise_distributions_tuple)}(table, precomputed_pointwise_distributions_tuple)
    end
end
PlanetRelAstromLikelihood(observations::NamedTuple...) = PlanetRelAstromLikelihood(observations)
export PlanetRelAstromLikelihood

function likeobj_from_epoch_subset(obs::PlanetRelAstromLikelihood, obs_inds)
    return PlanetRelAstromLikelihood(obs.table[obs_inds,:,1]...)
end

# Plot recipe for astrometry data
using LinearAlgebra

# PlanetRelAstromLikelihood likelihood function
function ln_like(astrom::PlanetRelAstromLikelihood, θ_system, θ_planet, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)

    # Note: since astrometry data is stored in a typed-table, the column name
    # checks using `hasproperty` ought to be compiled out completely.

    this_orbit = orbits[i_planet]
    T = _system_number_type(θ_system)
    ll = zero(T)
    for i_epoch in eachindex(astrom.table.epoch)

        # Epicycle approximation: consider outer planets
        # orbit the centre of mass of the star and any inner planets.
        # Account for the shifting photocentre (assume the star dominates)
        # due to inner orbits.
        # Note! this is essentially free to calculate, since we already pre-solved
        # all orbits at all epochs.

        ra_host_perturbation = zero(T)
        dec_host_perturbation = zero(T)
        for i_other_planet in eachindex(orbits)
            orbit_other = orbits[i_other_planet]
            # Only account for inner planets
            if semimajoraxis(orbit_other) < semimajoraxis(this_orbit)
                θ_planet′ = θ_system.planets[i_other_planet]
                mass_other = θ_planet′.mass*Octofitter.mjup2msol
                sol = orbit_solutions[i_other_planet][i_epoch + orbit_solutions_i_epoch_start]
                # Note about `total mass`: for this to be correct, user will have to specify
                # `M` at the planet level such that it doesn't include the outer planets.
                ra_host_perturbation += raoff(sol)*mass_other/totalmass(orbit_other)   
                dec_host_perturbation += decoff(sol)*mass_other/totalmass(orbit_other)
            end
        end

        # Take the measurement, and *add* the Delta, to get what we compare to the model
        sol = orbit_solutions[i_planet][i_epoch + orbit_solutions_i_epoch_start]

        ra_model = (raoff(sol) + ra_host_perturbation )
        dec_model = (decoff(sol) + dec_host_perturbation )

        # PA and Sep specified
        if hasproperty(astrom.table, :pa) && hasproperty(astrom.table, :sep)
            ρ = hypot(ra_model,dec_model)
            pa = atan(ra_model,dec_model)

            pa_diff = ( astrom.table.pa[i_epoch] - pa + π) % 2π - π;
            pa_diff = pa_diff < -π ? pa_diff + 2π : pa_diff;
            resid1 = pa_diff
            resid2 = astrom.table.sep[i_epoch] - ρ

        # RA and DEC specified
        else
            resid1 = astrom.table.ra[i_epoch]  - ra_model
            resid2 = astrom.table.dec[i_epoch]  -  dec_model
        end

        ll += logpdf(astrom.precomputed_pointwise_distributions[i_epoch], @SVector[resid1, resid2])
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
