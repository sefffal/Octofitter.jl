
# PlanetRelAstromLikelihood Data type
const astrom_cols1 = (:epoch, :ra, :dec, :σ_ra, :σ_dec)
const astrom_cols3 = (:epoch, :pa, :sep, :σ_pa, :σ_sep)

"""
    data = Table(
        (epoch = 5000, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 5050, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
        (epoch = 5100, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
    )
    PlanetRelAstromLikelihood(data)

Represents a likelihood function of relative astometry between a host star and a secondary body.
`:epoch` is a required column, in addition to either `:ra`, `:dec`, `:σ_ra`, `:σ_dec` or `:pa`, `:sep`, `:σ_pa`, `:σ_sep`.
All units are in **milliarcseconds** or **radians** as appropriate.

In addition to the example above, any Tables.jl compatible source can be provided.
"""
struct PlanetRelAstromLikelihood{TTable<:Table,TDistTuple} <: AbstractLikelihood
    table::TTable
    priors::Priors
    derived::Derived
    precomputed_pointwise_distributions::TDistTuple
    name::String
    function PlanetRelAstromLikelihood(
            observations;
            variables::Tuple{Priors,Derived}=(@variables begin;end),
            name
        )
        (priors,derived)=variables
        table = Table(observations)
        if !equal_length_cols(table)
            error("The columns in the input data do not all have the same length")
        end
        if !issubset(astrom_cols1, Tables.columnnames(table)) && 
           !issubset(astrom_cols3, Tables.columnnames(table))
            error("Expected columns $astrom_cols1 or $astrom_cols3")
        end

        if any(>=(mjd("2050")),  table.epoch) || any(<=(mjd("1950")),  table.epoch)
            @warn "The data you entered fell outside the range year 1950 to year 2050. The expected input format is MJD (modified julian date). We suggest you double check your input data!"
        end


        ii = sortperm(vec(table.epoch))
        table = table[ii]

        # For data points with a non-zero correlation, we can speed things up slightly
        # by pre-computing the appropriate 2x2 matrix factorizations. 
        # Distributions.jl handles this all for us--we just create N MvNormal distributions
        # PA and Sep specified
        if hasproperty(table, :pa) && hasproperty(table, :sep)
            σ₁ = table.σ_pa
            σ₂ = table.σ_sep

            if any(>=(2pi),  table.pa) || any(<=(-2pi),  table.pa)
                @warn "The data you entered fell outside the range [-2pi, +2pi]. The expected input format is radians (you can use `deg2rad` to convert). We suggest you double check your input data!"
            end
        # RA and DEC specified
        else
            σ₁ = table.σ_ra
            σ₂ = table.σ_dec
        end

        # Add non-zero correlation if present
        if hasproperty(table, :cor)
            cor = table.cor

            if any(abs.(cor) .> 1 - 1e-5)
                error("Correlation values may not be well-specified: $cor")
            end

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
        return new{typeof(table),typeof(precomputed_pointwise_distributions_tuple),}(
            table, priors, derived, precomputed_pointwise_distributions_tuple, name
        )
    end
end
# PlanetRelAstromLikelihood(observations::NamedTuple..., (priors,derived)::Tuple{Priors,Derived}; kwargs...) = PlanetRelAstromLikelihood(observations, (priors,derived); kwargs...)
export PlanetRelAstromLikelihood


# In-place simulation logic for PlanetRelAstromLikelihood (performance-critical)
function simulate!(ra_model_buf, dec_model_buf, astrom::PlanetRelAstromLikelihood, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
    T = _system_number_type(θ_system)
    this_orbit = orbits[i_planet]
    
    # Compute model astrometry for each epoch
    for i_epoch in eachindex(astrom.table.epoch)
        # Get the orbit solution for this planet at this epoch
        sol = orbit_solutions[i_planet][i_epoch + orbit_solutions_i_epoch_start]
        @assert isapprox(astrom.table.epoch[i_epoch], PlanetOrbits.soltime(sol), rtol=1e-2)
        
        # Calculate perturbations from inner planets
        ra_host_perturbation = zero(T)
        dec_host_perturbation = zero(T)
        for (i_other_planet, key) in enumerate(keys(θ_system.planets))
            orbit_other = orbits[i_other_planet]
            # Only account for inner planets with non-zero mass
            if semimajoraxis(orbit_other) < semimajoraxis(this_orbit)
                θ_planet′ = θ_system.planets[key]
                if !hasproperty(θ_planet′, :mass)
                    continue
                end
                mass_other = θ_planet′.mass*Octofitter.mjup2msol
                sol′ = orbit_solutions[i_other_planet][i_epoch + orbit_solutions_i_epoch_start]
                
                ra_host_perturbation += raoff(sol′, mass_other)
                dec_host_perturbation += decoff(sol′, mass_other)
                
                @assert isapprox(astrom.table.epoch[i_epoch], PlanetOrbits.soltime(sol′), rtol=1e-2)
            end
        end

        # Compute the model astrometry values
        # This is distance(outer planet, inner barycentre) + distance(inner barycentre, star)
        ra_model_buf[i_epoch] = raoff(sol) - ra_host_perturbation
        dec_model_buf[i_epoch] = decoff(sol) - dec_host_perturbation
    end
    
    return (ra_model = ra_model_buf, dec_model = dec_model_buf, epochs = astrom.table.epoch)
end

# Allocating simulation logic for PlanetRelAstromLikelihood (convenience method)
function simulate(astrom::PlanetRelAstromLikelihood, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
    T = _system_number_type(θ_system)
    L = length(astrom.table.epoch)
    ra_model_buf = Vector{T}(undef, L)
    dec_model_buf = Vector{T}(undef, L)
    return simulate!(ra_model_buf, dec_model_buf, astrom, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
end


function likeobj_from_epoch_subset(obs::PlanetRelAstromLikelihood, obs_inds)
    return PlanetRelAstromLikelihood(
        obs.table[obs_inds,:,1];
        obs.name,
        variables=(obs.priors, obs.derived,)
    )
end

# Plot recipe for astrometry data
using LinearAlgebra

# PlanetRelAstromLikelihood likelihood function
function ln_like(astrom::PlanetRelAstromLikelihood, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
    T = Octofitter._system_number_type(θ_system)
   
    jitter = hasproperty(θ_obs, :jitter) ? getproperty(θ_obs, :jitter) : zero(T)
    platescale = hasproperty(θ_obs, :platescale) ? getproperty(θ_obs, :platescale) : one(T)
    northangle = hasproperty(θ_obs, :northangle) ? getproperty(θ_obs, :northangle) : zero(T)

    L = length(astrom.table.epoch)
    ll = zero(T)
    
    @no_escape begin
        # Allocate buffers using bump allocator
        ra_model_buf = @alloc(T, L)
        dec_model_buf = @alloc(T, L)
        
        # Use in-place simulation method to get model values
        sim = Octofitter.simulate!(ra_model_buf, dec_model_buf, astrom, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)

        # Process each epoch to compute residuals and likelihood
        for i_epoch in eachindex(astrom.table.epoch)
            ra_model = ra_model_buf[i_epoch]
            dec_model = dec_model_buf[i_epoch]

            # PA and Sep specified
            if hasproperty(astrom.table, :pa) && hasproperty(astrom.table, :sep)
                ρ = hypot(ra_model, dec_model)
                pa = atan(ra_model, dec_model)

                pa_dat = astrom.table.pa[i_epoch] + northangle
                pa_diff = (pa_dat - pa + π) % 2π - π
                pa_diff = pa_diff < -π ? pa_diff + 2π : pa_diff
                resid1 = pa_diff
                resid2 = astrom.table.sep[i_epoch] * platescale - ρ

            # RA and DEC specified
            else
                pa_dat = atan(astrom.table.dec[i_epoch], astrom.table.ra[i_epoch]) + northangle
                sep_dat = hypot(astrom.table.dec[i_epoch], astrom.table.ra[i_epoch]) * platescale
                ra_dat = sep_dat * cos(pa_dat)
                dec_dat = sep_dat * sin(pa_dat)
                resid1 = ra_dat - ra_model
                resid2 = dec_dat - dec_model
            end

            if jitter == 0.
                ll += logpdf(astrom.precomputed_pointwise_distributions[i_epoch], @SVector[resid1, resid2])
            else
                # For data points with a non-zero correlation, we can speed things up slightly
                # by pre-computing the appropriate 2x2 matrix factorizations. 
                # Distributions.jl handles this all for us--we just create N MvNormal distributions
                # PA and Sep specified
                if hasproperty(astrom.table, :pa) && hasproperty(astrom.table, :sep)
                    σ₁ = astrom.table.σ_pa[i_epoch]
                    σ₂ = astrom.table.σ_sep[i_epoch]
                # RA and DEC specified
                else
                    σ₁ = astrom.table.σ_ra[i_epoch]
                    σ₂ = astrom.table.σ_dec[i_epoch]
                end
                # Add jitter in quadrature
                σ₁ = hypot(σ₁, jitter)
                σ₂ = hypot(σ₂, jitter)
                # we have to compute the factorization on the fly
                if hasproperty(astrom.table, :cor)
                    cor = astrom.table.cor[i_epoch]
                    Σ = @SArray[
                        σ₁^2        cor*σ₁*σ₂
                        cor*σ₁*σ₂   σ₂^2
                    ]
                    dist = MvNormal(Σ)
                else
                    Σ = Diagonal(@SArray[σ₁^2, σ₂^2])
                    dist = MvNormal(Σ)
                end
                ll += logpdf(dist, @SVector[resid1, resid2])
            end
        end
    end
    return ll
end

# Generate new astrometry observations
function generate_from_params(like::PlanetRelAstromLikelihood, θ_system,  θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start; add_noise)
    # Get epochs and uncertainties from observations
    epoch = like.table.epoch

    # Use the same simulation method as ln_like to generate model astrometry values  
    sim = Octofitter.simulate(like, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
    ra_model = sim.ra_model
    dec_model = sim.dec_model

    # Apply platescale and northangle corrections (which ln_like also handles)
    platescale = hasproperty(θ_obs, :platescale) ? θ_obs.platescale : 1.0
    northangle = hasproperty(θ_obs, :northangle) ? θ_obs.northangle : 0.0
    jitter = hasproperty(θ_obs, :jitter) ? getproperty(θ_obs, :jitter) : 0.0

    if hasproperty(like.table, :pa) && hasproperty(like.table, :sep)
        σ_sep = like.table.σ_sep 
        σ_pa = like.table.σ_pa

        # Convert ra/dec model to pa/sep
        sep = hypot.(ra_model, dec_model) ./ platescale
        pa = atan.(ra_model, dec_model) .- northangle
        
        if hasproperty(like.table, :cov)
            astrometry_table = Table(;epoch, sep, pa, σ_sep, σ_pa, like.table.cov)
        else
            astrometry_table = Table(;epoch, sep, pa, σ_sep, σ_pa)
        end

        if add_noise
            astrometry_table.sep .+= randn.() .* σ_sep
            astrometry_table.pa .+= randn.() .* σ_pa
        end
    else
        σ_ra = like.table.σ_ra 
        σ_dec = like.table.σ_dec

        # Convert ra/dec model to observed ra/dec with corrections
        sep_model = hypot.(ra_model, dec_model)
        pa_model = atan.(ra_model, dec_model)
        
        # Apply corrections
        pa_corrected = pa_model .- northangle
        sep_corrected = sep_model ./ platescale
        
        ra = sep_corrected .* sin.(pa_corrected)
        dec = sep_corrected .* cos.(pa_corrected)
        
        if hasproperty(like.table, :cov)
            astrometry_table = Table(;epoch, ra, dec, σ_ra, σ_dec, like.table.cov)
        else
            astrometry_table = Table(;epoch, ra, dec, σ_ra, σ_dec)
        end

        if add_noise
            astrometry_table.ra .+= randn.() .* hypot.(σ_ra, jitter)
            astrometry_table.dec .+= randn.() .* hypot.(σ_dec, jitter)
        end
        display(astrometry_table)
    end

    return PlanetRelAstromLikelihood(astrometry_table; like.name)
end
