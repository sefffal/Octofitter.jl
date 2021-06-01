module DirectDetections
using ComponentArrays
using Distributions: mode, logpdf

import KissMCMC

using Statistics
# using MCMCChains
using NamedTupleTools
using DirectImages: lookup_coord
using DirectOrbits
using Base.Threads: @threads
using StaticArrays

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# We use a sonora model grid to tie fluxes to physical properties
include("sonora.jl")


# This is a straight forward implementation that unfortunately is not type stable.
# This is because we are looping over a heterogeneous container
# function make_ln_prior(priors)
#     return function ln_prior(params)
#         lp = zero(first(params))
#         for i in eachindex(params)
#             pd = priors[i]
#             param = params[i]
#             lp += logpdf(pd, param)
#         end
#         return lp 
#     end
# end

# This implementation is ~5x faster for the same result.
# It uses metaprogramming to unroll the loop over the different
# prior types. Note that this also ensures the priors can't be modified
# after building the ln_prior function.
function make_ln_prior(θ)

    body = Expr[]
    for i in eachindex(θ)
        pd = θ[i]
        ex = :(lp += logpdf($pd, params[$i]))
        push!(body, ex)
    end

    ex = :(function(params)
        lp = zero(first(params))
        $(body...)
        return lp
    end)
    
    ln_prior =  @RuntimeGeneratedFunction(ex)
    return ln_prior
end


function ln_photometric_likelihood(θ_planet, θ_band, phot_observations, model_interpolator)
    ll = 0.0
    elements = KeplerianElements((;θ_planet.μ, θ_planet.plx, θ_planet.i, θ_planet.Ω, θ_planet.ω, θ_planet.e, θ_planet.τ, θ_planet.a))
    # for (θ_epoch_f, obs) in zip(θ_band.epochs, phot_observations)
    for obs in phot_observations

        # Calculate position at this epoch
        ra, dec = kep2cart(elements, obs.epoch)
        x = -ra
        y = dec

        # Get the photometry in this image at that location
        # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
        f̃ₓ = lookup_coord(obs.image, x,y, obs.platescale)
        
        # Find the uncertainty in that photometry value (i.e. the contrast)
        r = √(x^2 + y^2)
        σₓ = obs.contrast(r/obs.platescale)

        # When we get a position that falls outside of our available
        # data (e.g. under the coronagraph) we cannot say anything
        # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
        # of zero.
        if !isfinite(σₓ) || !isfinite(f̃ₓ)
            continue
        end

        # Ruffio et al 2017, eqn 31
        # ll += -1/(2σₓ^2) * (θ_epoch_f^2 - 2θ_epoch_f*f̃ₓ)
        ll += -1/(2σₓ^2) * (θ_band.f^2 - 2θ_band.f*f̃ₓ)
    end

    # Connect the flux at each epoch to an overall flux in this band for this planet
    # fᵢ = θ_band.epochs
    # ll += -1/2 * sum(
    #     (fᵢ .- θ_band.f).^2
    # ) / (θ_band.σ_f² * mean(fᵢ)^2)

    # And connect that flux to a modelled Teff and mass
    # @show θ_planet.Teff θ_planet.mass
    f_model = model_interpolator(θ_planet.Teff, θ_planet.mass)
    # @show f_model
    ll += -1/2 * (f_model - θ_band.f)^2 /  (θ_band.σ_f_model² * f_model^2)

    return ll
end



function make_ln_photometric_likelihood(phot_observations, model_interpolator)
    
    :(function ln_photometric_likelihood(θ_planet, θ_band)
        ll = 0.0
        elements = KeplerianElements((;θ_planet.μ, θ_planet.plx, θ_planet.i, θ_planet.Ω, θ_planet.ω, θ_planet.e, θ_planet.τ, θ_planet.a))
        # for (θ_epoch_f, obs) in zip(θ_band.epochs, phot_observations)
        for obs in phot_observations

            # Calculate position at this epoch
            ra, dec = kep2cart(elements, obs.epoch)
            x = -ra
            y = dec

            # Get the photometry in this image at that location
            # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
            f̃ₓ = lookup_coord(obs.image, x,y, obs.platescale)
            
            # Find the uncertainty in that photometry value (i.e. the contrast)
            r = √(x^2 + y^2)
            σₓ = obs.contrast(r/obs.platescale)

            # When we get a position that falls outside of our available
            # data (e.g. under the coronagraph) we cannot say anything
            # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
            # of zero.
            if !isfinite(σₓ) || !isfinite(f̃ₓ)
                continue
            end

            # Ruffio et al 2017, eqn 31
            # ll += -1/(2σₓ^2) * (θ_epoch_f^2 - 2θ_epoch_f*f̃ₓ)
            ll += -1/(2σₓ^2) * (θ_band.f^2 - 2θ_band.f*f̃ₓ)
        end

        # Connect the flux at each epoch to an overall flux in this band for this planet
        # fᵢ = θ_band.epochs
        # ll += -1/2 * sum(
        #     (fᵢ .- θ_band.f).^2
        # ) / (θ_band.σ_f² * mean(fᵢ)^2)

        # And connect that flux to a modelled Teff and mass
        f_model = model_interpolator(θ_planet.Teff, θ_planet.mass)
        ll += -1/2 * (f_model - θ_band.f)^2 /  (θ_band.σ_f_model² * f_model^2)

        return ll
    end)

end


function ln_astrometric_likelihood(observations, θ_planet)

    elements = KeplerianElements((;θ_planet.μ, θ_planet.plx, θ_planet.i, θ_planet.Ω, θ_planet.ω, θ_planet.e, θ_planet.τ, θ_planet.a))

    ll = 0.0
    for obs in observations
        x, y = kep2cart(elements, obs.epoch)
        residx = obs.ra - x
        residy = obs.dec - y
        σ²x = obs.σ_ra^2
        σ²y = obs.σ_dec^2
        χ²x = -0.5residx^2 / σ²x - log(sqrt(2π*σ²x))
        χ²y = -0.5residy^2 / σ²y - log(sqrt(2π*σ²y))
        ll += χ²x + χ²y
    end

    return ll
end

    
function ln_likelihood(data, interpolators, θ)

    # The ln likelihood:
    ll = 0.0

    # The model can support multiple planets
    for θ_planet in θ.planets
        # We can have observations from multiple bands
        for band in keys(data.phot)
            ll += ln_photometric_likelihood(θ_planet, θ_planet.phot[band], data.phot[band], interpolators[band])
        end

        # TODO: RV likelihood
        # TODO: Astrom. likelihood

        ll += ln_astrometric_likelihood(data.astrom, θ_planet)

    end


    # At this point, a NaN or Inf log-likelihood implies
    # an error in preparing the datas or in this code.
    # if !isfinite(ll)
    #     error("Non-finite log-likelihood encountered")
    # end
    return ll
end

function mcmc(
    priors, data;
    burnin,
    numwalkers,
    numsamples_perwalker,
    thinning=1,
    squash=true
)
    # column_names = ComponentArrays.labels(priors)

    # Prepare interpolators for any different bands we want to model
    bands = unique(reduce(vcat, (keys(planet.phot) for planet in priors.planets)))
    
    interpolators = namedtuple(bands, [sonora_interpolator_grid(band) for band in bands])

    ln_prior = make_ln_prior(priors)
    ln_post(θ) = ln_prior(θ) + ln_likelihood(data, interpolators, θ)

    @info "Finding starting point"
    initial_walkers = find_starting_walkers(ln_post, priors, numwalkers)

    # Convert the initial walkers into static arrays for stack allocation.
    # This messy line should have no impact on the semantics of the code.
    initial_walkers_static = [
        ComponentVector{SVector{length(cv)}}(;NamedTuple(cv)...)
        for cv in initial_walkers
    ]

    # Run the MCMC
    thetase, _accept_ratioe = KissMCMC.emcee(
        ln_post,
        initial_walkers_static;
        nburnin=burnin*numwalkers,
        use_progress_meter=true,
        nthin=thinning,
        niter=numsamples_perwalker*numwalkers
    )

    # Convert the output into an MCMCChains.Chain.
    # Use reinterpret to avoid re-allocating all that memory
    if squash
        thetase′, _ = KissMCMC.squash_walkers(thetase, _accept_ratioe)
        reinterptted = reinterpret(reshape, eltype(first(thetase′)), thetase′);
        chains = ComponentArray(collect(eachrow(reinterptted)), getaxes(thetase′[1]));
    else
        # We can reinterpret the vector of SVectors as a matrix directly without copying!
        # This can save massive amounts of memory and time on large changes
        reinterptted = cat(
            [reinterpret(reshape, eltype(first(θ)), θ) for θ in thetase]...,
            dims=3
        )
        chains = ComponentArray(collect(eachslice(reinterptted,dims=1)), getaxes(thetase[1][1]))
    end

    # return Chains(reinterptted, column_names)
    return chains
end



function find_starting_point(ln_post, priors)
    initial_guess = rand.(priors)
    i = 0
    while !isfinite(ln_post(initial_guess))
        i+=1
        initial_guess = rand.(priors)
        if i > 1000
            error("Could not find a starting point in the posterior that is finite by drawing from the priors after 1000 attempts")
        end
    end
    return initial_guess
end



# Start walkers in a gaussian ball around the MLE, while ensuring we don't
# step outside the ranges defined by the priors
function find_starting_walkers(ln_post, priors, numwalkers)
    # initial_position = find_starting_point(ln_post, priors)
    initial_walkers = map(1:numwalkers) do i
        initial_position = find_starting_point(ln_post, priors)
        # # This used to intiialize the walkers in a Gaussian ball around the MAP.
        # # But now we just draw starting points randomly from the priors, this isn't needed.
        # map(eachindex(initial_position)) do i
        #     p = NaN
        #     while !(minimum(priors[i]) < p < maximum(priors[i]))
        #         p = initial_position[i] + 0.01randn()*initial_position[i]
        #     end
        #     p
        # end
        # initial_position .* (1 .+ randn(length(initial_position)))
    end
    return initial_walkers
end



# include("analysis.jl")
end
