module DirectDetections
using ComponentArrays
using Distributions: mode, logpdf

import KissMCMC

using MCMCChains
using NamedTupleTools
using DirectImages: lookup_coord
using DirectOrbits
using Base.Threads: @threads
using StaticArrays


# This is a straight forward implementation that unfortunately is not type stable.
# This is because we are looping over a heterogeneous tuple
function make_ln_prior(priors)
    function ln_prior(params)
        lp = zero(first(params))
        for i in eachindex(params)
            pd = priors[i]
            param = params[i]
            lp += logpdf(pd, param)
        end
        return lp 
    end
    return ln_prior
end

    
function make_ln_like(priors, images, contrasts, times, platescale)
    # props is tuple of symbols
    # static is named tuple of values


    if !(all(size(images)  == size(times) == size(contrasts) == size(planet.epochs) for planet in priors.planets))
        error("All values must have the same length")
    end

    fi = MVector{length(images)}(zeros(length(images)))

    # Return our
    return function ln_like(θ)

        # The ln likelihood:
        ll = 0.0
        for (priors_planet, θ_planet) in zip(priors.planets, θ.planets)
            for i in eachindex(priors_planet.epochs)

                # TODO: this merging needs to be worked out at compile time, or at least when building the function!
                # TODO: see ComponentArrays.label2index
                # We already know the layout, can even just look up by index.                # Merge the three levels together. This gives us the deepest nested value for any given variable.
                θ_planet_epoch = θ_planet.epochs[i]

                fi[i] = θ_planet_epoch.f

                # @time elements = KeplerianElements(merge(NamedTuple(θ), NamedTuple(θ_planet), NamedTuple(θ_planet_epoch)))
                elements = KeplerianElements((;θ.μ, θ.plx, θ_planet.i, θ_planet.Ω, θ_planet.ω, θ_planet.e, θ_planet.τ, θ_planet.a))
                f = θ_planet_epoch.f

                image = images[i]
                contrast = contrasts[i]
                t = times[i]

                # Calculate position at this epoch
                ra, dec = kep2cart(elements, t)
                x = -ra
                y = dec

                # Get the photometry in this image at that location
                # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
                f̃ₓ = lookup_coord(image, x,y, platescale)
                
                # Find the uncertainty in that photometry value (i.e. the contrast)
                r = √(x^2 + y^2)
                σₓ = contrast(r/platescale)

                # When we get a position that falls outside of our available
                # data (e.g. under the coronagraph) we cannot say anything
                # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
                # of zero.
                if !isfinite(σₓ) || !isfinite(f̃ₓ)
                    continue
                end

                # Ruffio et al 2017, eqn 31
                ll += -1/(2σₓ^2) * (f^2 - 2f*f̃ₓ)

            end

            # Spread in flux between epochs
            # ll += 
            # *exp(-0.5*(K0/L0 - Model_ratio)^2/model_spread^2)

            # exp(-0.5*sum[(Ki-K0)^2]/(0.1*<Ki>)^2)*exp(-0.5*sum[(Li-L0)^2]/(0.1*<Li>)^2)*exp(-0.5*(K0/L0 - Model_ratio)^2/model_spread^2)

            # if mean(fi) != 0
                ll += -1/2 * sum(
                    fi .- θ_planet.f
                ) / (0.1 * mean(fi)^2)
            # end
            # @show fi θ_planet.f
        end

        # At this point, a NaN or Inf log-likelihood implies
        # an error in preparing the inputs or in this code.
        if !isfinite(ll)
            # @show r σₓ f f̃ₓ l
            error("Non-finite log-likelihood encountered")
        end
        return ll
    end
end



function make_ln_post(priors, images, contrasts, times, platescale)
    ln_prior = make_ln_prior(priors)
    ln_like = make_ln_like(priors, images, contrasts, times, platescale)
    ln_post(params) = ln_prior(params) + ln_like(params)
    return ln_post
end


function mcmc(
    priors, images, contrasts, times;
    platescale,
    burnin,
    numwalkers=10,
    thinning = 1,
    numsamples_perwalker,
    squash=true
    )
    # column_names = string.(collect(keys(priors)))
    column_names = ComponentArrays.labels(priors)

    ln_post = make_ln_post(priors, images, contrasts, times, platescale)

    
    @info "Finding starting point"
    # TODO: kissmcmc has a better method for creating the ball and rejecting some starting points
    initial_walkers = find_starting_walkers(ln_post, priors, numwalkers)

    # Convert the initial walkers into static arrays for stack allocation.
    # This messy line should have no impact on the semantics of the code.
    initial_walkers_static = [
        ComponentVector{SVector{length(cv)}}(;NamedTuple(cv)...)
        for cv in initial_walkers
    ]
    # initial_walkers = SVector{length(priors),Float64}.(initial_walkers)

    thetase, _accept_ratioe = KissMCMC.emcee(
        ln_post,
        initial_walkers_static;
        nburnin=burnin*numwalkers,
        use_progress_meter=true,
        nthin=thinning,
        niter=numsamples_perwalker*numwalkers
    )

    if squash
        thetase′, _ = KissMCMC.squash_walkers(thetase, _accept_ratioe)
        reinterptted = reinterpret(reshape, eltype(thetase′), thetase′)
    else
        # We can reinterpret the vector of SVectors as a matrix directly without copying!
        # This can save massive amounts of memory and time on large changes
        reinterptted = cat(
            [transpose(reinterpret(reshape, eltype(first(θ)), θ)) for θ in thetase]...,
            dims=3
        )
    end

    return Chains(reinterptted, column_names)
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
    # initial_walkers = mapreduce(hcat, 1:numwalkers) do _
    initial_walkers = map(1:numwalkers) do _
        return initial_position = find_starting_point(ln_post, priors)
        # This used to intiialize the walkers in a Gaussian ball around the MAP.
        # But now we just draw starting points randomly from the priors, this isn't needed.
        # @showprogress "Finding initial positions" map(eachindex(initial_position)) do i
        #     p = NaN
        #     while !(minimum(priors[i]) < p < maximum(priors[i]))
        #         p = initial_position[i] + 0.01randn()*initial_position[i]
        #     end
        #     p
        # end
    end
    return initial_walkers
end



include("analysis.jl")
end
