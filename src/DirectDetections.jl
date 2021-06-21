module DirectDetections
using Base: Symbol
using StaticArrays: getproperty
using ComponentArrays: haskey
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

include("types.jl")


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

    ex = :(function (params)
        lp = zero(first(params))
        $(body...)
        return lp
    end)

    ln_prior = @RuntimeGeneratedFunction(ex)
    return ln_prior
end

function ln_like_phot(phot_observations, model_interpolator, elements, θ_planet, θ_band)
    ll = 0.0
    for obs in phot_observations


        # Calculate position at this epoch
        ra, dec = kep2cart(elements, obs.epoch)
        x = -ra
        y = dec

        # Get the photometry in this image at that location
        # Note in the following equations, subscript x (ₓ) represents the current position (both x and y)
        f̃ₓ = lookup_coord(obs.image, (x, y), obs.platescale)

        # Find the uncertainty in that photometry value (i.e. the contrast)
        r = √(x^2 + y^2)
        σₓ = obs.contrast(r / obs.platescale)

        # When we get a position that falls outside of our available
        # data (e.g. under the coronagraph) we cannot say anything
        # about the likelihood. This is equivalent to σₓ→∞ or log likelihood 
        # of zero.
        if !isfinite(σₓ) || !isfinite(f̃ₓ)
            continue
        end

        # Ruffio et al 2017, eqn 31
        # ll += -1/(2σₓ^2) * (θ_epoch_f^2 - 2θ_epoch_f*f̃ₓ)
        # ll += -1/(2σₓ^2) * (θ_band^2 - 2θ_band*f̃ₓ)

        σₓ² = σₓ^2
        ll += -1 / (2σₓ²) * (θ_band^2 - 2θ_band * f̃ₓ) - log(sqrt(2π * σₓ²))
    end

    # # Connect the flux at each epoch to an overall flux in this band for this planet
    # # fᵢ = θ_band.epochs
    # # ll += -1/2 * sum(
    # #     (fᵢ .- θ_band.f).^2
    # # ) / (θ_band.σ_f² * mean(fᵢ)^2)

    # And connect that flux to a modelled Teff and mass
    # f_model = model_interpolator(θ_planet.Teff, θ_planet.mass)
    # ll += -1/2 * (f_model - θ_band)^2 /  (θ_planet.σ_f_model² * f_model^2)

    return ll
end


function ln_like_astrom(elements, observations)

    ll = 0.0
    for obs in observations
        x, y = kep2cart(elements, obs.epoch)
        residx = obs.ra - x
        residy = obs.dec - y
        σ²x = obs.σ_ra^2
        σ²y = obs.σ_dec^2
        χ²x = -0.5residx^2 / σ²x - log(sqrt(2π * σ²x))
        χ²y = -0.5residy^2 / σ²y - log(sqrt(2π * σ²y))
        ll += χ²x + χ²y
    end

    return ll
end


function ln_like_propmotionanom(θ, elements, observations)
    # raoffel = t -> raoff(elements, t)
    # decoffel = t -> decoff(elements, t)

    ll = 0.0
    for obs in observations

        pmra_star = 0.0
        pmdec_star = 0.0
        
        # The model can support multiple planets
        for θ_planet in θ.planets

            # pmra_planet = -ForwardDiff.derivative(raoffel, obs.ra_epoch) * 365.25 # mas/day -> mas/yr
            # pmdec_planet = -ForwardDiff.derivative(decoffel, obs.dec_epoch) * 365.25 # mas/day -> mas/yr
            # pmra_star = pmra_planet .* planet_mass ./ elements.μ
            # pmdec_star = pmdec_planet .* planet_mass ./ elements.μ

            # RA and dec epochs are usually slightly different
            pmra_star += propmotionanom(elements, obs.ra_epoch, elements.μ, θ_planet.mass)[1]
            pmdec_star += propmotionanom(elements, obs.dec_epoch, elements.μ, θ_planet.mass)[2]
        end

        residx = pmra_star - obs.pmra
        residy = pmdec_star - obs.pmdec
        σ²x = obs.σ_pmra^2
        σ²y = obs.σ_pmdec^2
        χ²x = -0.5residx^2 / σ²x - log(sqrt(2π * σ²x))
        χ²y = -0.5residy^2 / σ²y - log(sqrt(2π * σ²y))

        ll += χ²x + χ²y
    end

    # return 0.0
    return ll
end

function make_ln_like(data, interpolators)

    # Prepare photometric likelihood functions for each band

    # We use RuntimeGeneratedFunctions to unroll the loop
    # over the different bands directly

    phot_likes = Expr[]
    if haskey(data, :phot)
        for band in keys(data.phot)
            ex = :(
                ll += ln_like_phot(
                    data.phot.$band,
                    interpolators.$band,
                    elements,
                    θ_planet,
                    θ_planet.phot.$band,
                )
            )
            push!(phot_likes, ex)
        end
    end

    if haskey(data, :astrom) && length(data.astrom) > 0
        astrom_like = :(ll += ln_like_astrom(elements, $(data.astrom)))
    else
        astrom_like = nothing
    end


    if haskey(data, :propmotionanom)
        propmotionanom_like = :(ll += ln_like_propmotionanom(θ, elements, $(data.propmotionanom)))
    else
        propmotionanom_like = nothing
    end


    ex = :(
        function (data, interpolators, θ)

            # The ln likelihood:
            ll = 0.0

            # The model can support multiple planets
            for θ_planet in θ.planets

                elements = KeplerianElements((;
                    θ_planet.μ,
                    θ_planet.plx,
                    θ_planet.i,
                    θ_planet.Ω,
                    θ_planet.ω,
                    θ_planet.e,
                    θ_planet.τ,
                    θ_planet.a,
                ))

                # We can have observations from multiple bands
                $(phot_likes...)

                # TODO: RV likelihood

                $astrom_like
                $propmotionanom_like

            end

            # At this point, a NaN or Inf log-likelihood implies
            # an error in preparing the datas or in this code.
            # if !isfinite(ll)
            #     error("Non-finite log-likelihood encountered")
            # end
            return ll
        end
    )

    return @RuntimeGeneratedFunction(ex)
end


function mcmc(
    priors,
    data;
    burnin,
    numwalkers,
    numsamples_perwalker,
    thinning = 1,
    squash = true,
)
    # column_names = ComponentArrays.labels(priors)

    # Prepare interpolators for any different bands we want to model
    if haskey(data, :phot)
        bands = unique(reduce(vcat, [collect(keys(planet.phot)) for planet in priors.planets]))
        interpolators = namedtuple(bands, [sonora_interpolator_grid(band) for band in bands])
    else
        interpolators = nothing
    end


    ln_prior = make_ln_prior(priors)

    ln_like_0 = make_ln_like(data, interpolators)
    ln_like(θ) = ln_like_0(data, interpolators, θ)


    ln_post(θ) = ln_prior(θ) + ln_like(θ)

    @info "Finding starting point"
    initial_walkers = find_starting_walkers(ln_post, priors, numwalkers)

    # Convert the initial walkers into static arrays for stack allocation.
    # This messy line should have no impact on the semantics of the code.
    initial_walkers_static = [
        ComponentVector{SVector{length(cv)}}(; NamedTuple(cv)...) for cv in initial_walkers
    ]

    # Run the MCMC
    thetase, _accept_ratioe = KissMCMC.emcee(
        ln_post,
        initial_walkers_static;
        nburnin = burnin * numwalkers,
        use_progress_meter = true,
        nthin = thinning,
        niter = numsamples_perwalker * numwalkers,
    )

    # Convert the output into an MCMCChains.Chain.
    # Use reinterpret to avoid re-allocating all that memory
    if squash
        thetase′, _ = KissMCMC.squash_walkers(thetase, _accept_ratioe)
        reinterptted = reinterpret(reshape, eltype(first(thetase′)), thetase′)
        chains = ComponentArray(collect(eachrow(reinterptted)), getaxes(thetase′[1]))
    else
        # We can reinterpret the vector of SVectors as a matrix directly without copying!
        # This can save massive amounts of memory and time on large changes
        reinterptted =
            cat([reinterpret(reshape, eltype(first(θ)), θ) for θ in thetase]..., dims = 3)
        chains = ComponentArray(
            collect(eachslice(reinterptted, dims = 1)),
            getaxes(thetase[1][1]),
        )
    end

    # return Chains(reinterptted, column_names)
    return chains
end


using AdvancedHMC
using NamedTupleTools
using ForwardDiff
using Logging
function hmc(
    priors,
    data;
    numwalkers=1,
    burnin,
    numsamples_perwalker,
)
    
    # Choose parameter dimensionality and initial parameter value
    initial_θ_0 = rand.(priors)
    D = length(initial_θ_0)

    # Define the target distribution
    bands = unique(reduce(vcat, [collect(keys(planet.phot)) for planet in priors.planets]))
    interpolators = namedtuple(bands, [DirectDetections.sonora_interpolator_grid(band) for band in bands])
    ln_prior = DirectDetections.make_ln_prior(priors)
    ln_like_0 = DirectDetections.make_ln_like(data, interpolators)
    ln_like(θ) = ln_like_0(data.phot, interpolators, θ)

    ax = getaxes(initial_θ_0)
    ℓπ = let ax=ax
        function (θ)
            θ_cv = ComponentArray(θ, ax)
            ln_prior(θ_cv) + ln_like(θ_cv)
        end
    end

    chains = []
    stats = []
    Threads.@threads for _ in 1:numwalkers
        initial_θ = rand.(priors)

        # Define a Hamiltonian system
        metric = DiagEuclideanMetric(D)
        # metric = DenseEuclideanMetric(D)
        hamiltonian = Hamiltonian(metric, ℓπ, ForwardDiff)

        # Define a leapfrog solver, with initial step size chosen heuristically
        initial_ϵ = find_good_stepsize(hamiltonian, getdata(initial_θ))
        # initial_ϵ = 0.02
        integrator = Leapfrog(initial_ϵ)

        proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator, 10)
        adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator)) #0.08

        logger = SimpleLogger(stdout, Logging.Error)
        samples, stat = with_logger(logger) do
            sample(hamiltonian, proposal, getdata(initial_θ), numsamples_perwalker, adaptor, burnin; progress=(numwalkers==1), drop_warmup=true)
        end

        sample_grid = reduce(hcat, samples);
        chain = ComponentArray(collect(eachrow(sample_grid)), ax)
        
        push!(chains,chain)
        push!(stats,stat)
    end
    return chains, stats
end


# using Optim, ForwardDiff

function find_starting_point(ln_post, priors)
    θ₀ = rand.(priors)
    i = 0
    while !isfinite(ln_post(θ₀))
        i += 1
        θ₀ = rand.(priors)
        if i > 1000
            error(
                "Could not find a starting point in the posterior that is finite by drawing from the priors after 1000 attempts",
            )
        end
    end
    return θ₀

    # goal(θ) = -ln_post(θ)

    # m = optimize(goal, θ₀, BFGS(), Optim.Options(show_trace=true,x_tol=-1,g_tol=1e-1,f_tol=-1); autodiff = :forward)
    # # m = optimize(goal, θ₀, Newton(), Optim.Options(show_trace=true,); autodiff = :forward)
    # # m = optimize(goal, θ₀, NelderMead(), Optim.Options(show_trace=true,); autodiff = :forward)
    # # m = optimize(goal, θ₀, Optim.SimulatedAnnealing(), Optim.Options(show_trace=true,iterations=1_000_000); autodiff = :forward)
    # # m = optimize(goal, θ₀, Optim.ParticleSwarm(), Optim.Options(show_trace=true,iterations=100_000); autodiff = :forward)

    # display(m)
    # return Optim.minimizer(m)

end



# Start walkers in a gaussian ball around the MLE, while ensuring we don't
# step outside the ranges defined by the priors
function find_starting_walkers(ln_post, priors, numwalkers)
    initial_walkers = map(1:numwalkers) do i
        return initial_position = find_starting_point(ln_post, priors)
    end
    #     initial_position
    # end
    #     # # This used to intiialize the walkers in a Gaussian ball around the MAP.
    #     # # But now we just draw starting points randomly from the priors, this isn't needed.
    #     # map(eachindex(initial_position)) do i
    #     #     p = NaN
    #     #     while !(minimum(priors[i]) < p < maximum(priors[i]))
    #     #         p = initial_position[i] + 0.01randn()*initial_position[i]
    #     #     end
    #     #     p
    #     # end
    #     # initial_position .* (1 .+ randn(length(initial_position)))
    # end
    return initial_walkers
end


include("analysis.jl")

function __init__()
    init_models()
    init_plots()
    return
end

end
