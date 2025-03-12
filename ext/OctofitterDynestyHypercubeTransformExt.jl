module OctofitterDynestyHypercubeTransformExt

using Octofitter
using Dynesty
using HypercubeTransform
using Random
using StatsBase
using MCMCChains


function Dynesty.dysample(model::Octofitter.LogDensityModel, args...; kwargs...)
    
    priors = Octofitter._list_priors(model.system)
    priors_tuple = Tuple(_ for _ in priors)
    hc = ascube(priors_tuple)
    function pri_tr(x)
        θ_natural = HypercubeTransform.transform(hc, x)
        return θ_natural
    end
    
    # Generate a new likelihood only function to call.
    θ_natural_sample_1 = model.sample_priors(Random.default_rng())
    ln_lik′ = Octofitter.make_ln_like(model.system,model.arr2nt(θ_natural_sample_1))
    # ln_like can now be called just with model.system and θ_natural
    function log_lik(θ_natural)
        # Map from flat array into nested named tuple structure, resolving deterministic variables etc.
        θ_nt = model.arr2nt(θ_natural)
        return ln_lik′(model.system, θ_nt)
    end

    start_time = time()
    res = dysample(log_lik, pri_tr, model.D, args...; kwargs...)
    stop_time = time()

    samples = Dynesty.PythonCall.pyconvert(Array, res["samples"])
    weights = exp.(Dynesty.PythonCall.pyconvert(Vector, res["logwt"] - res["logz"][-1]))
    stats = (logl = Dynesty.PythonCall.pyconvert(Vector, res["logl"]), weights = weights,)
    logz = Dynesty.PythonCall.pyconvert(Float64, res["logz"][-1])
    logzerr = Dynesty.PythonCall.pyconvert(Float64, res["logzerr"][-1])
    inds = sample(1:size(samples)[1], Weights(stats.weights), min(length(samples)÷50, 50_000))
    esamples = @views samples[inds, :]

    ln_prior = Octofitter.make_ln_prior(model.system)
    chain_res = map(eachrow(esamples)) do θ
        # Reconstruct the parameter named tuple structure.
        resolved_namedtuple = model.arr2nt(θ)

        # Add log posterior, tree depth, and numerical error reported by
        # the sampler.
        logprior = ln_prior( θ)
        # Also recompute the log-likelihood and add that too.
        loglike = log_lik(θ)
        logpost = loglike + logprior

        return (;loglike, logpost, resolved_namedtuple...,)
    end
    # Then finally flatten and convert into an MCMCChain object / table.
    # Mark the posterior, likelihood, numerical error flag, and tree depth as internal
    mcmcchains = Octofitter.result2mcmcchain(
        chain_res,
    )
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
            logz,
            logzerr,
            model_name=model.system.name
        )
    )
    return mcmcchains_with_info
    
end

end