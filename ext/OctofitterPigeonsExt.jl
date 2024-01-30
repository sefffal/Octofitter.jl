module OctofitterPigeonsExt
using Random
using Octofitter
using Pigeons
using MCMCChains

function (model::Octofitter.LogDensityModel)(θ)
    return model.ℓπcallback(θ)
end
function Pigeons.initialization(model::Octofitter.LogDensityModel, rng::AbstractRNG, ::Int)
    initial_θ, mapv = Octofitter.guess_starting_position(rng,model.system,1)
    initial_θ_t = model.link(initial_θ)
    return initial_θ_t
end


"""
```octofit_pigeons(
    target=model_target,
    reference=model_reference;
    pigeons_kw...
)`

Use Pigeons.jl to sample from intractable posterior distributions.
Requires two model: a target and reference. 
The reference model is typically the prior only model. You can construct this 
by removing all data from the system, or by passing inverse_temp=0.0 when constructing
the log density model.

```julia
model_target = Octofitter.LogDensityModel(System, autodiff=:Enzyme, verbosity=4)
model_reference = Octofitter.LogDensityModel(System, inverse_temp=0.0, autodiff=:Enzyme, verbosity=4)
chain, pt = octofit_pigeons(target=model_target, reference=model_reference)
```
"""
function Octofitter.octofit_pigeons(;
    target,
    reference,
    n_chains=cld(16,Threads.nthreads())*Threads.nthreads(),
    pigeons_kw...
)
    start_time = time()
    pt = pigeons(;
            target,
            reference,
            # explorer = AutoMALA(default_autodiff_backend = :ForwardDiff),
            record = [traces; round_trip; record_default()],
            multithreaded=true,
            show_report=true,
            n_rounds=12,
            n_chains,
            pigeons_kw...
    )
    stop_time = time()

    ln_like = Octofitter.make_ln_like(target.system, target.arr2nt(Octofitter.sample_priors(target.system)))

    # Resolve the array back into the nested named tuple structure used internally.
    # Augment with some internal fields
    chain_res = map(get_sample(pt)) do θ_t
        # Map the variables back to the constrained domain and reconstruct the parameter
        # named tuple structure.
        θ = target.invlink(θ_t)
        resolved_namedtuple = target.arr2nt(θ)
        # Add log posterior, tree depth, and numerical error reported by
        # the sampler.
        # Also recompute the log-likelihood and add that too.
        loglike = ln_like(target.system, resolved_namedtuple)
        logpost = target.ℓπcallback(θ_t)
        return merge((;
            loglike,
            logpost,
        ), resolved_namedtuple)
    end
    # Then finally flatten and convert into an MCMCChain object / table.
    # Mark the posterior, likelihood, numerical error flag, and tree depth as internal
    mcmcchains = Octofitter.result2mcmcchain(
        chain_res,
        Dict(:internals => [
            :loglike
            :logpost
        ])
    )
    # Concatenate the log posteriors and make them the same shape as the chains (N_iters,N_vars,N_chains)
    # logposts_mat = reduce(hcat, logposts)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
        )
    )
    return (;chain=mcmcchains_with_info, pt)
end


end