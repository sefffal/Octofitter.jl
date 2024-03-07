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

# Valid for reference model only
function Pigeons.sample_iid!(model_reference::Octofitter.LogDensityModel, replica, shared)
    # This could in theory be done without any array allocations
    θ = sample_priors(replica.rng, model_reference.system)
    θ_t = model_reference.link(θ)
    replica.state .= θ_t
end

function Pigeons.default_reference(target::Octofitter.LogDensityModel)
    reference_sys = prior_only_model(target.system)
    # Note we could run into issues if their priors aren't well handled by the default
    # autodiff backend
    reference = Octofitter.LogDensityModel(reference_sys)
    return reference
end


"""
octofit_pigeons(model; nrounds, n_chains=[auto])

Use Pigeons.jl to sample from intractable posterior distributions.

```julia
model = Octofitter.LogDensityModel(System, autodiff=:Enzyme, verbosity=4)
chain, pt = octofit_pigeons(model)
```
"""
Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    target::Octofitter.LogDensityModel;
    n_rounds::Int,
    n_chains::Int=cld(8,Threads.nthreads())*Threads.nthreads(),
    pigeons_kw...
)
    @nospecialize
    inputs = Pigeons.Inputs(;
        target,
        # explorer = AutoMALA(default_autodiff_backend = :ForwardDiff),
        record = [traces; round_trip; record_default()],
        multithreaded=true,
        show_report=true,
        n_rounds,
        n_chains,
        pigeons_kw...
    )
    return octofit_pigeons(inputs)
end

Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    pt::Pigeons.PT
)
    @nospecialize

    start_time = time()
    pt = pigeons(pt)
    stop_time = time()

    mcmcchains = Chains(pt.inputs.target, pt)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
        )
    )
    return (;chain=mcmcchains_with_info, pt)
end
Base.@nospecializeinfer function Octofitter.octofit_pigeons(
    inputs::Pigeons.Inputs
)
    @nospecialize

    start_time = time()
    pt = pigeons(inputs)
    stop_time = time()

    mcmcchains = Chains(inputs.target, pt)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
        )
    )
    return (;chain=mcmcchains_with_info, pt)
end
function MCMCChains.Chains(model::Octofitter.LogDensityModel, pt::Pigeons.PT)
    ln_like = Octofitter.make_ln_like(model.system, model.arr2nt(Octofitter.sample_priors(model.system)))

    # Resolve the array back into the nested named tuple structure used internally.
    # Augment with some internal fields
    chain_res = map(get_sample(pt)) do sample 
        θ_t = @view(sample[begin:begin+model.D-1])
        logpost = sample[model.D+1]
        # Map the variables back to the constrained domain and reconstruct the parameter
        # named tuple structure.
        θ = model.invlink(θ_t)
        resolved_namedtuple = model.arr2nt(θ)
        # Add log posterior, tree depth, and numerical error reported by
        # the sampler.
        # Also recompute the log-likelihood and add that too.
        loglike = ln_like(model.system, resolved_namedtuple)
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
            :loglike,
            :logpost
        ])
    )
    return mcmcchains
end

end