using DiffResults
using LinearAlgebra
using Preferences
using Pathfinder
using CovarianceEstimation



# (Prior sampling, chain <-> NamedTuple conversion, and orbit reconstruction
#  live in chains.jl, rewritten for the flat (bodies, observations) model.)

# Fallback when no random number generator is provided (as is usually the case)
Base.@nospecializeinfer function advancedhmc(model::LogDensityModel, target_accept::Number=0.8; kwargs...)
    return advancedhmc(Random.default_rng(), model, target_accept; kwargs...)
end

"""
    octofit(
        [rng::Random.AbstractRNG],
        model::Octofitter.LogDensityModel,
        target_accept::Number = 0.8,
        ensemble::AbstractMCMC.AbstractMCMCEnsemble = MCMCSerial();
        adaptation = 1000,
        iterations = 1000,
        drop_warmup = true,
        max_depth = 12,
        initial_samples = pathfinder ? 500 : 250_000,  # deprecated
        initial_parameters = nothing,                  # deprecated
        step_size = nothing,
        verbosity = 2,
    )

Sample from the posterior defined by `model` using Hamiltonian Monte Carlo with
the No U-Turn Sampler from AdvancedHMC.jl.

!!! warning "`target_accept` is positional, not a keyword"
    It is the third **positional** argument, so

        octofit(model, 0.6, iterations=2000, adaptation=2000)   # correct
        octofit(model, iterations=2000, target_accept=0.6)      # MethodError

    Lower it (towards ~0.5) when the sampler is taking very small steps or
    reporting many divergences on a difficult posterior; raise it (towards
    ~0.95) when divergences persist at the default.

For posteriors with widely separated modes, or with a discrete variable, reach
for [`octofit_pigeons`](@ref) instead — NUTS cannot cross a low-density gap and
cannot move a discrete parameter at all.

See also [`initialize!`](@ref), [`octofit_rejection`](@ref),
[`octofit_pigeons`](@ref).
"""
Base.@nospecializeinfer function octofit(args...; kwargs...)
    return advancedhmc(args...; kwargs...)
end
export octofit


"""
    octofit_rejection(
        [rng::Random.AbstractRNG],
        model::Octofitter.LogDensityModel;
        draws=100_000,
        verbosity=2,
    )

Sample from the posterior defined by `model` using rejection sampling with the
prior as the proposal distribution.

This sampler draws `draws` samples from the prior, evaluates the likelihood at
each point, and accepts each sample with probability proportional to its
likelihood. The accepted samples are independent (no autocorrelation), but the
method can be very inefficient for high-dimensional problems or when the
posterior is much narrower than the prior.

Returns an `MCMCChains.Chains` object, consistent with `octofit` and
`octofit_pigeons`.
"""
function octofit_rejection end
octofit_rejection(model::LogDensityModel; kwargs...) = octofit_rejection(Random.default_rng(), model; kwargs...)
function octofit_rejection(
    rng::AbstractRNG,
    model::LogDensityModel;
    draws::Int=100_000,
    verbosity::Int=2,
)
    start_time = fill(time(), 1)

    # Sample from the prior (proposal distribution)
    verbosity >= 1 && @info "Drawing $draws samples from prior..."
    prior_samples = [model.sample_priors(rng) for _ in 1:draws]

    # Build the likelihood function
    θ_test = model.arr2nt(first(prior_samples))
    # Priors included, deliberately: the proposal is `sample_priors`, which
    # draws only from the declared distributions. A prior-shaped term (a `~`
    # line in a `@variables` block, `UnitLengthPrior`, a stability prior) is not
    # in the proposal, so dropping it from the acceptance weight would drop the
    # constraint entirely.
    ln_like = make_ln_like(model.system, θ_test)
    # …but the reported `loglike` column is the data terms only, matching every
    # other sampler. The rest is reported as `logprior`.
    ln_like_data = make_ln_like(model.system, θ_test; include_priors=false)

    # Evaluate log-likelihood for each prior sample.
    # Use a function barrier so the compiler can specialize on the concrete
    # types of arr2nt, ln_like, and system (these are abstract when accessed
    # through model fields in the outer closure).
    verbosity >= 1 && @info "Evaluating likelihoods..."
    log_likes = _rejection_evaluate_likelihoods(
        model.arr2nt, ln_like, model.system, prior_samples
    )

    # Find the maximum log-likelihood for numerical stability
    max_ll = maximum(log_likes)

    if !isfinite(max_ll)
        error("All $(draws) prior samples produced non-finite log-likelihoods. Check your model and priors.")
    end

    # Accept/reject: accept with probability exp(loglike - max_loglike)
    accepted_indices = Int[]
    for i in 1:draws
        if log_likes[i] == -Inf
            continue
        end
        acceptance_prob = exp(log_likes[i] - max_ll)
        if rand(rng) < acceptance_prob
            push!(accepted_indices, i)
        end
    end

    n_accepted = length(accepted_indices)
    if n_accepted == 0
        error("No samples were accepted out of $(draws) draws. The posterior may be extremely concentrated relative to the prior. Consider increasing `draws` or using a different sampler.")
    end

    acceptance_rate = n_accepted / draws
    if verbosity >= 1
        @info "Rejection sampling complete" draws n_accepted acceptance_rate
    end
    if acceptance_rate < 0.001 && verbosity >= 1
        @warn "Very low acceptance rate ($(round(acceptance_rate*100, sigdigits=2))%). Consider using `octofit` (HMC) for more efficient sampling."
    end

    # Build the chain results in the same format as octofit.
    # Again use a function barrier for the hot path.
    chain_res = _rejection_build_chain(
        model.arr2nt, model.link, model.ℓπcallback, ln_like_data, model.system,
        prior_samples, accepted_indices,
    )

    mcmcchains = Octofitter.result2mcmcchain(
        chain_res,
        Dict(:internals => [
            :loglike,
            :logprior,
            :logpost,
        ])
    )

    stop_time = fill(time(), 1)

    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
            model_name=model.system.name,
            corrections=model.system.corrections,
            observing_geometry=model.system.observing_geometry,
            barycentric_lighttime=model.system.barycentric_lighttime,
            sampler="rejection",
            draws,
            n_accepted,
            acceptance_rate,
        )
    )
    return mcmcchains_with_info
end
export octofit_rejection

# Function barriers for rejection sampling so the compiler can specialize
# on the concrete types of the closures stored in LogDensityModel.
function _rejection_evaluate_likelihoods(arr2nt, ln_like, system, prior_samples)
    log_likes = Vector{Float64}(undef, length(prior_samples))
    for (i, θ) in enumerate(prior_samples)
        resolved = arr2nt(θ)
        ll = ln_like(system, resolved)
        log_likes[i] = isfinite(ll) ? ll : -Inf
    end
    return log_likes
end

function _rejection_build_chain(arr2nt, link, ℓπcallback, ln_like_data, system,
                               prior_samples, accepted_indices)
    return map(accepted_indices) do i
        θ = prior_samples[i]
        resolved_namedtuple = arr2nt(θ)
        loglike = ln_like_data(system, resolved_namedtuple)
        θ_t = link(θ)
        logpost = ℓπcallback(θ_t)
        logprior = logpost - loglike
        return merge((;loglike, logprior, logpost), resolved_namedtuple)
    end
end


# Define some wrapper functions that hide type information
# so that we don't have to recompile pathfinder() with each 
# new model.
# It's worth it for the sampler, but not pathfinder.
struct CallableAny
    func::Function
end
(ca::CallableAny)(args...;kwargs...) = ca.func(args...;kwargs...)

struct LogDensityModelAny
    ldm::LogDensityModel
end
LogDensityProblems.logdensity(ldm_any::LogDensityModelAny, θ) = LogDensityProblems.logdensity(ldm_any.ldm, θ)
LogDensityProblems.logdensity_and_gradient(ldm_any::LogDensityModelAny, θ) = LogDensityProblems.logdensity_and_gradient(ldm_any.ldm, θ)
LogDensityProblems.dimension(ldm_any::LogDensityModelAny) = LogDensityProblems.dimension(ldm_any.ldm)
LogDensityProblems.capabilities(ldm_any::Type{<:LogDensityModelAny}) = LogDensityProblems.capabilities(ldm_any.ldm)



"""
The method signature of Octofitter.hmc is as follows:

    chain = advancedhmc(
        [rng::Random.AbstractRNG],
        model::Octofitter.LogDensityModel
        target_accept::Number=0.8,
        adaptation=1000,
        iterations=1000,
        drop_warmup=true,
        max_depth=12,
    )

Sample from the posterior defined by `model` using Hamiltonian Monte Carlo with 
the No U-Turn Sampler from AdvancedHMC.jl.
"""
Base.@nospecializeinfer function advancedhmc(
    rng::Union{AbstractRNG, AbstractVector{<:AbstractRNG}},
    model::LogDensityModel,
    target_accept::Number=0.8;
    adaptation::Int=1000,
    iterations::Int=1000,
    drop_warmup::Bool=true,
    max_depth::Int=12,
    verbosity::Int=2,
)
    @nospecialize
    if adaptation < 1000
        @warn "At least 1000 steps of adapation are recomended for good sampling"
    end

    # inialize if not already done or set by user
    get_starting_point!!(model)

    # Precondition warmup with the spread of the starting points, when they have
    # one. Often they do not: `initialize!` sets every starting point to a single
    # value whenever it falls back to the optimizer's result — pathfinder erroring
    # out, or its draws being rejected as far worse than the mode — so `cov` comes
    # back exactly zero.
    #
    # There is no information to conditionon in that case, so use a plain identity
    # and let warmup adapt from scratch. The old code reached the same place by a
    # more confusing route: it regularised the zero matrix up a ladder until
    # `1e-8*I` was positive definite, which reads as if it preserved something but
    # is only an identity metric scaled by 1e-8 — and a uniform rescaling of M⁻¹ is
    # absorbed exactly by the step size (for M⁻¹ = c*I the position step goes as
    # ε√c), so it was identity with a misleading constant attached.
    # Identical points do NOT give an exactly zero covariance: `cov` subtracts
    # `mean = sum/n`, and summing n copies of x then dividing does not round back
    # to exactly x, so the deviations come out at the 1e-16 level rather than 0.
    # Testing `iszero` here silently never fires. Ask instead whether the spread
    # is large enough to mean anything: starting points live in the unconstrained
    # space, where O(1) is the natural scale, so a variance below `1e-12` is
    # numerical noise about a single point.
    S = cov(SimpleCovariance(), stack(model.starting_points)')
    local metric = nothing
    if all(isfinite, S) && maximum(diag(S)) > 1e-12
        # There is real spread here. A single rank-deficient direction should not
        # discard what the other directions know, so regularise by the smallest
        # amount that makes the matrix positive definite rather than giving up.
        for diag_eps in [0; 10.0 .^ range(-8, 0)]
            try
                metric = DenseEuclideanMetric(S .+ Diagonal(diag_eps .* ones(model.D)))
                break
            catch err
                continue
            end
        end
    end
    if isnothing(metric)
        verbosity > 1 && @info(
            "Starting points carry no usable covariance (identical points, or a " *
            "non-finite estimate); starting from an identity mass matrix and " *
            "letting warmup adapt it.")
        metric = DenseEuclideanMetric(model.D)
    end
    if verbosity >= 3
        print("Initial mass matrix M⁻¹\n")
        display(metric.M⁻¹)
    end

    initial_θ_t = rand(rng, model.starting_points)
    initial_θ = model.invlink(initial_θ_t)

    if verbosity >= 4
        @info "flat starting point" initial_θ
        @info "transformed flat starting point" initial_θ_t
    end

    verbosity >= 3 && @info "Creating hamiltonian"
    hamiltonian = Hamiltonian(metric, model)
    verbosity >= 3 && @info "Finding good stepsize"
    ϵ = find_good_stepsize(rng, hamiltonian, initial_θ_t)
    verbosity >= 3 && @info "Found initial stepsize" ϵ 

    # Create integrator
    integrator = Leapfrog(ϵ)
    # integrator = JitteredLeapfrog(ϵ, 0.05) # 5% normal distribution on step size to help in areas of high curvature. 

    verbosity >= 3 && @info "Creating kernel"
    # Qualified: PlanetOrbits v1 also exports `Trajectory`, and Octofitter
    # re-exports it.
    kernel = HMCKernel(AdvancedHMC.Trajectory{MultinomialTS}(integrator, GeneralisedNoUTurn(;max_depth)))

    
    verbosity >= 3 && @info "Creating adaptor"
    # if isnothing(result_pf)
    # if metric isa Pathfinder.RankUpdateEuclideanMetric
        # adaptor = StepSizeAdaptor(target_accept, integrator)
    # else
        mma = MassMatrixAdaptor(metric)
        ssa = StepSizeAdaptor(target_accept, integrator)
        adaptor = StanHMCAdaptor(mma, ssa) 
    # end


    # Thread the trajectory solve inside the likelihood (see
    # `_kepsolve_use_threads` in model/codegen.jl): HMC runs one chain, so for
    # a many-epoch model spare threads have nothing better to do — design
    # §10.5's `rv-8000` sits at ~6.0 ms per gradient with the other cores
    # idle. Switched on *here*, after initialization, deliberately:
    # `guess_starting_position` reads the flag to keep its own outer-loop
    # threading from nesting inside a threaded solve, and its 500k prior
    # evaluations parallelize better over draws than within one solve (v1 set
    # the flag before initialization and paid for exactly that). The threshold
    # is 2 chunks' worth of epochs, below which `orbitsolve!` would run serial
    # regardless.
    Octofitter._kepsolve_use_threads[] =
        Threads.nthreads() > 1 &&
        _count_epochs(model.system) >= 2 * PlanetOrbits.MIN_EPOCHS_PER_TASK
    if verbosity >= 2 && Octofitter._kepsolve_use_threads[]
        @info "Likelihood trajectory solve will use $(Threads.nthreads()) threads"
    end

    initial_parameters = initial_θ_t

    verbosity >= 1 && @info "Sampling, beginning with adaptation phase..."
    start_time = fill(time(), 1)
    mc_samples, stats = AdvancedHMC.sample(
        rng,
        hamiltonian,
        kernel,
        initial_parameters,
        iterations+adaptation,
        adaptor,
        adaptation;
        verbose=verbosity>=4,
        progress=verbosity>=1,
        drop_warmup,
    )
    stop_time = fill(time(), 1)
    
    verbosity >= 1 && @info "Sampling compete. Building chains."

    # Rebuild just the likelihood function (should already be compiled anyways).
    # `include_priors=false`: the chain's `loglike` column is the *data* terms
    # only. Prior-shaped observations (a `~` line in a `@variables` block, the
    # `UnitLengthPrior` behind a `UniformCircular`, the dynamical stability
    # priors) reshape the prior, and are reported under `logprior` with the
    # parameter priors where they belong.
    ln_like_data = make_ln_like(model.system, model.arr2nt(initial_θ); include_priors=false)

    # Go through chain and repackage results
    numerical_error = getproperty.(stats, :numerical_error)
    actual_tree_depth = getproperty.(stats, :tree_depth)
     
    mean_accept = mean(getproperty.(stats, :acceptance_rate))
    ratio_divergent_transitions = mean(numerical_error)
    mean_tree_depth = mean(actual_tree_depth)
    max_tree_depth_frac = mean(==(max_depth), actual_tree_depth)

    total_steps = round(Int,sum(getproperty.(stats, :n_steps))*(iterations+adaptation)/iterations)

    elapsed = only(stop_time) - only(start_time)
    verbosity >= 1 && println("""
    Sampling report for chain:
    mean_accept                 = $mean_accept
    ratio_divergent_transitions = $ratio_divergent_transitions
    mean_tree_depth             = $mean_tree_depth
    max_tree_depth_frac         = $max_tree_depth_frac
    total_steps                 = $total_steps
    μs/step (approx.)           = $(round(elapsed/total_steps*1e6,sigdigits=3))\
    """)

    # Report some warnings if sampling did not work well
    if ratio_divergent_transitions == 1.0
        @error "Numerical errors encountered in ALL iterations. Check model and priors."
    elseif ratio_divergent_transitions > 0.05
        @warn "Numerical errors encountered in more than 5% of iterations. Results are likely incorrect. Your model might have high curvature, and could be improved. Otherwise, increasing the target acceptance rate (second argument to `octofit`) might help" ratio_divergent_transitions
    end
    # if max_tree_depth_frac > 0.1
    #     @warn "Maximum tree depth hit in more than 10% of iterations (reduced efficiency)" max_tree_depth_frac
    # end

    # Resolve the array back into the nested named tuple structure used internally.
    # Augment with some internal fields
    chain_res = map(zip(stats, mc_samples)) do (stat, θ_t)
        # Map the variables back to the constrained domain and reconstruct the parameter
        # named tuple structure.
        θ = model.invlink(θ_t)
        resolved_namedtuple = model.arr2nt(θ)
        # Add log posterior, tree depth, and numerical error reported by
        # the sampler.
        # Also recompute the log-likelihood and add that too.
        loglike = ln_like_data(model.system, resolved_namedtuple)
        # Everything in the log posterior that is not data: the parameter
        # priors plus the prior-shaped terms `ln_like_data` skipped. Computed
        # as a difference so the two columns always reconstruct `logpost`
        # exactly, including the link's log-Jacobian.
        logprior = stat.log_density - loglike
        return merge((;
            stat.n_steps,
            stat.is_accept,
            stat.acceptance_rate,
            stat.hamiltonian_energy,
            stat.hamiltonian_energy_error,
            stat.max_hamiltonian_energy_error,
            stat.tree_depth,
            stat.numerical_error,
            stat.step_size,
            stat.nom_step_size,
            stat.is_adapt,
            loglike = loglike,
            logprior = logprior,
            logpost = stat.log_density,
        ), resolved_namedtuple)
    end
    # Then finally flatten and convert into an MCMCChain object / table.
    # Mark the posterior, likelihood, numerical error flag, and tree depth as internal
    mcmcchains = Octofitter.result2mcmcchain(
        chain_res,
        Dict(:internals => [
            :n_steps
            :is_accept
            :acceptance_rate
            :hamiltonian_energy
            :hamiltonian_energy_error
            :max_hamiltonian_energy_error
            :tree_depth
            :numerical_error
            :step_size
            :nom_step_size
            :is_adapt
            :loglike
            :logprior
            :logpost
            :tree_depth
            :numerical_error
        ])
    )

    # Concatenate the log posteriors and make them the same shape as the chains (N_iters,N_vars,N_chains)
    # logposts_mat = reduce(hcat, logposts)
    mcmcchains_with_info = MCMCChains.setinfo(
        mcmcchains,
        (;
            start_time,
            stop_time,
            samples_transformed=mc_samples,
            adaptor,
            initial_metric=metric,
            model_name=model.system.name,
            corrections=model.system.corrections,
            observing_geometry=model.system.observing_geometry,
            barycentric_lighttime=model.system.barycentric_lighttime,
            sampler="nuts"
        )
    )
    return mcmcchains_with_info
end

