using Optimization
using OptimizationOptimJL

"""
    mapchain = Octofitter.findmap(model::LogDensityModel)

Given an Octofitter model, find the maximum a-posteriori
point using optimization. Returns a Chains object with a 
single row.
Returning a Chains object is a bit weird, but this way 
it can be handled the same as our posteriors, plotted, etc.
"""
function findmap(model::LogDensityModel;starting_position=nothing,N=10_000,verbosity=0)
    if isnothing(starting_position)
        starting_position, _ = guess_starting_position(model,N)
    end

    # logpost = model.ℓπcallback(model.link(starting_position))

    θ_t′ = _findmap(model,model.link(starting_position);verbosity)
    θ′ = model.invlink(θ_t′)
    # Evaluate log post and log like
    logpost = model.ℓπcallback(θ_t′)
    resolved_namedtuple = model.arr2nt(θ′)
    # Add log posterior, tree depth, and numerical error reported by
    # the sampler.
    # Also recompute the log-likelihood and add that too.
    ln_like = make_ln_like(model.system, resolved_namedtuple)
    loglike = ln_like(model.system, resolved_namedtuple)
    nt = (; logpost, loglike, model.arr2nt(θ′)...)
    return result2mcmcchain(
        [nt], 
        Dict(:internals => [:logpost, :loglike])
    )
end

# Returns the raw parameter vector
function _findmap(model::LogDensityModel,initial_θ_t;verbosity=0)
    func = OptimizationFunction(
        (θ,model)->-model.ℓπcallback(θ),
        grad=(G,θ,model)->G.=.-model.∇ℓπcallback(θ)[2],
    )
    verbosity > 1 && @info "Guessing starting position" N
    # # Start with Simulated Annealing
    prob = OptimizationProblem(func, initial_θ_t, model)
    verbosity > 1 && @info "Simualted annealing optimization" N
    sol = solve(prob, SimulatedAnnealing(), iterations=1_000_000, x_tol=0)

    # Then iterate with qusi-Newton
    prob = OptimizationProblem(func, sol.u, model)
    verbosity > 1 && @info "LBFGS optimization" N
    sol = solve(prob, 
        Optim.LBFGS(;
            m=6,
            linesearch=Pathfinder.Optim.LineSearches.BackTracking(),
            alphaguess=Pathfinder.Optim.LineSearches.InitialHagerZhang()
        )
    , g_tol=1e-12, iterations=100000, allow_f_increases=true)
    θ_map2 = sol.u

    return θ_map2

    #     logpost = model.ℓπcallback(model.link(θ′))
    #     if sol.retcode == ReturnCode.Success && isfinite(logpost)
    #         return θ′
    #     end
    # end
    # error("Solution did not converge after 10 attempts")
end



# Returns the raw parameter vector
function _refine(model::LogDensityModel,initial_θ_t;verbosity=0)
    func = OptimizationFunction(
        (θ,model)->-model.ℓπcallback(θ),
        grad=(G,θ,model)->G.=.-model.∇ℓπcallback(θ)[2],
    )
    
    # Then iterate with qusi-Newton
    prob = OptimizationProblem(func, initial_θ_t, model)
    verbosity > 1 && @info "LBFGS optimization"
    sol = solve(prob, 
        Optim.LBFGS(;
            m=6,
            linesearch=Pathfinder.Optim.LineSearches.BackTracking(),
            alphaguess=Pathfinder.Optim.LineSearches.InitialHagerZhang()
        )
    , g_tol=1e-12, show_trace=true, show_every=100, iterations=100000, allow_f_increases=true)
    display(sol)
    θ_map2 = sol.u

    return θ_map2
end