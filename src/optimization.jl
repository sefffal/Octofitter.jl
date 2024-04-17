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
function findmap(model::LogDensityModel,N=100_000;verbosity=0)
    θ′ = _findmap(model,N;verbosity)
    logpost = model.ℓπcallback(model.link(θ′))
    nt = (; logpost, model.arr2nt(θ′)...)
    return result2mcmcchain(
        [nt], 
        Dict(:internals => [:logpost])
    )
end

# Returns the raw parameter vector
function _findmap(model::LogDensityModel,N=100_000;verbosity=0)
    func = OptimizationFunction(
        (θ,model)->-model.ℓπcallback(θ),
        grad=(G,θ,model)->G.=.-model.∇ℓπcallback(θ)[2],
    )
    verbosity > 1 && @info "Guessing starting position" N
    θ0, _ = guess_starting_position(model.system,N)

    # Start with Simulated Annealing
    prob = OptimizationProblem(func, θ0, model)
    verbosity > 1 && @info "Simualted annealing optimization" N
    sol = solve(prob, SimulatedAnnealing(), iterations=1_00_000, x_tol=0)
    θ0 = sol.u

    # Then iterate with qusi-Newton
    prob = OptimizationProblem(func, sol.u, model)
    verbosity > 1 && @info "LBFGS optimization" N
    sol = solve(prob, LBFGS(), g_tol=1e-12, iterations=10000, allow_f_increases=true)
    θ0 = sol.u

    θ′ = model.invlink(θ0)
    return θ′

    #     logpost = model.ℓπcallback(model.link(θ′))
    #     if sol.retcode == ReturnCode.Success && isfinite(logpost)
    #         return θ′
    #     end
    # end
    # error("Solution did not converge after 10 attempts")
end