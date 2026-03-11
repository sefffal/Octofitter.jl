#!/usr/bin/env julia
#=
Benchmark comparison script for per-obs-diff vs main.
Tests relative astrometry and RV model performance (logdensity + gradient).

Usage:
  julia --project=test test/benchmark_branch_comparison.jl
=#

using Octofitter
using OctofitterRadialVelocity
using Distributions
using TypedTables
using Random
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10

# ═══════════════════════════════════════════════════════════════════════════
# Model 1: Relative Astrometry
# ═══════════════════════════════════════════════════════════════════════════
function make_astrom_system()
    astrom_table = Table(
        epoch = collect(range(50_000.0, 50_500.0, length=10)),
        ra    = [10.0, -15.0, 20.0, -5.0, 8.0, -12.0, 18.0, -3.0, 14.0, -9.0],
        dec   = [5.0, 12.0, -8.0, 15.0, -3.0, 10.0, -6.0, 11.0, -2.0, 7.0],
        σ_ra  = fill(5.0, 10),
        σ_dec = fill(5.0, 10),
    )
    planet_b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[PlanetRelAstromLikelihood(astrom_table; name="astrom1")],
        variables=@variables begin
            a  ~ LogUniform(0.5, 50)
            e  ~ Uniform(0, 0.5)
            ω  ~ Uniform(0, 2π)
            i  ~ Sine()
            Ω  ~ Uniform(0, 2π)
            θ  ~ Uniform(0, 2π)
            tp = θ_at_epoch_to_tperi(θ, 50_000.0; M=system.M, e, a, i, ω, Ω)
            mass ~ LogUniform(1, 1000)
        end
    )
    System(
        name="AstromTest",
        companions=[planet_b],
        observations=[],
        variables=@variables begin
            M = 1.0
            plx ~ Uniform(10, 100)
        end
    )
end

# ═══════════════════════════════════════════════════════════════════════════
# Model 2: Radial Velocity
# ═══════════════════════════════════════════════════════════════════════════
function make_rv_system()
    rv_table = Table(
        epoch = [50000.0, 50100.0, 50200.0, 50300.0, 50400.0],
        rv = [100.0, 110.0, 105.0, 95.0, 108.0],
        σ_rv = [10.0, 10.0, 10.0, 10.0, 10.0],
    )
    rvlike = StarAbsoluteRVLikelihood(
        rv_table;
        name="rv1",
        variables=@variables begin
            jitter ~ LogUniform(0.1, 100)
            offset ~ Normal(0, 100)
        end
    )
    planet_b = Planet(
        name="b",
        basis=RadialVelocityOrbit,
        observations=[],
        variables=@variables begin
            e ~ Uniform(0.0, 0.5)
            ω ~ UniformCircular()
            τ ~ Uniform(0, 1)
            P ~ Uniform(0.001, 10)
            a = ∛(P^2 * system.M)
            tp = τ * P * 365.25 + 50000
            mass ~ Uniform(0, 100)
        end
    )
    System(
        name="RVTest",
        companions=[planet_b],
        observations=[rvlike],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        end
    )
end

# ═══════════════════════════════════════════════════════════════════════════
# Run benchmarks
# ═══════════════════════════════════════════════════════════════════════════

function run_benchmarks()
    results = Dict{String, Any}()

    for (label, make_sys) in [("RelAstrom", make_astrom_system), ("RV", make_rv_system)]
        println("\n", "="^60)
        println("Model: $label")
        println("="^60)

        sys = make_sys()
        model = Base.invokelatest(Octofitter.LogDensityModel, sys; verbosity=2)

        # Get a valid parameter vector
        Random.seed!(42)
        θ_t = Base.invokelatest(model.link, Base.invokelatest(model.sample_priors, Random.default_rng()))

        # Report AD backend
        ad_type = typeof(model).parameters[end]
        println("  AD backend (type param): $ad_type")

        # Warmup + correctness via invokelatest
        logpost, grad, b_primal, b_grad = Base.invokelatest() do
            logpost = model.ℓπcallback(θ_t)
            println("  logposterior = $logpost")

            grad_result = model.∇ℓπcallback(θ_t)
            if grad_result isa Tuple
                logpost_g, grad = grad_result
            else
                logpost_g = grad_result
                grad = nothing
            end
            println("  logposterior (from gradient) = $logpost_g")
            if grad !== nothing
                println("  gradient = ", repr(grad))
                println("  gradient norm = $(sqrt(sum(abs2, grad)))")
            end

            # Benchmark primal
            println("\n  Benchmarking primal evaluation...")
            b_primal = @benchmark $(model.ℓπcallback)($θ_t)
            display(b_primal)

            # Benchmark gradient
            println("\n  Benchmarking gradient evaluation...")
            b_grad = @benchmark $(model.∇ℓπcallback)($θ_t)
            display(b_grad)

            logpost, grad, b_primal, b_grad
        end

        results["$(label)_primal"] = b_primal
        results["$(label)_gradient"] = b_grad
    end

    # Print summary
    println("\n\n", "="^60)
    println("SUMMARY")
    println("="^60)
    for key in sort(collect(keys(results)))
        b = results[key]
        println("  $key:")
        println("    median: $(BenchmarkTools.prettytime(median(b).time))")
        println("    allocs: $(median(b).allocs)")
        println("    memory: $(BenchmarkTools.prettymemory(median(b).memory))")
    end

    return results
end

run_benchmarks()
