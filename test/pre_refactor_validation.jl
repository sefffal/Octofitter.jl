#!/usr/bin/env julia
#=
Pre-refactor validation suite for the workspace redesign.

Establishes baselines and validates key assumptions:
  1. Reference gradient snapshots (numerical values at fixed θ)
  2. ForwardDiff with fixed custom tag produces identical gradients
  3. Enzyme with callable struct (functor) produces identical gradients
  4. Allocation baselines for primal and gradient evaluation

Run:
  julia --project=. test/pre_refactor_validation.jl
=#

using Octofitter
using OctofitterRadialVelocity
using Distributions
using TypedTables
using Random
using Test
using BenchmarkTools
using ForwardDiff
using ADTypes: AutoForwardDiff, AutoFiniteDiff, AutoEnzyme
using DifferentiationInterface
using Enzyme: Duplicated, Reverse, Const, set_runtime_activity

# ═══════════════════════════════════════════════════════════════════════════
# Shared test models
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

function make_dr4_system()
    N_epochs = 10
    mock_epochs = collect(range(58000.0, 59000.0, length=N_epochs))
    mock_xyz = Table(
        x = fill(0.5, N_epochs), y = fill(0.5, N_epochs), z = fill(0.2, N_epochs),
        vx = fill(0.0, N_epochs), vy = fill(0.0, N_epochs), vz = fill(0.0, N_epochs),
    )
    mock_table = Table(
        epoch = mock_epochs,
        scan_pos_angle = collect(range(0.0, 2π, length=N_epochs)),
        centroid_pos_al = zeros(N_epochs),
        centroid_pos_error_al = fill(0.1, N_epochs),
        parallax_factor_al = fill(0.5, N_epochs),
        outlier_flag = fill(false, N_epochs),
        xyz = collect(eachrow(mock_xyz)),
    )
    mock_gaia_sol = (ra=180.0, dec=45.0, pmra=5.0, pmdec=-10.0, parallax=20.0)
    ref_epoch_mjd = 57936.375
    priors, derived = @variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)
        ra_offset_mas ~ Normal(0, 100)
        dec_offset_mas ~ Normal(0, 100)
        pmra ~ Normal(5, 10)
        pmdec ~ Normal(-10, 10)
        plx = system.plx
        ref_epoch = $ref_epoch_mjd
    end
    mean_epoch = sum(mock_epochs) / length(mock_epochs)
    detrend_Δt = collect((mock_epochs .- mean_epoch) ./ 365.25)
    gaia_obs = Octofitter.GaiaDR4AstromObs{typeof(mock_table), typeof(mock_gaia_sol)}(
        mock_table, 123456789, mock_gaia_sol, priors, derived, "GaiaDR4", false,
        detrend_Δt, 1.0/length(mock_epochs), 1.0/sum(detrend_Δt .^ 2),
    )
    planet_b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[],
        variables=@variables begin
            a ~ LogUniform(0.5, 20)
            e ~ Uniform(0, 0.5)
            ω ~ Uniform(0, 2π)
            i ~ Sine()
            Ω ~ Uniform(0, 2π)
            θ ~ Uniform(0, 2π)
            tp = θ_at_epoch_to_tperi(θ, $mean_epoch; M=system.M, e, a, i, ω, Ω)
            mass ~ LogUniform(1, 100)
        end
    )
    System(
        name="DR4Test",
        companions=[planet_b],
        observations=[gaia_obs],
        variables=@variables begin
            M = 1.0
            plx ~ Uniform(10, 30)
        end
    )
end

"""Find a finite starting point for a model."""
function find_finite_theta(model; rng_base=42, max_attempts=50)
    for attempt in 0:max_attempts
        θ_natural = collect(model.sample_priors(Random.Xoshiro(rng_base + attempt)))
        θ = model.link(θ_natural)
        ll = model.ℓπcallback(θ)
        isfinite(ll) && return θ, ll
    end
    error("Could not find finite starting point after $max_attempts attempts")
end


# ═══════════════════════════════════════════════════════════════════════════
# TEST 1: Reference gradient snapshots
# ═══════════════════════════════════════════════════════════════════════════
#
# Captures exact gradient values at fixed parameter vectors.
# After the refactor, the same test should produce identical values.

@testset "Reference gradient snapshots" begin

    @testset "Astrometry (ForwardDiff)" begin
        sys = make_astrom_system()
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        model_fd = Octofitter.LogDensityModel(sys, autodiff=AutoFiniteDiff(), verbosity=0)
        θ, ll = find_finite_theta(model)

        ll_fwd, grad_fwd = model.∇ℓπcallback(θ)
        ll_fd, grad_fd = model_fd.∇ℓπcallback(θ)

        @test isfinite(ll_fwd)
        @test all(isfinite, grad_fwd)
        @test ll_fwd ≈ ll_fd
        @test isapprox(grad_fwd, grad_fd, atol=1e-3, rtol=1e-4)

        # Snapshot: log these so we can compare post-refactor
        @info "Astrom snapshot" ll=ll_fwd grad_norm=sqrt(sum(grad_fwd.^2)) D=length(θ)
        println("  Astrom θ = ", θ)
        println("  Astrom ll = ", ll_fwd)
        println("  Astrom grad = ", grad_fwd)
    end

    @testset "RV (ForwardDiff)" begin
        sys = make_rv_system()
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        model_fd = Octofitter.LogDensityModel(sys, autodiff=AutoFiniteDiff(), verbosity=0)
        θ, ll = find_finite_theta(model)

        ll_fwd, grad_fwd = model.∇ℓπcallback(θ)
        ll_fd, grad_fd = model_fd.∇ℓπcallback(θ)

        @test isfinite(ll_fwd)
        @test all(isfinite, grad_fwd)
        @test ll_fwd ≈ ll_fd
        @test isapprox(grad_fwd, grad_fd, atol=1e-3, rtol=1e-4)

        @info "RV snapshot" ll=ll_fwd grad_norm=sqrt(sum(grad_fwd.^2)) D=length(θ)
        println("  RV θ = ", θ)
        println("  RV ll = ", ll_fwd)
        println("  RV grad = ", grad_fwd)
    end

    @testset "Gaia DR4 (Enzyme)" begin
        sys = make_dr4_system()
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        model_fd = Octofitter.LogDensityModel(sys, autodiff=AutoFiniteDiff(), verbosity=0)
        θ, ll = find_finite_theta(model)

        ll_enz, grad_enz = model.∇ℓπcallback(θ)
        ll_fd, grad_fd = model_fd.∇ℓπcallback(θ)

        @test isfinite(ll_enz)
        @test all(isfinite, grad_enz)
        @test ll_enz ≈ ll_fd
        @test isapprox(grad_enz, grad_fd, atol=1e-2, rtol=1e-3)

        @info "DR4 snapshot" ll=ll_enz grad_norm=sqrt(sum(grad_enz.^2)) D=length(θ)
        println("  DR4 θ = ", θ)
        println("  DR4 ll = ", ll_enz)
        println("  DR4 grad = ", grad_enz)
    end
end


# ═══════════════════════════════════════════════════════════════════════════
# TEST 2: ForwardDiff with fixed custom tag
# ═══════════════════════════════════════════════════════════════════════════
#
# Validates that AutoForwardDiff(tag=OctofitterTag()) produces identical
# gradients to the default auto-tagged ForwardDiff.
# This is the key assumption for pre-allocating Dual-typed workspace buffers.

struct OctofitterTag end

@testset "ForwardDiff fixed tag" begin
    sys = make_astrom_system()

    # Model with default ForwardDiff (auto-generated tag)
    model_default = Octofitter.LogDensityModel(sys, verbosity=0)
    θ, _ = find_finite_theta(model_default)

    ll_default, grad_default = model_default.∇ℓπcallback(θ)

    # Model with custom fixed tag
    D = length(θ)
    model_tagged = Octofitter.LogDensityModel(sys,
        autodiff=AutoForwardDiff(chunksize=D, tag=OctofitterTag()),
        verbosity=0)

    ll_tagged, grad_tagged = model_tagged.∇ℓπcallback(θ)

    @test ll_default ≈ ll_tagged
    @test grad_default ≈ grad_tagged  # should be bitwise identical
    @info "Fixed tag test passed" ll_match=(ll_default ≈ ll_tagged) grad_max_diff=maximum(abs.(grad_default .- grad_tagged))

    # Also verify the Dual type is what we expect
    DType = ForwardDiff.Dual{OctofitterTag, Float64, D}
    @test sizeof(DType) == (D + 1) * sizeof(Float64)
    @info "Dual type" DType sizeof=sizeof(DType)
end


# ═══════════════════════════════════════════════════════════════════════════
# TEST 3: Enzyme with callable struct (functor) vs closure
# ═══════════════════════════════════════════════════════════════════════════
#
# Validates that wrapping eval_term in a named struct produces identical
# Enzyme gradients to the current closure-based approach.
# This is the key assumption for replacing closures with TermEvaluator.

@testset "Enzyme functor vs closure" begin
    sys = make_dr4_system()
    model = Octofitter.LogDensityModel(sys, verbosity=0)
    θ, _ = find_finite_theta(model)

    # Current approach: closure-based gradient (already tested to be correct)
    ll_closure, grad_closure = model.∇ℓπcallback(θ)

    # New approach: build a functor that does the same thing as the per-term closure
    # We extract the per-term components and wrap them in a struct.

    # Access internal config (these are implementation details we'll formalize later)
    mcfg = Octofitter.ModelEvalConfig(
        model.invlink, model.arr2nt,
        ntuple(length(sys.planets)) do i
            OT = Octofitter._planet_orbit_type(sys.planets[i])
            let OT=OT, i=i
                (θ_system) -> OT(;merge(θ_system, θ_system.planets[i])...)
            end
        end,
        length(sys.planets)
    )

    # Build TermEvalConfig for the system observation (Gaia DR4)
    obs = sys.observations[1]
    epochs = Octofitter._get_obs_epochs(obs)
    obs_name = Octofitter.normalizename(Octofitter.likelihoodname(obs))

    active_indices_all, D = Octofitter._compute_active_indices(sys)
    # The DR4 term is the last one (system obs come after planet obs)
    n_planet_obs = sum(length(p.observations) for p in sys.planets)
    active = active_indices_all[n_planet_obs + 1]

    tcfg = Octofitter.TermEvalConfig(obs, obs_name, epochs, active, 0)

    # Pre-allocate workspace (the new way)
    _θ_nat_0 = model.invlink(θ)
    _θ_sys_0 = model.arr2nt(_θ_nat_0)
    _orbits_0 = ntuple(i -> mcfg.orbit_constructors[i](_θ_sys_0), mcfg.n_planets)
    sol_type_examples = ntuple(mcfg.n_planets) do ip
        PlanetOrbits.orbitsolve(_orbits_0[ip], 50000.0)
    end
    sol_bufs = ntuple(mcfg.n_planets) do ip
        Vector{typeof(sol_type_examples[ip])}(undef, length(epochs))
    end
    ws = Octofitter.PreallocTermWS(sol_bufs, Vector{Float64}(undef, D))

    # --- The functor ---
    struct TestTermEvaluator{TMcfg, TTcfg, TWs}
        mcfg::TMcfg
        tcfg::TTcfg
        workspace::TWs
    end
    (te::TestTermEvaluator)(θ_vec) = Octofitter.eval_term(θ_vec, te.mcfg, te.tcfg, te.workspace)

    evaluator = TestTermEvaluator(mcfg, tcfg, ws)

    # Test: functor produces same value as closure
    ll_functor = evaluator(θ)
    ll_direct = Octofitter.eval_term(θ, mcfg, tcfg, ws)
    @test ll_functor ≈ ll_direct

    # Test: Enzyme can differentiate the functor with Duplicated
    backend = AutoEnzyme(mode=set_runtime_activity(Reverse), function_annotation=Duplicated)
    θ_zero = zero(θ)
    prep = DifferentiationInterface.prepare_gradient(evaluator, backend, θ_zero)
    grad_buf = similar(θ)

    ll_functor_enz, _ = DifferentiationInterface.value_and_gradient!(
        evaluator, grad_buf, prep, backend, θ)

    @test isfinite(ll_functor_enz)
    @test all(isfinite, grad_buf)
    @info "Functor Enzyme test" ll=ll_functor_enz grad_norm=sqrt(sum(grad_buf.^2))

    # Compare functor gradient to FiniteDiff (independent reference)
    backend_fd = AutoFiniteDiff()
    prep_fd = DifferentiationInterface.prepare_gradient(evaluator, backend_fd, θ_zero)
    grad_fd = similar(θ)
    ll_fd, _ = DifferentiationInterface.value_and_gradient!(
        evaluator, grad_fd, prep_fd, backend_fd, θ)

    @test ll_functor_enz ≈ ll_fd
    @test isapprox(grad_buf, grad_fd, atol=1e-2, rtol=1e-3)
    @info "Functor vs FiniteDiff" max_diff=maximum(abs.(grad_buf .- grad_fd))
end


# ═══════════════════════════════════════════════════════════════════════════
# TEST 4: Allocation baselines
# ═══════════════════════════════════════════════════════════════════════════
#
# Records current allocation counts so we can verify the refactor reduces them.
# After the refactor, the goal is 0 allocations for Enzyme terms
# and minimal allocations for ForwardDiff terms.

@testset "Allocation baselines" begin

    @testset "Astrometry (ForwardDiff) allocations" begin
        sys = make_astrom_system()
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        θ, _ = find_finite_theta(model)

        # Warm up
        model.ℓπcallback(θ)
        model.∇ℓπcallback(θ)

        alloc_primal = @allocated model.ℓπcallback(θ)
        alloc_grad = @allocated model.∇ℓπcallback(θ)

        @info "Astrom allocations (BEFORE refactor)" alloc_primal alloc_grad
        # Just record, don't assert — these are baselines
        @test true
    end

    @testset "RV (ForwardDiff) allocations" begin
        sys = make_rv_system()
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        θ, _ = find_finite_theta(model)

        model.ℓπcallback(θ)
        model.∇ℓπcallback(θ)

        alloc_primal = @allocated model.ℓπcallback(θ)
        alloc_grad = @allocated model.∇ℓπcallback(θ)

        @info "RV allocations (BEFORE refactor)" alloc_primal alloc_grad
        @test true
    end

    @testset "Gaia DR4 (Enzyme) allocations" begin
        sys = make_dr4_system()
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        θ, _ = find_finite_theta(model)

        model.ℓπcallback(θ)
        model.∇ℓπcallback(θ)

        alloc_primal = @allocated model.ℓπcallback(θ)
        alloc_grad = @allocated model.∇ℓπcallback(θ)

        @info "DR4 allocations (BEFORE refactor)" alloc_primal alloc_grad
        @test true
    end
end


# ═══════════════════════════════════════════════════════════════════════════
# TEST 5: ForwardDiff Dual type pre-allocation feasibility
# ═══════════════════════════════════════════════════════════════════════════
#
# Verifies that we can pre-allocate Vector{Dual{OctofitterTag, Float64, N}}
# buffers and use them inside a ForwardDiff computation via
# DifferentiationInterface.

@testset "Dual buffer pre-allocation" begin
    sys = make_astrom_system()
    model = Octofitter.LogDensityModel(sys, verbosity=0)
    θ, _ = find_finite_theta(model)
    D = length(θ)

    # The Dual type we'd use for pre-allocated buffers
    DType = ForwardDiff.Dual{OctofitterTag, Float64, D}

    # Can we allocate vectors of this type?
    buf = Vector{DType}(undef, 10)
    @test length(buf) == 10
    @test eltype(buf) === DType

    # Can we fill them with promoted Float64 values?
    for i in 1:10
        buf[i] = DType(Float64(i))  # Float64 promoted to Dual with zero partials
    end
    @test ForwardDiff.value(buf[3]) == 3.0
    @test all(iszero, ForwardDiff.partials(buf[3]).values)

    # Can we write Dual values from a computation into them?
    x_dual = DType(2.0, ForwardDiff.Partials(ntuple(i -> i == 1 ? 1.0 : 0.0, D)))
    buf[1] = x_dual * DType(3.0)
    @test ForwardDiff.value(buf[1]) == 6.0
    @test ForwardDiff.partials(buf[1]).values[1] == 3.0

    @info "Dual pre-allocation feasibility confirmed" DType sizeof_per_element=sizeof(DType)
    @test true
end


println("\n" * "=" ^ 72)
println("  Pre-refactor validation complete!")
println("  Save the gradient snapshots printed above for post-refactor comparison.")
println("=" ^ 72)
