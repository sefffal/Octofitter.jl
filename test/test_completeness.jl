using Octofitter, Distributions

# Simple test system with relative astrometry
astrom = PlanetRelAstromObs(
    Table(
        epoch=[50000.0, 50120.0, 50240.0, 50360.0],
        ra=[500.0, 400.0, 200.0, -100.0],
        dec=[300.0, 600.0, 700.0, 500.0],
        σ_ra=[10.0, 10.0, 10.0, 10.0],
        σ_dec=[10.0, 10.0, 10.0, 10.0],
    ),
    name="astrom",
)

planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[astrom],
    variables=Octofitter.@variables(begin
        e ~ Uniform(0, 0.5)
        mass ~ LogUniform(0.1, 100)
        a ~ LogUniform(1, 50)
        i ~ Sine()
        Ω ~ UniformCircular()
        ω ~ UniformCircular()
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ, 50000.0; M=system.M, e, a, i, ω, Ω)
    end),
)

sys = System(
    name="test",
    companions=[planet_b],
    observations=Octofitter.AbstractObs[],
    variables=Octofitter.@variables(begin
        M ~ truncated(Normal(1.0, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.5), lower=1.0)
    end),
)

# ── Test 1: _apply_overrides! ──
println("Test 1: _apply_overrides!")
arr2nt = Octofitter.make_arr2nt(sys)
θ_flat = Octofitter.sample_priors(sys)
θ_nt = arr2nt(θ_flat)
println("  Before override, mass = ", θ_nt.planets.b.mass)

overrides = (; planets=(; b=(; mass=42.0, a=10.0)))
Octofitter._apply_overrides!(θ_flat, arr2nt, overrides)
θ_nt2 = arr2nt(θ_flat)
println("  After override, mass = ", θ_nt2.planets.b.mass, ", a = ", θ_nt2.planets.b.a)
@assert θ_nt2.planets.b.mass ≈ 42.0
@assert θ_nt2.planets.b.a ≈ 10.0
println("  PASSED")

# ── Test 2: generate_from_params with overrides ──
println("\nTest 2: generate_from_params with overrides")
sim_sys = generate_from_params(sys, θ_nt2; add_noise=true)
println("  Simulated system has $(length(sim_sys.planets[1].observations)) observations")
println("  Simulated data ra[1] = ", sim_sys.planets[1].observations[1].table.ra[1])
println("  PASSED")

# ── Test 3: _initialize_at_truth! ──
println("\nTest 3: _initialize_at_truth!")
model = Octofitter.LogDensityModel(sim_sys; verbosity=0)
Octofitter._initialize_at_truth!(model, θ_flat)
println("  Starting points set: ", length(model.starting_points), " points")
lp = model.ℓπcallback(model.starting_points[1])
println("  Log posterior at truth: ", lp)
@assert isfinite(lp)
println("  PASSED")

# ── Test 4: Full single trial (no detection decision) ──
println("\nTest 4: run_completeness_trial")
job = CompletenessJob(1, 1, 1, 5.0, 8.0, UInt64(12345))

inject = (mass, sep) -> (; planets=(; b=(; mass=mass, a=sep)))

result = run_completeness_trial(
    job, sys,
    model -> octofit(model, iterations=200, adaptation=200, verbosity=0);
    inject=inject,
    add_noise=true,
    verbosity=1,
)
println("  chain type: ", typeof(result.chain))
println("  θ_true mass: ", result.θ_true.planets.b.mass)
@assert result.θ_true.planets.b.mass ≈ 5.0
@assert result.θ_true.planets.b.a ≈ 8.0
@assert length(vec(result.chain["b_mass"])) > 0
println("  PASSED")

# ── Test 5: assemble_completeness with detection criterion ──
println("\nTest 5: assemble_completeness with detection criterion")

# Detection applied here, not during the trial
detection = function(chain, θ_true)
    mass_samples = vec(chain["b_mass"])
    return 0.33 * θ_true.planets.b.mass < median(mass_samples) < 3.0 * θ_true.planets.b.mass
end

cmap = assemble_completeness(
    [result], detection;
    masses=[5.0], separations=[8.0],
)
println("  Completeness: ", cmap.completeness)
println("  PASSED")

# ── Test 6: re-apply different detection criterion ──
println("\nTest 6: re-apply different threshold")

# Stricter criterion on same result
cmap_strict = assemble_completeness(
    [result],
    (chain, θ) -> quantile(vec(chain["b_mass"]), 0.05) > 1.0;
    masses=[5.0], separations=[8.0],
)
println("  Strict completeness: ", cmap_strict.completeness)

# Very lenient criterion
cmap_lenient = assemble_completeness(
    [result],
    (chain, θ) -> median(vec(chain["b_mass"])) > 0.01;
    masses=[5.0], separations=[8.0],
)
println("  Lenient completeness: ", cmap_lenient.completeness)
println("  PASSED")

println("\n✓ All completeness tests passed!")
