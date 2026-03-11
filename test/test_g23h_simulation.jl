#=
Test that generate_from_params works for G23HObs by running a single
injection-recovery trial: inject a 150 Mjup companion at 3 AU and check
that the sampler recovers the injected parameters.

Uses Pigeons (octofit_pigeons) as in the production G23H workflow.
=#

using Octofitter, Distributions, Statistics, Random, Pigeons

println("Setting up G23H observation...")
absastrom = G23HObs(; hip_id=18512, freeze_epochs=true, include_rv=false)

ref_epoch = Octofitter.meta_gaia_DR3.ref_epoch_mjd

planet = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[],
    variables=@variables begin
        a    ~ LogUniform(0.1, 100)
        e    ~ Uniform(0, 0.99)
        ω    ~ Uniform(0, 2pi)
        i    ~ Sine()
        Ω    ~ Uniform(0, 2pi)
        θ    ~ Uniform(0, 2pi)
        tp   = θ_at_epoch_to_tperi(θ, $ref_epoch; M=system.M, e, a, i, ω, Ω)
        mass ~ LogUniform(0.01, 1000)
    end
)

sys = System(
    name="test_g23h_sim",
    companions=[planet],
    observations=[absastrom],
    variables=@variables begin
        M   ~ truncated(Normal(0.71, 0.04), lower=0.1)
        plx ~ truncated(Normal(absastrom.catalog.parallax, absastrom.catalog.parallax_error),
                         lower=max(0.1, absastrom.catalog.parallax - 10*absastrom.catalog.parallax_error))
        pmra  ~ Uniform(absastrom.catalog.pmra_dr3 - 100, absastrom.catalog.pmra_dr3 + 100)
        pmdec ~ Uniform(absastrom.catalog.pmdec_dr3 - 100, absastrom.catalog.pmdec_dr3 + 100)
        dec = $(absastrom.catalog.dec)
        ra  = $(absastrom.catalog.ra)
        ref_epoch = $ref_epoch
    end
)

# ── Injection parameters ──
inject_mass = 150.0  # Mjup — very large, should be easily detected
inject_sep  = 3.0    # AU

inject = (mass, sep) -> (; planets=(; b=(; mass=mass, a=sep)))

println("Running completeness trial with mass=$inject_mass Mjup, sep=$inject_sep AU...")

job = CompletenessJob(1, 1, 1, inject_mass, inject_sep, UInt64(42))

# Use Pigeons parallel-tempered sampling (as in G23H production workflow)
sampler = function(model)
    chain, pt = octofit_pigeons(
        model;
        n_rounds=6,
        n_chains=8,
        n_chains_variational=0,
        variational=nothing,
    )
    return chain
end

result = run_completeness_trial(
    job, sys, sampler;
    inject=inject,
    add_noise=true,
    verbosity=2,
)

# ── Check recovery ──
chain = result.chain
θ_true = result.θ_true

mass_samples = vec(chain["b_mass"])
a_samples = vec(chain["b_a"])

mass_med = median(mass_samples)
mass_lo  = quantile(mass_samples, 0.05)
mass_hi  = quantile(mass_samples, 0.95)

a_med = median(a_samples)
a_lo  = quantile(a_samples, 0.05)
a_hi  = quantile(a_samples, 0.95)

true_mass = θ_true.planets.b.mass
true_a    = θ_true.planets.b.a

println("\n=== Recovery Results ===")
println("Injected mass:  $true_mass Mjup")
println("Recovered mass: $mass_med Mjup  (90% CI: $mass_lo — $mass_hi)")
println("Mass within 10x factor: $(mass_med > true_mass/10 && mass_med < true_mass*10)")
println()
println("Injected a:     $true_a AU")
println("Recovered a:    $a_med AU  (90% CI: $a_lo — $a_hi)")
println("a within 3x factor:     $(a_med > true_a/3 && a_med < true_a*3)")
println()

# Basic sanity checks
@assert isfinite(mass_med) "Mass posterior is not finite"
@assert mass_med > 0 "Mass posterior median is non-positive"
@assert mass_med > true_mass / 10 "Mass underestimated by >10x: $mass_med vs $true_mass"
@assert mass_med < true_mass * 10 "Mass overestimated by >10x: $mass_med vs $true_mass"
@assert a_med > true_a / 3 "Semi-major axis underestimated by >3x: $a_med vs $true_a"
@assert a_med < true_a * 3 "Semi-major axis overestimated by >3x: $a_med vs $true_a"

println("✓ All recovery checks passed!")
