# OFTI-style Semi-Linear Orbit Fitting with Octofitter
#
# This example demonstrates how to use the Thiele-Innes semi-linear approach
# to efficiently fit orbits. The key idea: for fixed nonlinear parameters
# (e, a, tp, M, plx), the sky-plane positions are LINEAR in the Thiele-Innes
# constants (A, B, F, G). We analytically marginalize over A, B, F, G,
# reducing the problem from ~11 dimensions to just 5.
#
# We show two approaches:
#   1. Rejection sampling (octofit_rejection) -- independent samples, no tuning
#   2. HMC (octofit) -- more efficient for this dimensionality
#
# Both use the `LL +=` syntax to inject the marginalized log-likelihood
# computed by `ofti_linear_solve`.

using Octofitter
using Distributions
using Random
using PlanetOrbits
using Statistics

# ── Generate synthetic data from a known orbit ──────────────────────────────

rng = Random.Xoshiro(42)

true_orbit = Visual{KepOrbit}(;
    a   = 10.0,     # AU
    e   = 0.3,
    i   = 1.0,      # rad
    ω   = 0.5,      # rad
    Ω   = 2.0,      # rad
    tp  = 50000.0,   # MJD
    M   = 1.2,       # solar masses
    plx = 50.0,      # mas
)

epochs = collect(range(50000, 50840, length=8))
σ_astrom = 10.0  # mas per axis

ra_true  = [raoff(orbitsolve(true_orbit, ep))  for ep in epochs]
dec_true = [decoff(orbitsolve(true_orbit, ep)) for ep in epochs]

# Add Gaussian noise
ra_obs  = ra_true  .+ σ_astrom .* randn(rng, length(epochs))
dec_obs = dec_true .+ σ_astrom .* randn(rng, length(epochs))

σ_ra  = fill(σ_astrom, length(epochs))
σ_dec = fill(σ_astrom, length(epochs))
cor   = fill(0.0, length(epochs))

println("True parameters: a=$(semimajoraxis(true_orbit)), e=$(eccentricity(true_orbit)), M=$(totalmass(true_orbit)), plx=$(true_orbit.plx)")
println("Generated $(length(epochs)) astrometry epochs with σ=$(σ_astrom) mas")

# ── Define the OFTI model ───────────────────────────────────────────────────
#
# Free parameters: e, a, τ, M, plx  (5 dimensions)
# Marginalized:    A, B, F, G        (analytically via ofti_linear_solve)
#
# The `LL +=` syntax injects the marginal log-likelihood directly,
# and derived variables recover the best-fit A, B, F, G for post-processing.

sys = System(
    name = "OFTI_Demo",
    companions = [],
    observations = [],
    variables = @variables begin
        M   ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.5), lower=0.1)
        e   ~ Uniform(0, 0.99)
        a   ~ LogUniform(1, 100)
        τ   ~ Uniform(0, 1)

        # Derive tp from the phase fraction τ
        _period_days = √(a^3 / M) * $(PlanetOrbits.kepler_year_to_julian_day_conversion_factor)
        tp = $(epochs[1]) + τ * _period_days

        # Solve for A, B, F, G analytically and get the marginal likelihood
        _result = Octofitter.ofti_linear_solve(
            $epochs, $ra_obs, $dec_obs, $σ_ra, $σ_dec, $cor,
            1000.0,   # σ_ABFG: prior width on Thiele-Innes constants (mas)
            e, a, tp, M, plx,
        )

        # Inject the marginalized log-likelihood
        LL += _result.log_marginal_likelihood

        # Store best-fit Thiele-Innes constants as derived variables
        A = _result.A
        B = _result.B
        F = _result.F
        G = _result.G
    end
)

model = Octofitter.LogDensityModel(sys, verbosity=1)
println("Model has $(model.D) free parameters")

# ── Approach 1: Rejection sampling ──────────────────────────────────────────
#
# Draw from the prior, accept/reject based on the marginalized likelihood.
# No tuning needed, samples are independent (zero autocorrelation).
# Works well here because we only have 5 free parameters.

println("\n" * "="^60)
println("Approach 1: Rejection Sampling (octofit_rejection)")
println("="^60)

chain_rej = octofit_rejection(rng, model, draws=1_000_000)

println("\nRejection sampling results ($(size(chain_rej, 1)) accepted):")
println("  e:   median=$(round(median(chain_rej[:e]),   digits=2)),  true=$(eccentricity(true_orbit))")
println("  a:   median=$(round(median(chain_rej[:a]),   digits=1)),  true=$(semimajoraxis(true_orbit))")
println("  M:   median=$(round(median(chain_rej[:M]),   digits=2)),  true=$(totalmass(true_orbit))")
println("  plx: median=$(round(median(chain_rej[:plx]), digits=1)),  true=$(true_orbit.plx)")

# ── Approach 2: HMC ─────────────────────────────────────────────────────────
#
# Hamiltonian Monte Carlo via AdvancedHMC (NUTS). More sample-efficient than
# rejection sampling, but requires gradient computation and adaptation.
# The `LL +=` term is automatically differentiated by ForwardDiff.

println("\n" * "="^60)
println("Approach 2: HMC (octofit)")
println("="^60)

chain_hmc = octofit(rng, model, iterations=5000, adaptation=1000)

println("\nHMC results ($(size(chain_hmc, 1)) samples):")
println("  e:   median=$(round(median(chain_hmc[:e]),   digits=2)),  true=$(eccentricity(true_orbit))")
println("  a:   median=$(round(median(chain_hmc[:a]),   digits=1)),  true=$(semimajoraxis(true_orbit))")
println("  M:   median=$(round(median(chain_hmc[:M]),   digits=2)),  true=$(totalmass(true_orbit))")
println("  plx: median=$(round(median(chain_hmc[:plx]), digits=1)),  true=$(true_orbit.plx)")

# ── Reconstruct full orbits from the chain ──────────────────────────────────
#
# The chain contains the derived A, B, F, G values (posterior means from the
# linear solve at each accepted sample). These can be used to construct
# ThieleInnesOrbit objects for plotting and further analysis.

println("\n" * "="^60)
println("Reconstructing orbits from HMC chain")
println("="^60)

# Build a ThieleInnesOrbit from the first sample
i = 1
orb = ThieleInnesOrbit(;
    e  = chain_hmc[:e][i],
    tp = chain_hmc[:tp][i],
    M  = chain_hmc[:M][i],
    plx = chain_hmc[:plx][i],
    A  = chain_hmc[:A][i],
    B  = chain_hmc[:B][i],
    F  = chain_hmc[:F][i],
    G  = chain_hmc[:G][i],
)
println("Sample 1 orbit: a=$(round(semimajoraxis(orb), digits=1)) AU, e=$(round(eccentricity(orb), digits=2))")
