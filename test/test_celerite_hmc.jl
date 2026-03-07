#!/usr/bin/env julia
"""
Test HMC sampling with Celerite GP — the ultimate test!
Based on the K2-131 tutorial from docs/src/rv-gp.md.
Previously required octofit_pigeons (no-gradient slice sampler) or AutoFiniteDiff.
Now with Enzyme, we can use octofit (HMC/NUTS) directly.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Octofitter
import Enzyme
Enzyme.API.looseTypeAnalysis!(true)
using OctofitterRadialVelocity
using OctofitterRadialVelocity.Celerite
using PlanetOrbits
using CSV
using Distributions

# ──────────────────────────────────────────────
# 1. Load the K2-131 data (same as tutorial)
# ──────────────────────────────────────────────

println("═══ Loading K2-131 RV data ═══")
rv_file = download("https://raw.githubusercontent.com/California-Planet-Search/radvel/master/example_data/k2-131.txt")
using DelimitedFiles
# Parse the space-delimited file manually
raw_lines = readlines(rv_file)
header = split(raw_lines[1])
data_lines = raw_lines[2:end]

times = Float64[]
mnvels = Float64[]
errvels = Float64[]
tels = String[]
for line in data_lines
    parts = split(line)
    length(parts) >= 4 || continue
    push!(times, parse(Float64, parts[1]))
    push!(mnvels, parse(Float64, parts[2]))
    push!(errvels, parse(Float64, parts[3]))
    push!(tels, parts[4])
end

# Convert JD to MJD
jd2mjd(jd) = jd - 2400000.5
epochs = jd2mjd.(times)

println("  Loaded $(length(epochs)) RV observations")
println("  Instruments: $(sort(unique(tels)))")

# Split by instrument
harps_mask = tels .== "harps-n"
pfs_mask = tels .== "pfs"

harps_data = Table(
    epoch=epochs[harps_mask],
    rv=mnvels[harps_mask],
    σ_rv=errvels[harps_mask]
)
pfs_data = Table(
    epoch=epochs[pfs_mask],
    rv=mnvels[pfs_mask],
    σ_rv=errvels[pfs_mask]
)

println("  HARPS-N: $(length(harps_data)) obs")
println("  PFS: $(length(pfs_data)) obs")

# ──────────────────────────────────────────────
# 2. Build model (same as tutorial, but no AutoFiniteDiff!)
# ──────────────────────────────────────────────

println("\n═══ Building model ═══")

rvlike_harps = StarAbsoluteRVObs(
    harps_data,
    name="harps_n",
    gaussian_process = θ_obs -> Celerite.CeleriteGP(
        Celerite.RealTerm(
            log(θ_obs.B*(1+θ_obs.C)/(2+θ_obs.C)),
            log(1/θ_obs.L)
        ) + Celerite.ComplexTerm(
            log(θ_obs.B/(2+θ_obs.C)),
            -Inf,
            log(1/θ_obs.L),
            log(2pi/θ_obs.Prot)
        )
    ),
    variables=@variables begin
        offset ~ Normal(-6693, 100)
        jitter ~ LogUniform(0.1, 100)
        B ~ Uniform(0.00001, 2000000)
        C ~ Uniform(0.00001, 200)
        L ~ Uniform(2, 200)
        Prot ~ Uniform(8.5, 20)
    end
)

rvlike_pfs = StarAbsoluteRVObs(
    pfs_data,
    name="pfs",
    gaussian_process = θ_obs -> Celerite.CeleriteGP(
        Celerite.RealTerm(
            log(θ_obs.B*(1+θ_obs.C)/(2+θ_obs.C)),
            log(1/θ_obs.L)
        ) + Celerite.ComplexTerm(
            log(θ_obs.B/(2+θ_obs.C)),
            -Inf,
            log(1/θ_obs.L),
            log(2pi/θ_obs.Prot)
        )
    ),
    variables=@variables begin
        offset ~ Normal(0, 100)
        jitter ~ LogUniform(0.1, 100)
        B ~ Uniform(0.00001, 2000000)
        C ~ Uniform(0.00001, 200)
        L ~ Uniform(2, 200)
        Prot ~ Uniform(8.5, 20)
    end
)

planet_1 = Planet(
    name="b",
    basis=RadialVelocityOrbit,
    observations=[],
    variables=@variables begin
        e = 0
        ω = 0.0
        P ~ truncated(
            Normal(0.3693038/365.256360417, 0.0000091/365.256360417),
            lower=0.0001
        )
        M = system.M
        a = cbrt(M * P^2)
        τ ~ UniformCircular(1.0)
        tp = τ*P*365.256360417 + 57782
        mass ~ LogUniform(0.001, 10)
    end
)

sys = System(
    name = "k2_131",
    companions=[planet_1],
    observations=[rvlike_harps, rvlike_pfs],
    variables=@variables begin
        M ~ truncated(Normal(0.82, 0.02), lower=0.1)
    end
)

# THE KEY: no autodiff= override → uses Enzyme by default!
model = Octofitter.LogDensityModel(sys)

# ──────────────────────────────────────────────
# 3. Verify gradient works
# ──────────────────────────────────────────────

println("\n═══ Gradient verification ═══")
θ = Octofitter.sample_priors(sys)
ll = model.ℓπcallback(θ)
println("  Primal ll = $ll")

ll_grad, grad = model.∇ℓπcallback(θ)
println("  Gradient ll = $ll_grad")
println("  Gradient norm = $(sqrt(sum(grad.^2)))")
println("  All finite? = $(all(isfinite, grad))")

# ──────────────────────────────────────────────
# 4. Run HMC!
# ──────────────────────────────────────────────

println("\n═══ Running octofit (HMC/NUTS) with Celerite GP! ═══")
println("  This was previously impossible — Celerite only worked with")
println("  finite differences or gradient-free samplers.")
println()

chain = octofit(
    model,
    adaptation=200,
    iterations=200,
)

println("\n═══ Results ═══")
display(chain)

println("\n\n═══ Key parameters ═══")
for param in [:b_mass, :b_P, :M]
    if hasproperty(chain, param)
        vals = vec(Array(chain[param]))
        println("  $param: median=$(round(median(vals), sigdigits=4)), std=$(round(std(vals), sigdigits=4))")
    end
end

println("\n✅ HMC sampling with Celerite GP succeeded!")
