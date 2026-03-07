using Octofitter, Distributions, CSV, PlanetOrbits, BenchmarkTools

df = CSV.read(joinpath(@__DIR__, "..", "docs", "src", "target_1.csv"), FlexTable)
ref_epoch_mjd = 57936.375
gaia_dr4_obs = GaiaDR4AstromObs(df, gaia_id=4373465352415301632,
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)
        ra_offset_mas  ~ Normal(0, 10000)
        dec_offset_mas ~ Normal(0, 10000)
        pmra ~ Uniform(-1000, 1000)
        pmdec ~ Uniform(-1000, 1000)
        plx = system.plx
        ref_epoch = $ref_epoch_mjd
    end)
mjup2msol = Octofitter.mjup2msol
ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd
orbit_ref_epoch = mean(gaia_dr4_obs.table.epoch)
b = Planet(name="b", basis=Visual{KepOrbit}, observations=[],
    solver=PlanetOrbits.Markley(),
    variables=@variables begin
        a ~ LogUniform(0.01, 100); e ~ Uniform(0, 0.99)
        ω ~ Uniform(0,2pi); i ~ Sine(); Ω ~ Uniform(0,2pi)
        θ ~ Uniform(0,2pi)
        tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
        mass ~ LogUniform(0.01, 1000)
    end)
sys = System(name="target_1", companions=[b], observations=[gaia_dr4_obs],
    variables=@variables begin; M = 1.0; plx ~ Uniform(0.01,100); end)

t = [-1.6495356151026312, -6.491790601708638, 0.029738626170779867,
     0.0832180348920664, 0.010814871781942643, -0.04826697440973689,
     -3.357219015634567, 0.1154498724014276, -0.012059033319817353,
     -0.4030642907378578, -2.391709856749899, 0.21726752142249556, -4.866247384968242]

# Test primal (no autodiff)
println("=== Primal (autodiff=false) ===")
model_no_ad = Octofitter.LogDensityModel(sys, verbosity=4, autodiff=false)
result = model_no_ad.ℓπcallback(t)
println("Result: $result")
println("Expected ≈ 132.50998916711984")
@assert isapprox(result, 132.50998916711984, atol=1e-6) "Primal result mismatch!"
println("Primal benchmark:")
display(@benchmark $model_no_ad.ℓπcallback($t))
println()

# Test gradient
println("\n=== Gradient (autodiff=default) ===")
model = Octofitter.LogDensityModel(sys, verbosity=4)
primal, grad = model.∇ℓπcallback(t)
println("Primal: $primal")
println("Expected ≈ 132.50998916711984")
@assert isapprox(primal, 132.50998916711984, atol=1e-6) "Gradient primal mismatch!"

expected_grad = [-0.0005554893442992759, 1.592991672477338e-7, -4.298554218188531e-5, 6.223104889555368e-5, -0.020372812771732196, 0.02627009994154831, 8.459059177186312e-6, -4.59962866458663e-6, 6.786563689957503e-6, 1.1863688740509915e-5, 9.518254022400363e-6, -2.811051135115372e-6, -6.8209395289661745e-6]
println("Gradient: $grad")
println("Expected: $expected_grad")
@assert isapprox(grad, expected_grad, rtol=1e-4) "Gradient mismatch!"
println("Gradient benchmark:")
display(@benchmark $model.∇ℓπcallback($t))
