# Minimal reproduction of the Enzyme IllegalTypeAnalysisException
# on Julia 1.12 for GaiaDR4AstromObs.
#
# Run with: julia +1.12 --project=. test/repro_gaia_enzyme.jl

using Octofitter
using Distributions
using Statistics

# --- Build minimal GaiaDR4 model ---
N_epochs = 10
mock_epochs = collect(range(58000.0, 59000.0, length=N_epochs))

# Mock Earth position data
mock_xyz = Octofitter.TypedTables.Table(
    x = fill(0.5, N_epochs),
    y = fill(0.5, N_epochs),
    z = fill(0.2, N_epochs),
    vx = fill(0.0, N_epochs),
    vy = fill(0.0, N_epochs),
    vz = fill(0.0, N_epochs),
)

mock_table = Octofitter.TypedTables.Table(
    epoch = mock_epochs,
    scan_pos_angle = collect(range(0.0, 2π, length=N_epochs)),
    centroid_pos_al = zeros(N_epochs),
    centroid_pos_error_al = fill(0.1, N_epochs),
    parallax_factor_al = fill(0.5, N_epochs),
    outlier_flag = fill(false, N_epochs),
    xyz = collect(eachrow(mock_xyz)),
)

mock_gaia_sol = (
    ra = 180.0, dec = 45.0,
    pmra = 5.0, pmdec = -10.0, parallax = 20.0,
)

# Use the keyword-based outer constructor, which handles inner field computation
ref_epoch_mjd = 57936.375
gaia_obs = GaiaDR4AstromObs(
    mock_table;
    gaia_id=123456789,
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)
        ra_offset_mas ~ Normal(0, 100)
        dec_offset_mas ~ Normal(0, 100)
        pmra ~ Normal(5, 10)
        pmdec ~ Normal(-10, 10)
        plx = system.plx
        ref_epoch = $ref_epoch_mjd
    end,
    name="GaiaDR4",
    primary_star_perturbation=false,
)

orbit_ref_epoch = mean(mock_epochs)
b = Planet(
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
        tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
        mass ~ LogUniform(1, 100)
    end
)

sys = System(
    name="gaia_test",
    companions=[b],
    observations=[gaia_obs],
    variables=@variables begin
        M = 1.0
        plx ~ Uniform(10, 30)
    end
)

println("Julia version: ", VERSION)
println("Building LogDensityModel...")
try
    model = Octofitter.LogDensityModel(sys, verbosity=2)
    println("SUCCESS: model built without errors")
catch e
    println("ERROR: ", typeof(e))
    println(sprint(showerror, e; context=:limit=>true))

    # If it's an Enzyme error, try with strictAliasing=false
    if e isa Exception
        println("\n\nRetrying with Enzyme.API.strictAliasing!(false)...")
        import Enzyme
        Enzyme.API.strictAliasing!(false)
        try
            model2 = Octofitter.LogDensityModel(sys, verbosity=2)
            println("SUCCESS with strictAliasing=false!")
        catch e2
            println("STILL FAILED: ", typeof(e2))
            println(sprint(showerror, e2; context=:limit=>true))
        end
    end
end
