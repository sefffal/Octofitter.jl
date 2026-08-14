#=
Setup for the Gaia DR4 completeness study. Run this once, on a node that has
network access, before submitting the array job.

It queries the Gaia DR3 archive and GOST, reads the target's measured
along-scan noise budget, and caches all of it to `dr4_target.csv` +
`dr4_transits.csv` beside these scripts — so the array tasks need neither the
network nor the ~14 GB G23H catalog. It then builds the template model once,
which precompiles the code path every trial takes.

Usage:
    julia --project=. setup.jl

Set `OCTO_G23H_CATALOG=1` to read the real G23H catalog rather than the
recorded row for this source (same numbers, ~70 s).
=#

include(joinpath(@__DIR__, "common.jl"))

use_catalog = lowercase(get(ENV, "OCTO_G23H_CATALOG", "0")) ∉ ("0", "", "false", "no")

println("Querying Gaia DR3, the G23H noise calibration, and the GOST scan forecast...")
target = load_target(@__DIR__; refresh=true,
    catalog=use_catalog ? nothing : CATALOG_STANDIN)
describe_target(target)

println()
println("Forecast DR4 transits: $(target.n_transits) over MJD " *
        "$(round(minimum(target.transits.epoch), digits=1))–" *
        "$(round(maximum(target.transits.epoch), digits=1)) " *
        "($(round((maximum(target.transits.epoch) - minimum(target.transits.epoch)) / 365.25, digits=2)) yr)")
println("Cached to $(target_csv(@__DIR__)) and $(transits_csv(@__DIR__))")
println()
println("NOTE: GOST forecasts every *scheduled* transit. Real DR4 loses some to")
println("dead time and more to AGIS's outlier rejection, which is harshest for")
println("bright stars — Gaia-4 goes 122 forecast → 109 published → 93 used. Drop")
println("rows from dr4_transits.csv if that matters for your question.")

# Build the template model once: this precompiles the likelihood and gradient
# code path, and fails here rather than in 720 array tasks if anything is off.
println()
println("Building the template model...")
sys = build_system(target)
model = Octofitter.LogDensityModel(sys; verbosity=1)
println(model)

jobs = completeness_jobs(; masses=MASSES, separations=SEPARATIONS, n_trials=N_TRIALS)
println()
println("Grid: $(length(MASSES)) masses × $(length(SEPARATIONS)) separations × " *
        "$(N_TRIALS) trials = $(length(jobs)) jobs" * (QUICK ? "  [OCTO_COMPLETENESS_QUICK]" : ""))
println("Set `#SBATCH --array=1-$(length(jobs))` in submit.sh, then `sbatch submit.sh`.")
