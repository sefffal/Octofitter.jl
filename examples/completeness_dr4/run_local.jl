#=
The whole Gaia DR4 completeness workflow in one process, for a laptop.

    OCTO_COMPLETENESS_QUICK=1 julia --project=. --threads=auto run_local.jl

With `OCTO_COMPLETENESS_QUICK=1` this is a 2 × 2 grid, one trial per cell and
short chains — a few minutes, and the right way to check the model and the
detection criterion before spending a cluster allocation on the full 12 × 12 ×
5 grid. Without it, it runs the production grid serially, which you do not
want; use `submit.sh` for that.

`completeness_map` is just `completeness_jobs` → `run_completeness_trial` →
`assemble_completeness` in sequence, and returns the raw results as well as
the map, so a different detection threshold costs nothing.
=#

using Serialization
using CairoMakie

include(joinpath(@__DIR__, "common.jl"))

target = load_target(@__DIR__)
describe_target(target)

sys = build_system(target)

n_jobs = length(MASSES) * length(SEPARATIONS) * N_TRIALS
@info "Running the completeness grid locally" n_masses = length(MASSES) n_separations = length(SEPARATIONS) N_TRIALS n_jobs quick = QUICK

cmap, results = completeness_map(sys, sampler, detection;
    inject, masses=MASSES, separations=SEPARATIONS, n_trials=N_TRIALS,
    add_noise=true, verbosity=2)   # 2 so each trial announces itself

println("\nCompleteness matrix (rows = mass [MJup], cols = separation [AU]):")
println("  masses      = ", round.(MASSES ./ mjup, digits=2))
println("  separations = ", round.(SEPARATIONS, digits=2))
display(round.(cmap.completeness, digits=2))
println()

# What the detection criterion actually saw, trial by trial.
println("Per-trial summary:")
for r in results
    m = vec(r.chain["b_mass"]) ./ mjup
    println("  mass = ", rpad(round(r.job.mass / mjup, digits=2), 7),
        "MJup  sep = ", rpad(round(r.job.separation, digits=2), 6),
        "AU  ->  median ", rpad(round(median(m), digits=2), 8),
        " 5th pct ", rpad(round(quantile(m, 0.05), digits=3), 8),
        " MJup   detected = ", detection(r.chain, r.θ_true))
end

serialize(joinpath(@__DIR__, "completeness_map.jls"), cmap)
completenessplot(cmap, joinpath(@__DIR__, "dr4_completeness_map.png");
    title="Gaia DR4 completeness (source $(target.gaia_id))", show_counts=true)
@info "Saved completeness_map.jls and dr4_completeness_map.png"
