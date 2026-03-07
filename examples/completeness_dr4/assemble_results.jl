#=
Assemble completeness trial results and generate the completeness map plot.

Usage:
    julia --project=/scratch/wthompso/completeness_dr4 assemble_results.jl

Reads serialized CompletenessResult files from results/ directory,
assembles them into a CompletenessMap, saves the plot, and serializes
the map for use in documentation.
=#

using Octofitter
using CairoMakie
using Serialization

# Grid definition — must match completeness_trial.jl
const MASSES      = 10.0 .^ range(-1, 2, length=12)
const SEPARATIONS = 10.0 .^ range(-0.3, 1.7, length=12)

# ── Load results ──
result_dir = "results"
result_files = filter(f -> endswith(f, ".jls"), readdir(result_dir, join=true))
@info "Found $(length(result_files)) result files"

results = CompletenessResult[]
n_errors = 0
for f in result_files
    try
        r = deserialize(f)
        push!(results, r)
    catch e
        n_errors += 1
        @warn "Failed to load $f" exception=e
    end
end

@info "Loaded $(length(results)) results ($(n_errors) errors)"

# ── Assemble ──
cmap = assemble_completeness(results; masses=MASSES, separations=SEPARATIONS)

@info "Completeness summary" total_trials=sum(cmap.n_total) total_detected=sum(cmap.n_detected) overall_rate=round(sum(cmap.n_detected)/max(sum(cmap.n_total),1), digits=3)

# ── Save completeness map data ──
serialize("completeness_map.jls", cmap)
@info "Saved completeness map to completeness_map.jls"

# ── Plot ──
fig = Octofitter.completenessplot(
    cmap,
    "dr4_completeness_map.png";
    title="Gaia DR4 Completeness (source $5064625130502952704)",
    show_counts=true,
)
@info "Saved plot to dr4_completeness_map.png"

# Also save without counts overlay for docs
Octofitter.completenessplot(
    cmap,
    "dr4_completeness_map_clean.png";
    title="Gaia DR4 Detection Completeness",
)
@info "Saved clean plot to dr4_completeness_map_clean.png"

println("\nCompleteness matrix (rows=mass, cols=separation):")
display(round.(cmap.completeness, digits=2))
