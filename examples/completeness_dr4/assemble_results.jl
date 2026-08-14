#=
Assemble the completeness trials and draw the completeness map.

Usage:
    julia --project=. assemble_results.jl

Reads the serialized `CompletenessResult`s from `results/`, applies the
detection criterion from `common.jl`, and writes the map and its plots.

Because each result carries its whole chain, a different criterion can be
applied here without re-running a single trial — see the second map below.
=#

using Serialization
using CairoMakie

include(joinpath(@__DIR__, "common.jl"))

# ── Load results ──────────────────────────────────────────────────────
result_dir = joinpath(@__DIR__, "results")
isdir(result_dir) || error("No `results/` directory beside $(@__DIR__) — run the trials first (submit.sh, or completeness_trial.jl).")
result_files = filter(f -> endswith(f, ".jls"), readdir(result_dir, join=true))
isempty(result_files) && error("`$result_dir` holds no .jls results.")
@info "Found $(length(result_files)) result files"

results = CompletenessResult[]
failed = String[]
for f in result_files
    try
        push!(results, deserialize(f))
    catch e
        push!(failed, f)
        @warn "Failed to load $f" exception = e
    end
end
@info "Loaded $(length(results)) results ($(length(failed)) failed to deserialize)"

# ── Assemble ──────────────────────────────────────────────────────────
# Detection is applied here, not in the trial.
cmap = assemble_completeness(results, detection; masses=MASSES, separations=SEPARATIONS)

@info "Completeness summary" total_trials = sum(cmap.n_total) total_detected = sum(cmap.n_detected) overall_rate = round(sum(cmap.n_detected) / max(sum(cmap.n_total), 1), digits=3)

# Missing (e.g. timed out) trials are simply absent, and `assemble_completeness`
# reports `completeness = 0` for a cell with no trials — indistinguishable from
# a cell where nothing was recovered. `n_total` carries the distinction.
n_empty = count(iszero, cmap.n_total)
n_empty > 0 && @warn "$n_empty of $(length(cmap.n_total)) grid cells have no trials; their completeness reads 0 but is undefined."

serialize(joinpath(@__DIR__, "completeness_map.jls"), cmap)
@info "Saved completeness map to completeness_map.jls"

# ── Plot ──────────────────────────────────────────────────────────────
target = load_target(@__DIR__)

completenessplot(cmap, joinpath(@__DIR__, "dr4_completeness_map.png");
    title="Gaia DR4 completeness (source $(target.gaia_id))",
    show_counts=true)
@info "Saved plot to dr4_completeness_map.png"

# Without the counts overlay, for the documentation.
completenessplot(cmap, joinpath(@__DIR__, "dr4_completeness_map_clean.png");
    title="Gaia DR4 detection completeness")
@info "Saved clean plot to dr4_completeness_map_clean.png"

# ── Re-thresholding, for free ─────────────────────────────────────────
# The same trials under a stricter criterion: the 5th percentile of the mass
# posterior must clear 1 M_Jup, with no requirement that the median be near
# the truth.
cmap_strict = assemble_completeness(results,
    (chain, θ) -> quantile(vec(chain["b_mass"]), 0.05) > 1mjup;
    masses=MASSES, separations=SEPARATIONS)
@info "Strict criterion" overall_rate = round(sum(cmap_strict.n_detected) / max(sum(cmap_strict.n_total), 1), digits=3)
completenessplot(cmap_strict, joinpath(@__DIR__, "dr4_completeness_map_strict.png");
    title="Gaia DR4 completeness — 5th percentile > 1 MJup")

println("\nCompleteness matrix (rows = mass, cols = separation):")
display(round.(cmap.completeness, digits=2))
println()
