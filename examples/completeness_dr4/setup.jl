#=
Setup script for Gaia DR4 completeness mapping demo.
Run this on the login node to pre-cache GOST scan law and precompile.

Usage:
    julia --project=/scratch/wthompso/completeness_dr4 setup.jl
=#

using Octofitter, Distributions
using DataFrames, CSV

# ── Target star ──
# Use a bright nearby star with good Gaia DR3 solution.
# This determines the scan law and noise properties.
gaia_id = 5064625130502952704

# ── Query DR3 catalog and GOST scan law ──
# These results are cached locally, so subsequent runs skip the HTTP queries.
println("Querying Gaia DR3 catalog...")
dr3 = Octofitter._query_gaia_dr3(; gaia_id)
println("  ra=$(dr3.ra), dec=$(dr3.dec), plx=$(dr3.parallax)")

println("Querying GOST scan law (DR4 baseline)...")
gost = DataFrame(Octofitter.GOST_forecast(dr3.ra, dr3.dec; baseline=:dr4))
println("  $(size(gost, 1)) visibility windows")

# ── Noise parameters ──
# For a representative simulation, use UEVA_single approach:
# Set centroid_pos_error_al = σ_true and astrometric_jitter ≈ 0.
# See Kiefer et al 2025 / Thompson et al in-prep for details.
σ_att = 0.04   # attitude noise [mas]
σ_AL = 0.04    # along-scan CCD noise [mas]
σ_cal = 0.04   # calibration noise [mas]
σ_true = sqrt(σ_att^2 + σ_AL^2 + σ_cal^2)
println("  σ_true = $(round(σ_true, sigdigits=3)) mas per scan")

println("\nSetup complete. GOST cache saved locally.")
println("You can now submit the array job.")
