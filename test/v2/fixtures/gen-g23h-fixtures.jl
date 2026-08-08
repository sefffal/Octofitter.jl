# Generate the checked-in G23H fixtures, and the Octofitter v1 reference
# values the v2 port is gated against.
#
# Run under an environment that dev-links the **v1** clones (Octofitter v8.2.5
# + PlanetOrbits v0.11.4), e.g.
#
#     julia --project=benchenv-v1 test/v2/fixtures/gen-g23h-fixtures.jl
#
# It needs, and is the only thing that needs, data that is not in this repo:
#
#   CATALOG   the G23H catalog (Arrow, ~14 GB) — one row is extracted
#   SIDECAR   the DR2 matched-transit sidecar (Arrow)
#   GOST      a cached GOST scan forecast for the target
#   Hipparcos IAD, via the `Hipparcos_IAD` data dependency
#   DE440,    via the `DE440_Ephemeris` data dependency (Earth ephemeris)
#
# Everything it writes into this directory is small, offline, and is what the
# test suite reads. Regenerate only when the reference model changes; the
# outputs are committed.
#
# Target: HIP 1475 / Gaia DR3 385334230892516480 (GX And). Chosen because it
# is a real production-shaped row — Hipparcos present, RV present, the
# σ_AL/σ_att/σ_calib calibration present — with a locally cached scan forecast
# whose transit pool comfortably covers the catalog's matched-transit count.

using Octofitter, PlanetOrbits, Distributions, StaticArrays, Random
using Octofitter: Arrow, Tables, TypedTables, CSV
using Octofitter.TypedTables: Table, FlexTable
using DataFrames: DataFrame

const OUT = @__DIR__
const CATALOG = get(ENV, "G23H_CATALOG",
    "/arc/projects/planet-astrometry/gaia/GaiaHipparcosUEVA-joint-catalog-v6.1.0.arrow")
const SIDECAR = get(ENV, "G23H_DR2_SIDECAR",
    "/arc/projects/planet-astrometry/gaia/G23H-v1.0.dr2_matched_observations.feather")
const GOST = get(ENV, "G23H_GOST",
    "/arc/projects/planet-astrometry/gaia/15pc-sample-1/GOST-4.613226257557736-44.02478674398518.csv")
const GAIA_ID = 385334230892516480
const HIP_ID = 1475

# Columns G23HObs reads. Anything else in the catalog is the production
# model's business, not this observation's.
const CAT_COLS = (
    :gaia_source_id, :hip_id, :ra, :dec, :parallax,
    :epoch_ra_hip, :epoch_dec_hip, :epoch_ra_hg, :epoch_dec_hg,
    :epoch_ra_dr2, :epoch_dec_dr2, :epoch_ra_dr32, :epoch_dec_dr32,
    :epoch_ra_dr3, :epoch_dec_dr3,
    :pmra_hip, :pmdec_hip, :pmra_hg, :pmdec_hg, :pmra_dr2, :pmdec_dr2,
    :pmra_dr32, :pmdec_dr32, :pmra_dr3, :pmdec_dr3,
    :pmra_hip_error, :pmdec_hip_error, :pmra_hg_error, :pmdec_hg_error,
    :pmra_dr2_error, :pmdec_dr2_error, :pmra_dr32_error, :pmdec_dr32_error,
    :pmra_dr3_error, :pmdec_dr3_error,
    :pmra_pmdec_hip, :pmra_pmdec_hg, :pmra_pmdec_dr2, :pmra_pmdec_dr32, :pmra_pmdec_dr3,
    :rho_dr2_dr3,
    :ra_error_central_dr3, :dec_error_central_dr3,
    :ra_error_central_dr2, :dec_error_central_dr2,
    :ra_dec_corr_central_dr3, :ra_dec_corr_central_dr2,
    :astrometric_excess_noise_dr3, :astrometric_matched_transits_dr3,
    :astrometric_n_good_obs_al_dr3, :astrometric_params_solved_dr3,
    :astrometric_chi2_al_dr3, :phot_g_mean_mag_dr3, :ruwe_dr3,
    :nonlinear_dpmra, :nonlinear_dpmdec,
    :sig_AL, :sig_AL_sigma, :sig_att_radec, :sig_att_radec_sigma,
    :sig_cal, :sig_cal_sigma,
    :rv_ln_uncert_dr3, :rv_ln_uncert_err_dr3, :radial_velocity_error, :rv_nb_transits,
)

# ── 1. the catalog row ────────────────────────────────────────────────────
@info "Reading the catalog row" CATALOG
cat_row = let t = Arrow.Table(CATALOG)
    i = findfirst(==(GAIA_ID), Tables.getcolumn(t, :gaia_source_id))
    isnothing(i) && error("Gaia source $GAIA_ID not found in $CATALOG")
    full = NamedTuple(Table(t)[i])
    vals = map(CAT_COLS) do c
        v = getproperty(full, c)
        (ismissing(v) || (v isa AbstractFloat && isnan(v))) &&
            error("catalog column $c is missing/NaN for this source; pick another target")
        # Integer-valued catalog fields stay integers: `gaia_source_id` does
        # not survive a round trip through Float64.
        v isa Integer ? Int64(v) : Float64(v)
    end
    NamedTuple{CAT_COLS}(vals)
end
CSV.write(joinpath(OUT, "g23h-catalog-row.csv"), Table(map(v -> [v], cat_row)))

# ── 2. the DR2 matched-transit count ──────────────────────────────────────
n_dr2_matched = let t = Arrow.Table(SIDECAR)
    col = hasproperty(t, :astrometric_matched_observations_dr2) ?
          :astrometric_matched_observations_dr2 : :astrometric_matched_observations
    i = findfirst(==(GAIA_ID), Tables.getcolumn(t, :gaia_source_id))
    isnothing(i) && error("Gaia source $GAIA_ID not found in the DR2 sidecar")
    Int(Tables.getcolumn(t, col)[i])
end
@info "DR2 matched transits" n_dr2_matched
CSV.write(joinpath(OUT, "g23h-dr2-transits.csv"),
    Table(gaia_source_id=[GAIA_ID], astrometric_matched_observations_dr2=[n_dr2_matched]))

# ── 3. the Hipparcos IAD, exactly as v1 reconstructs it ───────────────────
# `G23HObs` mutates this table in place (the +0.140 mas recalibration), so
# load a separate copy for the fixture *before* building the observation.
J2000_mjd = 51544.5
_yr2mjd(y) = (y - 2000) * Octofitter.julian_year + J2000_mjd
hip_like = Octofitter.HipparcosIADLikelihood(;
    hip_id=HIP_ID,
    ref_epoch_ra=_yr2mjd(cat_row.epoch_ra_hip),
    ref_epoch_dec=_yr2mjd(cat_row.epoch_dec_hip),
    variables=(Octofitter.Priors(), Octofitter.Derived()))
hip_cols = (:iorb, :epoch_yrs, :parf, :cosϕ, :sinϕ, :res, :sres, :reject,
    :sres_renorm, :epoch, :x, :y, :z, :rv_kms, :plx_vs_time, :Δα✱, :Δδ,
    :scanAngle_rad, :parallaxFactorAlongScan, :proj_meas_alongscan)
CSV.write(joinpath(OUT, "g23h-hip-iad.csv"),
    Table(NamedTuple{hip_cols}(map(c -> collect(getproperty(hip_like.table, c)), hip_cols))))
hip_sol_cols = (:plx, :pm_ra, :pm_de, :radeg, :dedeg, :isol_n, :f2, :hp)
CSV.write(joinpath(OUT, "g23h-hip-sol.csv"),
    Table(NamedTuple{hip_sol_cols}(map(c -> [Float64(getproperty(hip_like.hip_sol, c))], hip_sol_cols))))
@info "Hipparcos IAD" n_transits = length(hip_like.table.epoch)

# ── 4. build the v1 observation ───────────────────────────────────────────
# v1 fetches its own scan forecast; point its on-disk cache at the local copy
# by running in a directory holding the file under the name it looks for.
workdir = mktempdir()
cp(GOST, joinpath(workdir, "GOST-$(cat_row.ra)-$(cat_row.dec)-dr3.csv"))
cd(workdir)

# v1 selects its row with DataFrame-style `catalog[idx, :]`.
cat_table = DataFrame(map(v -> [v], cat_row))
sidecar_table = Table(gaia_source_id=[GAIA_ID],
    astrometric_matched_observations_dr2=[n_dr2_matched])

like_v1 = Octofitter.G23HObs(;
    gaia_id=GAIA_ID,
    catalog=cat_table,
    dr2_transits_catalog=sidecar_table,
    ueva_mode=:RUWE,
    include_rv=true,
)
cd(OUT)

n_pool = length(like_v1.gaia_table.epoch)
@info "Gaia forecast pool" n_pool matched = Int(cat_row.astrometric_matched_transits_dr3)

# The forecast pool, in the shape `G23HObs(; forecast_table=…)` takes. v2
# re-derives cosϕ/sinϕ and re-applies the span/gap masks, so this round-trips
# to exactly the same pool without needing GOST or an ephemeris.
CSV.write(joinpath(OUT, "g23h-forecast.csv"),
    Table(epoch=collect(like_v1.gaia_table.epoch),
        scanAngle_rad=collect(like_v1.gaia_table.scanAngle_rad),
        parallaxFactorAlongScan=collect(like_v1.gaia_table.parallaxFactorAlongScan)))

# ── 5. the evaluation point ───────────────────────────────────────────────
# Fixed, deterministic, and written out so the v2 test drives exactly these
# numbers. The epoch selections are supplied explicitly rather than sampled:
# the gate is on the sky-path and reduction arithmetic, not on the epoch
# marginalization prior, and a shared selection makes the two runs comparable.
const N_MATCHED = Int(cat_row.astrometric_matched_transits_dr3)
const N_RV = Int(cat_row.rv_nb_transits)
transits = collect(1:min(N_MATCHED, n_pool))
transits_dr2 = collect(1:min(n_dr2_matched, n_pool))
transits_rv = collect(1:min(N_RV, length(transits)))
transit_priorities = collect(range(-1.0, 1.0, length=n_pool))

const θ_OBS = (;
    σ_AL=0.0721, σ_att=0.0854, σ_calib=0.5031,
    hip_iad_jitter=0.35,
    iad_Δra=0.12, iad_Δdec=-0.08, iad_Δplx=0.03,
    iad_pmra=cat_row.pmra_hip + 0.4, iad_pmdec=cat_row.pmdec_hip - 0.25,
    σ_rv_per_transit=0.31,
    transit_priorities, transits, transits_dr2, transits_rv,
)

# Single luminous companion, non-absolute frame — the configuration in which
# v1's per-companion linear photocentre coefficient is algebraically exact, so
# the two implementations must agree to roundoff.
const M_TOT = 0.486                      # M⊙, host + companion
const M_COMP = 0.030                     # M⊙
const F_G = 0.037                        # G-band contrast
const F_HP = 0.021                       # Hp-band contrast
const ELEMENTS = (a=4.7, e=0.23, i=0.94, ω=1.31, Ω=2.44, tp=52000.0)

orbit_v1 = PlanetOrbits.orbit(; M=M_TOT, plx=cat_row.parallax, ELEMENTS...)
θ_system = (; M=M_TOT, plx=cat_row.parallax,
    pmra=cat_row.pmra_dr3, pmdec=cat_row.pmdec_dr3,
    planets=(; b=(; mass=M_COMP / Octofitter.mjup2msol)))
θ_obs_v1 = merge(θ_OBS, (; fluxratio=[F_G], fluxratio_hip=[F_HP]))

ctx = Octofitter.SystemObservationContext(θ_system, θ_obs_v1, (orbit_v1,),
    (PlanetOrbits.orbitsolve.(Ref(orbit_v1), [0.0]),), 1)

sim = Octofitter.simulate(like_v1, θ_system, θ_obs_v1, (orbit_v1,), ctx.orbit_solutions, 1)
ll = Octofitter.ln_like(like_v1, ctx)
@info "v1 reference" ll

CSV.write(joinpath(OUT, "g23h-v1-reference.csv"), Table(
    name=["ll", "mu_h_ra", "mu_h_dec", "mu_hg_ra", "mu_hg_dec",
        "mu_dr2_ra", "mu_dr2_dec", "mu_dr32_ra", "mu_dr32_dec",
        "mu_dr3_ra", "mu_dr3_dec", "UEVA_model", "UEVA_unc", "mu_1_3",
        "sample_variance", "deflation_factor_dr3", "hip_bias_pm_sq",
        "delta_alpha_dr3", "delta_delta_dr3", "delta_pmra_dr3", "delta_pmdec_dr3",
        "n_pool", "n_hip"],
    value=[ll, sim.μ_h[1], sim.μ_h[2], sim.μ_hg[1], sim.μ_hg[2],
        sim.μ_dr2[1], sim.μ_dr2[2], sim.μ_dr32[1], sim.μ_dr32[2],
        sim.μ_dr3[1], sim.μ_dr3[2], sim.UEVA_model, sim.UEVA_unc, sim.μ_1_3,
        sim.sample_variance, sim.deflation_factor_dr3, sim.hip_bias_pm_sq,
        sim.Δα_dr3, sim.Δδ_dr3, sim.Δpmra_dr3, sim.Δpmdec_dr3,
        Float64(n_pool), Float64(length(hip_like.table.epoch))]))

# The evaluation point, so the v2 test does not have to restate it.
CSV.write(joinpath(OUT, "g23h-eval-point.csv"), Table(
    name=["M_tot", "M_comp", "fluxratio_G", "fluxratio_Hp",
        "a", "e", "i", "omega", "Omega", "tp",
        "sigma_AL", "sigma_att", "sigma_calib", "hip_iad_jitter",
        "iad_dra", "iad_ddec", "iad_dplx", "iad_pmra", "iad_pmdec",
        "sigma_rv_per_transit", "n_transits", "n_transits_dr2", "n_transits_rv"],
    value=[M_TOT, M_COMP, F_G, F_HP,
        ELEMENTS.a, ELEMENTS.e, ELEMENTS.i, ELEMENTS.ω, ELEMENTS.Ω, ELEMENTS.tp,
        θ_OBS.σ_AL, θ_OBS.σ_att, θ_OBS.σ_calib, θ_OBS.hip_iad_jitter,
        θ_OBS.iad_Δra, θ_OBS.iad_Δdec, θ_OBS.iad_Δplx, θ_OBS.iad_pmra, θ_OBS.iad_pmdec,
        θ_OBS.σ_rv_per_transit,
        Float64(length(transits)), Float64(length(transits_dr2)), Float64(length(transits_rv))]))

@info "wrote fixtures" OUT
