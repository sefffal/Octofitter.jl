# Regenerate `test/G23H-test-subset.feather`, the offline Hipparcos-Gaia catalog
# subset that the documentation and the test suite read in place of the ~14 GB
# G23H catalog DataDep.
#
#     julia --project=. test/v2/fixtures/gen-g23h-docs-subset.jl
#
# It needs data that is not in this repo, and is the only thing that does:
#
#   G23H_CATALOG      the G23H catalog (Arrow, ~14 GB)
#   G23H_DR2_SIDECAR  the DR2 matched-transit sidecar (Arrow, ~300 MB)
#
# The output is a three-row Arrow file of about 20 kB, and it is committed.
#
# Rows:
#   Gaia DR3 2738776816458107136 (HIP 384)     — `g23h.md`
#   Gaia DR3 756291174721509376  (HD 91312)    — `pma.md`, `astrom-pma-rv.md`,
#                                                `limits.md`, `data-simulation.md`
#   Gaia DR3 5164707970261890560 (eps Eridani) — `rv.md`
#
# The DR2 matched-transit count is folded in as a column, so a page reading this
# file needs neither the catalog DataDep nor the sidecar DataDep. `G23HObs`
# returns the catalog untouched when that column is already present.
#
# The column set is taken from the previous subset rather than restated, so a
# regeneration cannot silently widen or narrow what the docs exercise.

using Octofitter
using Octofitter: Arrow, Tables
using Octofitter.TypedTables: Table

const OUT = joinpath(@__DIR__, "..", "..", "G23H-test-subset.feather")
const CATALOG = get(ENV, "G23H_CATALOG",
    "/arc/projects/planet-astrometry/gaia/GaiaHipparcosUEVA-joint-catalog-v6.1.0.arrow")
const SIDECAR = get(ENV, "G23H_DR2_SIDECAR",
    "/arc/projects/planet-astrometry/gaia/G23H-v1.0.dr2_matched_observations.feather")

const IDS = [
    2738776816458107136,   # HIP 384
    756291174721509376,    # HD 91312
    5164707970261890560,   # eps Eridani
]

# The columns the previous subset carried. Regenerating must not change the
# schema out from under the pages that read it.
const COLS = let t = Arrow.Table(OUT)
    filter(!=(:astrometric_matched_observations_dr2), collect(Tables.columnnames(t)))
end
@info "Existing subset" ncols = length(COLS) rows = collect(Arrow.Table(OUT).gaia_source_id)

# ── the catalog rows ──────────────────────────────────────────────────────
rows = let t = Arrow.Table(CATALOG)
    ids = Tables.getcolumn(t, :gaia_source_id)
    map(IDS) do id
        i = findfirst(==(id), ids)
        isnothing(i) && error("Gaia source $id not found in $CATALOG")
        full = NamedTuple(Table(t)[i])
        vals = map(COLS) do c
            haskey(full, c) || error("catalog has no column $c; the subset schema has drifted")
            getproperty(full, c)
        end
        NamedTuple{Tuple(COLS)}(Tuple(vals))
    end
end

# Everything `G23HObs` reads must actually be there for every target. A missing
# value here surfaces as a confusing failure deep inside the likelihood, so
# check it while we still know which column it was.
const REQUIRED = (
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
)
for (id, row) in zip(IDS, rows), c in REQUIRED
    v = getproperty(row, c)
    (ismissing(v) || (v isa AbstractFloat && isnan(v))) &&
        error("catalog column $c is missing/NaN for Gaia source $id")
end

# ── the DR2 matched-transit counts ────────────────────────────────────────
n_dr2 = let t = Arrow.Table(SIDECAR)
    col = hasproperty(t, :astrometric_matched_observations_dr2) ?
          :astrometric_matched_observations_dr2 : :astrometric_matched_observations
    ids = Tables.getcolumn(t, :gaia_source_id)
    map(IDS) do id
        i = findfirst(==(id), ids)
        isnothing(i) && error("Gaia source $id not found in the DR2 sidecar $SIDECAR")
        Int32(Tables.getcolumn(t, col)[i])
    end
end
@info "DR2 matched transits" pairs = collect(zip(IDS, n_dr2))

# ── write ─────────────────────────────────────────────────────────────────
# Via a temporary file: Arrow memory-maps what it reads, and `COLS` above was
# taken from the very file being replaced.
outcols = NamedTuple{(Tuple(COLS)..., :astrometric_matched_observations_dr2)}((
    (map(r -> getproperty(r, c), rows) for c in COLS)..., n_dr2))
tmp = tempname() * ".feather"
Arrow.write(tmp, outcols)
mv(tmp, OUT; force=true)

@info "wrote" OUT size_kB = round(filesize(OUT) / 1024, digits=1) nrows = length(IDS)
for (id, row) in zip(IDS, rows)
    @info "  row" gaia_id = id hip_id = row.hip_id parallax = row.parallax G = row.phot_g_mean_mag_dr3
end
