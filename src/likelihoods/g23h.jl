using Arrow

"""
    G23HObs(; gaia_id=nothing, hip_id=nothing, include_rv=true, ...)
    G23HObs(; gaia_id=nothing, hip_id=nothing, include_rv=true, ...)

A likelihood for joint Gaia-Hipparcos astrometry using the G23H catalog, including
proper motion accelerations and unit weight error variance analysis (UEVA).

# Arguments
- `gaia_id`: Gaia DR3 source ID (provide either this or `hip_id`)
- `hip_id`: Hipparcos ID (provide either this or `gaia_id`)
- `catalog`: Path to G23H catalog file, or a loaded DataFrame/Table.
  Defaults to the automatically downloaded G23H catalog via DataDeps.
- `include_rv`: Whether to include Gaia RV variability constraints (default: true)
- `rv_ln_uncert_err_floor`: Minimum adopted uncertainty on the catalog per-transit
  RV ln-σ calibration `rv_ln_uncert_err_dr3` (default: 0.30). The GP calibration is
  Malmquist-biased for cool dwarfs — their calibration bins are dominated by distant
  giants with sharper lines — so the quoted uncertainty can be unrealistically small
  and over-constrain the σ_rv_per_transit prior. Set to `nothing` for the raw value.
- `ueva_mode`: `:RUWE` (default), `:EAN`, or `:none` for astrometric excess noise modeling.
  `:none` opts the star out of the UEVA channel entirely: the `:ueva_dr3` row is not
  added to the observation table, the σ_AL/σ_att/σ_calib nuisance parameters become
  fixed constants rather than catalog-calibrated priors (so no `variables` block is
  required), and the UEVA-driven deflation of the published DR3/DR32 covariances is
  disabled. Use it for sources where the G23H σ_AL/σ_att_radec/σ_cal calibration is
  absent — in the released catalog these are predominantly the very brightest stars
  (median G ≈ 5.7), where the calibration was not extrapolated. Everything else
  (Hipparcos, the DR2/DR3 positions and proper motions, the PM anomaly and the Gaia
  RV channel) is modelled exactly as usual. `:none` also does not read the
  `ruwe_dr3` / `astrometric_excess_noise_dr3` columns at all, so it works on
  catalog rows where those are absent; and `generate_from_params` writes
  `missing` for `ruwe_dr3`, `astrometric_chi2_al_dr3` and
  `astrometric_excess_noise_dr3` in simulated catalogs, since under `:none`
  there is no calibrated σ to back a meaningful value out of.
- `freeze_epochs`: If true, fix Gaia observation epochs for faster sampling (default: false)
- `variables`: Optional custom priors (defaults are set from catalog values)

# Examples
```julia
# Using Gaia DR3 source ID
absastrom = G23HObs(gaia_id=756291174721509376)

# Using Hipparcos ID (automatically resolved to Gaia ID)
absastrom = G23HObs(hip_id=21547)
```

The G23H catalog (~14 GB) is automatically downloaded on first use.

# Variable Priors
Default priors are set automatically from the catalog. Custom priors can be specified:
```julia
variables=@variables begin
    # G-band flux ratio per companion — used for the Gaia DR2/DR3 photocentre
    fluxratio     ~ Product([LogUniform(1e-6, 1e-1) for _ in 1:N_planets])
    # Hp-band flux ratio per companion — used for the Hipparcos abscissa branch.
    # Required for any star with `has_hipparcos`. The Hipparcos likelihood uses
    # this value with the BINARYS atan2 photocentre formula and σ inflation
    # (Leclerc et al. 2023, A&A 672 A82).
    fluxratio_hip ~ Product([LogUniform(1e-6, 1e-1) for _ in 1:N_planets])
    σ_att   ~ LogUniform(0.01, 1.0)   # Attitude error (mas)
    σ_AL    ~ LogUniform(0.01, 1.0)   # Along-scan error (mas)
    σ_calib ~ LogUniform(0.01, 1.0)   # Calibration error (mas)
end
```
"""
struct G23HObs{TTable,TTableH,TTableG,TCat,THip} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    hip_table::TTableH
    gaia_table::TTableG
    catalog::TCat
    hip_sol::THip
    A_prepared_5_hip::Matrix{Float64}
    A_prepared_5_dr2::Matrix{Float64}
    A_prepared_5_dr3::Matrix{Float64}
    include_iad::Bool
    ueva_mode::Symbol
    # Cached weighted pseudo-inverse for the Hipparcos 5-param fit:
    # Q = pinv(A_prepared_5_hip ./ sres) ./ sres', so the per-call weighted
    # LSQ reduces to a single 5×N matrix-vector multiply.  Both A and σ are
    # fixed at construction, so this is initialised eagerly by the
    # constructors (the lazy `_ensure_pinv_5_hip!` path remains only as a
    # fallback and to keep `simulate!` self-contained).  The empty 0×0
    # sentinel indicates "uninitialised".
    _pinv_5_hip::Base.RefValue{Matrix{Float64}}
    # Cached Hipparcos catalog-residual 5-param fit: x_const = Q_hip * residuals
    # where residuals = like.hip_table.res (when include_iad=true) or zero.
    # Used in the all-inactive (n=0) simulate! fast path to skip a per-call
    # matrix-vector multiply.  Initialised eagerly alongside _pinv_5_hip.
    _hip_x_const::Base.RefValue{NTuple{4, Float64}}
    _hip_x_const_initialised::Base.RefValue{Bool}
end

# Diagnostic hook (off by default; zero hot-path cost when nothing).  When set to
# a Vector, ln_like pushes a NamedTuple of the per-channel astrometry residuals
# (catalog vs model μ and marginal σ) for the PM/UEVA MvNormal block — used to
# numerically verify simulator/likelihood agreement (pulls ~ N(0,1) at truth).
const _G23H_DEBUG_PULLS = Ref{Any}(nothing)

# G23H short-circuits on mass == 0 everywhere it consumes per-planet orbit
# solutions (DR2/DR3 skypath perturbations, the Hippacentre combined loop,
# and the Gaia RV simulation all `continue` on zero-mass companions), so
# `make_ln_like` may safely elide the per-epoch Kepler solves for absent
# companions.
requires_solutions_for_zero_mass(::G23HObs) = false


function likelihoodname(like::G23HObs)
    return "G23H"
end

# ──────────────────────────────────────────────────────────────────────
# DR2/DR3 matched-transit sidecar plumbing
#
# The Gaia DR2 astrometric solution shares its input-span start with DR3 but
# stops earlier. The two transit sets overlap but are NOT nested: DR3's
# rebuilt crossmatch commonly used more transits inside the DR2 window than
# DR2 itself matched (audit 2026-07-11: 23 of 153 sampled targets are
# outright infeasible under a nesting assumption, and it is near-binding for
# half the rest), and DR2 occasionally used transits DR3 later dropped
# (n_dr2 > n_dr3 for ~10%). The epoch simulator therefore selects the two
# sets separately — DR3's from the whole AGIS span, DR2's from the DR2
# window — with the DR2 selection sized by the DR2 matched-transit count.
# For bright stars that count additionally includes doubly-downlinked
# transits (two windows per FoV crossing, both tallied by DR2's crossmatch;
# resolved to one in (E)DR3): the DR2 count exceeds the geometric transit
# pool by up to ~1.6× for G ≲ 6, monotonically in G, while DR3's never does.
# These are modeled as repeated epochs with the distinct-crossing count
# marginalized (see the DR2 epoch-selection block in the constructor).
# The DR2/DR3 catalog-covariance correlation is NOT derived from these
# counts; it is adopted directly from the published `rho_dr2_dr3`, against
# which the catalog uncertainties were calibrated.
#
# That count is NOT carried in the published G23H catalog, so it is supplied by
# a mandatory companion sidecar (the `G23H_DR2Transits` DataDep, or the
# `dr2_transits_catalog` keyword) via the column
# `astrometric_matched_observations_dr2`. In Gaia DR2 this field is the number
# of *FoV transits* matched to the source — it was renamed
# `astrometric_matched_transits` in (E)DR3 — so it is directly comparable to
# the DR3 count and needs no per-transit conversion. There is no fallback: if
# the count cannot be resolved for a source, construction errors.
# ──────────────────────────────────────────────────────────────────────

const _G23H_DR2_SIDECAR = Ref{Any}(nothing)
const _G23H_DR2_SIDECAR_LOCK = ReentrantLock()

# Lazily load (and memoize) the sidecar table registered as the
# `G23H_DR2Transits` DataDep. Returns a Tables.jl table.
function _g23h_dr2_sidecar_datadep()
    lock(_G23H_DR2_SIDECAR_LOCK) do
        if isnothing(_G23H_DR2_SIDECAR[])
            dir = datadep"G23H_DR2Transits"
            feathers = filter(f -> endswith(f, ".feather"), readdir(dir; join=true))
            isempty(feathers) && error("No .feather file found in the G23H_DR2Transits DataDep at $dir")
            _G23H_DR2_SIDECAR[] = Arrow.Table(first(feathers))
        end
        return _G23H_DR2_SIDECAR[]
    end
end

# Merge the DR2 matched-transit count for `gaia_id` into the catalog row from
# the mandatory sidecar. No-op if the catalog already carries the column.
# Errors (no fallback) if the sidecar or the source's row cannot be resolved.
function _g23h_merge_dr2_sidecar(catalog, gaia_id, dr2_transits_catalog)
    hasproperty(catalog, :astrometric_matched_observations_dr2) && return catalog
    tbl = if !isnothing(dr2_transits_catalog)
        Tables.istable(dr2_transits_catalog) ? dr2_transits_catalog : Arrow.Table(dr2_transits_catalog)
    else
        _g23h_dr2_sidecar_datadep()
    end
    idx = findfirst(==(gaia_id), Tables.getcolumn(tbl, :gaia_source_id))
    isnothing(idx) && error(
        "Gaia source $gaia_id was not found in the G23H_DR2Transits sidecar. The DR2 " *
        "matched-transit count (`astrometric_matched_observations_dr2`) is required to " *
        "set the DR2/DR3 correlation; there is no fallback.")
    srow = NamedTuple(Table(tbl)[idx])
    # The published G23H_DR2Transits sidecar (CANFAR DOI 26.0016) carries the DR2
    # count under its native Gaia column name `astrometric_matched_observations`
    # (renamed `astrometric_matched_transits` in (E)DR3); accept either that or the
    # `_dr2`-suffixed internal name, and store it internally as the latter.
    n_dr2 = if haskey(srow, :astrometric_matched_observations_dr2)
        srow.astrometric_matched_observations_dr2
    elseif haskey(srow, :astrometric_matched_observations)
        srow.astrometric_matched_observations
    else
        error("The G23H_DR2Transits sidecar is missing the DR2 matched-observation " *
              "count (`astrometric_matched_observations_dr2` or " *
              "`astrometric_matched_observations`).")
    end
    return merge(catalog, (; astrometric_matched_observations_dr2 = n_dr2))
end

# Total size of the DR2 epoch selection (drawn from the DR2-window forecast
# pool), taken from the mandatory DR2 matched-transit count. NOTE: this total
# includes DR2's doubly-downlinked bright-star transits, so it can exceed the
# number of distinct crossings — see `_g23h_select_dr2_epochs`. Errors if it
# is unavailable or non-finite (no fallback by design).
function _g23h_dr2_target_transits(catalog, n_dr3::Integer)
    hasproperty(catalog, :astrometric_matched_observations_dr2) || error(
        "G23HObs requires the Gaia DR2 matched-transit count " *
        "(`astrometric_matched_observations_dr2`) from the G23H_DR2Transits sidecar " *
        "or the `dr2_transits_catalog` keyword; it was not found for this source.")
    v = catalog.astrometric_matched_observations_dr2
    (ismissing(v) || !isfinite(v)) && error(
        "Gaia DR2 `astrometric_matched_observations_dr2` is missing/non-finite for this source.")
    n_dr2 = round(Int, v)
    if n_dr2 > 3 * n_dr3
        # n_dr2 modestly above n_dr3 is genuine (DR2's set is not nested in
        # DR3's, and DR3 dropped some DR2 transits), but a several-fold excess
        # suggests the column holds CCD-level AL observations (~9 per FoV
        # transit) rather than FoV transits. Flag it; downstream clamps to the
        # window pool size anyway.
        @warn "DR2 matched transits ($n_dr2) far exceed DR3 matched transits ($n_dr3); the " *
              "sidecar `astrometric_matched_observations_dr2` may be CCD-level AL " *
              "observations rather than FoV transits."
    end
    return max(n_dr2, 0)
end

# Materialize the DR2 epoch set from the sampled transit priorities:
# `n_distinct` distinct crossings (the top-priority epochs of the DR2-window
# pool) padded with repeats up to `n_total` entries. A repeated index puts
# that epoch into the DR2 least-squares fit twice — the model of a
# doubly-downlinked bright-star transit that DR2's matched-observation count
# tallied twice. Which crossings carry the repeats is taken from the top of
# the priority ordering; by exchangeability of the iid priorities this is
# distributed as a uniform choice among the selected crossings. The returned
# length is always `n_total` (fixed across samples, so chain storage and
# likelihood buffers stay rectangular), except in the pathological case of an
# empty pool.
function _g23h_select_dr2_epochs(priorities::AbstractVector, pool::AbstractVector{<:Integer}, n_distinct::Integer, n_total::Integer)
    (n_distinct <= 0 || n_total <= 0) && return Int[]
    sel = pool[partialsortperm(priorities[pool], 1:n_distinct, rev=true)]  # priority order
    n_rep = n_total - n_distinct
    # mod1 wrap covers the extreme case n_rep > n_distinct (multiplicity ≥ 3).
    reps = n_rep > 0 ? sel[mod1.(1:n_rep, n_distinct)] : Int[]
    return sort!(vcat(sel, reps))
end

function G23HObs(;
        gaia_id=nothing,
        hip_id=nothing,
        scanlaw_table=nothing,
        catalog=joinpath(datadep"G23H_Catalog", "G23H-v1.0.feather"),
        variables::Union{Nothing,Tuple{Priors,Derived}}=nothing,
        include_rv=true,
        rv_ln_uncert_err_floor::Union{Nothing,Real}=0.30,
        ueva_mode::Symbol=:RUWE,
        freeze_epochs=false,
        dr2_transits_catalog=nothing,
        # G-band threshold below which the DR2 duplicate-transit count is
        # marginalized (doubly-downlinked bright-star transits; see the DR2
        # epoch-selection comment). Set to -Inf to disable.
        dr2_dup_gmag_threshold::Real=6.5,
    )
    include_iad=false

    if ueva_mode ∉ (:RUWE, :EAN, :none)
        error("ueva_mode should be :RUWE, :EAN, or :none, was $(ueva_mode)")
    end

    # Validate that exactly one of gaia_id or hip_id is provided
    if isnothing(gaia_id) && isnothing(hip_id)
        error("Either gaia_id or hip_id must be specified")
    end
    if !isnothing(gaia_id) && !isnothing(hip_id)
        error("Specify either gaia_id or hip_id, not both")
    end

    # allow passing in table directly
    if Tables.istable(catalog)
        idx = findfirst(==(gaia_id), catalog.gaia_source_id)
        if isnothing(idx)
            error("The requested gaia source ID $gaia_id was not found in the catlog file $catalog.")
        end
        catalog = NamedTuple(catalog[idx,:])
    else
        # Load the catalog row for this system
        catalog =let t = Arrow.Table(catalog)
            # If hip_id provided, look up gaia_id
            if !isnothing(hip_id)
                hip_matches = findall(==(hip_id), t.hip_id)
                if isempty(hip_matches)
                    error("The requested Hipparcos ID $hip_id was not found in the catalog file $catalog.")
                end
                gaia_id = t.gaia_source_id[hip_matches[1]]
                @info "Resolved HIP $hip_id to Gaia DR3 source ID $gaia_id"
            end
            idx = findfirst(==(gaia_id), t.gaia_source_id)
            if isnothing(idx)
                error("The requested gaia source ID $gaia_id was not found in the catalog file $catalog.")
            end
            NamedTuple(Table(t)[idx])
        end
    end

    # Convert measurement epochs to MJD.
    # Careful: these are Julian years, not decimal years (T. Brant., private communications)
    J2000_mjd = 51544.5 # year J2000 in MJD
    catalog = (;
        catalog...,
        epoch_ra_hip_mjd=(catalog.epoch_ra_hip - 2000) * julian_year + J2000_mjd,
        epoch_dec_hip_mjd=(catalog.epoch_dec_hip - 2000) * julian_year + J2000_mjd,
        epoch_ra_dr2_mjd=(catalog.epoch_ra_dr2 - 2000) * julian_year + J2000_mjd,
        epoch_dec_dr2_mjd=(catalog.epoch_dec_dr2 - 2000) * julian_year + J2000_mjd,
        epoch_ra_dr3_mjd=(catalog.epoch_ra_dr3 - 2000) * julian_year + J2000_mjd,
        epoch_dec_dr3_mjd=(catalog.epoch_dec_dr3 - 2000) * julian_year + J2000_mjd,
    )


    # Fire the DR3 top-up when the columns are absent *or* present-but-null.
    # G23H carries astrometric_chi2_al_dr3 and parallax_error as nullable
    # columns, so a source whose row simply has no value passed the old
    # `hasproperty` gate and reached the model with `missing` in place of a
    # parallax uncertainty (Proxima Centauri is the one such source in
    # G23H v1.1 / v6.1.0 — the very target this matters most for).
    needs_dr3_topup =
        !hasproperty(catalog, :astrometric_chi2_al_dr3) ||
        !hasproperty(catalog, :rv_nb_transits) ||
        !hasproperty(catalog, :parallax_error) ||
        ismissing(catalog.astrometric_chi2_al_dr3) ||
        ismissing(catalog.parallax_error)
    if needs_dr3_topup
        @warn "Column missing or null in catalog, querying Gaia DR3 TAP server (or using cached value)"

        dr3 = Octofitter._query_gaia_dr3(;gaia_id)
        catalog = (;
            catalog...,
            astrometric_chi2_al_dr3=dr3.astrometric_chi2_al,
            parallax_error=dr3.parallax_error,
            rv_nb_transits=dr3.rv_nb_transits,
            radial_velocity_error=dr3.radial_velocity_error,
        )
    end

    # Resolve the Gaia DR2 matched-transit count (used to size the DR2 epoch
    # selection). Not in the published G23H catalog;
    # pulled from the `dr2_transits_catalog` keyword when given, else the
    # mandatory `G23H_DR2Transits` sidecar DataDep, and merged into the catalog
    # row. Errors if it cannot be resolved (no fallback).
    # See `_g23h_merge_dr2_sidecar` / `_g23h_dr2_target_transits`.
    if !isnothing(gaia_id)
        catalog = _g23h_merge_dr2_sidecar(catalog, gaia_id, dr2_transits_catalog)
    end

    # Floor the uncertainty of the per-transit RV ln-σ calibration. The GP
    # calibration behind rv_ln_uncert_err_dr3 is Malmquist-biased for cool
    # dwarfs: their (colour, magnitude) calibration bins are dominated by
    # distant giants with sharper lines, so the quoted uncertainty on ln σ_rv
    # is unrealistically small and over-constrains the σ_rv_per_transit prior.
    if !isnothing(rv_ln_uncert_err_floor) &&
            hasproperty(catalog, :rv_ln_uncert_err_dr3) &&
            !ismissing(catalog.rv_ln_uncert_err_dr3) &&
            isfinite(catalog.rv_ln_uncert_err_dr3)
        catalog = merge(catalog, (;
            rv_ln_uncert_err_dr3 = max(catalog.rv_ln_uncert_err_dr3, rv_ln_uncert_err_floor)))
    end

    if isnan(catalog.hip_id)
        @warn "No Hipparcos data found; will skip HGCA and IAD modelling"
        hip_like = nothing
        hip_sol = nothing
        dist_hip = nothing
        dist_hg  = nothing
        hip_table = Table(
            @NamedTuple{iorb::Int64, epoch_yrs::Float64, parf::Float64, cosϕ::Float64, sinϕ::Float64, res::Float64, sres::Float64, reject::Bool, sres_renorm::Float64, epoch::Float64, x::Float64, y::Float64, z::Float64, vx::Float64, vy::Float64, vz::Float64, rv_kms::Float64, Δα✱::Float64, Δδ::Float64, plx_vs_time::Float64, α✱ₐ::Float64, δₐ::Float64, α✱ₘ::SubArray{Float64, 1, Matrix{Float64}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}, δₘ::SubArray{Float64, 1, Matrix{Float64}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}, scanAngle_rad::Float64, parallaxFactorAlongScan::Float64}[]
        )
        A_prepared_5_hip = fill(0.0, 0,0)
    else
        # Load the Hipparcos IAD data for epochs and scan angles
        hip_like = HipparcosIADLikelihood(;
            catalog.hip_id,
            ref_epoch_ra=catalog.epoch_ra_hip_mjd,
            ref_epoch_dec=catalog.epoch_dec_hip_mjd,
            variables=(Priors(),Derived())
        )
        A_prepared_5_hip = hip_like.A_prepared_5
        hip_table = hip_like.table
        hip_sol = hip_like.hip_sol

        # Following "Statistical properties of Hipparcos 2, caveats on its use, and a recalibration of the intermediate astrometric data"
        # by G Mirek Brandt,  Daniel Michalik,  Timothy D Brandt; 
        # we add 0.140 mas to the residuals and 2.25 mas additional dispersion to the unceratinties
        # this mitigates overfitting.
        hip_table.res .+= 0.140
        hip_table.sres_renorm .= hypot.(hip_table.sres_renorm, 2.25)

        # The HipparcosIADLikelihood constructor baked res-derived columns
        # (α✱ₐ, δₐ, α✱ₘ, δₘ, proj_meas_alongscan) BEFORE this in-place
        # recalibration, so recompute them here — otherwise the +0.140 mas
        # shift never reaches the IAD channel, which reads
        # proj_meas_alongscan (audit 2026-07-02).
        hip_table.α✱ₐ .= hip_table.res .* hip_table.cosϕ .+ hip_table.Δα✱
        hip_table.δₐ  .= hip_table.res .* hip_table.sinϕ .+ hip_table.Δδ
        for i in eachindex(hip_table.α✱ₘ)
            hip_table.α✱ₘ[i] .= [-1, 1] .* hip_table.sinϕ[i] .+ hip_table.α✱ₐ[i]
            hip_table.δₘ[i]  .= [1, -1] .* hip_table.cosϕ[i] .+ hip_table.δₐ[i]
        end
        hip_table.proj_meas_alongscan .= hip_table.res .+
            hip_table.Δα✱ .* hip_table.cosϕ .+ hip_table.Δδ .* hip_table.sinϕ

        # Precompute MvNormal distributions for correlation between ra and dec
        # Hipparcos epoch
        c = catalog.pmra_pmdec_hip[1] * catalog.pmra_hip_error[1] * catalog.pmdec_hip_error[1]
        dist_hip = MvNormal(
            @SVector([catalog.pmra_hip, catalog.pmdec_hip]),
            @SArray[
                catalog.pmra_hip_error[1]^2 c
                c catalog.pmdec_hip_error[1]^2
            ]
        )

        # Hipparcos - GAIA epoch
        c = catalog.pmra_pmdec_hg[1] * catalog.pmra_hg_error[1] * catalog.pmdec_hg_error[1]
        dist_hg = MvNormal(
            @SVector([catalog.pmra_hg, catalog.pmdec_hg]),
            @SArray [
                catalog.pmra_hg_error[1]^2 c
                c catalog.pmdec_hg_error[1]^2
            ]
        )
    end

    # Load the Gaia scanlaw etc
    # gaia_like = GaiaCatalogFitLikelihood(; gaia_id_dr3=gaia_id)

    # Besides epoch and catalog, I'm not sure we will really use this data table
    # except maybe for plotting
    # table = Table(;
    #     epoch=[hipparcos_catalog_epoch_mjd, meta_gaia_DR2.ref_epoch_mjd],
    #     catalog=[:hipparcos, :gaia],
    #     ra=[hip_like.hip_sol.radeg, gaia_like.gaia_sol.ra],
    #     dec=[hip_like.hip_sol.dedeg, gaia_like.gaia_sol.dec],
    #     plx=[hip_like.hip_sol.plx, gaia_like.gaia_sol.parallax],
    #     pmra=[hip_like.hip_sol.pm_ra, gaia_like.gaia_sol.pmra],
    #     pmdec=[hip_like.hip_sol.pm_de, gaia_like.gaia_sol.pmdec],
    # )

   

    # GAIA DR2 epoch
    c = catalog.pmra_pmdec_dr2[1] * catalog.pmra_dr2_error[1] * catalog.pmdec_dr2_error[1]
    dist_dr2 = MvNormal(
        @SVector([catalog.pmra_dr2, catalog.pmdec_dr2]),
        @SArray [
            catalog.pmra_dr2_error[1]^2 c
            c catalog.pmdec_dr2_error[1]^2
        ]
    )


    # GAIA DR3-DR2 epoch
    c = catalog.pmra_pmdec_dr32[1] * catalog.pmra_dr32_error[1] * catalog.pmdec_dr32_error[1]
    dist_dr32 = MvNormal(
        @SVector([catalog.pmra_dr32, catalog.pmdec_dr32]),
        @SArray [
            catalog.pmra_dr32_error[1]^2 c
            c catalog.pmdec_dr32_error[1]^2
        ]
    )

    # GAIA DR3 epoch
    c = catalog.pmra_pmdec_dr3[1] * catalog.pmra_dr3_error[1] * catalog.pmdec_dr3_error[1]
    dist_dr3 = MvNormal(
        @SVector([catalog.pmra_dr3, catalog.pmdec_dr3]),
        @SArray [
            catalog.pmra_dr3_error[1]^2 c
            c catalog.pmdec_dr3_error[1]^2
        ]
    )

    catalog = (; catalog..., dist_hip, dist_hg, dist_dr2, dist_dr32, dist_dr3)

    if isnothing(scanlaw_table)
        # @warn "No scan law table provided. We will fetch an approximate solution from the GOST webservice, but for best results please use the `scanninglaw` python package, installable via pip, to query the RA and Dec of this target and supply it as `scanlaw_table`. Run: `import astropy.coordinates, scanninglaw, pandas; o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg');t = scanninglaw.times.Times(version='dr3_nominal'); t.query(o,return_angles=True)`"
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GOST_forecast(catalog.ra, catalog.dec))
        forecast_table.epoch = jd2mjd.(forecast_table.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_)
        forecast_table.scanAngle_rad = forecast_table.scanAngle_rad_
    else
        @info "Scanlaw table from the `scanninglaw` python package was provided, will not use GOST."
        forecast_table = FlexTable(scanlaw_table)
        forecast_table.epoch = tcb_at_gaia_2mjd.(forecast_table.times)
        forecast_table.scanAngle_rad = deg2rad.(forecast_table.angles)

        earth_pos_vel = FlexTable(geocentre_position_query.(forecast_table.epoch))

        f = @. earth_pos_vel.x * sind(catalog.ra)-earth_pos_vel.y*cosd(catalog.ra)
        g = @. earth_pos_vel.x * cosd(catalog.ra) * sind(catalog.dec) + 
            earth_pos_vel.y * sind(catalog.ra) * sind(catalog.dec) -
            earth_pos_vel.z * cosd(catalog.dec)
        forecast_table.parallaxFactorAlongScan = @. f*sin(forecast_table.scanAngle_rad) + g*cos(forecast_table.scanAngle_rad)

    end
    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    forecast_table.cosϕ = cos.(π/2 .+ forecast_table.scanAngle_rad)
    forecast_table.sinϕ = sin.(π/2 .+ forecast_table.scanAngle_rad)

    # Get the Earth's position at those epochs
    earth_pos_vel = geocentre_position_query.(forecast_table.epoch)

    # merge the Gaia scan prediction and geocentre position results into one table
    gaia_table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)



    # Trim to the (E)DR3 AGIS input span: GOST forecasts start at the beginning
    # of science operations and can extend well past the DR3 data span (often
    # to 2018+), but the DR2 and DR3 astrometric solutions both used only data
    # from OBMT 1192.13 rev onward — excluding the first month of EPSL
    # scanning — up to their respective stops (Lindegren et al. 2018, 2021).
    # Out-of-span epochs can never contribute to any modeled channel, but left
    # in the table they consume transit_priorities selection slots:
    # chi_squared_astro is then summed over only the in-window survivors
    # while its UEVA normalization uses the catalog astrometric_matched_transits
    # count — understating the predicted companion-induced UEVA excess by ~the
    # out-of-window fraction (median ~29% over the Hipparcos sample). Single
    # stars are unaffected (zero signal either way); companion hosts need a
    # spuriously larger companion to reproduce the observed RUWE/UEVA.
    gaia_table = Table(gaia_table[gaia_agis_span_dr3.start_mjd .<= vec(gaia_table.epoch) .<= gaia_agis_span_dr3.stop_mjd, :])

    # Known data gaps, applied PER RELEASE (audit 2026-07-11): the gap lists
    # differ inside the DR2 window by ~13 days of DR2-valid time (VPU reset,
    # PAA anomalies, lunar eclipses, ...) that only the (E)DR3 processing
    # excluded, so removing the union from one shared pool undercounts the
    # epochs DR2 could use. Each release gets its own usability mask instead:
    #   DR2: its own hard (persistent=TRUE) gap list, within the DR2 span;
    #   DR3: the EDR3 gap list.
    # Rows dead for both releases are dropped. Non-persistent DR2 entries
    # (e.g. lunar eclipses) stay DR2-selectable: empirically (HIP 56379,
    # whose DR2 matched-transit count equals the raw forecast count minus
    # only the hard gaps) DR2's AGIS did use those transits.
    # Gap tables sourced from HTOF.py; authors G.M. Brandt et al.
    gaps_dr2 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiadr2_08252020.csv"), FlexTable)
    gaps_edr23 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiaedr3_12232020.csv"), FlexTable)
    _persistent(x) = x === true || (x isa AbstractString && uppercase(strip(x)) == "TRUE")
    dr2_hard = _persistent.(gaps_dr2.persistent)
    gap_starts_dr2 = obmt2mjd.(gaps_dr2.start[dr2_hard])
    gap_stops_dr2  = obmt2mjd.(gaps_dr2.end[dr2_hard])
    gap_starts_dr3 = obmt2mjd.(gaps_edr23.start)
    gap_stops_dr3  = obmt2mjd.(gaps_edr23.end)
    _in_gaps(e, starts, stops) = any(s <= e <= p for (s, p) in zip(starts, stops))
    epochs_all = vec(gaia_table.epoch)
    dr2_ok_mask = (epochs_all .<= gaia_agis_span_dr2.stop_mjd) .&
        .!_in_gaps.(epochs_all, Ref(gap_starts_dr2), Ref(gap_stops_dr2))
    dr3_ok_mask = .!_in_gaps.(epochs_all, Ref(gap_starts_dr3), Ref(gap_stops_dr3))
    keep = dr2_ok_mask .| dr3_ok_mask
    if !all(keep)
        @info "Removed forecast transits falling in data gaps of every applicable release." n_removed=count(.!keep)
    end
    n_dr2_only = count(dr2_ok_mask .& .!dr3_ok_mask)
    if n_dr2_only > 0
        @info "Keeping forecast transits usable by DR2 but excluded by the EDR3 gap list (DR2-only pool)." n=n_dr2_only
    end
    gaia_table = Table(gaia_table[keep, :])
    dr2_ok_mask = dr2_ok_mask[keep]
    dr3_ok_mask = dr3_ok_mask[keep]

    # DR2
    A_prepared_5_dr2 = prepare_A_5param(gaia_table, catalog.epoch_ra_dr2_mjd,  catalog.epoch_dec_dr2_mjd)
    
    # DR3
    A_prepared_5_dr3 = prepare_A_5param(gaia_table, catalog.epoch_ra_dr3_mjd,  catalog.epoch_dec_dr3_mjd)

    # This table serves only for plotting and keeping track of subsetting for cross-validation.
    # The actual likelihood calculations happen against the prepared linear system matrices above.
    table = Table(
        epoch=[
            isnothing(hip_like) ? NaN : mean(hip_like.table.epoch),
            isnothing(hip_like) ? NaN : years2mjd(catalog.epoch_ra_hip),
            isnothing(hip_like) ? NaN : years2mjd(catalog.epoch_dec_hip),
            isnothing(hip_like) ? NaN : years2mjd(catalog.epoch_ra_hg),
            isnothing(hip_like) ? NaN : years2mjd(catalog.epoch_dec_hg),
            years2mjd(catalog.epoch_ra_dr2),
            years2mjd(catalog.epoch_dec_dr2),
            years2mjd(catalog.epoch_ra_dr32),
            years2mjd(catalog.epoch_dec_dr32),
            years2mjd(catalog.epoch_ra_dr3),
            years2mjd(catalog.epoch_dec_dr3),
            years2mjd((catalog.epoch_dec_dr3+catalog.epoch_dec_dr2)/2),
        ],
        start_epoch=[
            isnothing(hip_like) ? 0 : minimum(hip_table.epoch),
            isnothing(hip_like) ? 0 : minimum(hip_table.epoch),
            isnothing(hip_like) ? 0 : minimum(hip_table.epoch),
            isnothing(hip_like) ? 0 : years2mjd(catalog.epoch_ra_hip),
            isnothing(hip_like) ? 0 : years2mjd(catalog.epoch_dec_hip),
            first(gaia_table.epoch[gaia_table.epoch.>=(gaia_agis_span_dr2.start_mjd)]),
            first(gaia_table.epoch[gaia_table.epoch.>=(gaia_agis_span_dr2.start_mjd)]),
            years2mjd(catalog.epoch_ra_dr2),
            years2mjd(catalog.epoch_dec_dr2),
            first(gaia_table.epoch[gaia_table.epoch.>=(gaia_agis_span_dr3.start_mjd)]),
            first(gaia_table.epoch[gaia_table.epoch.>=(gaia_agis_span_dr3.start_mjd)]),
            first(gaia_table.epoch[gaia_table.epoch.>=(gaia_agis_span_dr3.start_mjd)]),
        ],
        stop_epoch=[
            isnothing(hip_like) ? 0 : maximum(hip_table.epoch),
            isnothing(hip_like) ? 0 : maximum(hip_table.epoch),
            isnothing(hip_like) ? 0 : maximum(hip_table.epoch),
            isnothing(hip_like) ? 0 : years2mjd(catalog.epoch_ra_dr3),
            isnothing(hip_like) ? 0 : years2mjd(catalog.epoch_dec_dr3),
            last(gaia_table.epoch[gaia_table.epoch.<=(gaia_agis_span_dr2.stop_mjd)]),
            last(gaia_table.epoch[gaia_table.epoch.<=(gaia_agis_span_dr2.stop_mjd)]),
            years2mjd(catalog.epoch_ra_dr3),
            years2mjd(catalog.epoch_dec_dr3),
            last(gaia_table.epoch[gaia_table.epoch.<=(gaia_agis_span_dr3.stop_mjd)]),
            last(gaia_table.epoch[gaia_table.epoch.<=(gaia_agis_span_dr3.stop_mjd)]),
            last(gaia_table.epoch[gaia_table.epoch.<=(gaia_agis_span_dr3.stop_mjd)]),
        ],
        pm = [
            NaN,
            catalog.pmra_hip,
            catalog.pmdec_hip,
            catalog.pmra_hg,
            catalog.pmdec_hg, 
            catalog.pmra_dr2,
            catalog.pmdec_dr2,
            catalog.pmra_dr32,
            catalog.pmdec_dr32,
            catalog.pmra_dr3,
            catalog.pmdec_dr3,
            ueva_mode == :EAN  ? catalog.astrometric_excess_noise_dr3 :
            ueva_mode == :RUWE ? catalog.ruwe_dr3                     :
                                 NaN
        ],
        σ_pm = [
            NaN,
            catalog.pmra_hip_error,
            catalog.pmdec_hip_error,
            catalog.pmra_hg_error,
            catalog.pmdec_hg_error,
            catalog.pmra_dr2_error,
            catalog.pmdec_dr2_error,
            catalog.pmra_dr32_error,
            catalog.pmdec_dr32_error,
            catalog.pmra_dr3_error,
            catalog.pmdec_dr3_error,
            NaN
        ],
        kind=[
            :iad_hip,
            :ra_hip,
            :dec_hip,
            :ra_hg,
            :dec_hg,
            :ra_dr2,
            :dec_dr2,
            :ra_dr32,
            :dec_dr32,
            :ra_dr3,
            :dec_dr3,
            :ueva_dr3,
        ],

    )

    # `ueva_mode == :none`: drop the UEVA row before anything else touches the
    # table, so the row indices the rest of the constructor (and any external
    # `obs_subset` range) sees are simply the remaining channels in order.
    # Every downstream consumer selects channels by `kind ∈ table.kind`, so no
    # further bookkeeping is needed to remove UEVA from the likelihood.
    if ueva_mode == :none
        splice!(table, 12:12)
    end

    has_rv = include_rv && (
        hasproperty(catalog, :rv_ln_uncert_dr3) && !ismissing(catalog.rv_ln_uncert_dr3) && !ismissing(catalog.rv_ln_uncert_err_dr3) &&
        isfinite(catalog.rv_ln_uncert_dr3) && isfinite(catalog.rv_ln_uncert_err_dr3)
    )
    if has_rv
        push!(table.epoch, mean(gaia_table.epoch))  # or use RV-specific epoch if available
        push!(table.start_epoch, first(gaia_table.epoch))
        push!(table.stop_epoch, last(gaia_table.epoch))
        push!(table.pm, NaN)  # RV doesn't have a "pm" equivalent
        push!(table.σ_pm, NaN)
        push!(table.kind, :rv_dr3)
    end
    if isempty(hip_table)
        splice!(table, 1:5)
    end



    if isnothing(variables)

        len_epochs = length(gaia_table.epoch)   # union pool; transit_priorities spans it
        astrometric_matched_transits_dr3 = catalog.astrometric_matched_transits_dr3
        dec = catalog.dec
        ra = catalog.ra

        variables = if ueva_mode == :none
            # INERT PLACEHOLDERS.  With no UEVA channel and no deflation these
            # three cannot influence the likelihood at all; they are declared
            # only so that `σ_formal = √(σ_att² + σ_AL²)` stays finite and
            # positive for the DR3 5-parameter refit, and as constants rather
            # than sampled nuisances so the model dimension does not grow by
            # three pure prior draws.
            #
            # Why they cannot matter: `σ_formal` reaches the model solely as a
            # *scalar* weight in `fit_5param_prepared`, and a scalar weight
            # cancels out of a linear least-squares solution (see the comment
            # there) — it only rescales `chi_squared_astro`, which is consumed
            # exclusively by the now-absent UEVA channel.  Every other use
            # (μ_UEVA_single, σ_UEVA_single, the change-of-variables Jacobian)
            # lives behind the `:ueva_dr3 ∈ table.kind` mask, and the deflation
            # factor is pinned to 1.  Verified numerically: ln_like is
            # bit-identical over a 200x range in these values under :none, and
            # moves by ~1e7 over the same range under :RUWE.
            #
            # Values are the G23H population medians, so anything that reads
            # them back off a chain sees a plausible number rather than a made-
            # up one — but nothing depends on the choice.
            @variables begin
                σ_AL = 0.132
                σ_att = 0.0779
                σ_calib = 0.0795
                # G-band flux ratio (used by Gaia DR2/DR3 photocentre branch).
                fluxratio = hasproperty(sys, :fluxratio) ? sys.fluxratio : 0.0
                fluxratio_hip = hasproperty(sys, :fluxratio_hip) ? sys.fluxratio_hip : 0.0
            end
        else
            @variables begin
                σ_AL ~ truncated(Normal(catalog.sig_AL, catalog.sig_AL_sigma), lower=eps(), upper=10.0)
                σ_att ~ truncated(Normal(catalog.sig_att_radec, catalog.sig_att_radec_sigma), lower=eps(), upper=10.0)
                σ_calib ~ truncated(Normal(catalog.sig_cal, catalog.sig_cal_sigma), lower=eps(), upper=10.0)
                # G-band flux ratio (used by Gaia DR2/DR3 photocentre branch).
                fluxratio = hasproperty(sys, :fluxratio) ? sys.fluxratio : 0.0
                fluxratio_hip = hasproperty(sys, :fluxratio_hip) ? sys.fluxratio_hip : 0.0
            end
        end

        # Per-release selection pools (see the gap-mask construction above).
        # The table is already trimmed to the common AGIS start (DR2 and DR3
        # both begin at OBMT 1192.13 rev), so the DR3-only "tail" pool is
        # exactly the DR3-usable epochs after the DR2 stop.
        dr2_stop_mjd = gaia_agis_span_dr2.stop_mjd
        epochs_mjd = vec(gaia_table.epoch)
        dr2_pool = findall(dr2_ok_mask)
        dr3_win  = findall(dr3_ok_mask .& (epochs_mjd .<= dr2_stop_mjd))
        dr3_tail = findall(dr3_ok_mask .& (epochs_mjd .> dr2_stop_mjd))
        n_dr3_pool = length(dr3_win) + length(dr3_tail)
        @info "Count of missed or rejected transits:" dr3=max(0, n_dr3_pool - astrometric_matched_transits_dr3)

        # ---- DR3 epoch selection ----
        # Sample the *indices* of the GOST-forecast epochs Gaia actually used
        # via continuous `transit_priorities` (highest values win):
        # `astrometric_matched_transits_dr3` epochs, split between the DR2
        # window and the DR3-only tail in proportion to the pool sizes,
        # clamped into the hard feasibility bounds (n_tail cannot exceed the
        # tail pool, n2_win the window pool).
        degenerate_dr3 = n_dr3_pool < astrometric_matched_transits_dr3
        if degenerate_dr3
            # Every DR3-usable forecast epoch is selected; the count shortfall
            # is unmodelable from GOST.
            @warn "Fewer usable epochs in GOST forecast than `astrometric_matched_transits` reported by Gaia. Results may be inaccurate."
            n2_win = length(dr3_win)
            n_tail = length(dr3_tail)
        else
            # The clamp interval is non-empty because n_dr3_pool ≥ n_dr3 here.
            n2_win = clamp(
                round(Int, astrometric_matched_transits_dr3 * length(dr3_win) / n_dr3_pool),
                max(astrometric_matched_transits_dr3 - length(dr3_tail), 0),
                min(length(dr3_win), astrometric_matched_transits_dr3))
            n_tail = astrometric_matched_transits_dr3 - n2_win
        end

        # ---- DR2 epoch selection ----
        # Selected separately from the DR3 set (the two are NOT nested; see
        # the sidecar-plumbing comment above). `n_dr2_total` counts DR2's
        # matched observations INCLUDING doubly-downlinked bright-star
        # transits (audit 2026-07-11: DR2's count exceeds the geometric
        # transit pool by up to ~1.6× for G ≲ 6, monotonically in G, while
        # DR3's never does — DR2 tallied both windows of a doubled crossing
        # where (E)DR3's rebuilt crossmatch resolves them to one). The DR2
        # set is therefore `n_dr2_total` entries drawn as `n_dr2_distinct`
        # distinct crossings plus repeats — a repeated index enters the DR2
        # LSQ twice, exactly like a doubly-downlinked transit. For bright
        # stars (G < dr2_dup_gmag_threshold) the distinct count is a latent
        # marginalized between "every crossing doubled" (n_dr2_total/2) and
        # "no duplicates"; for fainter stars duplicates are empirically rare,
        # so only geometrically-forced repeats (count exceeding the pool)
        # occur. The DR2/DR3 catalog correlation is NOT derived from these
        # counts; it is adopted from the published `rho_dr2_dr3`.
        n_dr2_total = _g23h_dr2_target_transits(catalog, astrometric_matched_transits_dr3)
        n_dr2_hi = min(n_dr2_total, length(dr2_pool))    # max distinct crossings
        if length(dr2_pool) < n_dr2_total
            @warn "Gaia DR2 matched-transit count exceeds the geometric DR2-window pool; the excess must be duplicated (doubly-downlinked) transits and is modeled as repeated epochs." n_pool=length(dr2_pool) n_dr2_total
        end
        gmag = hasproperty(catalog, :phot_g_mean_mag_dr3) ? catalog.phot_g_mean_mag_dr3 : NaN
        bright = !ismissing(gmag) && isfinite(gmag) && gmag < dr2_dup_gmag_threshold
        n_dr2_lo = bright ? clamp(cld(n_dr2_total, 2), min(1, n_dr2_hi), n_dr2_hi) : n_dr2_hi
        marginalize_dup = n_dr2_lo < n_dr2_hi
        @info "DR2/DR3 epoch selection" n2_win n_tail n_dr2_total n_dr2_distinct_range=(n_dr2_lo, n_dr2_hi)

        # sort: downstream window logic requires chronological order;
        # partialsortperm returns priority order. Both selections read the
        # SAME priorities vector, so within the DR2 window the smaller
        # selection is automatically the top-k subset of the larger —
        # maximal epoch overlap, on the grounds that a transit usable by
        # DR2's pipeline was almost certainly reused by DR3.
        if freeze_epochs
            # Optional speed-up: draw the epoch sets once and fix them.
            transit_priorities = (randn(len_epochs)...,)
            n_dr2_distinct = rand(n_dr2_lo:max(n_dr2_lo, n_dr2_hi))
            transits = sort(vcat(
                dr3_win[ partialsortperm(SVector(transit_priorities)[dr3_win],  1:n2_win, rev=true)],
                dr3_tail[partialsortperm(SVector(transit_priorities)[dr3_tail], 1:n_tail, rev=true)]))
            transits_dr2 = _g23h_select_dr2_epochs(SVector(transit_priorities), dr2_pool, n_dr2_distinct, n_dr2_total)
            variables = vcat(variables, @variables begin
                transit_priorities = $transit_priorities
                transits = $transits
                transits_dr2 = $transits_dr2
            end)
        else
            variables = vcat(variables, @variables begin
                transit_priorities ~ MvNormal(zeros(len_epochs), I)
                transits = sort(vcat(
                    $(dr3_win)[ partialsortperm(SVector(transit_priorities)[$(dr3_win)],  1:$(n2_win), rev=true)],
                    $(dr3_tail)[partialsortperm(SVector(transit_priorities)[$(dr3_tail)], 1:$(n_tail), rev=true)]))
            end)
            if marginalize_dup
                # u_dup_dr2 maps uniformly onto the integer range of distinct
                # crossing counts [n_dr2_lo, n_dr2_hi].
                variables = vcat(variables, @variables begin
                    u_dup_dr2 ~ Uniform(0, 1)
                    transits_dr2 = _g23h_select_dr2_epochs(SVector(transit_priorities), $(dr2_pool),
                        min($(n_dr2_lo) + floor(Int, u_dup_dr2 * $(n_dr2_hi - n_dr2_lo + 1)), $(n_dr2_hi)),
                        $(n_dr2_total))
                end)
            else
                variables = vcat(variables, @variables begin
                    transits_dr2 = _g23h_select_dr2_epochs(SVector(transit_priorities), $(dr2_pool), $(n_dr2_hi), $(n_dr2_total))
                end)
            end
        end

        # The IAD nuisance parameters are only meaningful while the :iad_hip
        # likelihood row exists (likeobj_from_epoch_subset strips them again
        # if a subset later removes that row).
        if !isnothing(hip_like) && :iad_hip ∈ table.kind
            variables_iad = @variables begin
                hip_iad_jitter ~ LogUniform(0.001, 100)
                iad_Δra     ~ Uniform(-1000, 1000)
                iad_Δdec    ~ Uniform(-1000, 1000)
                iad_Δplx    ~ Uniform(-10, 10)
                iad_Δpmra   ~ Uniform(-1000, 1000) 
                iad_Δpmdec  ~ Uniform(-1000, 1000)
                iad_pmra = $(hip_sol.pm_ra) + iad_Δpmra
                iad_pmdec = $(hip_sol.pm_de) + iad_Δpmdec
            end
            variables = vcat(variables, variables_iad)
        end

        if has_rv

            # The paired GP calibration reports the per-transit RV uncertainty in
            # log space: `rv_ln_uncert_dr3` is the GP posterior mean of ln σ, and
            # `rv_ln_uncert_err_dr3` is its posterior std dev. That makes σ itself
            # LogNormal(μ_ln, σ_ln); we sample it directly to preserve the paired
            # pipeline's uncertainty faithfully.
            variables_rv = @variables begin
                σ_rv_per_transit ~ LogNormal(catalog.rv_ln_uncert_dr3, catalog.rv_ln_uncert_err_dr3)
            end
            variables = vcat(variables, variables_rv)

            n_rv = Int(catalog.rv_nb_transits)
            # Number of transits the astrometric selection actually contains
            # (shorter than the catalog count only in the degenerate case).
            n_astro_sel = min(Int(astrometric_matched_transits_dr3), n_dr3_pool)
            @info "Count of RV transits:" n_rv n_astro_sel

            # RV transits are modeled as a subset of the astrometric-used
            # transits: entirely-missed visits (e.g. gaps in Gaia coverage)
            # are assumed missed for both astrometry and RV. The astrometric
            # set is not a plain global top-k, so select the top-`n_rv`
            # priorities from *within* `transits` to keep transits_rv ⊆
            # transits.
            if 0 < n_rv < n_astro_sel
                rv_vars = @variables begin
                    transits_rv = sort(transits[partialsortperm(SVector(transit_priorities)[transits], 1:$n_rv, rev=true)])
                end
                variables = vcat(variables, rv_vars)
            elseif n_rv >= n_astro_sel > 0
                # Gaia took an RV measurement on at least every modeled
                # astrometric visit; under the RV ⊆ astrometry assumption the
                # RV epochs are exactly the astrometric epochs. The modeled
                # sample variance stays unbiased (it is count-normalized,
                # unlike the astrometric chi² sum) but samples the window
                # more coarsely than Gaia did.
                if n_rv > n_astro_sel
                    @warn "More Gaia RV transits than modeled astrometric transits; using all astrometric transits for RV." n_rv n_astro_sel
                end
                rv_vars = @variables begin
                    transits_rv = transits
                end
                variables = vcat(variables, rv_vars)
            end

        end

        @info "Added the following observation variables:"
        display(variables[1])
        display(variables[2])
    end
    (priors,derived)=variables


    obj = G23HObs{
        typeof(table),
        typeof(hip_table),
        typeof(gaia_table),
        typeof(catalog),
        typeof(hip_sol),
    }(
        table,
        priors,
        derived,
        hip_table,
        gaia_table,
        catalog,
        hip_sol,
        A_prepared_5_hip,
        A_prepared_5_dr2,
        A_prepared_5_dr3,
        include_iad,
        ueva_mode,
        Ref{Matrix{Float64}}(zeros(0, 0)),
        Ref{NTuple{4, Float64}}((0.0, 0.0, 0.0, 0.0)),
        Ref{Bool}(false),
    )
    # Eagerly initialise the cached Hipparcos weighted pseudo-inverse and
    # catalog-residual projection so no lazy mutation happens during
    # (possibly multi-threaded) sampling.
    _ensure_pinv_5_hip!(obj)
    _ensure_hip_x_const!(obj)
    return obj

end

function Octofitter.likeobj_from_epoch_subset(like::G23HObs, obs_inds)
    (;  table,
        priors,
        derived,
        hip_table,
        gaia_table,
        catalog,
        hip_sol,
        A_prepared_5_hip,
        A_prepared_5_dr2,
        A_prepared_5_dr3,
        include_iad,
        ueva_mode ) = like

    table = table[obs_inds,:]
    # NB: test the POST-subset table (`table`), not `like.table` — the whole
    # point is to react to rows the subset removed. (Previously checked the
    # original table, so the dist_hip/dist_hg nulling could never fire.)
    if  (
            :iad_hip ∉ table.kind &&
            :ra_hip ∉ table.kind &&
            :dec_hip ∉ table.kind &&
            :ra_hg ∉ table.kind &&
            :dec_hg ∉ table.kind
        )
        catalog = (;catalog..., dist_hip = nothing, dist_hg=nothing)
    end
    # When the subset removes the :iad_hip row, the six IAD nuisance
    # parameters (hip_iad_jitter, iad_Δra/Δdec/Δplx/Δpmra/Δpmdec) and their
    # two derived quantities (iad_pmra/iad_pmdec) no longer touch the
    # likelihood anywhere — they would be sampled as pure prior draws,
    # inflating the model dimension by 6 for nothing. Strip them so the
    # model dimension matches the data channels actually in use.
    if :iad_hip ∉ table.kind
        iad_keys = (:hip_iad_jitter, :iad_Δra, :iad_Δdec, :iad_Δplx,
                    :iad_Δpmra, :iad_Δpmdec, :iad_pmra, :iad_pmdec)
        priors = Priors(OrderedDict(k => v for (k, v) in priors.priors if k ∉ iad_keys))
        derived = Derived(
            OrderedDict(k => v for (k, v) in derived.variables if k ∉ iad_keys),
            derived.captured_names, derived.captured_vals)
    end
    obj = G23HObs{
        typeof(table),
        typeof(hip_table),
        typeof(gaia_table),
        typeof(catalog),
        typeof(hip_sol),
    }(
        table,
        priors,
        derived,
        hip_table,
        gaia_table,
        catalog,
        hip_sol,
        A_prepared_5_hip,
        A_prepared_5_dr2,
        A_prepared_5_dr3,
        include_iad,
        ueva_mode,
        Ref{Matrix{Float64}}(zeros(0, 0)),
        Ref{NTuple{4, Float64}}((0.0, 0.0, 0.0, 0.0)),
        Ref{Bool}(false),
    )
    _ensure_pinv_5_hip!(obj)
    _ensure_hip_x_const!(obj)
    return obj
end

function _ensure_pinv_5_hip!(like::G23HObs)
    # Hip fit uses *per-epoch* σ from like.hip_table.sres, so the LSQ is
    # weighted: x = (Aᵀ W A)⁻¹ Aᵀ W b with W = diag(1/σ_i²). Both A and σ
    # are fixed at construction, so we can cache `Q = pinv(A/σ) ./ σ'` —
    # one matvec `Q · b_orig` then reproduces the weighted-LSQ x exactly.
    if isempty(like._pinv_5_hip[]) && !isempty(like.A_prepared_5_hip)
        # sres ≤ 0 encodes scans the Hipparcos reduction REJECTED from its
        # solution; give them zero weight (σ = Inf) so this catalog-reproduction
        # LSQ matches what the catalog actually fit.  (A raw negative sres would
        # sneak the scan back in at weight 1/|sres|, and sres == 0 would produce
        # an Inf-weighted row.)  The IAD residual channel skips these rows
        # separately (audit 2026-07-02).
        σ_hip = map(s -> s > 0 ? float(s) : Inf, like.hip_table.sres)
        A = like.A_prepared_5_hip
        # A_scaled[i, j] = A[i, j] / σ[i] (row-wise scale).  Then
        # Q = pinv(A_scaled) is 5×n; weight back to x = Q·(b/σ) form by
        # dividing each column by σ once at construction.
        A_scaled = A ./ σ_hip
        Q = pinv(A_scaled) ./ permutedims(σ_hip)
        like._pinv_5_hip[] = Q
    end
    return like._pinv_5_hip[]
end

# Cache the Hipparcos catalog-residual LSQ projection x_const = Q_hip * b_const
# where b_const = (Δα·cosϕ + Δδ·sinϕ + residuals) evaluated at zero
# perturbations (Δα = Δδ = 0).  When include_iad=true, residuals = hip_table.res
# (constant catalog residuals); otherwise residuals = 0 and x_const is zero.
# Used in the all-inactive simulate! fast path.  Returned in
# fit_5param_pinv's reordered convention: (Δα_h, Δδ_h, Δpmra_h, Δpmdec_h).
function _ensure_hip_x_const!(like::G23HObs)
    if like._hip_x_const_initialised[]
        return like._hip_x_const[]
    end
    if isempty(like.A_prepared_5_hip) || isnothing(like.catalog.dist_hip)
        like._hip_x_const[] = (0.0, 0.0, 0.0, 0.0)
        like._hip_x_const_initialised[] = true
        return like._hip_x_const[]
    end
    Q_hip = _ensure_pinv_5_hip!(like)
    n_obs = size(like.hip_table, 1)
    if like.include_iad
        residuals = like.hip_table.res
        # b = 0 + 0 + residuals; x_buf = Q_hip * residuals (5-vector)
        x_buf = Q_hip * residuals
        # fit_5param_pinv reorders parameters as (x[1], x[2], x[4], x[5], x[3])
        # but only returns the first 4 in (Δα, Δδ, Δpmra, Δpmdec) form
        like._hip_x_const[] = (x_buf[1], x_buf[2], x_buf[4], x_buf[5])
    else
        like._hip_x_const[] = (0.0, 0.0, 0.0, 0.0)
    end
    like._hip_x_const_initialised[] = true
    return like._hip_x_const[]
end

function ln_like(like::G23HObs, ctx::SystemObservationContext)
    (; θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = ctx

    T = _system_number_type(θ_system)
    ll = zero(T)

    # TODO: optimize this, we only need to grab the epochs here -- it'll be faster
    if hasproperty(θ_obs, :transits)
        (;transits) = θ_obs 
        if eltype(transits) <: AbstractFloat
            transits = Int.(transits)
        end
        # The list of missed transits must be unique
        if length(unique(transits)) < length(transits)
            return nothing
        end
        ii = transits
        gaia_table = like.gaia_table[ii,:]
    else
        gaia_table = like.gaia_table
    end

    # The DR2 solution is simulated from its own epoch selection; every
    # constructor branch emits this variable, so its absence means a custom
    # `variables` set omitted it.
    hasproperty(θ_obs, :transits_dr2) || error(
        "G23HObs requires a `transits_dr2` observation variable (the DR2-used " *
        "epoch selection). It is generated automatically unless custom " *
        "`variables` are supplied — custom variable sets must define it.")

    istart_dr3 = findfirst(>=(gaia_agis_span_dr3.start_mjd), vec(gaia_table.epoch))
    iend_dr3 = findlast(<=(gaia_agis_span_dr3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(gaia_table.epoch)
    end

    @no_escape begin


        iad_resid  = @alloc(T, size(like.hip_table,1)); fill!(iad_resid, 0)
        Δα_mas_hip = @alloc(T, size(like.hip_table,1)); fill!(Δα_mas_hip, 0)
        Δδ_mas_hip = @alloc(T, size(like.hip_table,1)); fill!(Δδ_mas_hip, 0)
        # DR2 buffers match the DR2 epoch selection (see simulate!).
        n_dr2_buf = length(θ_obs.transits_dr2)
        Δα_mas_dr2 = @alloc(T, n_dr2_buf); fill!(Δα_mas_dr2, 0)
        Δδ_mas_dr2 = @alloc(T, n_dr2_buf); fill!(Δδ_mas_dr2, 0)
        Δα_mas_dr3 = @alloc(T, iend_dr3-istart_dr3+1); fill!(Δα_mas_dr3, 0)
        Δδ_mas_dr3 = @alloc(T, iend_dr3-istart_dr3+1); fill!(Δδ_mas_dr3, 0)
        # σ-inflation buffer for the BINARYS first-harmonic correction
        # (Leclerc et al. 2023, Eq. 15). Initialised to 1 (no inflation);
        # populated in simulate! per transit from the *combined* multi-companion
        # modulated signal. Consumed in two places: (a) the IAD-residual
        # likelihood when include_iad=true, scaling the per-transit residual
        # σ; (b) the catalog covariance Σ_h below, scaling it by the
        # transit-averaged f_σ² to reflect binary-induced uncertainty
        # inflation on the published Hipparcos catalog point. Must NOT be
        # used in the catalog-bias LSQ inside simulate!, because that LSQ
        # reproduces the *catalog construction* under the Hipparcos
        # pipeline's point-source σ.
        σ_inflation_hip = @alloc(T, size(like.hip_table,1)); fill!(σ_inflation_hip, 1)
        buffers = (;iad_resid, Δα_mas_hip, Δδ_mas_hip, σ_inflation_hip, Δα_mas_dr2, Δδ_mas_dr2, Δα_mas_dr3, Δδ_mas_dr3)

        sim = simulate!(buffers, like, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

        if isnothing(sim)
            ll = convert(T,-Inf)
        else

            (; μ_h, μ_hg, μ_dr2, μ_dr32, μ_dr3, UEVA_model, UEVA_unc, μ_1_3, n_dr3, n_dr2) = sim
            (;deflation_factor_dr3) = sim
            (;hip_bias_pm_sq) = sim
            # Check if we have absolute orbits
            absolute_orbits = false
            for orbit in orbits
                absolute_orbits |= orbit isa AbsoluteVisual
            end

            # Get distribution objects from catalog
            dist_hip = like.catalog.dist_hip
            dist_hg = like.catalog.dist_hg
            dist_dr2 = like.catalog.dist_dr2
            dist_dr32 = like.catalog.dist_dr32
            dist_dr3 = like.catalog.dist_dr3

            # If the user is fitting jitters, we have to regenerate these distributions 
            if !isnothing(dist_hip) && (hasproperty(θ_obs, :σ_hip_pmra) || hasproperty(θ_obs, :σ_hip_pmdec))
                c = like.catalog.pmra_pmdec_hip[1] * like.catalog.pmra_hip_error[1] * like.catalog.pmdec_hip_error[1]
                pmra_var = like.catalog.pmra_hip_error[1]^2 + (hasproperty(θ_obs, :σ_hip_pmra) ? θ_obs.σ_hip_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_hip_error[1]^2 + (hasproperty(θ_obs, :σ_hip_pmdec) ? θ_obs.σ_hip_pmdec : zero(T))^2
                dist_hip = MvNormal(
                    @SVector([like.catalog.pmra_hip, like.catalog.pmdec_hip]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end

            if !isnothing(dist_hg) && (hasproperty(θ_obs, :σ_hg_pmra) || hasproperty(θ_obs, :σ_hg_pmdec))
                c = like.catalog.pmra_pmdec_hg[1] * like.catalog.pmra_hg_error[1] * like.catalog.pmdec_hg_error[1]
                pmra_var = like.catalog.pmra_hg_error[1]^2 + (hasproperty(θ_obs, :σ_hg_pmra) ? θ_obs.σ_hg_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_hg_error[1]^2 + (hasproperty(θ_obs, :σ_hg_pmdec) ? θ_obs.σ_hg_pmdec : zero(T))^2
                dist_hg = MvNormal(
                    @SVector([like.catalog.pmra_hg, like.catalog.pmdec_hg]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end

            if !isnothing(dist_dr2) && (hasproperty(θ_obs, :σ_dr2_pmra) || hasproperty(θ_obs, :σ_dr2_pmdec))
                c = like.catalog.pmra_pmdec_dr2[1] * like.catalog.pmra_dr2_error[1] * like.catalog.pmdec_dr2_error[1]
                pmra_var = like.catalog.pmra_dr2_error[1]^2 + (hasproperty(θ_obs, :σ_dr2_pmra) ? θ_obs.σ_dr2_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_dr2_error[1]^2 + (hasproperty(θ_obs, :σ_dr2_pmdec) ? θ_obs.σ_dr2_pmdec : zero(T))^2
                dist_dr2 = MvNormal(
                    @SVector([like.catalog.pmra_dr2, like.catalog.pmdec_dr2]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end

            if !isnothing(dist_dr32) && (hasproperty(θ_obs, :σ_dr32_pmra) || hasproperty(θ_obs, :σ_dr32_pmdec))
                c = like.catalog.pmra_pmdec_dr32[1] * like.catalog.pmra_dr32_error[1] * like.catalog.pmdec_dr32_error[1]
                pmra_var = like.catalog.pmra_dr32_error[1]^2 + (hasproperty(θ_obs, :σ_dr32_pmra) ? θ_obs.σ_dr32_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_dr32_error[1]^2 + (hasproperty(θ_obs, :σ_dr32_pmdec) ? θ_obs.σ_dr32_pmdec : zero(T))^2
                dist_dr32 = MvNormal(
                    @SVector([like.catalog.pmra_dr32, like.catalog.pmdec_dr32]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end

            if !isnothing(dist_dr3) && (hasproperty(θ_obs, :σ_dr3_pmra) || hasproperty(θ_obs, :σ_dr3_pmdec))
                c = like.catalog.pmra_pmdec_dr3[1] * like.catalog.pmra_dr3_error[1] * like.catalog.pmdec_dr3_error[1]
                pmra_var = like.catalog.pmra_dr3_error[1]^2 + (hasproperty(θ_obs, :σ_dr3_pmra) ? θ_obs.σ_dr3_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_dr3_error[1]^2 + (hasproperty(θ_obs, :σ_dr3_pmdec) ? θ_obs.σ_dr3_pmdec : zero(T))^2
                dist_dr3 = MvNormal(
                    @SVector([like.catalog.pmra_dr3, like.catalog.pmdec_dr3]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end


            # Apply nonlinear correction for absolute orbits
            if absolute_orbits && !isnothing(dist_hip)
                # Add nonlinear corrections to model values
                μ_hg += @SVector [
                    like.catalog.nonlinear_dpmra,
                    like.catalog.nonlinear_dpmdec,
                ]

                # Remove HGCA's nonlinear correction from Hipparcos epoch
                # Factor of two needed since dpmra is defined to the HG epoch
                μ_h += @SVector [
                    2like.catalog.nonlinear_dpmra,
                    2like.catalog.nonlinear_dpmdec,
                ]
            end

            # Define the order of all components in our unified system
            # Order: hip_iad, ra_hip, dec_hip, ra_hg, dec_hg, ra_dr2, dec_dr2, ra_dr32, dec_dr32, ra_dr3, dec_dr3, ueva_dr3
            component_flags = @SVector [
                :ra_hip, :dec_hip,
                :ra_hg, :dec_hg,
                :ra_dr2, :dec_dr2,
                :ra_dr32, :dec_dr32,
                :ra_dr3, :dec_dr3,
                :ueva_dr3
            ]

            # Build index mask for which components are present.
            # `like.table.kind` is fixed at construction time but we lack a
            # cached form here, so we compute mask/indices into bump-allocated
            # buffers (already inside @no_escape) instead of heap-allocating
            # a Vector{Bool} + Vector{Int} per call.
            mask = @alloc(Bool, 11)
            @inbounds for i in 1:11
                mask[i] = component_flags[i] ∈ like.table.kind
            end
            n_components = 0
            @inbounds for i in 1:11
                n_components += mask[i]
            end
            indices = @alloc(Int, n_components)
            let k = 0
                @inbounds for i in 1:11
                    if mask[i]
                        k += 1
                        indices[k] = i
                    end
                end
            end

            # Change-of-variables Jacobian for the UEVA component (index 11).
            # The MvNormal below treats t = UEVA_Gaia^(1/3) = μ_1_3 as the
            # datum, but t is a parameter-dependent transform of the raw
            # catalog datum:
            #   :EAN  mode: t = (EAN² + σ_att² + σ_AL²)^(1/3),        datum EAN
            #   :RUWE mode: t = (χ²_AL·(σ_att² + σ_AL²)/(N-dof))^(1/3), datum χ²_AL
            #     (u0 cancels: (ruwe·u0)² = χ²_AL/(N-dof), and σ_formal² =
            #      σ_att² + σ_AL² is parameter-dependent — see simulate!).
            # A proper density in the raw datum needs log|dt/d(datum)|; the
            # parameter-dependent parts (data-only constants omitted) are
            #   :EAN : -(2/3)·log(UEVA_Gaia) = -2·log(μ_1_3)   [for EAN > 0]
            #   :RUWE: +(1/3)·log(σ_att² + σ_AL²)
            # Same class of bug as the rv_dr3 ξ² term (SBC 2026-07-04).
            if mask[11]
                if like.ueva_mode == :EAN
                    # EAN == 0 is a boundary atom of the Gaia reduction (the
                    # excess-noise fit pinned at zero); the continuous
                    # change-of-variables does not apply there, so those stars
                    # keep the untransformed density (flagged for a future
                    # censored-likelihood treatment).
                    if like.catalog.astrometric_excess_noise_dr3 > 0
                        ll += -2 * log(μ_1_3)
                    end
                elseif like.ueva_mode == :RUWE
                    ll += (T(1) / 3) * log(θ_obs.σ_att^2 + θ_obs.σ_AL^2)
                end
            end

            # We handle the iad separately from the big covariance matrix below, since it's not correlated
            # and we don't want to factorize a massively bigger matrix than necessary
            if :iad_hip ∈ like.table.kind
                # Hoist `hip_iad_jitter` extraction and the Normal logpdf
                # constant `-½log(2π)` out of the per-epoch loop, and inline
                # the logpdf math to skip the `Normal(0, σ_iad)` constructor.
                # Work in σ_iad² instead of σ_iad to drop the hypot's sqrt
                # and use log(σ²)·½ in place of log(σ).  Math:
                #   -½(z² + log σ²) - ½log 2π   where z² = r² / σ²
                # is algebraically identical to the original
                #   -½ z² - log σ - ½log 2π     with z = r/σ.
                (;hip_iad_jitter) = θ_obs
                _half_log2π = T(0.5) * log(T(2π))
                jitter_sq = hip_iad_jitter * hip_iad_jitter
                @inbounds for i in eachindex(iad_resid)
                    like.hip_table.reject[i] && continue
                    # Inflate per-transit residual σ by the BINARYS first-harmonic
                    # factor (Leclerc et al. 2023, Eq. 15). σ_inflation_hip[i] = 1
                    # in the unresolved limit and grows where the binary modulation
                    # reduces the signal amplitude — the IAD residual noise scales
                    # the same way.
                    s = like.hip_table.sres_renorm[i] * σ_inflation_hip[i]
                    σ_iad_sq = s * s + jitter_sq
                    r = iad_resid[i]
                    ll += -T(0.5) * (r * r / σ_iad_sq + log(σ_iad_sq)) - _half_log2π
                end
            end

            # # We handle the RV unceratinties separatley for the same reasons
            # if :rv_dr3 ∈ like.table.kind
            #     # The likelihood is based on chi-squared distribution
            #     # From the paper: ξ² follows χ²(N-1) distribution
            #     rv_chi2_dist = Chisq(sim.rv_dof)
                
            #     # Calculate log-likelihood
            #     # We want P(observing the catalog sample variance | our model)
            #     # The catalog reports error on median, we need to convert back to sample variance
            #     ε_catalog = like.catalog.radial_velocity_error
            #     N_rv = like.catalog.rv_nb_transits
            #     s_catalog_squared = (2 * N_rv / π) * (ε_catalog^2 - 0.113^2)  # Equation 4 from paper
                
            #     ξ_catalog_squared = (N_rv - 1) * s_catalog_squared / θ_obs.σ_rv_per_transit^2 # these are all in km/s

            #     @show s_catalog_squared ξ_catalog_squared rv_chi2_dist
                
            #     # Compare model to catalog
            #     ll += @show logpdf(rv_chi2_dist, ξ_catalog_squared)
            # end
            if :rv_dr3 ∈ like.table.kind
                ε_catalog = like.catalog.radial_velocity_error
                N_rv = like.catalog.rv_nb_transits
                σ_rv_per_transit = θ_obs.σ_rv_per_transit  # per-transit uncertainty in km/s

                # Convert catalog error to sample variance (Eq. A4, Chance et al. 2022)
                s_catalog_squared = (2 * N_rv / π) * (ε_catalog^2 - 0.113^2)

                # Non-centrality parameter (Eq. C2). Since sim.sample_variance is
                # computed from rv_model at N_rv RV epochs, the identity
                #   (N_rv - 1) · sample_variance = Σ_n (μ(t_n) - μ̄)²
                # reproduces λ = Σ_n ((μ_n - μ̄)/σ)² exactly.
                ncp = (N_rv - 1) * sim.sample_variance / σ_rv_per_transit^2

                # The paper's Eq. C2 states the sampling distribution for ξ² is
                # non-central χ² with N_k degrees of freedom, but the null in Eq. A6
                # uses N_k - 1. The standard sampling-theory result for
                # Σ(v_n - v̄)²/σ² with v_n ~ N(μ_n, σ²) is non-central χ²_{N-1}(λ),
                # which reduces to central χ²_{N-1} under the null. We use N-1 to
                # be self-consistent with the null.
                try
                    rv_chi2_dist = NoncentralChisq(N_rv - 1, ncp)

                    # Catalog's chi-squared statistic (Eq. A5)
                    ξ_catalog_squared = (N_rv - 1) * s_catalog_squared / σ_rv_per_transit^2

                    # Change-of-variables Jacobian: the raw datum is the catalog
                    # radial_velocity_error ε, but the density above is evaluated
                    # in ξ² = (N-1)·(2N/π)(ε² - 0.113²)/σ², a transform that
                    # depends on the parameter σ_rv_per_transit. A proper density
                    # in ε needs log|dξ²/dε| = log(2ε·(N-1)·2N/π) - 2·log(σ);
                    # the data-only constant is omitted, the parameter-dependent
                    # -2·log(σ) is required (without it the posterior is biased
                    # high by 2·sd(ln σ)² in log space; caught by SBC 2026-07-04).
                    ll += logpdf(rv_chi2_dist, ξ_catalog_squared) - 2 * log(σ_rv_per_transit)
                catch
                    ll += -Inf
                end
            end


            μ_dr3_cat, Σ_dr3 = params(dist_dr3)

            #############################################
            # Account for the UEVA-based potential uncertainty delfation of Gaia DR3 positions
            # we know the astrometric excess noise they applied and we are assuming a-priori that the
            # EAN is well-accounted for by a planet

            # DR3 position covariance at central epoch (already inflated by Gaia)
            σ_ra_dr3 = like.catalog.ra_error_central_dr3
            σ_dec_dr3 = like.catalog.dec_error_central_dr3
            ρ_radec_dr3 = like.catalog.ra_dec_corr_central_dr3

            Σ_pos_dr3 = @SMatrix [
                σ_ra_dr3^2                    ρ_radec_dr3*σ_ra_dr3*σ_dec_dr3
                ρ_radec_dr3*σ_ra_dr3*σ_dec_dr3    σ_dec_dr3^2
            ]

            # DR2 position covariance at central epoch
            σ_ra_dr2 = like.catalog.ra_error_central_dr2
            σ_dec_dr2 = like.catalog.dec_error_central_dr2
            ρ_radec_dr2 = like.catalog.ra_dec_corr_central_dr2

            Σ_pos_dr2 = @SMatrix [
                σ_ra_dr2^2                    ρ_radec_dr2*σ_ra_dr2*σ_dec_dr2
                ρ_radec_dr2*σ_ra_dr2*σ_dec_dr2    σ_dec_dr2^2
            ]
            
            ρ_dr3_dr2 = like.catalog.rho_dr2_dr3
            # ρ_dr3_dr2 = √(min(n_dr2, n_dr3) / max(n_dr2, n_dr3))
            # ρ_dr3_dr2 = θ_obs.ρ_dr3_dr2
            
            Σ_cross = @SMatrix [
                ρ_dr3_dr2*σ_ra_dr3*σ_ra_dr2                        ρ_dr3_dr2*ρ_radec_dr3*σ_ra_dr3*σ_dec_dr2
                ρ_dr3_dr2*ρ_radec_dr2*σ_dec_dr3*σ_ra_dr2          ρ_dr3_dr2*σ_dec_dr3*σ_dec_dr2
            ]
            
            # Time baselines for DR3-DR2 scaled position difference
            Δt_ra = (like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_dr2_mjd) / julian_year
            Δt_dec = (like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_dr2_mjd) / julian_year
            
            # Deflation adjustment for DR32 proper motions
            # Only the DR3-contributed terms get deflated
            # deflation_factor_dr3 = 1.0
            d = deflation_factor_dr3

            # Position covariance adjustment (in mas²)
            # ΔΣ_pos = (d^2 - 1) * Σ_pos_dr3 - 2 * (d - 1) * Σ_cross
            ΔΣ_pos = (d^2 - 1) * Σ_pos_dr3 - (d - 1) * (Σ_cross + Σ_cross')

            # Transform to proper motion covariance (mas²/yr²)
            # Different time baselines for RA and Dec
            Tr = @SMatrix [
                1/Δt_ra    0.0
                0.0        1/Δt_dec
            ]

            ΔΣ_dr32 = Tr * ΔΣ_pos * Tr'

            # Extract catalog parameters
            μ_h_cat, Σ_h = isnothing(dist_hip) ? (@SVector[0.,0.], @SMatrix zeros(2,2)) : params(dist_hip) 
            μ_hg_cat, Σ_hg = isnothing(dist_hg) ? (@SVector[0.,0.], @SMatrix zeros(2,2)) : params(dist_hg) 
            μ_dr2_cat, Σ_dr2 = params(dist_dr2)
            μ_dr32_cat, Σ_dr32 = params(dist_dr32)
            T = promote_type(
                eltype(Σ_dr32),
                eltype(Σ_h),
                eltype(Σ_hg),
                eltype(Σ_dr2),
                eltype(Σ_dr32),
                eltype(ΔΣ_dr32),
                eltype(Σ_dr3),
                typeof(deflation_factor_dr3)
            )
            Σ_h = SMatrix{2, 2, T, 4}(Σ_h)
            Σ_hg = SMatrix{2, 2, T, 4}(Σ_hg)
            Σ_dr2 = SMatrix{2, 2, T, 4}(Σ_dr2)

            # BINARYS f_σ inflation of the Hipparcos catalog covariance.  In
            # the include_iad=false path the σ_inflation_hip buffer (the
            # combined first-harmonic amplitude reduction from
            # `_simulate_skypath_hippacentre_combined!`, Leclerc 2023 Eq. 15
            # generalised) is otherwise dormant — the catalog 5-param fit at
            # the Hipparcos pipeline used point-source σ, so the LSQ that
            # reproduces the catalog bias must too, but the *uncertainty* on
            # the catalog point itself should reflect the binary-induced
            # increase in per-transit noise.  Multiply Σ_h by the
            # transit-averaged f_σ² so that the catalog-likelihood
            # covariance inflates regime-appropriately for luminous-binary
            # configurations.  In the dark-companion / single-star /
            # wide-resolved limit f_σ → 1 and this is a no-op.  Σ_hg is left
            # uninflated for v1 (the long-baseline HG covariance has a Gaia
            # endpoint as well, so the right multiplier would be < f_σ²; we
            # take the conservative course and let HG self-calibrate via
            # HGCA's renormalisation).
            if !isnothing(dist_hip) && size(like.hip_table, 1) > 0
                n_used_hip = 0
                sumsq_inflate = zero(T)
                for i in eachindex(σ_inflation_hip)
                    like.hip_table.reject[i] && continue
                    n_used_hip += 1
                    sumsq_inflate += σ_inflation_hip[i]^2
                end
                if n_used_hip > 0
                    hip_inflate_sq = sumsq_inflate / n_used_hip
                    Σ_h = Σ_h * hip_inflate_sq
                end
            end

            # BINARYS epistemic uncertainty on the catalog-bias correction.
            # The model's predicted bias (Δpmra_h, Δpmdec_h) is the BINARYS
            # photocentre modulation absorbed by the published Hipparcos
            # catalog 5-param fit; it's correct only to the extent that our
            # likelihood matches what the Hipparcos pipeline actually did.
            # Known approximations:
            #   * H1+H2 composite catalog point modelled with a single-IAD
            #     basis matrix (one reduction's scan angles, weights, and
            #     parallax factors);
            #   * Hp-band MLR systematics in the per-companion fluxratio_hip
            #     prediction, especially across the tip of the M-dwarf MS;
            #   * partial resolution gate (Gaussian taper anchored to the
            #     grating step rather than empirical detection efficiency).
            # We absorb residual model error by adding ε² · |Δpm_h|² to Σ_h
            # isotropically.  At the dark-companion / single-star limit the
            # bias is zero and this is a no-op; at high f / sep ≈ s the bias
            # grows and the catalog likelihood loosens by a regime-
            # appropriate amount.  ε_binarys is the relative trust on the
            # bias prediction (0.3 ≡ 30%, conservative first-pass value).
            ε_binarys = T(0.3)
            if !isnothing(dist_hip) && hip_bias_pm_sq > zero(T)
                Σ_h = Σ_h + (ε_binarys^2 * hip_bias_pm_sq) * SMatrix{2, 2, T, 4}(I)
            end


            # Apply deflation adjustment to DR32 covariance
            Σ_dr32 = SMatrix{2, 2, T, 4}(Σ_dr32 .+ ΔΣ_dr32)
            # Σ_dr32 = SMatrix{2, 2, T, 4}(Σ_dr32)
            # Σ_dr32 =SMatrix{2, 2, T, 4}( [
            #     Σ_dr32[1,1] 0
            #     0           Σ_dr32[1,1]
            # ] .+ ΔΣ_dr32)
            # Apply deflation adjustment to DR3 proper motions
            Σ_dr3 = SMatrix{2, 2, T, 4}(Σ_dr3 .* deflation_factor_dr3^2)

            

            μ_catalog_full = @SVector [
                μ_h_cat[1], μ_h_cat[2],     # ra_hip, dec_hip
                μ_hg_cat[1], μ_hg_cat[2],   # ra_hg, dec_hg
                μ_dr2_cat[1], μ_dr2_cat[2], # ra_dr2, dec_dr2
                μ_dr32_cat[1], μ_dr32_cat[2], # ra_dr32, dec_dr32
                μ_dr3_cat[1], μ_dr3_cat[2], # ra_dr3, dec_dr3
                μ_1_3                       # ueva_dr3 catalog value
            ]

            μ_model_full = @SVector [
                μ_h[1], μ_h[2],           # ra_hip, dec_hip
                μ_hg[1], μ_hg[2],         # ra_hg, dec_hg
                μ_dr2[1], μ_dr2[2],       # ra_dr2, dec_dr2
                μ_dr32[1], μ_dr32[2],     # ra_dr32, dec_dr32
                μ_dr3[1], μ_dr3[2],       # ra_dr3, dec_dr3
                UEVA_model                # ueva_dr3 model value
            ]

            # Build the full covariance matrix into a Bumper-allocated buffer
            # (we're already inside @no_escape — no heap allocation).
            Σ_full = @alloc(T, 11, 11); fill!(Σ_full, zero(T))
            Σ_full[1:2, 1:2] .= Σ_h     # Hipparcos
            Σ_full[3:4, 3:4] .= Σ_hg    # HGCA
            Σ_full[5:6, 5:6] .= Σ_dr2   # DR2
            Σ_full[7:8, 7:8] .= Σ_dr32  # DR3-DR2
            Σ_full[9:10, 9:10] .= Σ_dr3 # DR3
            Σ_full[11, 11] = UEVA_unc^2
            # Cross-epoch correlations between DR2 and DR3
            K = ρ_dr3_dr2 * sqrt(Σ_dr2) * sqrt(Σ_dr3)'
            Σ_full[5:6, 9:10] .= K
            Σ_full[9:10, 5:6] .= K'

            # Extract selected components into Bumper-allocated scratch.
            μ_catalog_selected = @alloc(T, n_components)
            μ_model_selected = @alloc(T, n_components)
            @inbounds for k in 1:n_components
                μ_catalog_selected[k] = μ_catalog_full[indices[k]]
                μ_model_selected[k] = μ_model_full[indices[k]]
            end
            Σ_selected = @alloc(T, n_components, n_components)
            @inbounds for kj in 1:n_components, ki in 1:n_components
                Σ_selected[ki, kj] = Σ_full[indices[ki], indices[kj]]
            end

            # Diagnostic capture (no-op unless _G23H_DEBUG_PULLS holds a Vector).
            if _G23H_DEBUG_PULLS[] !== nothing
                lbls = [component_flags[indices[k]] for k in 1:n_components]
                cat  = [Float64(μ_catalog_selected[k]) for k in 1:n_components]
                mod  = [Float64(μ_model_selected[k]) for k in 1:n_components]
                sig  = [sqrt(Float64(Σ_selected[k,k])) for k in 1:n_components]
                Σcap = [Float64(Σ_selected[ki,kj]) for ki in 1:n_components, kj in 1:n_components]
                push!(_G23H_DEBUG_PULLS[],
                      (; labels=lbls, catalog=cat, model=mod, sigma=sig,
                         pull=(cat .- mod) ./ sig, Σ=Σcap))
            end

            # Compute likelihood
            if n_components == 1
                ll += logpdf(Normal(μ_catalog_selected[1], sqrt(Σ_selected[1,1])), μ_model_selected[1])
            else
                # Manual MvNormal logpdf via in-place cholesky on a Bumper
                # buffer — skips the heap-allocating PDMat wrapper inside
                # `Distributions.MvNormal`, which re-factorises Σ on every
                # call. Math: log p = -½(δ'Σ⁻¹δ + n·log(2π) + log|Σ|),
                # with δ'Σ⁻¹δ = ||L⁻¹δ||² and log|Σ| = 2·Σ log(L_ii).
                L_buf = @alloc(T, n_components, n_components)
                @inbounds for kj in 1:n_components, ki in 1:n_components
                    L_buf[ki, kj] = Σ_selected[ki, kj]
                end
                δ_buf = @alloc(T, n_components)
                @inbounds for k in 1:n_components
                    δ_buf[k] = μ_model_selected[k] - μ_catalog_selected[k]
                end
                try
                    chol_F = cholesky!(Hermitian(L_buf, :L))
                    ldiv!(chol_F.L, δ_buf)
                    quad = zero(T)
                    @inbounds for k in 1:n_components
                        quad += δ_buf[k] * δ_buf[k]
                    end
                    logdet_Σ = zero(T)
                    @inbounds for k in 1:n_components
                        logdet_Σ += log(L_buf[k, k])
                    end
                    logdet_Σ *= 2
                    ll += -T(0.5) * (quad + n_components * log(T(2π)) + logdet_Σ)
                catch err
                    if err isa PosDefException
                        ll = convert(T, -Inf)
                    else
                        rethrow(err)
                    end
                end
            end
        end
    end

    if isnan(ll)
        return convert(T, -Inf)
    end

    return ll
end

function simulate(like::G23HObs, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    # TODO: optimize this, we only need to grab the epochs here -- it'll be faster
    if hasproperty(θ_obs, :transits)
        (;transits) = θ_obs 
        if eltype(transits) <: AbstractFloat
            transits = Int.(transits)
        end
        # The list of missed transits must be unique
        if length(unique(transits)) < length(transits)
            return nothing
        end
        ii = transits
        gaia_table = like.gaia_table[ii,:]
    else
        gaia_table = like.gaia_table
    end

    # The DR2 solution is simulated from its own epoch selection; every
    # constructor branch emits this variable, so its absence means a custom
    # `variables` set omitted it.
    hasproperty(θ_obs, :transits_dr2) || error(
        "G23HObs requires a `transits_dr2` observation variable (the DR2-used " *
        "epoch selection). It is generated automatically unless custom " *
        "`variables` are supplied — custom variable sets must define it.")

    istart_dr3 = findfirst(>=(gaia_agis_span_dr3.start_mjd), vec(gaia_table.epoch))
    iend_dr3 = findlast(<=(gaia_agis_span_dr3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(gaia_table.epoch)
    end

    iad_resid  = zeros(size(like.hip_table,1)); fill!(iad_resid, 0)
    Δα_mas_hip = zeros(size(like.hip_table,1)); fill!(Δα_mas_hip, 0)
    Δδ_mas_hip = zeros(size(like.hip_table,1)); fill!(Δδ_mas_hip, 0)
    # DR2 buffers match the DR2 epoch selection (see simulate!).
    n_dr2_buf = length(θ_obs.transits_dr2)
    Δα_mas_dr2 = zeros(n_dr2_buf); fill!(Δα_mas_dr2, 0)
    Δδ_mas_dr2 = zeros(n_dr2_buf); fill!(Δδ_mas_dr2, 0)
    Δα_mas_dr3 = zeros(iend_dr3-istart_dr3+1); fill!(Δα_mas_dr3, 0)
    Δδ_mas_dr3 = zeros(iend_dr3-istart_dr3+1); fill!(Δδ_mas_dr3, 0)
    σ_inflation_hip = ones(size(like.hip_table,1))

    buffers = (;iad_resid, Δα_mas_hip, Δδ_mas_hip, σ_inflation_hip, Δα_mas_dr2, Δδ_mas_dr2, Δα_mas_dr3, Δδ_mas_dr3)

    out = simulate!(buffers, like, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    return out 
end

function simulate!(buffers, like::G23HObs, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    (;Δα_mas_hip, Δδ_mas_hip, σ_inflation_hip, Δα_mas_dr2, Δδ_mas_dr2, Δα_mas_dr3, Δδ_mas_dr3, iad_resid, ) = buffers


    T = _system_number_type(θ_system)

    # Generate simulated observations from this sample draw
    # Get Gaia noise parameters from observation variables
    (;σ_att, σ_AL, σ_calib,) = θ_obs
    σ_formal = sqrt(σ_att^2 + σ_AL^2)

    gaia_n_dof = like.catalog.astrometric_params_solved_dr3 == 31 ? 5 : 6

    # The gaia_table and A_prepared_5_dr3 include all available visibility
    # windows, not filtered to specifically be DR2 or DR3.
    # Here we may further reject some more to marginalize over
    # unknown missed/rejected transits.
    if hasproperty(θ_obs, :transits)
        (;transits) = θ_obs
        if eltype(transits) <: AbstractFloat
            transits = Int.(transits)
        end
        # The list of missed transits must be unique
        if length(unique(transits)) < length(transits)
            return nothing
        end
        ii = transits
        gaia_table = like.gaia_table[ii,:]
        A_prepared_5_dr3 = view(like.A_prepared_5_dr3, ii,:)
    else
        gaia_table = like.gaia_table
        A_prepared_5_dr3 = like.A_prepared_5_dr3
    end

    # DR2's own epoch selection (see the construction-time selection block):
    # the DR2-used transit set is not nested in the DR3-used set, so the DR2
    # solution is simulated from exactly these `like.gaia_table` rows. Always
    # present: every constructor branch emits it, and custom `variables` must
    # define it (checked with a friendly error in ln_like / simulate).
    # REPEATED indices are legitimate here (unlike `transits`): a repeat is a
    # doubly-downlinked bright-star transit entering the DR2 fit twice.
    (;transits_dr2) = θ_obs
    if eltype(transits_dr2) <: AbstractFloat
        transits_dr2 = Int.(transits_dr2)
    end
    kk_dr2 = collect(Int, transits_dr2)

    if hasproperty(θ_obs, :transits_rv)
        (;transits_rv) = θ_obs 
        if eltype(transits_rv) <: AbstractFloat
            transits_rv = Int.(transits_rv)
        end
        # The list of RV transit indices must be unique
        if length(unique(transits_rv)) < length(transits_rv)
            return nothing
        end
        # `transits_rv` holds the indices of the RV-used transits (a subset of the
        # astrometric-used set), so we select them directly.
        jj = collect(transits_rv)
        gaia_table_rv = like.gaia_table[jj,:]
    else
        gaia_table_rv = like.gaia_table
    end

    # Now we fit a no-planet (zero mass planet) sky path model to this data.
    # These should be fit using the appropriate catalog reference epoch so 
    # that they can be compared correctly.

    absolute_orbits = false
    for orbit in orbits
        absolute_orbits |= orbit isa AbsoluteVisual
        # TODO: could check in a more user-friendly way
        # that we don't have a mismatch of different orbit types
        # for different planets?
    end

    # Fast path: all companions inactive (mass == 0).  Triggered by the
    # n_planets prior's P(N=0)=0.5 + the active-flag pp_*=0 multiplier in the
    # @variables block, so ~50% of prior draws (and most of the chain volume
    # for sparse-companion stars) skip the full simulate! body.  In the
    # all-inactive limit:
    #   * every per-transit perturbation Δα_mas_*/Δδ_mas_* is exactly 0;
    #   * the DR3/DR2 5-param fits on zero data return parameters = 0 and
    #     chi_squared_astro = 0, so UEVA_model_raw = 0;
    #   * the Hippacentre combined loop early-exits with σ_inflation = 1;
    #   * the Hipparcos 5-param fit collapses to the catalog-residual
    #     projection x_const = Q_hip * residuals, cached at construction
    #     since it depends only on hip_table.res (constant);
    #   * RV simulation produces an all-zero rv_model ⇒ sample_variance = 0.
    # We restrict the fast path to non-AbsoluteVisual orbits — AbsoluteVisual
    # needs the per-epoch propagate_astrom (barycentric motion, differential
    # light-travel time), which doesn't simplify under all-inactive.
    all_inactive = !absolute_orbits
    if all_inactive
        @inbounds for i in eachindex(orbits)
            if θ_system.planets[i].mass != zero(θ_system.planets[i].mass)
                all_inactive = false
                break
            end
        end
    end

    if all_inactive
        # Hipparcos catalog-residual cached LSQ projection.
        x_const = _ensure_hip_x_const!(like)
        Δα_h     = T(x_const[1])
        Δδ_h     = T(x_const[2])
        Δpmra_h  = T(x_const[3])
        Δpmdec_h = T(x_const[4])

        # IAD residual loop — identical to the slow path's loop but reading
        # Δα_mas_hip = Δδ_mas_hip ≡ 0 and σ_inflation_hip ≡ 1.
        if !isnothing(like.catalog.dist_hip)
            # Only when the :iad_hip likelihood row is active: it is the sole
            # consumer of iad_resid, and when it is excluded (epoch subset)
            # the iad_* nuisance parameters are stripped from θ_obs entirely.
            if :iad_hip ∈ like.table.kind
            (;iad_Δra, iad_Δdec, iad_pmra, iad_pmdec, iad_Δplx) = θ_obs
            plx_at_epoch = like.hip_sol.plx + iad_Δplx
            inv_julian_year = inv(julian_year)
            iad_pmra_eff  = iad_pmra  - Δpmra_h
            iad_pmdec_eff = iad_pmdec - Δpmdec_h
            iad_Δra_eff   = iad_Δra   - Δα_h
            iad_Δdec_eff  = iad_Δdec  - Δδ_h
            hip_epoch    = like.hip_table.epoch
            hip_cosϕ     = like.hip_table.cosϕ
            hip_sinϕ     = like.hip_table.sinϕ
            hip_plxFact  = like.hip_table.parallaxFactorAlongScan
            hip_projMeas = like.hip_table.proj_meas_alongscan
            @inbounds @simd for i_epoch in eachindex(hip_epoch, iad_resid)
                cosϕ = hip_cosϕ[i_epoch]
                sinϕ = hip_sinϕ[i_epoch]
                Δt = (hip_epoch[i_epoch] - hipparcos_catalog_epoch_mjd) * inv_julian_year
                α_off = iad_Δra_eff  + Δt * iad_pmra_eff
                δ_off = iad_Δdec_eff + Δt * iad_pmdec_eff
                proj_model = α_off * cosϕ + δ_off * sinϕ + plx_at_epoch * hip_plxFact[i_epoch]
                iad_resid[i_epoch] = abs(hip_projMeas[i_epoch] - proj_model)
            end
            end
            μ_h_fast = @SVector [θ_system.pmra + Δpmra_h, θ_system.pmdec + Δpmdec_h]
            hip_bias_pm_sq = Δpmra_h*Δpmra_h + Δpmdec_h*Δpmdec_h
        else
            μ_h_fast = @SVector [zero(T), zero(T)]
            hip_bias_pm_sq = zero(T)
        end

        pmra_sys = θ_system.pmra
        pmdec_sys = θ_system.pmdec
        μ_zero = @SVector [pmra_sys, pmdec_sys]
        μ_dr3_fast  = μ_zero
        μ_dr2_fast  = μ_zero
        # Mirror the slow path's HG model PM exactly (non-absolute branch):
        # with all perturbations zero, Δα_dr3 = Δδ_dr3 = 0, but the Hipparcos
        # catalog-residual projection (Δα_h, Δδ_h) = x_const[1:2] can be
        # nonzero, and it enters μ_hg through the long HG baseline.  Keep the
        # exact operation order of the slow path so the two are bit-identical.
        μ_hg_fast = if isnothing(like.catalog.dist_hip)
            @SVector [zero(T), zero(T)]
        else
            @SVector [
                (zero(T) - Δα_h) / (
                    like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_hip_mjd
                )*julian_year + pmra_sys,
                (zero(T) - Δδ_h) / (
                    like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_hip_mjd
                )*julian_year + pmdec_sys,
            ]
        end
        μ_dr32_fast = μ_zero

        istart_dr3 = findfirst(>=(gaia_agis_span_dr3.start_mjd), vec(gaia_table.epoch))
        iend_dr3 = findlast(<=(gaia_agis_span_dr3.stop_mjd), vec(gaia_table.epoch))
        if isnothing(istart_dr3); istart_dr3 = 1; end
        if isnothing(iend_dr3); iend_dr3 = length(gaia_table.epoch); end

        (;astrometric_chi2_al_dr3, astrometric_n_good_obs_al_dr3,
           astrometric_matched_transits_dr3) = like.catalog
        N = astrometric_n_good_obs_al_dr3
        N_FoV = astrometric_matched_transits_dr3
        N_AL = N / N_FoV
        if like.ueva_mode == :EAN
            UEVA_Gaia = like.catalog.astrometric_excess_noise_dr3^2 + σ_att^2 + σ_AL^2
        elseif like.ueva_mode == :RUWE
            ruwe_dr3 = like.catalog.ruwe_dr3
            u0 = 1/ruwe_dr3 * sqrt(astrometric_chi2_al_dr3/(N - gaia_n_dof))
            UEVA_Gaia = (ruwe_dr3 * u0)^2 * σ_formal^2
        elseif like.ueva_mode == :none
            # Placeholder only: the UEVA channel is absent from table.kind, so
            # μ_1_3 is masked out of the likelihood and the deflation factor is
            # pinned to 1 below.  Kept finite and positive so nothing downstream
            # sees a NaN.
            UEVA_Gaia = σ_formal^2
        else
            error("Unsupported mode (should be :EAN, :RUWE or :none, was $(like.ueva_mode))")
        end
        μ_UEVA_single = (N_AL / (N - gaia_n_dof)) *
                    ((N_FoV - gaia_n_dof) * σ_calib^2 + N_FoV * σ_AL^2)
        σ_UEVA_single = sqrt(
            2 * N_AL / (N - gaia_n_dof)^2 * (
                N_AL * (N_FoV - gaia_n_dof) * σ_calib^4 +
                N_FoV * σ_AL^4 +
                2 * N_FoV * σ_AL^2 * σ_calib^2
            )
        )
        μ_1_3 = UEVA_Gaia^(1/3)
        UEVA_unc = σ_UEVA_single * μ_UEVA_single^(-2/3) / 3
        # chi_squared_astro = 0 ⇒ UEVA_model_1 = 0 ⇒ UEVA_model = cbrt(μ_UEVA_single)
        UEVA_model = cbrt(μ_UEVA_single)
        deflation_factor_raw = sqrt(μ_UEVA_single / UEVA_Gaia)
        # With ueva_mode == :none there is no calibrated σ_AL/σ_att/σ_calib to
        # compute μ_UEVA_single from, so the deflation is not identified — take
        # Gaia's published DR3/DR32 uncertainties at face value.
        deflation_factor_dr3 = like.ueva_mode == :none ? one(deflation_factor_raw) :
            (deflation_factor_raw > 1.0 ? 1.0 : deflation_factor_raw)

        if :rv_dr3 ∈ like.table.kind
            N_rv = like.catalog.rv_nb_transits
            rv_dof_fast = N_rv - 1
            ε_catalog = like.catalog.radial_velocity_error
            s_catalog_squared_fast = (2 * N_rv / π) * (ε_catalog^2 - 0.113^2)
            rv_mean_fast = zero(T)
            sample_variance_fast = zero(T)
        else
            rv_dof_fast = convert(T, NaN)
            s_catalog_squared_fast = convert(T, NaN)
            rv_mean_fast = convert(T, NaN)
            sample_variance_fast = convert(T, NaN)
        end

        return (;
            UEVA_model,
            UEVA_unc,
            μ_1_3,
            μ_h = μ_h_fast,
            μ_hg = μ_hg_fast,
            μ_dr2 = μ_dr2_fast,
            μ_dr32 = μ_dr32_fast,
            μ_dr3 = μ_dr3_fast,
            μ = (@SVector [μ_h_fast[1],μ_h_fast[2],μ_hg_fast[1],μ_hg_fast[2],μ_dr2_fast[1],μ_dr2_fast[2],μ_dr32_fast[1],μ_dr32_fast[2],μ_dr3_fast[1],μ_dr3_fast[2],UEVA_model,sample_variance_fast]),
            hip_bias_pm_sq,
            n_dr3 = iend_dr3 - istart_dr3 + 1,
            n_dr2 = length(kk_dr2),
            rv_dof = rv_dof_fast,
            rv_mean = rv_mean_fast,
            sample_variance = sample_variance_fast,
            s_catalog_squared = s_catalog_squared_fast,
            deflation_factor_dr3,
            pmra_hip_model = μ_h_fast[1], pmdec_hip_model = μ_h_fast[2],
            pmra_hg_model = μ_hg_fast[1], pmdec_hg_model = μ_hg_fast[2],
            pmra_dr2_model = μ_dr2_fast[1], pmdec_dr2_model = μ_dr2_fast[2],
            pmra_dr32_model = μ_dr32_fast[1], pmdec_dr32_model = μ_dr32_fast[2],
            pmra_dr3_model = μ_dr3_fast[1], pmdec_dr3_model = μ_dr3_fast[2],
            Δα_dr3 = zero(T), Δδ_dr3 = zero(T),
            Δpmra_dr3 = zero(T), Δpmdec_dr3 = zero(T),
        )
    end


    # Helper functions to either get the static pmra from the orbital elements,
    # or, if using an AbsoluteVisualOrbit, get the propagated pmra at the
    # current epoch accounting for barycentric motion.
    function propagate_astrom(orbits::NTuple{N,<:PlanetOrbits.AbsoluteVisualOrbit} where N, epoch_ra, epoch_dec)
        o = first(orbits)
        sol_ra = orbitsolve(o, epoch_ra)
        cmp_ra = sol_ra.compensated
        sol_dec = orbitsolve(o, epoch_dec)
        cmp_dec = sol_dec.compensated
        # Account for the instantaneous differential light travel time apparent acceleration.
        # Treat as linear for the duration of Gaia or Hipparcos
        t1 = max(epoch_ra, epoch_dec)
        Δt = 100
        t2 = t1 + Δt
        sol = epoch_ra >= epoch_dec ? sol_ra : sol_dec
        sol′ = orbitsolve(o,t2)
        # This isn't right! This is double counting the proper motion which already goes into ra/dec
        # Take change in delta_time and multiply it by pmra/pmdec
        diff_lt_app_pmra = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmra2
        diff_lt_app_pmdec = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmdec2
        return cmp_ra.ra2, cmp_dec.dec2, cmp_ra.pmra2+diff_lt_app_pmra, cmp_dec.pmdec2+diff_lt_app_pmdec
        # return (
        #     cmp_ra.ra2 - Δα_dr3/60/60/1000/cosd(cmp_dec.dec2),
        #     cmp_dec.dec2 - Δδ_dr3/60/60/1000,
        #     cmp_ra.pmra2+diff_lt_app_pmra - Δpmra_dr3,
        #     cmp_dec.pmdec2+diff_lt_app_pmdec - Δpmdec_dr3
        # )
    end
    function propagate_astrom(orbits::Tuple{}, _, _)
        return 0.0, 0.0, θ_system.pmra, θ_system.pmdec
    end
    function propagate_astrom(orbits::Any, _, _)
        return 0.0, 0.0, θ_system.pmra, θ_system.pmdec
    end




    # Pre-compute orbit solutions per active companion at every hip_table +
    # gaia_table epoch.  Without this cache, simulate! would call orbitsolve
    # three times per active planet per Hipparcos transit (Hippacentre combined,
    # ~100 transits) and twice per planet per Gaia transit (DR3 and DR2 paths
    # iterate over overlapping epoch sets, so they double-count the shared
    # DR2/DR3 region).  Layout per planet: indices 1..n_hip = hip epochs,
    # n_hip+1..n_hip+n_gaia = (DR3-selected) gaia epochs, and — only when the
    # model carries a separate DR2 selection — n_hip+n_gaia+1..end = the
    # DR2-selected epochs (which need not be rows of the DR3-selected
    # gaia_table).  Bumper buffer is the caller's @no_escape (ln_like).
    # Allocated for every planet (3 in production) but only populated for
    # active ones; downstream loops skip mass==0 planets without indexing
    # into the array.
    n_planets = length(orbits)
    n_hip_cache = size(like.hip_table, 1)
    n_gaia_cache = length(gaia_table.epoch)
    n_dr2_cache = length(kk_dr2)
    n_total_cache = n_hip_cache + n_gaia_cache + n_dr2_cache
    _ref_epoch_cache = n_hip_cache > 0 ? like.hip_table.epoch[1] : gaia_table.epoch[1]
    # NOTE: all companions are assumed to share a single orbit-solution type
    # (true whenever the planets use the same orbit basis, as in the lumcomp
    # pipeline). A heterogeneous mix would fail loudly on the setindex!
    # below, not silently corrupt.
    _sol0_cache = orbitsolve(first(orbits), _ref_epoch_cache)
    _SolType = typeof(_sol0_cache)
    _bumper = Bumper.default_buffer()
    planet_sols_cache = ntuple(n_planets) do p
        Bumper.alloc!(_bumper, _SolType, n_total_cache)
    end
    @inbounds for p in 1:n_planets
        if θ_system.planets[p].mass != zero(θ_system.planets[p].mass)
            sols = planet_sols_cache[p]
            o = orbits[p]
            for i in 1:n_hip_cache
                sols[i] = orbitsolve(o, like.hip_table.epoch[i])
            end
            for i in 1:n_gaia_cache
                sols[n_hip_cache + i] = orbitsolve(o, gaia_table.epoch[i])
            end
            for i in 1:n_dr2_cache
                sols[n_hip_cache + n_gaia_cache + i] = orbitsolve(o, like.gaia_table.epoch[kk_dr2[i]])
            end
        end
    end

    ################################
    # DR3
    istart_dr3 = findfirst(>=(gaia_agis_span_dr3.start_mjd), vec(gaia_table.epoch))
    iend_dr3 = findlast(<=(gaia_agis_span_dr3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(gaia_table.epoch)
    end
    gaia_table_dr3 = @views gaia_table[istart_dr3:iend_dr3]
    # `_simulate_skypath_perturbations!` indexes `orbit_solutions[i_start + i]`
    # with i in 1:length(table); offset so index lands in the gaia portion at
    # the istart_dr3 row.
    _dr3_sol_start = n_hip_cache + istart_dr3 - 1
    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        if planet_mass_msol == 0.0
            continue
        end
        if hasproperty(θ_obs, :fluxratio)
            if θ_obs.fluxratio isa Number
                fluxratio = θ_obs.fluxratio
            else
                fluxratio = θ_obs.fluxratio[i_planet]
            end
        else
            fluxratio = 0.0
        end
        _simulate_skypath_perturbations!(
            Δα_mas_dr3, Δδ_mas_dr3,
            gaia_table_dr3, orbit,
            planet_mass_msol, fluxratio,
            planet_sols_cache[i_planet],
            _dr3_sol_start, T,
        )
    end

    out_dr3 = fit_5param_prepared(
        view(A_prepared_5_dr3, istart_dr3:iend_dr3,:),
        view(gaia_table, istart_dr3:iend_dr3),
        Δα_mas_dr3, Δδ_mas_dr3, 0.0, σ_formal;
        include_chi2=Val(true)
    )
    Δα_dr3, Δδ_dr3, Δpmra_dr3, Δpmdec_dr3 = out_dr3.parameters
    # Rigorously propagate the linear proper motion component in spherical coordinates
    # Account for within-gaia differential light travel time 
    α_dr3₀, δ_dr3₀, pmra_dr3₀, pmdec_dr3₀ = propagate_astrom(orbits, like.catalog.epoch_ra_dr3_mjd, like.catalog.epoch_dec_dr3_mjd)
    μ_dr3 = @SVector [pmra_dr3₀ + Δpmra_dr3, pmdec_dr3₀ + Δpmdec_dr3]

    # Note: we shift the entire reference frame so that the proper motion is defined on the primary star
    # all proper motions derived below are shifted the perturbation in DR3 
    # This vastly improves sampling efficiency.
    # Leave Δpmdec_dr3 - Δpmdec_dr3 above as an explicit reminder about this ^

    ################################
    # DR2
    # Simulate the DR2 solution from exactly the `transits_dr2` forecast
    # epochs (cached in the third planet_sols_cache segment; they need not be
    # rows of the DR3-selected gaia_table).
    gaia_table_dr2 = like.gaia_table[kk_dr2,:]
    A_5_dr2_sel = view(like.A_prepared_5_dr2, kk_dr2,:)
    _dr2_sol_start = n_hip_cache + n_gaia_cache
    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        if planet_mass_msol == 0.0
            continue
        end
        if hasproperty(θ_obs, :fluxratio)
            if θ_obs.fluxratio isa Number
                fluxratio = θ_obs.fluxratio
            else
                fluxratio = θ_obs.fluxratio[i_planet]
            end
        else
            fluxratio = 0.0
        end
        _simulate_skypath_perturbations!(
            Δα_mas_dr2, Δδ_mas_dr2,
            gaia_table_dr2, orbit,
            planet_mass_msol, fluxratio,
            planet_sols_cache[i_planet],
            _dr2_sol_start, T
        )
    end

    out = fit_5param_prepared(A_5_dr2_sel, gaia_table_dr2, Δα_mas_dr2, Δδ_mas_dr2)
    # out = fit_4param_prepared(hgca_like.gaialike.A_prepared_4, gaia_table, Δα_mas_dr2, Δδ_mas_dr2)
    Δα_dr2, Δδ_dr2, Δpmra_dr2, Δpmdec_dr2 = out.parameters
    # Rigorously propagate the linear proper motion component in spherical coordinates
    # Account for within-gaia differential light travel time 
    α_dr2₀, δ_dr2₀, pmra_dr2₀, pmdec_dr2₀ = propagate_astrom(orbits, like.catalog.epoch_ra_dr2_mjd, like.catalog.epoch_dec_dr2_mjd)
    μ_dr2 = @SVector [pmra_dr2₀ + Δpmra_dr2, pmdec_dr2₀ + Δpmdec_dr2]

        

    ################################
    # Hipparcos
    # Track |bias|² in PM space for the catalog-likelihood epistemic inflation
    # downstream — see the ε_binarys block in ln_like.  Initialised here so the
    # value flows through both branches of the if/else below into the return
    # named tuple.
    hip_bias_pm_sq = zero(T)
    if isnothing(like.catalog.dist_hip)
        # type stable since dist_hip is part of the likelihood type parameters
        # ie. we statically know which of these branches will be taken.
        μ_h = @SVector [zero(T), zero(T)]
    else


        # σ_inflation_hip comes from the caller's buffer; reset it to 1 at the
        # start of this evaluation. It will be populated *jointly* across all
        # companions in a single call to `_simulate_skypath_hippacentre_combined!`
        # below — multiplicative per-companion accumulation is wrong because the
        # multi-source modulated signal is not the sum of per-companion signals.
        n_hip = length(like.hip_table.epoch)
        fill!(σ_inflation_hip, one(T))

        planet_masses_msol_hip = ntuple(i_planet -> θ_system.planets[i_planet].mass * Octofitter.mjup2msol, n_planets)
        # Companions with mass = 0 contribute nothing — zero out their flux ratio so
        # they fall out of both the modulated-signal sum and the host-reflex sum.
        flux_ratios_hip = ntuple(n_planets) do i_planet
            planet_masses_msol_hip[i_planet] == 0.0 ? zero(T) :
                (θ_obs.fluxratio_hip isa Number ? θ_obs.fluxratio_hip : θ_obs.fluxratio_hip[i_planet])
        end
        # Hippacentre indexes `orbit_solutions_per_planet[k][i_start + i]` for
        # i in 1:n_hip; with the planet_sols_cache layout (hip rows first),
        # i_start = 0 maps i directly to the hip portion.
        orbit_sol_starts_hip = ntuple(_ -> 0, n_planets)

        _simulate_skypath_hippacentre_combined!(
            Δα_mas_hip, Δδ_mas_hip, σ_inflation_hip,
            like.hip_table,
            orbits, planet_masses_msol_hip, flux_ratios_hip,
            planet_sols_cache, orbit_sol_starts_hip, T,
            HIPPARCOS_GRID_STEP_ARCSEC,
        )

        # NB: extract the catalog 5-parameter bias with the *uninflated* σ — the
        # Hipparcos pipeline that produced the catalog used point-source σ, so to
        # reproduce the bias it absorbed we must weight the LSQ the same way.
        # σ_inflation_hip is propagated separately into the IAD-residual likelihood
        # (`σ_iad = hypot(sres_renorm * σ_inflation_hip, hip_iad_jitter)` below).
        # Cached weighted-LSQ path: `_ensure_pinv_5_hip!` returns the
        # `Q = pinv(A./σ) ./ σ'` matrix, so the LSQ solution is a single
        # 5×N matrix-vector multiply against the *unscaled* RHS.
        Q_hip = _ensure_pinv_5_hip!(like)
        if like.include_iad
            out = fit_5param_pinv(Q_hip, like.hip_table, Δα_mas_hip, Δδ_mas_hip, like.hip_table.res)
        else
            out = fit_5param_pinv(Q_hip, like.hip_table, Δα_mas_hip, Δδ_mas_hip, 0.0)
        end
        Δα_h, Δδ_h, Δpmra_h, Δpmdec_h = out.parameters
        # Track magnitude of the BINARYS-predicted PM bias for the catalog-
        # likelihood epistemic inflation downstream.  Σ_h compares against the
        # Hipparcos catalog PMs (`dist_hip` is in mas/yr — see line 174), so
        # the relevant bias components are Δpmra_h and Δpmdec_h.
        hip_bias_pm_sq = Δpmra_h^2 + Δpmdec_h^2
        α_h₀, δ_h₀, pmra_h₀, pmdec_h₀ = propagate_astrom(orbits, like.catalog.epoch_ra_hip_mjd, like.catalog.epoch_dec_hip_mjd)
        μ_h = @SVector [pmra_h₀ + Δpmra_h, pmdec_h₀ + Δpmdec_h]


        ################################
        # Hipparcos IAD
            

        # We can include the Hipparcos IAD in the following way
        #=
            * The Hipaprcos IAD from the Java tool (what we have) is not necessarily consistent with the composite catalog point used in the HGCA, which is more accurate
            * We want to keep fitting the composite point with well-calibrated uncertainties compared to Gaia; however, we can still additionally fit the Hipparcos IAD
            * using a flexible offset and linear trend, just like RV data.
            * We need the following user-supplied parameters:
                * iad_Δra ("zero point offset for pmra [mas]")
                * iad_Δdec ("zero point offset for pmdec [mas]")
                * iad_pmra ("slope in pmra [mas/yr] ")
                * iad_pmdec ("slope in pmra [mas/yr] ")
                * hip_iad_jitter (astrometric excess jitter term)
            * This includes information from the IAD in a way that's still consistent with the calibrated composite Hipparcos catalog values for ra, dec, pmra, pmdec.
            * This doesn't count any information twice, because we here we are only fitting curvature, while above we are only fitting positions and proper motions.

            Probably, one could put informed priors on most of these properies but we leave them wide open.
        =#
        # Only when the :iad_hip likelihood row is active: it is the sole
        # consumer of iad_resid, and when it is excluded (epoch subset) the
        # iad_* nuisance parameters are stripped from θ_obs entirely, so the
        # destructure below would throw.
        if :iad_hip ∈ like.table.kind
            (;iad_Δra,
                iad_Δdec,
                iad_pmra,
                iad_pmdec,
                iad_Δplx) = θ_obs


            # like.hip_table.res, like.hip_table.sres
            # Δα_mas_hip, Δδ_mas_hip

            # Hoist loop-invariant terms.  `hip_sol.plx + iad_Δplx`, the
            # division by `julian_year`, and the PM-error effective values
            # are all constant across i_epoch.
            plx_at_epoch = like.hip_sol.plx + iad_Δplx
            inv_julian_year = inv(julian_year)
            iad_pmra_eff  = iad_pmra  - Δpmra_h
            iad_pmdec_eff = iad_pmdec - Δpmdec_h
            iad_Δra_eff   = iad_Δra   - Δα_h
            iad_Δdec_eff  = iad_Δdec  - Δδ_h

            # The Hipparcos IAD residual is the perpendicular distance from
            # the predicted model point (α✱_model_w_p, δ_model_w_p) to the
            # measured abscissa line, which `distance_point_to_line` evaluates
            # as |cross(r₂-r₁, r₁-r₀)| / ‖r₂-r₁‖.  By construction the line
            # endpoints are α✱ₘ = α✱ₐ ± sinϕ and δₘ = δₐ ∓ cosϕ, so r₂-r₁ has
            # length √(sin²ϕ + cos²ϕ)·2 = 2.  After substituting that and
            # cancelling the per-axis sinϕ·cosϕ cross-terms, the distance
            # reduces to the scalar scan-projection
            #     resid = |α✱ₐ·cosϕ + δₐ·sinϕ
            #            − (α✱_model_w_p·cosϕ + δ_model_w_p·sinϕ)|
            # which we evaluate without ever forming the 2D vectors.  The
            # parallax-along-scan factor is already precomputed in the table
            # (`parallaxFactorAlongScan`), so the per-epoch x/y/z sind/cosd
            # multiplies that this loop previously performed collapse to a
            # single multiply.  α✱ₐ·cosϕ + δₐ·sinϕ = res + Δα✱·cosϕ +
            # Δδ·sinϕ because α✱ₐ = res·cosϕ + Δα✱ and δₐ = res·sinϕ + Δδ
            # (cos²+sin² = 1).
            hip_epoch    = like.hip_table.epoch
            hip_cosϕ     = like.hip_table.cosϕ
            hip_sinϕ     = like.hip_table.sinϕ
            hip_plxFact  = like.hip_table.parallaxFactorAlongScan
            hip_projMeas = like.hip_table.proj_meas_alongscan
            @inbounds @simd for i_epoch in eachindex(hip_epoch, Δα_mas_hip, Δδ_mas_hip)
                cosϕ = hip_cosϕ[i_epoch]
                sinϕ = hip_sinϕ[i_epoch]
                Δt = (hip_epoch[i_epoch] - hipparcos_catalog_epoch_mjd) * inv_julian_year
                α_off = iad_Δra_eff  + Δt * iad_pmra_eff  + Δα_mas_hip[i_epoch]
                δ_off = iad_Δdec_eff + Δt * iad_pmdec_eff + Δδ_mas_hip[i_epoch]
                proj_model = α_off * cosϕ + δ_off * sinϕ + plx_at_epoch * hip_plxFact[i_epoch]
                iad_resid[i_epoch] = abs(hip_projMeas[i_epoch] - proj_model)
            end

            # @show θ_system.ra  θ_system.dec iad_Δra iad_Δdec iad_pmra iad_pmdec
            # f,a,p=Main.scatterlines(
            #     α✱_models, δ_models
            # )
            # Main.stem(f[2,1], α✱_models, iad_resid)
            # f|>display
        end

    end


    ################################
    # H-G (if H is available) and DR3-DR2

    # Simple linear approximation: don't deal with curvature & secular acceleration directly
    if absolute_orbits

        if isnothing(like.catalog.dist_hip)
            pmra_hg_model = zero(T)
            pmdec_hg_model = zero(T)
        else
            # HG
            Δα_hg_prop = (α_dr3₀ - α_h₀)*60*60*1000*cosd((δ_dr3₀ + δ_h₀)/2)
            Δδ_hg_prop = (δ_dr3₀ - δ_h₀)*60*60*1000
            pmra_hg_model = (Δα_dr3 - Δα_h + Δα_hg_prop) / (
                like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_hip_mjd
            )*julian_year
            pmdec_hg_model = (Δδ_dr3 - Δδ_h + Δδ_hg_prop) / (
                like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_hip_mjd
            )*julian_year
        end


        # DR3-DR2
        Δα_dr32_prop = (α_dr3₀ - α_dr2₀)*60*60*1000*cosd((δ_dr3₀ + δ_dr2₀)/2)
        Δδ_dr32_prop = (δ_dr3₀ - δ_dr2₀)*60*60*1000
        pmra_dr32_model = (Δα_dr3 - Δα_dr2 + Δα_dr32_prop) / (
            like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_dr2_mjd
        )*julian_year
        pmdec_dr32_model = (Δδ_dr3 - Δδ_dr2 + Δδ_dr32_prop) / (
            like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_dr2_mjd
        )*julian_year

    else
        if isnothing(like.catalog.dist_hip)
            pmra_hg_model = zero(T)
            pmdec_hg_model = zero(T)
        else
            pmra_hg_model = (Δα_dr3 - Δα_h) / (
                    like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_hip_mjd
            )*julian_year + θ_system.pmra
            pmdec_hg_model = (Δδ_dr3 - Δδ_h) / (
                like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_hip_mjd
            )*julian_year + θ_system.pmdec
        end

        pmra_dr32_model = (Δα_dr3 - Δα_dr2) / (
                like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_dr2_mjd
        )*julian_year + θ_system.pmra
        pmdec_dr32_model = (Δδ_dr3 - Δδ_dr2) / (
            like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_dr2_mjd
        )*julian_year + θ_system.pmdec


    end


    # μ_hg = @SVector [pmra_hg_model - Δpmra_dr3, pmdec_hg_model - Δpmdec_dr3]
    # μ_dr32 = @SVector [pmra_dr32_model - Δpmra_dr3, pmdec_dr32_model - Δpmdec_dr3]
    μ_hg = @SVector [pmra_hg_model, pmdec_hg_model]
    μ_dr32 = @SVector [pmra_dr32_model, pmdec_dr32_model]

    ##############################
    # DR3 UEVA calculation and uncertainty deflation
    # From Gaia catalog:
    (;
        astrometric_chi2_al_dr3,
        astrometric_n_good_obs_al_dr3,
        astrometric_matched_transits_dr3,
    ) = like.catalog

    N = astrometric_n_good_obs_al_dr3
    N_FoV = astrometric_matched_transits_dr3
    N_AL = N / N_FoV

    # Calculate Gaia's published UEVA (what they measured)
    if like.ueva_mode == :EAN
        UEVA_Gaia = like.catalog.astrometric_excess_noise_dr3^2 + σ_att^2 + σ_AL^2
    elseif like.ueva_mode == :RUWE
        ruwe_dr3 = like.catalog.ruwe_dr3
        u0 = 1/ruwe_dr3 * sqrt(astrometric_chi2_al_dr3/(N - gaia_n_dof))
        UEVA_Gaia = (ruwe_dr3 * u0)^2 * σ_formal^2
    elseif like.ueva_mode == :none
        # See the matching branch in the all-inactive fast path: placeholder
        # only, masked out of the likelihood, deflation pinned to 1 below.
        UEVA_Gaia = σ_formal^2
    else
        error("Unsupported mode (should be :EAN, :RUWE or :none, was $(like.ueva_mode))")
    end

    # Calculate expected UEVA for a single star (Eq. D.8 from paper)
    μ_UEVA_single = (N_AL / (N - gaia_n_dof)) * 
                ((N_FoV - gaia_n_dof) * σ_calib^2 + N_FoV * σ_AL^2)

    # And its standard deviation (Eq. D.9)
    σ_UEVA_single = sqrt(
        2 * N_AL / (N - gaia_n_dof)^2 * (
            N_AL * (N_FoV - gaia_n_dof) * σ_calib^4 + 
            N_FoV * σ_AL^4 + 
            2 * N_FoV * σ_AL^2 * σ_calib^2
        )
    )

    μ_1_3 = UEVA_Gaia^(1/3)
    UEVA_unc = σ_UEVA_single * μ_UEVA_single^(-2/3) / 3 # divide by 3 due to cube root transformation

    # Calculate model-predicted UEVA from our fit.
    # chi_squared_astro sums squared residuals over the transits actually
    # modeled (the DR3 window slice), while the UEVA normalizations below
    # assume N_FoV of them. The counts agree except in the GOST-shortfall
    # case (fewer in-window forecast epochs than catalog matched transits,
    # where `transits` falls back to all epochs); rescale so the predicted
    # excess per companion stays consistent with the N_FoV normalization.
    n_dr3_modeled = iend_dr3 - istart_dr3 + 1
    chi2_astro_scaled = out_dr3.chi_squared_astro * N_AL * (N_FoV / n_dr3_modeled)
    UEVA_model_raw = (chi2_astro_scaled * σ_formal^2) / (N - gaia_n_dof)

    # For the UEVA likelihood, use cube-root transformation (Eq. 27, Sect 5.1.1)
    UEVA_model_1 = (chi2_astro_scaled * σ_formal^2) / (N_AL * N_FoV - gaia_n_dof)
    UEVA_model = cbrt(UEVA_model_1 + μ_UEVA_single)

    # Calculate the "deflation factor" -- the amount of Gaia's inflated uncertainties
    # that come from our now-explained companion model

    # What a 5-param fit would measure with this companion model
    UEVA_predicted = UEVA_model_raw + μ_UEVA_single

    # Deflation factor
    deflation_factor_raw = sqrt(μ_UEVA_single / UEVA_Gaia)
    # equivalent to :
    # deflation_factor_raw = sqrt(1 - UEVA_orbital_perturb / UEVA_Gaia) 


    # @show deflation_factor_raw 

    # Clamp to valid range.  With ueva_mode == :none the deflation is not
    # identified (no calibrated σ_AL/σ_att/σ_calib behind μ_UEVA_single), so
    # Gaia's published DR3/DR32 uncertainties are used as-is.
    deflation_factor_dr3 = like.ueva_mode == :none ? one(deflation_factor_raw) :
        (deflation_factor_raw > 1.0 ? 1.0 : deflation_factor_raw)

    # # for data simulation purposes, here is an estimate of what these parameters would produce for RUWE
    # # given everything we know about the gaia uncertainties etc.
    # # Given: UEVA, u0, σ_formal, calculate ruwe
    # # u0 = 1/ruwe_dr3*sqrt(astrometric_chi2_al_dr3/(astrometric_n_good_obs_al_dr3-gaia_n_dof))
    # # UEVA = ruwe_dr3^2 * u0^2 * σ_formal^2
    # # UEVA/( u0^2 * σ_formal^2) = ruwe_dr3^2
    # ruwe_dr3 = sqrt(UEVA_model/( u0^2 * σ_formal^2))

    # Forward-model the Gaia RV uncertainty using the approach from the "paired"
    # tool (Chance et al. 2022, https://arxiv.org/abs/2206.11275).
    # σ_rv_per_transit is parameterized directly in linear km/s (matching the prior
    # and likelihood usage); no log-space conversion here.
    if :rv_dr3 ∈ like.table.kind
        σ_rv_per_transit = θ_obs.σ_rv_per_transit  # per-transit RV uncertainty in km/s

        # Simulate RV measurements at Gaia epochs — into a bump-allocated
        # buffer (caller's @no_escape provides the scope) so we don't
        # heap-alloc a fresh Vector per evaluation, and compute sample
        # variance with an explicit Welford-style loop to avoid the
        # `sum((rv_model .- rv_mean).^2)` broadcast temporary.
        n_rv_ep = length(gaia_table_rv.epoch)
        rv_model = Bumper.alloc!(Bumper.default_buffer(), T, n_rv_ep)
        fill!(rv_model, zero(T))

        # When `θ_obs.transits_rv` is a prefix of `θ_obs.transits` (guaranteed
        # for constructor-built variables — both come from partialsortperm of
        # the same transit_priorities, so the top n_rv selected for RV equal
        # the first n_rv of the top N_FoV selected for astrometry),
        # gaia_table_rv.epoch[i] == gaia_table.epoch[i] for i in 1..n_rv_ep.
        # The planet_sols_cache's gaia portion (rows n_hip+1..n_hip+n_gaia)
        # thus already holds the required orbit solutions; skip the per-epoch
        # orbitsolve.  Verify the full epoch prefix (cheap, O(n_rv) compares)
        # rather than assuming it, so user-supplied `variables` with an
        # arbitrary transits_rv subset fall back to the exact path.
        rv_use_cache = n_rv_ep <= n_gaia_cache
        if rv_use_cache
            @inbounds for i in 1:n_rv_ep
                if gaia_table_rv.epoch[i] != gaia_table.epoch[i]
                    rv_use_cache = false
                    break
                end
            end
        end
        if rv_use_cache
            @inbounds for i_planet in 1:n_planets
                planet_mass_msol = θ_system.planets[i_planet].mass*Octofitter.mjup2msol
                planet_mass_msol == 0.0 && continue
                sols = planet_sols_cache[i_planet]
                for i in 1:n_rv_ep
                    sol = sols[n_hip_cache + i]
                    rv_model[i] += radvel(sol, planet_mass_msol)/1e3
                end
            end
        else
            for (i_planet, (orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
                planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
                if planet_mass_msol == 0.0
                    continue
                end
                # Accumulate the RV contribution from this planet at each epoch.
                for (i, epoch) in enumerate(gaia_table_rv.epoch)
                    sol = orbitsolve(orbit, epoch)
                    rv_model[i] += radvel(sol, planet_mass_msol)/1e3 # barycentric rv in km/s
                end
            end
        end

        # Calculate sample variance (Eq. A2 of Chance et al. 2022)
        rv_sum = zero(T)
        @inbounds for i in 1:n_rv_ep
            rv_sum += rv_model[i]
        end
        rv_mean = rv_sum / n_rv_ep
        rv_sumsq = zero(T)
        @inbounds for i in 1:n_rv_ep
            d = rv_model[i] - rv_mean
            rv_sumsq += d * d
        end
        sample_variance = rv_sumsq / (n_rv_ep - 1)

        N_rv = like.catalog.rv_nb_transits
        rv_dof = N_rv - 1

        ε_catalog = like.catalog.radial_velocity_error

        # Convert catalog error back to sample variance (Eq. A4)
        s_catalog_squared = (2 * N_rv / π) * (ε_catalog^2 - 0.113^2)
    else
        # rv_chi2_stat = convert(T, NaN)
        rv_dof = convert(T, NaN)
        s_catalog_squared  = convert(T, NaN)
        rv_dof = convert(T, NaN)
        rv_mean = convert(T, NaN)
        sample_variance = convert(T, NaN)
        s_catalog_squared = convert(T, NaN)
    end

    # Adjust the reference frame such that, effectively, the pmra/pmdec system variables are referring to the primary
    # instead of the barycentre.
    # Specifically, the primary's proper motion at this epoch:
    μ_h    = μ_h     .- @SVector [Δpmra_dr3, Δpmdec_dr3,]
    μ_hg   = μ_hg    .- @SVector [Δpmra_dr3, Δpmdec_dr3,]
    μ_dr2  = μ_dr2   .- @SVector [Δpmra_dr3, Δpmdec_dr3,]
    μ_dr32 = μ_dr32  .- @SVector [Δpmra_dr3, Δpmdec_dr3,]
    μ_dr3  = μ_dr3   .- @SVector [Δpmra_dr3, Δpmdec_dr3,]


    return (;

        # UEVA: EAN/RUWE
        UEVA_model,
        UEVA_unc,
        μ_1_3,

        # Packaged up nicely
        μ_h,
        μ_hg,
        μ_dr2,
        μ_dr32,
        μ_dr3,
        μ = (@SVector [μ_h[1],μ_h[2],μ_hg[1],μ_hg[2],μ_dr2[1],μ_dr2[2],μ_dr32[1],μ_dr32[2],μ_dr3[1],μ_dr3[2],UEVA_model,sample_variance]),

        # Magnitude (squared) of the BINARYS-predicted PM bias on the
        # Hipparcos catalog point.  Consumed by the ε_binarys epistemic
        # inflation in ln_like.  Zero when no Hipparcos data are present.
        hip_bias_pm_sq,

        n_dr3 = iend_dr3 - istart_dr3 + 1,
        n_dr2 = length(kk_dr2),

        # rv_chi2_stat,
        rv_dof,
        rv_mean,
        sample_variance,
        s_catalog_squared,

        deflation_factor_dr3,

        # Individual
        # TODO: get rid of these
        pmra_hip_model=μ_h[1],
        pmdec_hip_model=μ_h[2],
        pmra_hg_model=μ_hg[1],
        pmdec_hg_model=μ_hg[2],
        pmra_dr2_model=μ_dr2[1],
        pmdec_dr2_model=μ_dr2[2],
        pmra_dr32_model=μ_dr32[1],
        pmdec_dr32_model=μ_dr32[2],
        pmra_dr3_model=μ_dr3[1],
        pmdec_dr3_model=μ_dr3[2],

        Δα_dr3, Δδ_dr3, Δpmra_dr3, Δpmdec_dr3


    )
end



# Generate new astrometry observations for completeness simulations
function Octofitter.generate_from_params(like::G23HObs, ctx::SystemObservationContext; add_noise)
    (; θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = ctx

    sim = simulate(like, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    if isnothing(sim)
        error("G23HObs simulate returned nothing during data generation (duplicate transit indices?)")
    end
    (; μ_h, μ_hg, μ_dr2, μ_dr32, μ_dr3, UEVA_model, UEVA_unc, μ_1_3, sample_variance, n_dr2, n_dr3) = sim

    catalog = like.catalog
    has_hip = !isnothing(catalog.dist_hip)

    # ── 0. HGCA nonlinear proper-motion correction ──
    # `ln_like` adds catalog.nonlinear_dpmra/dpmdec to the model μ_hg (and 2× to
    # μ_h) for absolute orbits when Hipparcos is present (see the "Apply nonlinear
    # correction for absolute orbits" block above).  The simulated catalog PMs are
    # built below from the bare μ returned by `simulate`, so without mirroring the
    # correction here the synthetic data sits on a different convention than the
    # likelihood's model — a phantom HG/Hip proper-motion anomaly of magnitude
    # nonlinear_dpm at the true parameters.  Apply the identical correction so the
    # forward simulator and the likelihood agree.
    absolute_orbits = false
    for orbit in orbits
        absolute_orbits |= orbit isa AbsoluteVisual
    end
    if absolute_orbits && has_hip
        μ_hg = μ_hg + @SVector [catalog.nonlinear_dpmra, catalog.nonlinear_dpmdec]
        μ_h  = μ_h  + @SVector [2catalog.nonlinear_dpmra, 2catalog.nonlinear_dpmdec]
    end

    # ── 1. UEVA / RUWE / EAN simulation and DR3 uncertainty inflation ──
    #
    # When a companion is present, its astrometric perturbation degrades Gaia's
    # 5-parameter fit, increasing the chi² and RUWE/EAN. Gaia then inflates all
    # formal parameter uncertainties by this excess. Our likelihood deflates them
    # back when the companion model explains the excess. For the simulation, we
    # must produce catalog errors that are inflated by the companion we inject,
    # so that the deflation mechanism works correctly during sampling.
    #
    # The inflation factor is √(UEVA_new / UEVA_original), which transforms
    # the original catalog errors (already inflated by the real star's excess noise)
    # to what Gaia would report with our injected companion instead.
    # Correlations are invariant since the inflation is a scalar scaling of the
    # covariance (all along-scan observations get the same noise increase).

    (;σ_att, σ_AL, σ_calib) = θ_obs
    σ_formal = sqrt(σ_att^2 + σ_AL^2)
    N = catalog.astrometric_n_good_obs_al_dr3
    N_FoV = catalog.astrometric_matched_transits_dr3
    N_AL = N / N_FoV
    gaia_n_dof = catalog.astrometric_params_solved_dr3 == 31 ? 5 : 6

    # Add noise to cube-root UEVA (the space in which the likelihood operates)
    new_UEVA_cuberoot = UEVA_model + (add_noise ? randn() * UEVA_unc : 0.0)
    new_UEVA = max(new_UEVA_cuberoot, 0.0)^3

    # Original catalog UEVA (what Gaia measured for the real star)
    UEVA_original = μ_1_3^3

    # Expected single-star UEVA (formal noise only, no companions)
    μ_UEVA_single = (N_AL / (N - gaia_n_dof)) *
        ((N_FoV - gaia_n_dof) * σ_calib^2 + N_FoV * σ_AL^2)

    # DR3 uncertainty inflation factor:
    # Deflate original errors to formal level, then re-inflate by the new companion.
    # new_err = original_err × √(μ_UEVA_single / UEVA_original) × √(new_UEVA / μ_UEVA_single)
    #         = original_err × √(new_UEVA / UEVA_original)
    # ueva_mode == :none: ln_like neither reads the UEVA datum nor deflates the
    # published covariances, so the simulator must not inflate them either —
    # otherwise simulated catalogs would carry an excess the fit never removes.
    inflation_dr3 = like.ueva_mode == :none ? 1.0 :
        sqrt(max(1.0, new_UEVA / max(eps(), UEVA_original)))

    # Back-calculate astrometric_chi2_al from the new total UEVA
    # UEVA_Gaia = chi2/(N-dof) * σ_formal²  ⟹  chi2 = UEVA * (N-dof) / σ_formal²
    new_chi2_al = max(Float64(N - gaia_n_dof), new_UEVA * (N - gaia_n_dof) / σ_formal^2)

    # Preserve the u0 calibration: new_ruwe = old_ruwe * √(new_chi2/old_chi2)
    new_ruwe = if like.ueva_mode == :none
        missing
    else
        old_chi2 = catalog.astrometric_chi2_al_dr3
        old_chi2 > 0 ? catalog.ruwe_dr3 * sqrt(new_chi2_al / old_chi2) : catalog.ruwe_dr3
    end

    # Back-calculate EAN: UEVA_Gaia = ean² + σ_att² + σ_AL²
    new_ean = sqrt(max(0.0, new_UEVA - σ_att^2 - σ_AL^2))

    # Inflated DR3 errors (PM and central-epoch position)
    new_pmra_dr3_error = catalog.pmra_dr3_error[1] * inflation_dr3
    new_pmdec_dr3_error = catalog.pmdec_dr3_error[1] * inflation_dr3
    new_ra_error_central_dr3 = catalog.ra_error_central_dr3 * inflation_dr3
    new_dec_error_central_dr3 = catalog.dec_error_central_dr3 * inflation_dr3
    # DR32 errors also inflate (they depend on DR3 positions).
    # This slightly overestimates since DR2 contributes too, but the ΔΣ_dr32
    # correction in ln_like handles the fine structure.
    new_pmra_dr32_error = catalog.pmra_dr32_error[1] * inflation_dr3
    new_pmdec_dr32_error = catalog.pmdec_dr32_error[1] * inflation_dr3

    # ── 2. Generate new PM values with correlated noise ──
    # Helper: draw correlated noise for a PM (ra,dec) pair
    function _draw_correlated_pm(μ_model, pmra_err, pmdec_err, pmra_pmdec_corr)
        if add_noise
            c = pmra_pmdec_corr * pmra_err * pmdec_err
            Σ = @SArray [
                pmra_err^2 c
                c pmdec_err^2
            ]
            noise = cholesky(Hermitian(Σ)).L * @SVector[randn(), randn()]
            return μ_model .+ noise
        else
            return SVector{2,Float64}(μ_model)
        end
    end

    # Hipparcos and HG: use original errors (not inflated by Gaia's DR3 pipeline)
    if has_hip
        new_pm_hip = _draw_correlated_pm(μ_h,
            catalog.pmra_hip_error[1], catalog.pmdec_hip_error[1], catalog.pmra_pmdec_hip[1])
        new_pm_hg = _draw_correlated_pm(μ_hg,
            catalog.pmra_hg_error[1], catalog.pmdec_hg_error[1], catalog.pmra_pmdec_hg[1])
    end
    # DR2/DR32/DR3: draw from the SAME covariance structure ln_like will assemble
    # at the truth parameters, otherwise the synthetic data violates the
    # likelihood's correlation model.  ln_like (a) couples the DR2 and DR3 PMs
    # with the cross block K = ρ·√Σ_dr2·√Σ_dr3' where ρ is the published
    # catalog `rho_dr2_dr3` (the catalog uncertainties were calibrated against
    # it), (b) deflates Σ_dr3 by deflation_factor_dr3² recomputed
    # at fit time from the *simulated* catalog, and (c) adds the deflation-driven
    # ΔΣ_dr32 adjustment to Σ_dr32.  Independent draws leave the two conditional
    # DR2/DR3 directions over-dispersed by 1/(1-ρ²) in whitened space — the fit
    # reads the excess as astrometric acceleration → spurious decades-period
    # companions.  Mirror ln_like exactly (verified by the joint-χ² MC test).
    ρ_dr3_dr2 = catalog.rho_dr2_dr3
    # deflation the fit will apply at truth: UEVA_Gaia reconstructed from the
    # NEW catalog is max(σ_formal², new_UEVA) in both :RUWE and :EAN modes
    # (the new_chi2_al / new_ean clamps above make the two expressions agree).
    UEVA_Gaia_fit = max(σ_formal^2, new_UEVA)
    deflation_truth = like.ueva_mode == :none ? 1.0 :
        min(1.0, sqrt(μ_UEVA_single / UEVA_Gaia_fit))

    c2 = catalog.pmra_pmdec_dr2[1] * catalog.pmra_dr2_error[1] * catalog.pmdec_dr2_error[1]
    Σ_dr2_gen = @SArray [
        catalog.pmra_dr2_error[1]^2 c2
        c2 catalog.pmdec_dr2_error[1]^2
    ]
    c3 = catalog.pmra_pmdec_dr3[1] * new_pmra_dr3_error * new_pmdec_dr3_error
    Σ_dr3_gen = (@SArray [
        new_pmra_dr3_error^2 c3
        c3 new_pmdec_dr3_error^2
    ]) .* deflation_truth^2

    if add_noise
        K_gen = ρ_dr3_dr2 * sqrt(Σ_dr2_gen) * sqrt(Σ_dr3_gen)'
        Σ_23 = SMatrix{4,4,Float64,16}([Σ_dr2_gen K_gen; K_gen' Σ_dr3_gen])
        noise_23 = cholesky(Hermitian(Σ_23)).L * @SVector[randn(), randn(), randn(), randn()]
        new_pm_dr2 = μ_dr2 .+ @SVector[noise_23[1], noise_23[2]]
        new_pm_dr3 = μ_dr3 .+ @SVector[noise_23[3], noise_23[4]]
    else
        new_pm_dr2 = SVector{2,Float64}(μ_dr2)
        new_pm_dr3 = SVector{2,Float64}(μ_dr3)
    end

    # DR32 with the fit-time ΔΣ_dr32 deflation adjustment (mirrors ln_like).
    Σ_pos_dr3_gen = @SMatrix [
        new_ra_error_central_dr3^2  catalog.ra_dec_corr_central_dr3*new_ra_error_central_dr3*new_dec_error_central_dr3
        catalog.ra_dec_corr_central_dr3*new_ra_error_central_dr3*new_dec_error_central_dr3  new_dec_error_central_dr3^2
    ]
    Σ_cross_gen = @SMatrix [
        ρ_dr3_dr2*new_ra_error_central_dr3*catalog.ra_error_central_dr2  ρ_dr3_dr2*catalog.ra_dec_corr_central_dr3*new_ra_error_central_dr3*catalog.dec_error_central_dr2
        ρ_dr3_dr2*catalog.ra_dec_corr_central_dr2*new_dec_error_central_dr3*catalog.ra_error_central_dr2  ρ_dr3_dr2*new_dec_error_central_dr3*catalog.dec_error_central_dr2
    ]
    ΔΣ_pos_gen = (deflation_truth^2 - 1) * Σ_pos_dr3_gen - (deflation_truth - 1) * (Σ_cross_gen + Σ_cross_gen')
    Δt_ra_gen = (catalog.epoch_ra_dr3_mjd - catalog.epoch_ra_dr2_mjd) / julian_year
    Δt_dec_gen = (catalog.epoch_dec_dr3_mjd - catalog.epoch_dec_dr2_mjd) / julian_year
    Tr_gen = @SMatrix [1/Δt_ra_gen 0.0; 0.0 1/Δt_dec_gen]
    ΔΣ_dr32_gen = Tr_gen * ΔΣ_pos_gen * Tr_gen'
    c32 = catalog.pmra_pmdec_dr32[1] * new_pmra_dr32_error * new_pmdec_dr32_error
    Σ_dr32_gen = (@SArray [
        new_pmra_dr32_error^2 c32
        c32 new_pmdec_dr32_error^2
    ]) + ΔΣ_dr32_gen
    if add_noise
        noise_32 = cholesky(Hermitian(Σ_dr32_gen)).L * @SVector[randn(), randn()]
        new_pm_dr32 = μ_dr32 .+ noise_32
    else
        new_pm_dr32 = SVector{2,Float64}(μ_dr32)
    end

    # ── 3. Hipparcos IAD residuals ──
    # The residuals capture the companion's curvature signal (non-linear sky path
    # after removing position, PM, and parallax via a 5-parameter fit) plus noise.
    n_hip = size(like.hip_table, 1)
    new_hip_res = zeros(n_hip)
    if has_hip && n_hip > 0
        Δα_hip = zeros(n_hip)
        Δδ_hip = zeros(n_hip)
        σ_inflation_hip = ones(n_hip)
        n_planets = length(orbits)
        planet_masses_msol = ntuple(i_planet -> θ_system.planets[i_planet].mass * Octofitter.mjup2msol, n_planets)
        # Companions with mass = 0 contribute nothing — zero out their flux ratio.
        flux_ratios_hip = ntuple(n_planets) do i_planet
            planet_masses_msol[i_planet] == 0.0 ? 0.0 :
                (θ_obs.fluxratio_hip isa Number ? θ_obs.fluxratio_hip : θ_obs.fluxratio_hip[i_planet])
        end
        orbit_sol_starts = ntuple(_ -> -1, n_planets)
        _simulate_skypath_hippacentre_combined!(
            Δα_hip, Δδ_hip, σ_inflation_hip,
            like.hip_table,
            orbits, planet_masses_msol, flux_ratios_hip,
            orbit_solutions, orbit_sol_starts, Float64,
            HIPPARCOS_GRID_STEP_ARCSEC,
        )

        # Project perturbation along scan direction
        b_hip = zeros(n_hip)
        for i in 1:n_hip
            b_hip[i] = Δα_hip[i] * like.hip_table.cosϕ[i] + Δδ_hip[i] * like.hip_table.sinϕ[i]
        end

        # 5-param fit absorbs position, PM, parallax (linear part)
        x_hip = like.A_prepared_5_hip \ b_hip
        model_hip = like.A_prepared_5_hip * x_hip

        # Curvature residual = perturbation - linear model
        new_hip_res .= b_hip .- model_hip
        if add_noise
            # Inflate per-transit residual σ by the BINARYS first-harmonic factor
            # (Leclerc et al. 2023, Eq. 15) when generating synthetic noise.
            # hip_iad_jitter only exists in θ_obs while the :iad_hip likelihood
            # row is active (the nuisance block is stripped otherwise); without
            # it, simulate the residual noise with zero excess jitter — the
            # residuals are not consumed by any active likelihood term then.
            hip_jitter = :iad_hip ∈ like.table.kind ? θ_obs.hip_iad_jitter : 0.0
            for i in 1:n_hip
                new_hip_res[i] += randn() * hypot(like.hip_table.sres_renorm[i] * σ_inflation_hip[i], hip_jitter)
            end
        end
    end

    # ── 4. Gaia RV simulation ──
    has_rv = :rv_dr3 ∈ like.table.kind
    new_rv_error = hasproperty(catalog, :radial_velocity_error) ? catalog.radial_velocity_error : NaN
    if has_rv && isfinite(sample_variance)
        σ_rv = θ_obs.σ_rv_per_transit  # per-transit RV uncertainty in km/s
        N_rv = catalog.rv_nb_transits
        # Non-centrality parameter from the companion's RV signal
        ncp = max(0.0, (N_rv - 1) * sample_variance / σ_rv^2)
        if add_noise
            ξ² = rand(NoncentralChisq(max(1, N_rv - 1), ncp))
        else
            ξ² = ncp + (N_rv - 1)  # E[noncentral χ²] = dof + ncp
        end
        S² = ξ² * σ_rv^2 / max(1, N_rv - 1)
        # Convert sample variance back to Gaia's reported radial_velocity_error
        # s² = (2N/π)(ε² - 0.113²)  ⟹  ε = √(s²π/(2N) + 0.113²)
        new_rv_error = sqrt(max(0.0, S² * π / (2 * N_rv) + 0.113^2))
    end

    # ── 5. Rebuild catalog with new values ──
    new_catalog = (; catalog...)

    if has_hip
        new_catalog = (; new_catalog...,
            pmra_hip = new_pm_hip[1],
            pmdec_hip = new_pm_hip[2],
            pmra_hg = new_pm_hg[1],
            pmdec_hg = new_pm_hg[2],
        )
    end
    new_catalog = (; new_catalog...,
        pmra_dr2 = new_pm_dr2[1],
        pmdec_dr2 = new_pm_dr2[2],
        pmra_dr32 = new_pm_dr32[1],
        pmdec_dr32 = new_pm_dr32[2],
        pmra_dr3 = new_pm_dr3[1],
        pmdec_dr3 = new_pm_dr3[2],
        # UEVA-related fields
        astrometric_chi2_al_dr3 = like.ueva_mode == :none ? missing : new_chi2_al,
        ruwe_dr3 = new_ruwe,
        astrometric_excess_noise_dr3 = like.ueva_mode == :none ? missing : new_ean,
        # Inflated DR3 uncertainties (as Gaia would report with companion present)
        pmra_dr3_error = new_pmra_dr3_error,
        pmdec_dr3_error = new_pmdec_dr3_error,
        ra_error_central_dr3 = new_ra_error_central_dr3,
        dec_error_central_dr3 = new_dec_error_central_dr3,
        pmra_dr32_error = new_pmra_dr32_error,
        pmdec_dr32_error = new_pmdec_dr32_error,
    )
    if has_rv
        new_catalog = (; new_catalog..., radial_velocity_error = new_rv_error)
    end

    # ── 6. Recompute MvNormal distributions from new PM values and inflated errors ──
    if has_hip
        c = catalog.pmra_pmdec_hip[1] * catalog.pmra_hip_error[1] * catalog.pmdec_hip_error[1]
        dist_hip = MvNormal(
            @SVector([new_catalog.pmra_hip, new_catalog.pmdec_hip]),
            @SArray[
                catalog.pmra_hip_error[1]^2 c
                c catalog.pmdec_hip_error[1]^2
            ]
        )
        c = catalog.pmra_pmdec_hg[1] * catalog.pmra_hg_error[1] * catalog.pmdec_hg_error[1]
        dist_hg = MvNormal(
            @SVector([new_catalog.pmra_hg, new_catalog.pmdec_hg]),
            @SArray[
                catalog.pmra_hg_error[1]^2 c
                c catalog.pmdec_hg_error[1]^2
            ]
        )
    else
        dist_hip = nothing
        dist_hg = nothing
    end

    c = catalog.pmra_pmdec_dr2[1] * catalog.pmra_dr2_error[1] * catalog.pmdec_dr2_error[1]
    dist_dr2 = MvNormal(
        @SVector([new_catalog.pmra_dr2, new_catalog.pmdec_dr2]),
        @SArray[
            catalog.pmra_dr2_error[1]^2 c
            c catalog.pmdec_dr2_error[1]^2
        ]
    )

    # DR32 and DR3: use inflated errors in the distributions
    c = catalog.pmra_pmdec_dr32[1] * new_pmra_dr32_error * new_pmdec_dr32_error
    dist_dr32 = MvNormal(
        @SVector([new_catalog.pmra_dr32, new_catalog.pmdec_dr32]),
        @SArray[
            new_pmra_dr32_error^2 c
            c new_pmdec_dr32_error^2
        ]
    )

    c = catalog.pmra_pmdec_dr3[1] * new_pmra_dr3_error * new_pmdec_dr3_error
    dist_dr3 = MvNormal(
        @SVector([new_catalog.pmra_dr3, new_catalog.pmdec_dr3]),
        @SArray[
            new_pmra_dr3_error^2 c
            c new_pmdec_dr3_error^2
        ]
    )

    new_catalog = (; new_catalog..., dist_hip, dist_hg, dist_dr2, dist_dr32, dist_dr3)

    # ── 7. Build new hip_table with simulated IAD residuals ──
    if n_hip > 0
        # Replace only the res column; all other columns (epochs, scan angles, etc.) stay
        hip_cols = Tables.columntable(like.hip_table)
        new_hip_table = Table(merge(hip_cols, (; res = new_hip_res)))
    else
        new_hip_table = like.hip_table
    end

    # ── 8. Update summary table pm values ──
    new_pm = copy(like.table.pm)
    new_σ_pm = copy(like.table.σ_pm)
    for (i, kind) in enumerate(like.table.kind)
        if kind == :ra_hip && has_hip
            new_pm[i] = new_catalog.pmra_hip
        elseif kind == :dec_hip && has_hip
            new_pm[i] = new_catalog.pmdec_hip
        elseif kind == :ra_hg && has_hip
            new_pm[i] = new_catalog.pmra_hg
        elseif kind == :dec_hg && has_hip
            new_pm[i] = new_catalog.pmdec_hg
        elseif kind == :ra_dr2
            new_pm[i] = new_catalog.pmra_dr2
        elseif kind == :dec_dr2
            new_pm[i] = new_catalog.pmdec_dr2
        elseif kind == :ra_dr32
            new_pm[i] = new_catalog.pmra_dr32
            new_σ_pm[i] = new_pmra_dr32_error
        elseif kind == :dec_dr32
            new_pm[i] = new_catalog.pmdec_dr32
            new_σ_pm[i] = new_pmdec_dr32_error
        elseif kind == :ra_dr3
            new_pm[i] = new_catalog.pmra_dr3
            new_σ_pm[i] = new_pmra_dr3_error
        elseif kind == :dec_dr3
            new_pm[i] = new_catalog.pmdec_dr3
            new_σ_pm[i] = new_pmdec_dr3_error
        elseif kind == :ueva_dr3
            new_pm[i] = like.ueva_mode == :RUWE ? new_catalog.ruwe_dr3 : new_catalog.astrometric_excess_noise_dr3
        end
    end
    table_cols = Tables.columntable(like.table)
    new_table = Table(merge(table_cols, (; pm = new_pm, σ_pm = new_σ_pm)))

    # ── 9. Construct new G23HObs ──
    obj = G23HObs{
        typeof(new_table),
        typeof(new_hip_table),
        typeof(like.gaia_table),
        typeof(new_catalog),
        typeof(like.hip_sol),
    }(
        new_table,
        like.priors,
        like.derived,
        new_hip_table,
        like.gaia_table,
        new_catalog,
        like.hip_sol,
        like.A_prepared_5_hip,
        like.A_prepared_5_dr2,
        like.A_prepared_5_dr3,
        like.include_iad,
        like.ueva_mode,
        Ref{Matrix{Float64}}(zeros(0, 0)),
        Ref{NTuple{4, Float64}}((0.0, 0.0, 0.0, 0.0)),
        Ref{Bool}(false),
    )
    # Note: the caches derive from hip_table.{sres,res}, which may differ in
    # the newly generated table — recompute eagerly for the new object.
    _ensure_pinv_5_hip!(obj)
    _ensure_hip_x_const!(obj)
    return obj
end

export G23HObs
