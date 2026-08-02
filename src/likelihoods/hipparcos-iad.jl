# ---------------------------------------------------
# Hipparcos intermediate astrometric data (IAD)
#
# Loading and reconstruction of the van Leeuwen (2007) per-transit abscissa
# residuals from the Java-tool release, plus the two Hipparcos instrument
# constants the G23H abscissa branch needs.
#
# This is *data loading*, not a likelihood: it produces the per-transit scan
# geometry table (`epoch`, `cosϕ`, `sinϕ`, `parallaxFactorAlongScan`,
# `res`, `sres`, …) that `G23HObs` reduces. There is deliberately no
# `HipparcosIADObs` observation type on this branch — the only consumer is
# G23H, and the Hipparcos grating response it applies is G23H-internal (see
# `_hippacentre!` in `g23h.jl`).
#
# Ported from `src/legacy/likelihoods/hipparcos.jl`. Two deliberate
# differences, both recorded in the port notes:
#
#   * v1 computed a "rigorous" (α, δ, ϖ(t)) sky path from an
#     `AbsoluteVisual` orbit, assigned `Δα✱`/`Δδ` from it, and then
#     *immediately overwrote both* with the tangent-plane version below. The
#     rigorous block was dead code in v1 and is not reproduced here.
#   * The `α✱ₘ`/`δₘ` two-point abscissa-line columns are dropped. They exist
#     for a `distance_point_to_line` residual formulation that G23H replaced
#     with the algebraically identical scalar projection
#     `proj_meas_alongscan`, and nothing reads them.
# ---------------------------------------------------

using Combinatorics: combinations

"""
Hipparcos main-grid step [arcsec] (Lindegren 1997, ESA SP-1200 vol. 3). This
is the spatial period of the modulating grid, and therefore the period of the
BINARYS atan2 photocentre response in projected separation.
"""
const HIPPARCOS_GRID_STEP_ARCSEC = 1.2074

"""
Resolution scale [arcsec] of the per-transit taper applied to the BINARYS
modulated photocentre. At projected separations ρ ≫ s the Hipparcos pipeline
either resolved the components into separate entries or rejected the row, so
the catalog point reflects the primary alone and the modulation must taper to
zero — which the atan2 formula does not do on its own. Anchored to the grid
step (Lindegren 1997 §3.3; van Leeuwen 2007 §6.4).
"""
const HIPPARCOS_RESOLUTION_ARCSEC = 1.207

"""
    α_resolve_hip(ρ_arcsec)

Gaussian resolution taper on the *full* projected separation (not the
along-scan projection: the resolution criterion is geometric, not
abscissa-projected).
"""
α_resolve_hip(ρ_arcsec) = exp(-(ρ_arcsec / HIPPARCOS_RESOLUTION_ARCSEC)^2)

# Table 1.2.3 of the 1997 catalog (vol. 1): stars for which the Hipparcos
# reduction modelled a time-varying parallax from the barycentric RV. Only
# applicable to the 1997 reduction; the van Leeuwen reduction did not use it.
const hipparcos_stellar_rvs_used_for_cat_solutions = Table([
    (hip_id=439, rv_kms=+22.9), (hip_id=3829, rv_kms=-38.0),
    (hip_id=5336, rv_kms=-98.1), (hip_id=15510, rv_kms=+86.7),
    (hip_id=19849, rv_kms=-42.7), (hip_id=24186, rv_kms=+245.5),
    (hip_id=26857, rv_kms=+105.6), (hip_id=54035, rv_kms=-84.3),
    (hip_id=54211, rv_kms=+68.8), (hip_id=57939, rv_kms=-99.1),
    (hip_id=70890, rv_kms=-16.0), (hip_id=71681, rv_kms=-18.1),
    (hip_id=71683, rv_kms=-26.2), (hip_id=74234, rv_kms=+308.0),
    (hip_id=74235, rv_kms=+294.3), (hip_id=86990, rv_kms=-115.0),
    (hip_id=87937, rv_kms=-111.0), (hip_id=99461, rv_kms=-129.8),
    (hip_id=104214, rv_kms=-64.8), (hip_id=104217, rv_kms=-64.3),
    (hip_id=108870, rv_kms=-40.4),
])

"""
    hipparcos_iad(; hip_id, catalog=datadep"Hipparcos_IAD", renormalize=true,
                  attempt_correction=true, is_van_leeuwen=true)

Load one star's Hipparcos intermediate astrometric data and reconstruct the
per-transit sky-path model the abscissa residuals are measured against.

Returns `(; table, hip_sol)`:

  - `hip_sol` — the catalog five-parameter solution and header quantities
    parsed from the Java-tool file.
  - `table` — one row per transit, carrying the scan geometry (`cosϕ`,
    `sinϕ`, `parallaxFactorAlongScan`), the abscissa residual `res` and its
    formal error `sres` (`reject` marks the scans the original reduction
    threw out, which it encodes as `sres ≤ 0`), the Nielsen-renormalized
    `sres_renorm`, the reconstructed reference sky path `Δα✱`/`Δδ` [mas], and
    `proj_meas_alongscan = res + Δα✱·cosϕ + Δδ·sinϕ` — the measured abscissa
    in the tangent plane, which is what an abscissa likelihood compares
    against.

`renormalize` applies Nielsen et al. (2020) Eq. 10; `attempt_correction`
applies the Appendix-A repair of G. M. Brandt et al. (2021) for the known
Java-tool duplicate-row corruption.
"""
function hipparcos_iad(; hip_id,
                       catalog=nothing,
                       renormalize::Bool=true,
                       attempt_correction::Bool=true,
                       is_van_leeuwen::Bool=true)
    if isnothing(catalog)
        catalog = datadep"Hipparcos_IAD"
    end
    file = @sprintf("H%06d.d", hip_id)
    fname = joinpath(catalog, "ResRec_JavaTool_2014", file[1:4], file)
    isfile(fname) || error(
        "no Hipparcos IAD file for HIP $hip_id at $fname. Pass `catalog=` " *
        "pointing at the directory containing `ResRec_JavaTool_2014/`.")
    lines = readlines(fname)

    # See Table 1 of arXiv:1108.4971 for the units.
    #   HIP    MCE    NRES NC isol_n SCE  F2     F1
    hip, mce, nres, nc, isol_n, sce, f2, f1 = parse.(Float64, split(lines[7])[2:end])
    #   Hp      B-V    VarAnn NOB NR
    hp, b_m_v, varann, nob, nr = parse.(Float64, split(lines[9])[2:end])
    (radeg, dedeg, plx, pm_ra, pm_de, e_ra, e_de, e_plx, e_pmra, e_pmde,
     dpmra, dpmde, e_dpmra, e_dpmde, ddpmra, ddpmde, e_ddpmra, e_ddpmde,
     upsra, upsde, e_upsra, e_upsde, var) = tryparse.(Float64, split(lines[11])[2:end])
    hip_sol = (; hip, mce, nres, nc, isol_n, sce, f2, f1,
        hp, b_m_v, varann, nob, nr,
        radeg, dedeg, plx, pm_ra, pm_de, e_ra, e_de, e_plx, e_pmra, e_pmde,
        dpmra, dpmde, e_dpmra, e_dpmde, ddpmra, ddpmde, e_ddpmra, e_ddpmde,
        upsra, upsde, e_upsra, e_upsde, var)

    if isol_n ∉ (5, 7, 9)
        @warn "Only Hipparcos solution types 5, 7 and 9 are supported; this is $isol_n."
    end

    rows = NamedTuple[]
    for line in lines[13:end]
        startswith(line, '#') && continue
        iorb, epoch, parf, cosϕ, sinϕ, res, sres = split(line)
        push!(rows, (; iorb=parse(Int, iorb), epoch_yrs=parse(Float64, epoch),
            parf=parse(Float64, parf), cosϕ=parse(Float64, cosϕ),
            sinϕ=parse(Float64, sinϕ), res=parse(Float64, res),
            sres=parse(Float64, sres)))
    end
    iad = FlexTable(rows)

    # `sres ≤ 0` is how the reduction encodes a scan it rejected.
    iad.reject = iad.sres .<= 0
    any(iad.reject) && @warn "rejected Hipparcos scans present" count(iad.reject)

    if renormalize
        # Nielsen et al. (2020) Eq. 10 — undo the catalog's uncertainty
        # renormalization, which otherwise hides real scatter.
        D = length(iad.sres) - isol_n
        f = (f2 * √(2 / 9D) + 1 - (2 / 9D))^(3 / 2)
        iad.sres_renorm = iad.sres .* f
    else
        iad.sres_renorm = iad.sres
    end

    if attempt_correction
        iad, was_corrected = _correct_iad_corruption(iad)
        was_corrected && @info "Corrected corrupted Hipparcos IAD rows " *
                               "(Appendix A of G. M. Brandt et al. 2021)."
    end

    iad.epoch = hipparcos_catalog_epoch_mjd .+ iad.epoch_yrs .* julian_year
    earth = geocentre_position_query.(iad.epoch)
    table = FlexTable(eachcol(iad)..., eachcol(earth)...)

    # Reconstruct the reference sky path the residuals were measured against
    # (Nielsen et al., β Pic, Eqs. 1–2): a tangent-plane model about the
    # catalog position with the catalog proper motion and parallax.
    μα✱ = hip_sol.pm_ra
    μδ = hip_sol.pm_de
    α₀ = hip_sol.radeg
    δ₀ = hip_sol.dedeg

    # The 1997 reduction applied a time-varying parallax for 21 high-RV
    # stars; the van Leeuwen reduction did not, and tests show it must not be
    # applied when undoing the latter.
    rv_kms = if is_van_leeuwen
        0.0
    else
        i = findfirst(==(hip_id), hipparcos_stellar_rvs_used_for_cat_solutions.hip_id)
        isnothing(i) ? 0.0 : hipparcos_stellar_rvs_used_for_cat_solutions.rv_kms[i]
    end
    table.rv_kms = fill(rv_kms, length(table.epoch_yrs))

    δdist_pc_δt_sec = rv_kms / IAU_pc2km / PlanetOrbits.sec2day
    Δdist_pc = δdist_pc_δt_sec .* (table.epoch .- hipparcos_catalog_epoch_mjd)
    table.plx_vs_time = 1000 ./ (1000 / hip_sol.plx .+ Δdist_pc)
    table.Δα✱ = @. table.plx_vs_time * (table.x * sind(α₀) - table.y * cosd(α₀)) +
                   (table.epoch - hipparcos_catalog_epoch_mjd) / julian_year * μα✱
    table.Δδ = @. table.plx_vs_time * (table.x * cosd(α₀) * sind(δ₀) +
                                       table.y * sind(α₀) * sind(δ₀) -
                                       table.z * cosd(δ₀)) +
                  (table.epoch - hipparcos_catalog_epoch_mjd) / julian_year * μδ

    table.scanAngle_rad = @. atan(table.sinϕ, table.cosϕ)
    table.parallaxFactorAlongScan = @. (
        (table.x * sind(α₀) - table.y * cosd(α₀)) * table.cosϕ +
        (table.x * cosd(α₀) * sind(δ₀) + table.y * sind(α₀) * sind(δ₀) -
         table.z * cosd(δ₀)) * table.sinϕ)

    # The measured abscissa, projected on scan. Follows from
    # α✱ₐ = res·cosϕ + Δα✱, δₐ = res·sinϕ + Δδ and cos²ϕ + sin²ϕ = 1, so a
    # likelihood never has to form the 2D abscissa line.
    table.proj_meas_alongscan = @. table.res + table.Δα✱ * table.cosϕ +
                                   table.Δδ * table.sinϕ

    return (; table=Table(table), hip_sol)
end

"""
    hipparcos_recalibrate!(table)

Apply the Brandt, Michalik & Brandt recalibration of the van Leeuwen IAD:
+0.140 mas on the residuals and 2.25 mas of extra dispersion, which
mitigates the overfitting their statistical study documents.

Mutates `table` and recomputes `proj_meas_alongscan`, which is derived from
`res` — the v1 code path baked that column *before* the recalibration, so the
shift never reached the abscissa channel (audit 2026-07-02).
"""
function hipparcos_recalibrate!(table)
    table.res .+= 0.140
    table.sres_renorm .= hypot.(table.sres_renorm, 2.25)
    table.proj_meas_alongscan .= table.res .+
                                 table.Δα✱ .* table.cosϕ .+ table.Δδ .* table.sinϕ
    return table
end

# --- Java-tool duplicate-row corruption (G. M. Brandt et al. 2021, App. A) ---

"""
    _detect_iad_corruption(iad) -> Int

Number of corrupted trailing rows: the Java-tool export sometimes repeats the
last orbit's abscissa errors. Signature is four trailing rows from one orbit
whose formal errors read as a duplicated pair.
"""
function _detect_iad_corruption(iad)
    length(iad.sres_renorm) < 4 && return 0
    last_idx = length(iad.sres_renorm)
    idx = last_idx-3:last_idx
    orbits = iad.iorb[idx]
    all(==(orbits[1]), orbits) || return 0
    s = iad.sres_renorm[idx]
    if abs(s[1] - s[4]) < 0.0001 && abs(s[2] - s[3]) < 0.1
        return 3
    end
    return 0
end

"""
    _find_best_iad_correction(iad, n_corrupt)

Translation of `htof.parse.find_epochs_to_reject_java` (G. M. Brandt): the
repeated residuals are always dropped, and the orbit combination whose χ²
partial-derivative vector is closest to stationary identifies which orbits to
remove.
"""
function _find_best_iad_correction(iad, n_corrupt)
    n = length(iad.epoch_yrs)
    resid_reject_idx = [n - i + 1 for i in 1:n_corrupt]
    orbits_to_keep = trues(n)
    residuals_to_keep = trues(n)
    residuals_to_keep[resid_reject_idx] .= false
    residual_factors = (iad.res ./ (iad.sres_renorm .^ 2))[residuals_to_keep]
    dt = iad.epoch_yrs
    _orbit_factors = [iad.parf iad.cosϕ iad.sinϕ dt .* iad.cosϕ dt .* iad.sinϕ]

    best_reject = Int[]
    best_partial = Inf
    for orbit_to_reject in combinations(1:n, n_corrupt)
        orbits_to_keep[orbit_to_reject] .= false
        orbit_factors = @view _orbit_factors[orbits_to_keep, :]
        chi2_vector = (2 .* residual_factors .* orbit_factors)'
        partial = sqrt(sum(sum(chi2_vector, dims=2) .^ 2))
        if partial < best_partial
            best_partial = partial
            best_reject = copy(orbit_to_reject)
        end
        orbits_to_keep .= true
    end
    return best_reject, best_partial
end

function _correct_iad_corruption(iad)
    n_corrupt = _detect_iad_corruption(iad)
    n_corrupt == 0 && return iad, false
    indices_to_remove, _ = _find_best_iad_correction(iad, n_corrupt)
    n = length(iad.epoch_yrs)
    # Drop the repeated residuals from the tail and the identified orbits'
    # geometry, then re-pair what is left.
    keep_resid = trues(n)
    keep_resid[[n - i + 1 for i in 1:n_corrupt]] .= false
    keep_orbit = trues(n)
    keep_orbit[indices_to_remove] .= false
    out = FlexTable(
        iorb=iad.iorb[keep_orbit],
        epoch_yrs=iad.epoch_yrs[keep_orbit],
        parf=iad.parf[keep_orbit],
        cosϕ=iad.cosϕ[keep_orbit],
        sinϕ=iad.sinϕ[keep_orbit],
        res=iad.res[keep_resid],
        sres=iad.sres[keep_resid],
        reject=iad.reject[keep_resid],
        sres_renorm=iad.sres_renorm[keep_resid],
    )
    return out, true
end
