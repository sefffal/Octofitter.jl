#=
Shared utilities for Gaia astrometric likelihoods.

This file contains constants, helper functions, and SPICE-based ephemeris
queries used across multiple Gaia likelihood implementations (g23h.jl,
gaia-dr4.jl, hipparcos.jl, hgca-linfit.jl).
=#

using SPICE
using DataDeps
using Dates
using HTTP

# ──────────────────────────────────────────────────────────────────────
# Gaia data release metadata constants
# ──────────────────────────────────────────────────────────────────────

const meta_gaia_DR1 = (;
    start_mjd=mjd("2014-07-25"),
    stop_mjd =mjd("2015-09-16"),
    ref_epoch_mjd=57023.25 # J2015.0
)

const meta_gaia_DR2 = (;
    start_mjd=mjd("2014-07-25"), # 10:30 UTC
    stop_mjd =mjd("2016-05-23"), # 11:35 UTC
    ref_epoch_mjd=57205.875 # J2015.5
)

const meta_gaia_DR3 = (;
    start_mjd=mjd("2014-07-25"), # 10:30 UTC
    stop_mjd =mjd("2017-05-28"), # 11:35 UTC
    ref_epoch_mjd=57388.5, # J2016.0
)

# ──────────────────────────────────────────────────────────────────────
# Time conversion utilities
# ──────────────────────────────────────────────────────────────────────

# OBMT (On board mission time) to MJD converter
tcb_at_gaia_2mjd(tcb_gaia) = jd2mjd(tcb_gaia+2455197.5)

"""
    obmt2mjd(obmt::Float64)::Float64

Convert Gaia On-Board Mission Timeline (OBMT) to Modified Julian Date.

# Arguments
- `obmt`: On-board mission timeline in units of six-hour revolutions since launch

Based on the relationship defined in Gaia Data Release documentation:
https://gea.esac.esa.int/archive/documentation/GDR2/Introduction/chap_cu0int/cu0int_sec_release_framework/cu0int_ssec_time_coverage.html
"""
function obmt2mjd(obmt::Float64)
    # First convert to TCB Julian Year as in the Python version
    tcbjy = 2015.0 + (obmt - 1717.6256) / 1461.0

    # Convert Julian Year to Julian Date
    # Julian Year 2015.0 corresponds to JD 2457023.5
    jd_at_2015 = 2457023.75

    # 365.25 days per Julian year
    days_since_2015 = (tcbjy - 2015.0) * 365.25

    # Calculate Julian Date
    jd = jd_at_2015 + days_since_2015

    # Return Modified Julian Date
    return jd - 2400000.5
end

# ──────────────────────────────────────────────────────────────────────
# SPICE-based Earth ephemeris
# ──────────────────────────────────────────────────────────────────────

# Global variable to track if SPICE kernels are loaded
const _SPICE_KERNELS_LOADED = Ref(false)

"""
    _ensure_spice_kernels_loaded()

Ensure that SPICE kernels are loaded for Earth barycentric position calculations.
Downloads DE440 ephemeris and leap seconds kernel if not already present.
"""
function _ensure_spice_kernels_loaded()
    if _SPICE_KERNELS_LOADED[]
        return
    end

    # Get data directory from DataDeps
    data_dir = @datadep_str "DE440_Ephemeris"

    # Load leap seconds kernel
    lsk_path = joinpath(data_dir, "naif0012.tls")
    if isfile(lsk_path)
        furnsh(lsk_path)
    else
        error("Leap seconds kernel not found at $lsk_path")
    end

    # Load planetary ephemeris kernel
    spk_path = joinpath(data_dir, "de440.bsp")
    if isfile(spk_path)
        furnsh(spk_path)
    else
        error("DE440 ephemeris kernel not found at $spk_path")
    end

    _SPICE_KERNELS_LOADED[] = true
    @info "SPICE kernels loaded successfully for Earth barycentric position calculations"
end

"""
    geocentre_position_query(epoch_MJD)

Given a date+time in MJD format, return a named tuple of Earth position and velocity in AU
on that date. Uses SPICE.jl with JPL DE440 ephemeris data for offline calculations.

The positions and velocities represent the Geocenter of the Earth relative to the solar system
barycenter in the J2000 reference frame.
"""
function geocentre_position_query(epoch_MJD::Number)

    # Ensure SPICE kernels are loaded
    _ensure_spice_kernels_loaded()

    # Convert MJD to Julian Date
    jd = epoch_MJD + 2400000.5  # MJD to JD conversion

    # Convert JD to DateTime
    # JD 2451545.0 = January 1, 2000, 12:00:00 TT (J2000.0 epoch)
    # Use Dates.julian2datetime for conversion
    dt = Dates.julian2datetime(jd)

    # Convert to ephemeris time string for SPICE
    et = utc2et(string(dt))

    # Get Earth's state relative to Solar System Barycenter
    # 399 = Earth geocenter, 0 = Solar System Barycenter
    state, _ = spkez(399, et, "J2000", "NONE", 0)

    # Extract position and velocity, convert to AU and AU/day
    pos_km = state[1:3]    # Position in km
    vel_km_s = state[4:6]  # Velocity in km/s

    # Convert to AU (1 AU = 149,597,870.7 km)
    AU_KM = Octofitter.PlanetOrbits.au2m / 1e3
    pos_au = pos_km ./ AU_KM
    vel_au_day = vel_km_s .* 86400 ./ AU_KM  # km/s to AU/day

    return (; x=pos_au[1], y=pos_au[2], z=pos_au[3],
              vx=vel_au_day[1], vy=vel_au_day[2], vz=vel_au_day[3])
end

# ──────────────────────────────────────────────────────────────────────
# Gaia archive queries
# ──────────────────────────────────────────────────────────────────────

function _query_gaia_dr3(;gaia_id)
    fname = "_gaia_dr3_final/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiadr3.gaia_source WHERE source_id=$gaia_id"
            ],
            cookies=false,
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr3_final")
            mkdir("_gaia_dr3_final")
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end
    header_line, body_line = split(buf,"\n")
    headers = Symbol.(split(header_line,','))
    data = tryparse.(Float64, split(body_line,','))
    if length(data) < length(headers)
        error("Could not query DR3 for source")
    end
    return namedtuple(headers, data)
end

function _query_gaia_dr2(;gaia_id)
    fname = "_gaia_dr2/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiadr2.gaia_source WHERE source_id=$gaia_id"
            ],
            cookies=false
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr2")
            mkdir("_gaia_dr2")
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end
    header_line, body_line = split(buf,"\n")
    headers = Symbol.(split(header_line,','))
    data = tryparse.(Float64, split(body_line,','))
    if length(data) <= 1
        error("Could not query DR2 for source")
    end
    return namedtuple(headers, data)
end

function _query_gaia_dr1(;gaia_id)
    fname = "_gaia_dr1/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiadr1.gaia_source WHERE source_id=$gaia_id"
            ],
            cookies=false
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr1")
            mkdir("_gaia_dr1")
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end
    header_line, body_line = split(buf,"\n")
    headers = Symbol.(split(header_line,','))
    data = tryparse.(Float64, split(body_line,','))
    if length(data) <= 1
        error("could not find source")
    end
    return namedtuple(headers, data)
end

# ──────────────────────────────────────────────────────────────────────
# Astrometric model fitting
# ──────────────────────────────────────────────────────────────────────

function prepare_A_4param(
    table,
    reference_epoch_mjd_ra,
    reference_epoch_mjd_dec,
)
    n_obs = size(table, 1)

    A = zeros(n_obs, 4)
    for i in 1:n_obs
        # Position terms
        A[i, 1] = table.cosϕ[i]  # α
        A[i, 2] = table.sinϕ[i]  # δ

        # Proper motion terms
        A[i, 3] = table.cosϕ[i] * (table.epoch[i] - reference_epoch_mjd_ra)/ julian_year  # μα*
        A[i, 4] = table.sinϕ[i] * (table.epoch[i] - reference_epoch_mjd_dec)/ julian_year # μδ
    end

    return A
end


function prepare_A_5param(
    table,
    reference_epoch_mjd_ra,
    reference_epoch_mjd_dec,
)
    n_obs = size(table, 1)

    A = zeros(n_obs, 5)
    for i in 1:n_obs
        # Position terms
        A[i, 1] = table.cosϕ[i]  # α
        A[i, 2] = table.sinϕ[i]  # δ

        # Parallax term
        A[i, 3] = -table.parallaxFactorAlongScan[i]

        # Proper motion terms
        A[i, 4] = table.cosϕ[i] * (table.epoch[i] - reference_epoch_mjd_ra)/ julian_year  # μα*
        A[i, 5] = table.sinϕ[i] * (table.epoch[i] - reference_epoch_mjd_dec)/ julian_year # μδ
    end

    return A
end


function fit_4param_prepared(
    A_factored,
    table,
    Δα_mas,
    Δδ_mas,
    σ_formal=0.0;
)
    n_obs = size(table, 1)

    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))

    b = zeros(T, n_obs)
    x = zeros(T, 4)

    for i in 1:n_obs
        b[i] = Δα_mas[i] * table.cosϕ[i] + Δδ_mas[i] * table.sinϕ[i]
    end

    if σ_formal != 0.
        @. b *= 1/σ_formal
    end

    x = A_factored \ b

    parameters = @SVector [x[1], x[2], x[3], x[4]]

    return (; parameters)
end

function fit_5param_prepared(
    A_prepared,
    table,
    Δα_mas,
    Δδ_mas,
    residuals=0.0,
    σ_formal=0.0;
    include_chi2=Val(false),
)
    n_obs = size(table, 1)

    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))

    # For a *scalar* weight the explicit σ_formal scaling is redundant: the
    # LSQ x = (A/σ) \ (b/σ) is identical to x = A \ b — the σ cancels — so we
    # solve unweighted and fold 1/σ² into the chi² at the end.  For a
    # *per-epoch* σ vector (e.g. Hipparcos `sres`) the weights do NOT cancel
    # and we perform the genuinely weighted solve.  Either way we solve
    # against a Bumper-allocated copy of A_prepared (some `\` overloads
    # require an owned matrix rather than a plain view).
    @no_escape begin
        b = @alloc(T, n_obs)
        A_dense = @alloc(T, n_obs, size(A_prepared, 2))

        @. b = Δα_mas * table.cosϕ + Δδ_mas * table.sinϕ + residuals
        if σ_formal isa Number
            @. A_dense = A_prepared
        else
            @. A_dense = A_prepared / σ_formal
            @. b = b / σ_formal
        end

        x = A_dense \ b

        parameters = @SVector [x[1], x[2], x[4], x[5], x[3]]

        if include_chi2 == Val(true)
            if σ_formal == 0
                error("Asked for `include_chi2=true` but `σ_formal==0`")
            end
            model_predictions = @alloc(T, n_obs)
            residuals_buf = @alloc(T, n_obs)
            mul!(model_predictions, A_dense, x)
            residuals_buf .= b .- model_predictions

            # Scalar σ: residuals are unweighted, apply 1/σ² here.
            # Vector σ: residuals are already whitened by the weighted solve.
            chi_squared_astro = σ_formal isa Number ?
                dot(residuals_buf, residuals_buf) / (σ_formal * σ_formal) :
                dot(residuals_buf, residuals_buf)

            n_parameters = 5
            dof = n_obs - n_parameters

            chi2_reduced = chi_squared_astro / dof
        end
    end

    if include_chi2 != Val(true)
        return (; parameters)
    end
    return (;
        parameters,
        chi_squared_astro,
        chi2_reduced,
        dof
    )

end

"""
    fit_5param_pinv(pinv_A, table, Δα_mas, Δδ_mas, residuals) -> nt

Same as `fit_5param_prepared(...; include_chi2=Val(false))` but uses a
pre-computed pseudo-inverse `pinv_A` (5×N) for the LSQ solve.  Returns
just the SVector of parameters in the same order as fit_5param_prepared.
"""
function fit_5param_pinv(pinv_A, table, Δα_mas, Δδ_mas, residuals=0.0)
    n_obs = size(pinv_A, 2)
    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))
    @no_escape begin
        b = @alloc(T, n_obs)
        @. b = Δα_mas * table.cosϕ + Δδ_mas * table.sinϕ + residuals
        x_buf = @alloc(T, 5)
        mul!(x_buf, pinv_A, b)
        parameters = @SVector [x_buf[1], x_buf[2], x_buf[4], x_buf[5], x_buf[3]]
    end
    return (; parameters)
end

# ──────────────────────────────────────────────────────────────────────
# Sky path perturbation simulation
# ──────────────────────────────────────────────────────────────────────

# Hipparcos main grid step (Lindegren 1997, ESA SP-1200 vol 3). Used as the
# argument `s` to `_simulate_skypath_hippacentre_combined!` for the Hipparcos
# abscissa branch.
const HIPPARCOS_GRID_STEP_ARCSEC = 1.2074

# Per-transit resolution gate for the BINARYS modulated photocentre.  At
# projected separations ρ ≫ s, the Hipparcos pipeline either resolved the
# components into separate Tycho/Hipparcos entries or rejected the row, so
# the catalog point reflects the primary alone; the BINARYS atan2 modulation
# must therefore taper to zero in this limit.  The atan2 formula does not do
# this on its own — `(s/2π)·φ` keeps oscillating ±s/2 with separation.  We
# gate the per-companion modulated-signal flux ratio by a Gaussian taper in
# the *full* projected separation (not the along-scan projection ρ_pk, since
# the resolution criterion is geometric, not abscissa-projected).  The
# resolution scale is anchored to s = 1.207″ (Lindegren 1997 §3.3, ESA
# SP-1200; van Leeuwen 2007 §6.4); refine empirically as needed.
const HIPPARCOS_RESOLUTION_ARCSEC = 1.207
α_resolve_hip(ρ_arcsec) = exp(-(ρ_arcsec / HIPPARCOS_RESOLUTION_ARCSEC)^2)

"""
Given scan epochs and angles, and an orbit describing perturbations from a
(possibly luminous) companion, accumulate the astrometric perturbations to
the sky path using the linear photocentre formula `(host + f·planet)/(1+f)`.

This is the Gaia / small-separation form. For the Hipparcos abscissa branch
(BINARYS atan2 Hippacentre) use `_simulate_skypath_hippacentre_combined!`
instead — that routine must be called once with all companions because the
modulated-signal Hippacentre is not the sum of per-companion Hippacentres.
"""
function _simulate_skypath_perturbations!(
    Δα_model, Δδ_model,
    table,
    orbit::AbstractOrbit,
    planet_mass_msol,
    flux_ratio,
    orbit_solutions, orbit_solutions_i_epoch_start, T=Float64,
)
    # The photocentre offset is `raoff(sol) · coeff` (and `decoff(sol) · coeff`)
    # for the *same* `coeff` per call — derive it from the standard
    # `(host_reflex + planet_bary · f) / (1 + f)` formula and hoist out of the
    # loop. raoff(sol, m) = -m/M·raoff(sol) where M = totalmass(orbit). Use
    # the `totalmass(orbit)` API so this works for Visual/AbsoluteVisual
    # wrappers (no .M field) as well as bare KepOrbit.
    M_tot = PlanetOrbits.totalmass(orbit)
    m_host_eff = M_tot - planet_mass_msol
    coeff = (-planet_mass_msol + flux_ratio * m_host_eff) / (M_tot * (1 + flux_ratio))
    if orbit_solutions_i_epoch_start >= 0
        @inbounds for i in eachindex(table.epoch)
            sol = orbit_solutions[orbit_solutions_i_epoch_start+i]
            Δα_model[i] += raoff(sol) * coeff
            Δδ_model[i] += decoff(sol) * coeff
        end
    else
        @inbounds for i in eachindex(table.epoch)
            sol = orbitsolve(orbit, table.epoch[i])
            Δα_model[i] += raoff(sol) * coeff
            Δδ_model[i] += decoff(sol) * coeff
        end
    end
    return
end

"""
    _simulate_skypath_hippacentre_combined!(
        Δα_model, Δδ_model, σ_inflation,
        table, orbits, planet_masses_msol, flux_ratios,
        orbit_solutions_per_planet, orbit_solutions_i_epoch_starts, T, s,
    )

Compute the BINARYS Hippacentre along-scan offset Δν_B from the *combined*
multi-companion modulated signal (Leclerc et al. 2023, A&A 672 A82, Eq. 13 +
Eq. 15) and accumulate it into (Δα_model, Δδ_model) such that the downstream
scan projection `b = Δα·cosϕ + Δδ·sinϕ` recovers Δν_B exactly. Cross-scan
components are zero (unobservable from a Hipparcos abscissa).

For N companions (k = 1..N) at scan-projected separations ρ_p^(k) from the
host, with Hp-band flux ratios f_k = L_k/L_host **gated by a per-transit
resolution taper** α_k = α_resolve_hip(|ρ^(k)|), so that f_k_eff = f_k · α_k,
the combined Hippacentre phase in the *host frame* is

    φ = atan2( Σ_k f_k_eff sin ζ_k ,  1 + Σ_k f_k_eff cos ζ_k ),     ζ_k = 2π ρ_p^(k)/s

so the Hippacentre offset from the system barycentre, projected on scan, is

    Δν_B = (s/2π) · φ + Σ_k host_along^(k)

where host_along^(k) is the host-vs-barycentre offset *due to companion k*
projected on scan (encodes the per-companion mass fraction via
`raoff(sol_k, m_k)`). Summing per-companion host-reflex contributions matches
Octofitter's existing 2-body convention for `raoff(sol, m)`.  The host-reflex
sum is *not* gated by α — the host's barycentric motion is physical
regardless of whether Hipparcos resolved the binary; only the modulated
photocentre contribution tapers off in the resolved limit.

The σ-inflation factor is the *combined* first-harmonic amplitude reduction
(Leclerc et al. Eq. 15, generalised to N companions, using f_k_eff):

    f_σ = (1 + Σ_k f_k_eff) / sqrt( (1 + Σ_k f_k_eff cos ζ_k)² + (Σ_k f_k_eff sin ζ_k)² )

Pass `σ_inflation` as a buffer initialised to 1; it is multiplied in place by
f_σ per transit. For all-dark companions (every f_k = 0) this routine reduces
to Δν_B = Σ_k host_along^(k) and f_σ = 1, identical to the linear photocentre.
For the wide-separation limit (every α_k → 0) it again reduces to
Δν_B = Σ_k host_along^(k) and f_σ = 1 — the "resolved binary, primary alone"
answer the Hipparcos pipeline would give.

NOTE: The σ_inflation buffer is the noise-model correction for the *likelihood*
comparison of the predicted abscissa to observed IAD residuals. It must NOT be
folded into the weighting of any LSQ that reproduces the published catalog
5-parameter solution, because the catalog fit was performed by the Hipparcos
pipeline with point-source σ.
"""
function _simulate_skypath_hippacentre_combined!(
    Δα_model, Δδ_model, σ_inflation,
    table,
    orbits,                          # iterable of <:AbstractOrbit
    planet_masses_msol,              # iterable of mass (Msol)
    flux_ratios,                     # iterable of Hp-band L_k/L_host
    orbit_solutions_per_planet,      # iterable of orbit_solutions arrays
    orbit_solutions_i_epoch_starts,  # iterable of Int
    T,
    s::Float64,
)
    n_planets = length(orbits)
    # Early-exit when no companion contributes (n_planets=0 has prior P=0.5,
    # and any inactive companion has mass=0 from the @variables block).
    # In that case the per-epoch loop reduces to: Re=1, Im=0, atan(0,1)=0,
    # Δν_B=0 — so Δα/Δδ_model and σ_inflation are unchanged from their
    # caller-initialised values (0 and 1 respectively).  Skip the
    # transit loop entirely.
    any_active = false
    @inbounds for k in 1:n_planets
        if planet_masses_msol[k] != zero(eltype(planet_masses_msol))
            any_active = true
            break
        end
    end
    any_active || return

    inv_2π_over_s = s / (2π)
    two_π_over_s = (2π) / s
    # Squared resolution scale in mas² — saves a sqrt per transit per active
    # companion. α_resolve_hip(ρ) = exp(-(ρ/RES_arcsec)²) and ρ_arcsec =
    # √(ra²+dec²)/1000, so (ρ_arcsec/RES_arcsec)² = (ra²+dec²)/(1000·RES_arcsec)².
    inv_res_mas2 = 1 / (1000 * HIPPARCOS_RESOLUTION_ARCSEC)^2
    for i in eachindex(table.epoch)
        cosϕ = table.cosϕ[i]
        sinϕ = table.sinϕ[i]

        # Combined modulated signal in the host frame:
        #   1 (host, by convention) + Σ_k f_k · exp(i·ζ_k)
        # Initial (Re, Im) is the host (k=0) contribution: 1 + 0i; f_total starts at 0.
        Re = one(T)
        Im = zero(T)
        f_total = zero(T)
        host_along_total = zero(T)

        for k in 1:n_planets
            # Absent companions (mass == 0 → host reflex 0; the calling site
            # also forces flux_ratios[k] = 0 → modulated signal 0). Mirrors
            # the DR2/DR3 _simulate_skypath_perturbations! short-circuit.
            if planet_masses_msol[k] == zero(eltype(planet_masses_msol))
                continue
            end

            i_start_k = orbit_solutions_i_epoch_starts[k]
            sol_k = i_start_k >= 0 ?
                orbit_solutions_per_planet[k][i_start_k + i] :
                orbitsolve(orbits[k], table.epoch[i])

            # Host reflex from companion k (projected) = −B_k · ρ_p^(k)
            ra_h  = raoff(sol_k, planet_masses_msol[k])
            dec_h = decoff(sol_k, planet_masses_msol[k])
            host_along_total += ra_h * cosϕ + dec_h * sinϕ

            # Companion-k separation from host, projected on scan
            ra_p  = raoff(sol_k)
            dec_p = decoff(sol_k)
            ρ_pk = ra_p * cosϕ + dec_p * sinϕ

            # Per-transit resolution gate: the BINARYS modulated-signal
            # contribution from companion k tapers to zero at ρ ≫ s, because
            # at those separations the Hipparcos pipeline detected the pair
            # as resolved (or rejected the row) and the catalog point is the
            # primary alone.  Gate is applied to the modulated-signal f_k
            # only — the host-reflex sum below is barycentric and remains
            # physical regardless of resolution state.  In the resolved
            # limit (α → 0) Re → 1, Im → 0, φ → 0, f_σ → 1, so Δν_B reduces
            # to host_along_total exactly and no spurious σ inflation
            # remains.
            ρ_full_sq_mas2 = ra_p*ra_p + dec_p*dec_p
            α_k = exp(-ρ_full_sq_mas2 * inv_res_mas2)

            ζ_k = two_π_over_s * ρ_pk
            f_k = flux_ratios[k] * α_k
            sin_ζk, cos_ζk = sincos(ζ_k)
            Re += f_k * cos_ζk
            Im += f_k * sin_ζk
            f_total += f_k
        end

        # Combined Hippacentre phase (host frame) and barycentre-frame offset
        φ = atan(Im, Re)
        Δν_B = inv_2π_over_s * φ + host_along_total

        Δα_model[i] += Δν_B * cosϕ
        Δδ_model[i] += Δν_B * sinϕ

        if σ_inflation !== nothing
            amp = sqrt(Re*Re + Im*Im)
            σ_inflation[i] *= (1 + f_total) / amp
        end
    end
    return
end


# ──────────────────────────────────────────────────────────────────────
# GOST scan forecast
# ──────────────────────────────────────────────────────────────────────

"""
    forecast_table = GOST_forecast(ra_deg,dec_deg;baseline=:dr3)

Given an Ra and Dec position, retreive a forecast of Gaia observations from the GOST tool automatically.
See tool URL here: https://gaia.esac.esa.int/gost/

Please be aware that others  might be able to discover the target coordinates you searched for
(though not who performed the search) via information leaked to the external service.

Baseline can be :dr3, :dr4, or :dr5.
"""
function GOST_forecast(ra_deg,dec_deg;baseline=:dr3)
    if baseline == :dr3
        to = "2017-06-28T00:00:00"
    elseif baseline == :dr4
        to = "2020-01-20T00:00:00"
    elseif baseline == :dr5
        to = "2025-01-15T06:16:00"
    end

    if haskey(ENV, "OCTO_GOST_CATALOG") && !isempty(ENV["OCTO_GOST_CATALOG"])
        fname = ENV["OCTO_GOST_CATALOG"]
        @info "Using provided Gaia scan forecast database $fname"
        forecast_table = CSV.read(fname, Table, normalizenames=true)
        themin, idx = findmin(hypot.(
            (forecast_table.ra_rad_ .- deg2rad(ra_deg)) .*60 .*60 .*1000 .* cos(dec_deg),
            (forecast_table.dec_rad_ .- deg2rad(dec_deg)).*60 .*60 .*1000
        ))
        if themin > 500
            error("Could not find this target within the provided Gaia scan forecast database file set through OCTO_GOST_CATALOG=$fname Closest target: $themin [mas]")
        end
        ra_rad = forecast_table.ra_rad_[idx]
        dec_rad = forecast_table.dec_rad_[idx]
        mask = isapprox.(forecast_table.ra_rad_, ra_rad) .& isapprox.(forecast_table.dec_rad_, dec_rad)
        @info "Found forecasted visibility windows" windows=count(mask)
        if isempty(mask)
            error("Invalid condition: no visibility windows.")
        end
        return forecast_table[mask,:]
    end

    fname = "GOST-$ra_deg-$dec_deg-$baseline.csv"
    if isfile(fname)
        @info "Using cached Gaia scan forecast $fname"
        forecast_table = CSV.read(fname, Table, normalizenames=true)
        return forecast_table
    end

    # Just pick a cookie ID to use.
    # Might be better to let the service create one for us.
    cookiejar = HTTP.CookieJar()
    @info "Contacting the GAIA scan forecast tool GOST: https://gaia.esac.esa.int/gost/"
    resp0 = HTTP.get(
        "https://gaia.esac.esa.int/gost/",
        cookiejar=cookiejar,
    )
    if resp0.status != 200
        println(String(resp0.body))
        error("Could not contact the GAIA scan forecast tool GOST https://gaia.esac.esa.int See above error message.")
    end
    formdata = Dict([
        "serviceCode"=>"1",
        "inputmode"=>"single",
        "srcname"=>"009",
        "srcra" => string(round(ra_deg,digits=7)),
        "srcdec" => string(round(dec_deg,digits=7)),
        "from" => "2014-07-25T10:31:26",
        "to" => to,
    ])


    @info "Retrieving forecasted GAIA scans from GOST: https://gaia.esac.esa.int/gost/"
    resp = HTTP.post(
        "https://gaia.esac.esa.int/gost/GostServlet",
        body=HTTP.Form(formdata),
        cookiejar=cookiejar
    )
    if resp.status != 200 || contains(String(collect(resp.body)),"error")
        println(String(resp.body))
        error("Could not fetch GAIA scan forecast from GOST. See above error message. Do you have an internet connection available?")
    end

    m = match(r"Submitted with id (\d+)", String(resp.body))
    response_id = m.captures[1]
    session_id = cookiejar.entries["gaia.esac.esa.int"]["gaia.esac.esa.int;/gost;JSESSIONID"].value
    url = "https://gaia.esac.esa.int/gost/export.jsp?id=$session_id/$response_id&format=csv"
    resp_dat = HTTP.get(
        url,
        cookiejar=cookiejar
    )
    @info "done"
    body = collect(resp_dat.body)

    # Save for offline use eg on clusters
    write(fname, body)

    io = IOBuffer(body)
    if bytesavailable(io) == 0
        error("Empty response from GOST service. Rate limited?")
    end
    forecast_table = CSV.read(io, Table, normalizenames=true)

    return forecast_table
end


# ──────────────────────────────────────────────────────────────────────
# GaiaCatalogFitObs struct (used by HGCAObs)
# ──────────────────────────────────────────────────────────────────────

struct GaiaCatalogFitObs{TTable,TCat,TDist,TFact} <: AbstractObs
    table::TTable
    source_id::Int
    gaia_sol::TCat
    dist::TDist
    A_prepared_4::TFact
    A_prepared_5::TFact
end
const GaiaCatalogFitLikelihood = GaiaCatalogFitObs

function GaiaCatalogFitObs(;
    gaia_id_dr2=nothing,
    gaia_id_dr3=nothing,
    scanlaw_table=nothing,
    ref_epoch_ra=nothing,
    ref_epoch_dec=nothing
)
    if !isnothing(gaia_id_dr2)
        source_id = gaia_id_dr2
        gaia_sol = Octofitter._query_gaia_dr2(; gaia_id=gaia_id_dr2)
        if isnothing(ref_epoch_ra)
            ref_epoch_ra = meta_gaia_DR2.ref_epoch_mjd
        end
        if isnothing(ref_epoch_dec)
            ref_epoch_dec = meta_gaia_DR2.ref_epoch_mjd
        end
    elseif !isnothing(gaia_id_dr3)
        source_id = gaia_id_dr3
        gaia_sol = Octofitter._query_gaia_dr3(; gaia_id=gaia_id_dr3)
        if isnothing(ref_epoch_ra)
            ref_epoch_ra = meta_gaia_DR3.ref_epoch_mjd
        end
        if isnothing(ref_epoch_dec)
            ref_epoch_dec = meta_gaia_DR3.ref_epoch_mjd
        end
    else
        throw(ArgumentError("Please provide at least one of `gaia_id_dr2` or `gaia_id_dr3`"))
    end

    ra_deg = gaia_sol.ra
    dec_deg = gaia_sol.dec
    μ = [
        gaia_sol.parallax,
        gaia_sol.ra,
        gaia_sol.dec,
        gaia_sol.pmra,
        gaia_sol.pmdec,
    ]
    σ = [
        gaia_sol.parallax_error,
        gaia_sol.ra_error / 60 / 60 / 1000 / cosd(gaia_sol.dec),
        gaia_sol.dec_error / 60 / 60 / 1000,
        gaia_sol.pmra_error,
        gaia_sol.pmdec_error,
    ]
    C = [
        1 gaia_sol.ra_parallax_corr gaia_sol.dec_parallax_corr gaia_sol.parallax_pmra_corr gaia_sol.parallax_pmdec_corr
        gaia_sol.ra_parallax_corr 1 gaia_sol.ra_dec_corr gaia_sol.ra_pmra_corr gaia_sol.ra_pmdec_corr
        gaia_sol.dec_parallax_corr gaia_sol.ra_dec_corr 1 gaia_sol.dec_pmra_corr gaia_sol.dec_pmdec_corr
        gaia_sol.parallax_pmra_corr gaia_sol.ra_pmra_corr gaia_sol.dec_pmra_corr 1 gaia_sol.pmra_pmdec_corr
        gaia_sol.parallax_pmdec_corr gaia_sol.ra_pmdec_corr gaia_sol.dec_pmdec_corr gaia_sol.pmra_pmdec_corr 1
    ]
    Σ = Diagonal(σ) * C * Diagonal(σ)
    dist = MvNormal(μ, Hermitian(Σ))

    if isnothing(scanlaw_table)
        @info "No scan law table provided. We will fetch an approximate solution from the GOST webservice."
        forecast_table = FlexTable(GOST_forecast(ra_deg, dec_deg))
        forecast_table.epoch = jd2mjd.(forecast_table.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_)
        forecast_table.scanAngle_rad = forecast_table.scanAngle_rad_
    else
        @info "Scanlaw table was provided, will not query GOST."
        forecast_table = FlexTable(scanlaw_table)
        forecast_table.epoch = tcb_at_gaia_2mjd.(forecast_table.times)
        forecast_table.scanAngle_rad = deg2rad.(forecast_table.angles)
    end

    forecast_table.cosϕ = cos.(π / 2 .+ forecast_table.scanAngle_rad)
    forecast_table.sinϕ = sin.(π / 2 .+ forecast_table.scanAngle_rad)

    earth_pos_vel = geocentre_position_query.(forecast_table.epoch)

    table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)

    gaps_dr2 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiadr2_08252020.csv"), FlexTable)
    gaps_edr23 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiaedr3_12232020.csv"), FlexTable)
    gaps = Table(
        start_mjd=obmt2mjd.(vcat(gaps_dr2.start,gaps_edr23.start)),
        stop_mjd=obmt2mjd.(vcat(gaps_dr2.end,gaps_edr23.end)),
        note=[gaps_dr2.comment; gaps_edr23.description]
    )
    table = filter(eachrow(table)) do row
        row = row[]
        for gap in eachrow(gaps)
            gap = gap[]
            if gap.start_mjd <= row.epoch <= gap.stop_mjd
                @info "Detected known gap in Gaia scans; skipping." window=row.epoch note=gap.note
                return false
            end
        end
        return true
    end
    table = Table(map(dat->dat[], table))

    A_prepared_4 = prepare_A_4param(table, ref_epoch_ra, ref_epoch_dec)
    A_prepared_5 = prepare_A_5param(table, ref_epoch_ra, ref_epoch_dec)

    return GaiaCatalogFitObs(
        table,
        source_id,
        gaia_sol,
        dist,
        A_prepared_4,
        A_prepared_5,
    )
end
