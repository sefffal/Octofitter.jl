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
# AGIS astrometric input spans
# ──────────────────────────────────────────────────────────────────────
#
# Time spans of the *input data actually used* by each release's astrometric
# solution (AGIS), as published in the release astrometry papers. These differ
# from the nominal data-collection spans in `meta_gaia_DR2`/`meta_gaia_DR3`
# above: both the DR2 and (E)DR3 solutions excluded the first month of the
# operational phase (the ecliptic-pole scanning-law period, OBMT < 1192.13
# rev), and the boundaries fall mid-day, not at midnight UTC.
#   DR2:   OBMT 1192.13–3750.56 rev (Lindegren et al. 2018, Sect. 2)
#   (E)DR3: OBMT 1192.13–5230.09 rev (Lindegren et al. 2021, Sect. 2)
# Transits forecast outside these spans can never have contributed to the
# corresponding astrometric solution, regardless of the matched-transit
# counts reported in the source catalogs.
const gaia_agis_span_dr2 = (;
    start_mjd = obmt2mjd(1192.13), # 2014-08-22
    stop_mjd  = obmt2mjd(3750.56), # 2016-05-23 ~11:35 UTC
)
const gaia_agis_span_dr3 = (;
    start_mjd = obmt2mjd(1192.13), # 2014-08-22
    stop_mjd  = obmt2mjd(5230.09), # 2017-05-28 ~08:45 UTC
)

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

"""
    gaia_plx(; gaia_id)

A truncated Normal prior on parallax [mas] read from that source's Gaia DR3
astrometric solution, for use as `plx ~ gaia_plx(gaia_id=…)`.

The truncation at ±10σ is deliberate and inherited: it keeps a sampler from
walking into the negative-parallax tail, where the AU→mas conversion is
undefined.

!!! note "Source changed in v9"
    v8 read this out of the HGCA `HGCA_vEDR3.fits` data dependency, since the
    only caller was `HGCAObs`. The HGCA modelling stack is retired
    ([`HGCAObs`](@ref) is now a helper over [`G23HObs`](@ref)), so this reads
    the DR3 catalog directly instead. For a source in both, `parallax_gaia`
    in the HGCA *is* the DR3 parallax, so the numbers agree; the difference
    is that this no longer downloads the 30 MB HGCA catalog to read one row,
    and it works for sources the HGCA does not contain.
"""
function gaia_plx(; gaia_id)
    dr3 = _query_gaia_dr3(; gaia_id)
    μ = Float64(dr3.parallax)
    σ = Float64(dr3.parallax_error)
    return truncated(Normal(μ, σ), lower=μ - 10σ, upper=μ + 10σ)
end
export gaia_plx

"""
    gaia_dr3_solution(; gaia_id) -> NamedTuple

That source's row of the Gaia DR3 `gaia_source` catalog, as a NamedTuple keyed
by the catalog's own column names (`parallax`, `pmra`, `pmdec`, `ra`, `dec`,
`ref_epoch`, `phot_g_mean_mag`, `ruwe`, …).

Queried from the ESA TAP service and cached at
`_gaia_dr3_final/source-<gaia_id>.csv` in the working directory, so a repeated
call — or a re-run on a cluster node with the directory copied across — is
offline.

Useful for seeding an absolute-frame system block from the published solution:

```julia
cat = gaia_dr3_solution(gaia_id=756291174721509376)
sys = System(name="s", bodies=(A, b), observations=(pma,), variables=@variables begin
    plx   ~ truncated(Normal(cat.parallax, cat.parallax_error), lower=0)
    ra    = \$(cat.ra)
    dec   = \$(cat.dec)
    pmra  ~ Normal(cat.pmra, 10)
    pmdec ~ Normal(cat.pmdec, 10)
    rv    = 0.0
    ref_epoch = Octofitter.jd2mjd(2457388.5)
end)
```

(Note the asymmetry: `\$` interpolation is needed on `=` lines, whose
right-hand sides are quoted and evaluated later inside the model, and is
rejected on `~` lines, which already see the enclosing scope.)

!!! note "v8 → v9"
    `GaiaDR4AstromObs` used to carry the published solution as `obs.gaia_sol`,
    and took a `gaia_id=`. It no longer does either — the observation models a
    sky path and no longer needs the catalog row — so this is the supported way
    to get at it. It is the same function [`gaia_plx`](@ref) reads.
"""
gaia_dr3_solution(; gaia_id) = _query_gaia_dr3(; gaia_id)
export gaia_dr3_solution

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
            (forecast_table.ra_rad_ .- deg2rad(ra_deg)) .*60 .*60 .*1000 .* cosd(dec_deg),
            (forecast_table.dec_rad_ .- deg2rad(dec_deg)).*60 .*60 .*1000
        ))
        if themin > 500
            error("Could not find this target within the provided Gaia scan forecast database file set through OCTO_GOST_CATALOG=$fname Closest target: $themin [mas]")
        end
        ra_rad = forecast_table.ra_rad_[idx]
        dec_rad = forecast_table.dec_rad_[idx]
        mask = isapprox.(forecast_table.ra_rad_, ra_rad) .& isapprox.(forecast_table.dec_rad_, dec_rad)
        @info "Found forecasted visibility windows" windows=count(mask)
        if !any(mask)
            error("Invalid condition: no visibility windows.")
        end
        return _sort_dedup_gost(forecast_table[mask,:])
    end

    fname = "GOST-$ra_deg-$dec_deg-$baseline.csv"
    if isfile(fname)
        @info "Using cached Gaia scan forecast $fname"
        forecast_table = CSV.read(fname, Table, normalizenames=true)
        _check_gost_nonempty(forecast_table, ra_deg, dec_deg, baseline,
            "cached file $(abspath(fname))"; cached=true)
        return _sort_dedup_gost(forecast_table)
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

    io = IOBuffer(body)
    if bytesavailable(io) == 0
        error("Empty response from GOST service. Rate limited?")
    end
    forecast_table = CSV.read(io, Table, normalizenames=true)

    # Validate BEFORE writing the cache. GOST has been observed to return a
    # header-only body on a transient failure — with a 200 status and no error
    # string, so nothing above catches it. Caching that poisons every
    # subsequent run for this target: the file exists, so the branch above
    # reads it and the zero-row table fails far downstream in
    # `_g23h_forecast_pool` as a `DimensionMismatch`, with nothing pointing at
    # the cache. Write the file only once it is known to be worth keeping.
    _check_gost_nonempty(forecast_table, ra_deg, dec_deg, baseline,
        "the GOST service (https://gaia.esac.esa.int/gost/)")

    # Save for offline use eg on clusters
    write(fname, body)

    return _sort_dedup_gost(forecast_table)
end

"""
    _check_gost_nonempty(tbl, ra_deg, dec_deg, baseline, source; cached=false)

Fail with a message naming the target and the source when a GOST forecast has
no scans, rather than letting a zero-row table propagate into the likelihood.
"""
function _check_gost_nonempty(tbl, ra_deg, dec_deg, baseline, source; cached=false)
    n = hasproperty(tbl, :ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_) ?
        length(vec(collect(tbl.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_))) : 0
    n > 0 && return nothing
    hint = cached ?
        "Delete that file and retry; it was cached from a failed request." :
        "This is usually transient — retry. If it persists, check the coordinates."
    error("GOST returned no scans for ra=$ra_deg, dec=$dec_deg (baseline=$baseline) " *
          "from $source. $hint")
end

# Defensively sort and deduplicate a GOST forecast table.  User-supplied
# catalogs have been observed to contain a target's entire forecast block
# duplicated (even triplicated) verbatim, and everything downstream (window
# slicing by findfirst/findlast, first/last epoch selection, missed-transit
# accounting) assumes a time-sorted table of UNIQUE scans.  Real field-of-view
# transits are ≥ 1.7 h apart (Gaia's 6 h spin, two FoVs), so rows closer than
# ~10 s in time are duplicates of the same scan.
function _sort_dedup_gost(tbl)
    # vec(): table slices like `forecast_table[mask,:]` carry n×1 Matrix columns.
    times = vec(collect(tbl.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_))
    order = sortperm(times)
    keep = Int[]
    sizehint!(keep, length(order))
    last_t = -Inf
    for i in order
        if times[i] - last_t > 1e-4  # days; ≈ 8.6 s
            push!(keep, i)
            last_t = times[i]
        end
    end
    if length(keep) == length(times) && issorted(times)
        return tbl  # already clean: return unchanged
    end
    @warn "GOST forecast contained duplicate and/or unsorted scan rows; sorted and deduplicated. Downstream results assume unique time-ordered scans." n_raw=length(times) n_unique=length(keep)
    names = propertynames(tbl)
    out = Table(NamedTuple{names}(map(p -> vec(getproperty(tbl, p))[keep], names)))
    @assert issorted(vec(out.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_))
    return out
end
