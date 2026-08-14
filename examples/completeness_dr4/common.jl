#=
Gaia DR4 completeness study — everything the entry points share.

Included by `setup.jl`, `run_local.jl`, `completeness_trial.jl` and
`assemble_results.jl`, so the grid, the noise budget, the template model and
the detection criterion are defined exactly once. (In the previous version of
this example the target setup and the noise numbers were copy-pasted between
`setup.jl` and `completeness_trial.jl`, which is precisely how the two came to
disagree.)

Three pieces:

  1. `load_target()` — the star: its DR3 astrometry, its *measured* along-scan
     noise budget, and its forecast DR4 transits. Cached to CSV beside these
     scripts, so a compute node needs neither the network nor the ~14 GB G23H
     catalog.
  2. the completeness grid, sampler, `inject` and detection criterion.
  3. `build_system(target)` — the template `System` that every trial injects
     into.

See `docs/src/gaia-dr4-simulation.md` for the tutorial version of the
simulation, and `docs/src/completeness.md` for the completeness framework.
=#

using Octofitter
using Distributions
using Statistics
using CSV

# ──────────────────────────────────────────────────────────────────────
# 1. The target star
# ──────────────────────────────────────────────────────────────────────

# A bright (G = 6.94) nearby star with a clean DR3 solution. The scan law and
# the noise budget both come from this choice: Gaia's sensitivity is far from
# uniform, so a completeness map is a statement about *this* star, not about
# Gaia in general.
const GAIA_ID = 5064625130502952704

# The columns of this source's G23H catalog row that `g23h_scan_uncertainty`
# reads, recorded so the example runs without the ~14 GB `G23H_Catalog`
# DataDep. These are the real catalog values (verified against the catalog on
# 2026-08-14); passing `catalog=nothing` below reads the catalog itself and
# reproduces them, at the cost of a ~70 s scan of the Arrow file.
const CATALOG_STANDIN = (;
    gaia_source_id                   = GAIA_ID,
    phot_g_mean_mag_dr3              = 6.941057,
    sig_AL                           = 0.05813228279803717,
    sig_att_radec                    = 0.07765105456802447,
    sig_cal                          = 0.14901214069992602,
    astrometric_n_good_obs_al_dr3    = 882,
    astrometric_matched_transits_dr3 = 100,
)

target_csv(dir) = joinpath(dir, "dr4_target.csv")
transits_csv(dir) = joinpath(dir, "dr4_transits.csv")

"""
    load_target(dir=@__DIR__; refresh=false, catalog=CATALOG_STANDIN)

The star, its noise budget and its forecast DR4 transit table.

Returns a NamedTuple with the DR3 astrometry, every field of
[`g23h_scan_uncertainty`](@ref), and `transits`: a transit-level table in Gaia
DR4's own format (epochs in MJD, `scan_pos_angle` in **degrees**), with
`centroid_pos_al` left at zero for `generate_from_params` to fill in and
`centroid_pos_error_al` set to the measured per-transit scatter.

The first call queries the DR3 archive and GOST and writes `dr4_target.csv` +
`dr4_transits.csv` next to these scripts; later calls read those. Pass
`refresh=true` to re-query, or `catalog=nothing` to read the real G23H catalog
instead of the recorded row above.

Note that `GOST_forecast` caches its own response in the *working directory*,
so run `setup.jl` from a directory you will come back to.
"""
function load_target(dir::AbstractString=@__DIR__; refresh::Bool=false,
                     catalog=CATALOG_STANDIN)
    if !refresh && isfile(target_csv(dir)) && isfile(transits_csv(dir))
        info = NamedTuple(Table(CSV.File(target_csv(dir)))[1])
        transits = Table(CSV.File(transits_csv(dir)))
        return (; info..., transits)
    end

    # α, δ, ϖ from the DR3 archive (cached in a depot scratchspace).
    dr3 = gaia_dr3_solution(; gaia_id=GAIA_ID)

    # The measured along-scan noise budget — the same σ_AL/σ_att/σ_calib that
    # `G23HObs` builds its priors from for this source.
    σ = try
        g23h_scan_uncertainty(; gaia_id=GAIA_ID, catalog)
    catch err
        catalog === CATALOG_STANDIN && rethrow()
        @warn "Could not read the G23H catalog; using the recorded row for this source instead." exception = (err, catch_backtrace())
        g23h_scan_uncertainty(; gaia_id=GAIA_ID, catalog=CATALOG_STANDIN)
    end

    # The scan geometry, from GOST, in DR4's format and units. We simulate and
    # fit the same numbers, so the table quotes the *true* per-transit scatter
    # (`σ_transit_true`, calibration term included) rather than the formal
    # error a real DR4 table would carry. See the "What a real DR4 table would
    # quote" note in docs/src/gaia-dr4-simulation.md.
    transits = gaia_dr4_transit_template(;
        ra=dr3.ra, dec=dr3.dec, σ_al=σ.σ_transit_true, baseline=:dr4)

    info = (;
        gaia_id=GAIA_ID, ra=dr3.ra, dec=dr3.dec, plx=dr3.parallax,
        phot_g_mean_mag=σ.phot_g_mean_mag,
        σ.σ_AL, σ.σ_att, σ.σ_calib, σ.σ_formal, σ.n_ccd,
        σ.σ_transit_formal, σ.σ_transit_true,
        n_transits=length(transits.epoch),
    )

    mkpath(dir)
    CSV.write(target_csv(dir), [info])
    CSV.write(transits_csv(dir), transits)

    return (; info..., transits)
end


# ──────────────────────────────────────────────────────────────────────
# 2. The grid, the sampler, and the detection criterion
# ──────────────────────────────────────────────────────────────────────

# `OCTO_COMPLETENESS_QUICK=1` swaps in a 2 × 2 grid with one trial per cell, so
# the whole workflow runs locally in minutes. It is only a knob: the production
# grid below is what `submit.sh` dispatches, and the array range there must
# match `length(completeness_jobs(...))`.
#
# Note what QUICK does *not* shorten: the chains. On the identical trial (same
# seed, same simulated data) of a 30 MJup companion at 20 AU, 500 warmup + 500
# sampling steps returned a mass posterior with median 38 MJup and a 5th
# percentile of 3.9 — a confident "detection" — while 1000 + 1000 returned
# median 10 MJup with a 5th percentile of 0.1, i.e. essentially unconstrained,
# which is the honest answer for a 89-year period seen over 5.3 years. A
# completeness map made with chains too short to converge measures the sampler,
# not the data.
const QUICK = lowercase(get(ENV, "OCTO_COMPLETENESS_QUICK", "0")) ∉ ("0", "", "false", "no")

# Masses are in SOLAR masses throughout v9 — `mjup` is a plain multiplicative
# constant — and separations are semi-major axes in AU.
#
# The QUICK separations are deliberately 3 and 10 AU rather than the ends of
# the production range: a 1 AU orbit around a 1 M⊙ star has P = 1 yr and is
# nearly degenerate with parallax, and a 20 AU orbit has P = 89 yr against a
# 5.3 yr baseline, so both corners are non-detections for uninteresting reasons
# and a smoke test that lands on them tells you nothing.
const MASSES      = QUICK ? [3.0mjup, 30.0mjup] : 10 .^ range(-1, 2, length=12) .* mjup
const SEPARATIONS = QUICK ? [3.0, 10.0] : 10 .^ range(-0.3, 1.7, length=12)
const N_TRIALS    = QUICK ? 1 : 5

const ADAPTATION = 1000
const ITERATIONS = QUICK ? 1000 : 2000

"""
    sampler(model) -> Chains

HMC at a raised target acceptance: a single along-scan time series makes for a
bumpy posterior. For a production run consider `octofit_pigeons` instead — see
the note in docs/src/gaia-dr4-simulation.md on why parallel tempering suits
epoch astrometry — at a few times the cost per trial.
"""
sampler(model) = octofit(model, 0.85;
    adaptation=ADAPTATION, iterations=ITERATIONS, verbosity=0)

"""
    inject(mass, separation) -> NamedTuple

Grid point → parameter overrides. Overrides nest under `bodies`, matching the
model's own parameter structure, and only *free* (`~`) variables can be set.
"""
inject(mass, separation) = (; bodies=(; b=(; mass=mass, a=separation)))

"""
    detection(chain, θ_true) -> Bool

Recovered if the posterior median mass lands within a factor of three of the
truth *and* the 5th percentile clears 0.1 M_Jup — i.e. the fit is both
unbiased and confident that something is there.

Detection is applied at assembly time, not during the trial, so a different
threshold can be tried on the same results without re-sampling.
"""
function detection(chain, θ_true)
    mass_samples = vec(chain["b_mass"])
    med = median(mass_samples)
    low = quantile(mass_samples, 0.05)
    true_mass = θ_true.bodies.b.mass
    return (med > true_mass / 3) && (med < true_mass * 3) && (low > 0.1mjup)
end


# ──────────────────────────────────────────────────────────────────────
# 3. The template model
# ──────────────────────────────────────────────────────────────────────

# The Gaia DR4 reference epoch, J2017.5.
const DR4_REF_EPOCH_MJD = 57936.375

"""
    build_system(target) -> System

The template system: a host star, a dark companion, and one simulated Gaia DR4
epoch-astrometry observation. Its `centroid_pos_al` column is all zeros —
`run_completeness_trial` replaces the data with a simulation from the injected
parameters, so only the epochs, the scan geometry and the uncertainties of
`target.transits` matter here.
"""
function build_system(tgt)
    orbit_ref_epoch = mean(tgt.transits.epoch)
    plx_catalog = Float64(tgt.plx)

    gaia_obs = GaiaDR4AstromObs(tgt.transits;
        target=Photocentre,
        ref=Barycentre,
        name="GaiaDR4",
        variables=@variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)  # mas
            ra_offset_mas ~ Normal(0, 10000)
            dec_offset_mas ~ Normal(0, 10000)
            pmra ~ Uniform(-1000, 1000)                   # mas/yr
            pmdec ~ Uniform(-1000, 1000)                  # mas/yr
            ref_epoch = $DR4_REF_EPOCH_MJD
        end
    )

    # The host star is a body like any other. It carries the flux scale the
    # photocentre is weighted by; with `b.flux = 0` the photocentre is exactly
    # the host, which is the right model for a dark companion.
    A = Body(
        name="A",
        variables=@variables begin
            mass = 1.0   # M⊙
            flux = 1.0
        end
    )

    b = Body(
        name="b",
        about=A,
        variables=@variables begin
            mass ~ LogUniform(0.01mjup, 1000mjup)   # M⊙
            flux = 0.0
            a ~ LogUniform(0.1, 100)                # AU
            e ~ Uniform(0, 0.99)
            ω ~ Uniform(0, 2pi)
            i ~ Sine()
            Ω ~ Uniform(0, 2pi)
            θ ~ Uniform(0, 2pi)
            epoch = $orbit_ref_epoch
        end
    )

    return System(
        name="DR4_completeness",
        bodies=[A, b],
        observations=[gaia_obs],
        variables=@variables begin
            # No `$`: the right-hand side of a `~` is evaluated where the
            # `@variables` block is written, so a local is already visible.
            plx ~ truncated(Normal(plx_catalog, 0.5), lower=0.1)   # mas
        end
    )
end

"""
    describe_target(target)

One-line summary of the noise budget, so every entry point reports the same
numbers.
"""
function describe_target(tgt)
    @info "Gaia DR4 completeness target" gaia_id = tgt.gaia_id G = tgt.phot_g_mean_mag plx_mas = tgt.plx n_transits = tgt.n_transits
    @info "Measured G23H along-scan noise budget [mas]" σ_AL = tgt.σ_AL σ_att = tgt.σ_att σ_calib = tgt.σ_calib σ_formal_per_ccd = tgt.σ_formal n_ccd_per_transit = tgt.n_ccd σ_transit_formal = tgt.σ_transit_formal σ_transit_true = tgt.σ_transit_true
    return nothing
end
