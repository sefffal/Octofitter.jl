#=
Gaia DR4 Completeness Trial — single array job element.

Runs one injection-recovery trial from a completeness grid.
The job index comes from SLURM_ARRAY_TASK_ID (1-based).

Usage:
    julia --project=/scratch/wthompso/completeness_dr4 --threads=auto completeness_trial.jl

Environment variables:
    SLURM_ARRAY_TASK_ID  — 1-based job index into the completeness grid
=#

using Octofitter, Distributions
using DataFrames, Serialization

# ──────────────────────────────────────────────────────────────────
# Grid definition (must match assemble_results.jl)
# ──────────────────────────────────────────────────────────────────
const MASSES      = 10.0 .^ range(-1, 2, length=12)   # 0.1 to 100 Mjup
const SEPARATIONS = 10.0 .^ range(-0.3, 1.7, length=12) # ~0.5 to 50 AU
const N_TRIALS    = 5

# ──────────────────────────────────────────────────────────────────
# Target star setup
# ──────────────────────────────────────────────────────────────────
gaia_id = 5064625130502952704

# Noise: UEVA_single approach
σ_att = 0.04
σ_AL = 0.04
σ_cal = 0.04
σ_true = sqrt(σ_att^2 + σ_AL^2 + σ_cal^2)

# Query DR3 (cached) and build scan law table
dr3 = Octofitter._query_gaia_dr3(; gaia_id)
gost = DataFrame(Octofitter.GOST_forecast(dr3.ra, dr3.dec; baseline=:dr4))
N_epochs = size(gost, 1)

df = DataFrame(
    epoch             = Octofitter.jd2mjd.(gost.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_),
    scan_pos_angle    = gost.scanAngle_rad_,
    centroid_pos_al   = fill(0.0, N_epochs),
    centroid_pos_error_al = fill(σ_true, N_epochs),
    outlier_flag      = fill(false, N_epochs),
)

# Earth positions for parallax factors
earth_pos_vel = DataFrame(Octofitter.geocentre_position_query.(df.epoch))
df = [df earth_pos_vel]
df.parallax_factor_al = @. (
    (df.x * sind(dr3.ra) - df.y * cosd(dr3.ra)) * cos(df.scan_pos_angle) +
    (df.x * cosd(dr3.ra) * sind(dr3.dec) + df.y * sind(dr3.ra) * sind(dr3.dec) - df.z * cosd(dr3.dec)) * sin(df.scan_pos_angle)
)

ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd

# ──────────────────────────────────────────────────────────────────
# Model definition
# ──────────────────────────────────────────────────────────────────
gaiaIADobs = GaiaDR4AstromObs(df;
    gaia_id=gaia_id,
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)
        ra_offset_mas  ~ Normal(0, 10000)
        dec_offset_mas ~ Normal(0, 10000)
        pmra  ~ Uniform(-1000, 1000)
        pmdec ~ Uniform(-1000, 1000)
        plx = system.plx
        ref_epoch = $ref_epoch_mjd
    end
)

orbit_ref_epoch = mean(gaiaIADobs.table.epoch)

b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[],
    variables=@variables begin
        a    ~ LogUniform(0.1, 100)
        e    ~ Uniform(0, 0.99)
        ω    ~ Uniform(0, 2pi)
        i    ~ Sine()
        Ω    ~ Uniform(0, 2pi)
        θ    ~ Uniform(0, 2pi)
        tp   = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
        mass ~ LogUniform(0.01, 1000)
    end
)

sys = System(
    name="DR4_completeness",
    companions=[b],
    observations=[gaiaIADobs],
    variables=@variables begin
        M   = 1.0
        plx ~ truncated(Normal($(Float64(dr3.parallax)), 0.5), lower=0.1)
    end
)

# ──────────────────────────────────────────────────────────────────
# Run the trial
# ──────────────────────────────────────────────────────────────────
job_id = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

jobs = completeness_jobs(; masses=MASSES, separations=SEPARATIONS, n_trials=N_TRIALS)

if job_id > length(jobs)
    @warn "Job index $job_id exceeds total jobs $(length(jobs)); exiting."
    exit(0)
end

job = jobs[job_id]

@info "Starting completeness trial" job_id mass=job.mass separation=job.separation trial=job.i_trial

# Detection criterion: mass recovered within a factor of 3
# and 5th percentile of mass posterior > 0.1 Mjup
detection = function(chain, θ_true)
    mass_samples = vec(chain["b_mass"])
    med = median(mass_samples)
    low = quantile(mass_samples, 0.05)
    true_mass = θ_true.planets.b.mass
    # Detected if: median within factor of 3, and 95% CI excludes near-zero
    return (med > true_mass / 3) && (med < true_mass * 3) && (low > 0.1)
end

# Injection function: override mass and semi-major axis
inject = (mass, sep) -> (; planets=(; b=(; mass=mass, a=sep)))

# Sampler: HMC with moderate iterations
sampler = function(model)
    octofit(model, 0.85;
        adaptation=1000,
        iterations=2000,
        verbosity=0,
    )
end

result = run_completeness_trial(
    job, sys, sampler, detection;
    inject=inject,
    add_noise=true,
    verbosity=1,
)

# ──────────────────────────────────────────────────────────────────
# Save result
# ──────────────────────────────────────────────────────────────────
outdir = "results"
mkpath(outdir)
outfile = joinpath(outdir, "trial_$(lpad(job_id, 4, '0')).jls")
serialize(outfile, result)

@info "Trial complete" detected=result.detected outfile
