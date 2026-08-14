# [Simulating and Fitting Gaia DR4 Data](@id data-simulation-dr4)

This tutorial demonstrates how to use Octofitter to simulate Gaia DR4 data for a particular target. 

**Before reading this, take a look at the general [Generating and Fitting Simulated Data](@ref data-simulation) tutorial.**


----

To be as realistic as possible, we will select a particular star as the basis of our simulation. This will allow us to give reasonable estimates of
* the measurement epochs and scan angles (retrieved from GOST)
* the measurement noise, representative of that coordinate and target photometry/colour

> *"But wait, I just want to simulate a hypothetical planet, why do I have to pick a real star?"*

Gaia's sensitivity is far from uniform! If you're not sure what star to use, just pick a the Gaia ID of your favourite star that has a representative magnitude and colour.

## Setup

```@example 1
using Octofitter, Distributions
using CairoMakie, PairPlots, Pigeons
using CSV, DataFrames
using Statistics
```

## Prepare Star Info and Noise

We now query the Gaia positions etc. from DR3, the scan law from GOST, and the uncertainties from Thompson et al (2026).

```@example 1
gaia_id = 5064625130502952704
# We will use the σ_att and σ_AL for this target's position, photometry, and colour as the 
# formal uncertainty on each scan. They can be retrieved or guessed from Keifer et al 2025,
# or Thompson et al in-prep.

# look these up, or download the tables and interpolate
σ_att = 0.04
σ_AL = 0.04
σ_cal = 0.04
σ_formal = sqrt(σ_att^2 + σ_AL^2)
σ_true = sqrt(σ_att^2 + σ_AL^2 + σ_cal^2)



dr3 = Octofitter._query_gaia_dr3(;gaia_id)
# TODO: we have on average ~8 scans per epoch
gost = DataFrame(Octofitter.GOST_forecast(dr3.ra,dr3.dec;baseline=:dr4))

nothing # hide
```


**If you want to simulate per CCD measurement (10x slower), run the following:**

```julia
"""
    expand_scanlaw_to_scans(scanlaw_table::DataFrame, star_row::DataFrameRow, 
                           avg_scans_per_window::Float64=9.0)

Expand visibility window-level scanlaw to individual CCD scans.

# Arguments
- `scanlaw_table`: DataFrame with columns `times` (OBMT) and `angles` (degrees)
- `star_row`: Row from star catalog containing stellar parameters
- `avg_scans_per_window`: Average number of CCD scans per visibility window

# Returns
- DataFrame with expanded scan-level data
"""
function expand_scanlaw_to_scans(scanlaw_table, dr3, avg_scans_per_window=9.0)
    # Calculate actual average if available
    if hasproperty(dr3, :astrometric_n_good_obs_al_dr3) && 
       hasproperty(dr3, :astrometric_matched_transits_dr3) &&
       dr3.astrometric_matched_transits_dr3 > 0
        avg_scans_per_window = dr3.astrometric_n_good_obs_al_dr3 / 
                              dr3.astrometric_matched_transits_dr3
    end
    
    expanded_data = []
    
    for (i, row) in enumerate(eachrow(scanlaw_table))
        # Number of CCDs for this window (randomize around average)
        n_ccds = rand(Poisson(avg_scans_per_window))
        n_ccds = max(1, min(n_ccds, 12))  # Gaia has max 9 CCDs in practice
        
        # Time spread within window (typically ~40 seconds total)
        window_duration_seconds = 40.0
        dt_seconds = window_duration_seconds / max(n_ccds - 1, 1)
        
        # Small angle variation within window
        angle_variation = 0.1  # degrees
        
        for j in 1:n_ccds
            # Time offset from window center
            time_offset_days = (j - (n_ccds + 1) / 2) * dt_seconds / 86400
            obs_time_obmt = row.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_ + time_offset_days * 1461.0 / 365.25  # Convert days to OBMT
            
            # Small angle perturbation
            angle_perturbation = angle_variation * (j - (n_ccds + 1) / 2) / n_ccds
            scan_angle = row.scanAngle_rad_ + deg2rad(angle_perturbation)
            
            push!(expanded_data, (
                transit_id = i,
                ccd_id = j,
                ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_ = obs_time_obmt,
                scanAngle_rad_ = scan_angle,
            ))
        end
    end
    
    return DataFrame(expanded_data)
end

gost_ex = expand_scanlaw_to_scans(gost, dr3)
r = size(gost_ex,1)/size(gost,1)
σ_formal /= r
σ_true /= r
gost = gost_ex
```

```@example 1
N_epochs = size(gost,1)


# σ_cal is *not* part of the Gaia formal uncertainties, it's noise above and beyond the stated uncertainties.
# For an accurate simulation, we should randomize by the true uncertainty (both sigma_formal + sigma_cal) but report to the downstream fitter that the uncertainty is just σ_formal.

df = DataFrame(
    # Epoch of measurements in MJD
    epoch = jd2mjd.(gost.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_),
    # Scan angle in DEGREES. GOST reports radians (`scanAngle_rad_`), while
    # `GaiaDR4AstromObs` ingests degrees — the unit the Gaia archive publishes
    # the DR4 scan angle in — and converts to radians internally.
    scan_pos_angle = rad2deg.(gost.scanAngle_rad_),
    # In theory you could just populate this vector with whatever measurements you want,
    # BUT we can also leave it blank, and leverage Octofitter's built in data simulation 
    # capabilities below...
    centroid_pos_al       = fill(0.0, N_epochs), 
    centroid_pos_error_al = fill(σ_true, N_epochs), 
    outlier_flag          = fill(false, N_epochs), 
)


# The trickiest part: we have to calculate the parllax factors for each epoch
earth_pos_vel = DataFrame(Octofitter.geocentre_position_query.(df.epoch))
df = [df earth_pos_vel]
# we now have columns x, y, z of the earth in AU, and can calculate the parallax factor
# for each row

# Calculate parallax factors for each epoch

# For each epoch, calculate the parallax factor in along-scan direction
# The parallax displacement has components in RA and Dec:
# Δα* = plx * (x*sin(α) - y*cos(α))
# Δδ  = plx * (x*cos(α)*sin(δ) + y*sin(α)*sin(δ) - z*cos(δ))

# Then project onto the scan direction using the scan angle. `scan_pos_angle`
# is in degrees (see above), so this uses `sind`/`cosd` like the α/δ terms.
df.parallax_factor_al = @. (
    (df.x * sind(dr3.ra) - df.y * cosd(dr3.ra)) * cosd(df.scan_pos_angle) +
    (df.x * cosd(dr3.ra) * sind(dr3.dec) + df.y * sind(dr3.ra) * sind(dr3.dec) - df.z * cosd(dr3.dec)) * sind(df.scan_pos_angle)
)

# now construct the observation template
ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd

gaiaIADobs = GaiaDR4AstromObs(df;
    target = Photocentre,
    ref = Barycentre,
    name = "GaiaDR4",
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10) # mas
        ra_offset_mas  ~ Normal(0, 10000)
        dec_offset_mas ~ Normal(0, 10000)
        pmra ~ Uniform(-1000, 1000) # mas/yr
        pmdec ~  Uniform(-1000, 1000) # mas/yr
        ref_epoch = $ref_epoch_mjd
    end
)

nothing # hide
```

!!! note "What the observation carries"
    `GaiaDR4AstromObs` takes data, references and variables and nothing else — no
    `gaia_id=`, and no archive query. It declares `target=Photocentre` (the system's
    flux-weighted point, the default) measured against `ref=Barycentre`, and reads its
    parallax from the system's own `plx`.

## Define the Model

Now, we define a model that incorporates this data. Bodies are model nodes in their own
right — the host star included — and every mass is in **solar masses**:

```@example 1
orbit_ref_epoch = mean(gaiaIADobs.table.epoch)

A = Body(
    name="A",
    variables=@variables begin
        mass = 1.0     # Msol
        flux = 1.0     # sets the flux scale the photocentre is weighted by
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass ~ LogUniform(0.01mjup, 1000mjup)   # Msol
        flux = 0.0                              # dark companion
        a ~ LogUniform(0.01, 100)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0,2pi)
        i ~ Sine()
        Ω ~ Uniform(0,2pi)
        θ ~ Uniform(0,2pi)
        epoch = $orbit_ref_epoch
    end
)

sys = System(
    name="target_1",
    bodies=[A, b],
    observations=[gaiaIADobs],
    variables=@variables begin
        # Note: keep these physically plausible to prevent numerical errors
        plx ~ Uniform(0.01,100) # mas
    end
)
model = Octofitter.LogDensityModel(sys, verbosity=4)
```

A `Photocentre` target needs at least one body to declare a flux, so `A` carries
`flux = 1.0`; every other body's `flux` is then a contrast ratio against it. With
`b.flux = 0.0` the photocentre is exactly the host, which is the right model for a dark
companion.

## Simulating a Planet

Now we have defined our template model! Normally, we would just fit this...
Here, we will hook into Octofitter's simulation capabilities to generate a new model for us with a simulated planet.



## Generate Synthetic Data

We have three choices for generating simulated data:
1. Draw values from the priors
3. Specifying values manually

We will look at each.

### 1. Draw values from the priors

We can draw a value from the priors like so:
```@example 1
params_to_simulate = Octofitter.drawfrompriors(model.system)
```


### 3. Specifying values manually
We can also specify all values for the simulation manually. This process is a bit more involved. 


!!! warning
    Note that the output below is just an example, you must generate your own template from your model and modify it as needed. The exact structure is not garuanteed to be stable between versions of Octofitter.

Note the shape: system variables at the top level, then **`bodies`** — the host star
included — then `observations`.

```@example 1
params_to_simulate = (
    plx = 16.138978209522527,
    bodies = (
        A = (
            mass = 1.0,
            flux = 1.0,
        ),
        b = (
            mass = 10.0 * Octofitter.mjup,
            a = 2.0,
            e = 0.01,
            ω = 0.01,
            i = 0.01,
            Ω = 0.6,
            θ = 3.5226970272017826,
            epoch = orbit_ref_epoch,
            flux = 0.0,
        ),
    ),
    observations = (
        GaiaDR4 = (
            astrometric_jitter = 0.015590772157762368,
            ra_offset_mas = 0.05413791355838311,
            dec_offset_mas = 0.0889816167388366,
            pmra = 5.301074192615374,
            pmdec = -24.188882826919325,
            ref_epoch = 57936.375,
        ),
    ),
)
nothing # hide
```


## Generate synthetic system with simulated data

We will call `Octofitter.generate_from_params` to generate a model with new synthetic observations.
If you set `add_noise = true`, the generated data points will have scatter according to the `centroid_pos_error_al` specified above. 


```@example 1
sim_system = Octofitter.generate_from_params(model.system, params_to_simulate; add_noise=true)
sim_model = Octofitter.LogDensityModel(sim_system)
```


## Fit simulated data

Follow the usual flow documented elsewhere.

Find starting point for MCMC (via variational approximation)
```@example 1
init_chain = initialize!(sim_model)
octoplot(sim_model, init_chain)
```


Fit:
```@example 1
chain, pt = octofit_pigeons(sim_model, n_rounds=9)
```

!!! note "Why parallel tempering here"
    A single epoch-astrometry time series constrains the orbit only through the
    along-scan abscissa, and that leaves genuinely separated solutions — most
    conspicuously the ±180° ambiguity in Ω and the degeneracy between period and
    the fraction of an orbit the scan window covers. [`octofit_pigeons`](@ref) runs
    a ladder of tempered chains down to the prior, where those solutions merge, so
    replicas can move between them and each stays represented in the posterior.

    [`octofit`](@ref) will sample this model too and is considerably cheaper. If you
    use it, seed it with [`initialize!`](@ref) as above and run a few chains from
    different starting points, so a mode one chain missed shows up as disagreement
    rather than as false confidence.


```@example 1
octoplot(sim_model, chain)
```

`GaiaDR4AstromObs` declares one plot channel (the along-scan abscissa), so `octoplot` draws
the data, the modelled abscissae, a residual strip and a residual histogram automatically.
Restrict the [`PosteriorSeries`](@ref) to get the "one particular draw" version:

```@example 1
octoplot(Octofitter.PosteriorSeries(sim_model, chain; ii=[1]))
```

[`gaiastarplot`](@ref) shows that single draw in the sky plane instead — the reflex track
with each transit's along-scan residual re-projected along its own scan angle, which is
where the one-dimensional nature of a Gaia measurement becomes visible:

```@example 1
Octofitter.gaiastarplot(sim_model, chain, 1)
```

Finally, compare the recovered orbit against the truth: rebuild the
`PlanetOrbits.System` for each draw with `construct_system` and plot the tracks directly,
which is what the sky panel does internally:

```@example 1
ts = range(minimum(df.epoch), minimum(df.epoch) + 4000, length=300)

fig = Figure(size=(500,500))
ax = Axis(fig[1,1], xlabel="Δα⋆ [mas]", ylabel="Δδ [mas]",
          aspect=DataAspect(), xreversed=true)

for i in rand(1:size(chain,1), 50)
    posys = construct_system(sim_model, chain, i)
    traj = orbitsolve(posys, ts)
    lines!(ax, [raoff(traj[k], :b, :A) for k in eachindex(ts)],
               [decoff(traj[k], :b, :A) for k in eachindex(ts)],
           color=(:black, 0.1))
end

# the true orbit, in red
true_sys = construct_system(sim_model, params_to_simulate)
true_traj = orbitsolve(true_sys, ts)
lines!(ax, [raoff(true_traj[k], :b, :A) for k in eachindex(ts)],
           [decoff(true_traj[k], :b, :A) for k in eachindex(ts)],
       color=:red, linewidth=3, label="truth")
scatter!(ax, [0], [0], color=:orange, markersize=12)
axislegend(ax)
fig
```

```@example 1
octocorner(sim_model, chain, small=true)
```
