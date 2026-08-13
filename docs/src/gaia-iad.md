# Fitting Gaia DR4 IAD

This tutorial shows how you can use Octofitter to fit preliminary RV and absolute astrometry data from DR4, using various already published data, including:
* Gaia-BH3 black hole (the only real DR4 data released so far)
* Data from https://dace.unige.ch/openData/?record=10.82180/dace-gaia-ohp

The final format of the Gaia IAD data may change, in which case this tutorial will be updated.

```@example 1
using CairoMakie
using Octofitter
using Pigeons
using Distributions
using CSV, DataFrames
using Statistics
```

## OHP Gaia Splinter Session

Download and plot the data:
```julia
fname = download("https://dace.unige.ch/downloads/open_data/dace-gaia-ohp/files/target_1.csv")
df = CSV.read(fname, DataFrame);
```

```@example 1
df = CSV.read(joinpath(@__DIR__, "target_1.csv"), DataFrame)
nothing # hide
```

We can now construct a likelihood object for this data.

```@example 1
ref_epoch_mjd = 57936.375

gaia_dr4_obs = GaiaDR4AstromObs(
    df,
    target = Photocentre, # <-- this is how we handle blended sources
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
```

You will need to add a flux variable to at least the primary star. The units are up to you. You can just put `flux = 1.0` for `A`, and `flux = 0.0` for `b` (or leave it off) for a dark companion.
Use `flux_G` and `target=Photocentre(:G)` if your has more than one band.

```@example 1
orbit_ref_epoch = mean(gaia_dr4_obs.table.epoch)

A = Body(
    name="A",
    variables=@variables begin
        mass = 1.0     # Msol
        flux = 1.0     # the host defines the flux scale for the photocentre
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        flux = 0.0     # dark companion: the photocentre is then just the host
        a ~ LogUniform(0.01, 100)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0,2pi)
        i ~ Sine()
        Ω ~ Uniform(0,2pi)
        θ ~ Uniform(0,2pi)
        epoch = $orbit_ref_epoch
        mass ~ LogUniform(0.01mjup, 1000mjup)   # Msol
    end
)

sys = System(
    name="target_1",
    bodies=[A, b],
    observations=[gaia_dr4_obs],
    variables=@variables begin
        # Note: keep these physically plausible to prevent numerical errors
        plx ~ Uniform(0.01,100) # mas
    end
)
model = Octofitter.LogDensityModel(sys, verbosity=4)
```

Find a starting point via global optimization and variational approximation, and plot the initial points against the data:
```@example 1
init_chain = initialize!(model)
# optionally plot it
# octoplot(model, init_chain)
```

Sample using parallel tempering (could also use HMC for these unimodal distributions):
```@example 1
chain, pt = octofit_pigeons(model, n_rounds=10)
```

!!! note "Sampling faster on a multi-core machine"
    For a long fit on a machine with many cores, add `cores=8` (or however many
    you have free) to any `octofit_pigeons` call on this page to run the
    sampler in separate worker processes. Costs a minute or two to startup per run, then goes faster.
    See [`octofit_pigeons`](@ref).


Plot the results:
```@example 1
octoplot(model, chain)
```

This shows a sample of posterior draws. Lines show the locations in the sky plane of the companion relative to the primary. 
It's not valid to show Gaia DR4 datapoints on this plot, since they are 1D only and would move around with every draw. 


Gaia constrains one direction per transit, so a measurement is a *line*, not a point. That means any 2D view has to use our assumed model and a particular posterior draw to decide where to draw the datapoint, and the result varies per draw.
If we do pick an individual draw, we can plot the star's orbit against the Gaia data more directly.
Like RV, this only works with individual draws, because the Gaia points are "detrended" in a
sense from the values of a particular draw (they move around the plot depending on the
draw).


```@example 1
Octofitter.gaiastarplot(model, chain) # plots one draw only
```


A good practice is to plot a few different values from the posterior, or e.g. plot draws from 5th, 50th, and 95th percentile in a key orbit parameter like
*  semi-major axis
*  eccentricity 
*  inclination 

Here, we show semi-major axis/period
```@example 1
percentile_positions = round.(Int, [0.05, 0.50, 0.95] .* size(chain,1))
indices = [partialsortperm(chain["b_a"][:], k) for k in percentile_positions]

# hint! Try `"b_e", "b_a", and "b_i"

fig = Figure(size=(920,345))
for (col, idx) in enumerate(indices)
    fig = Octofitter.gaiastarplot(model, chain, idx)
    display(fig)
end
fig
```

[`skytrackplot`](@ref) zooms out to the source's *whole* path: parallactic loops, proper
motion, and the orbital wobble superimposed. If you want, pass `keplerian_mult` to exaggerate the orbital motion term so you can see it against the parallax.

The parallax ellipse has to be projected onto a sky direction. The observation does not
carry a catalog solution, so give the direction explicitly, or let
[`gaia_dr3_solution`](@ref) fetch the published one:

```julia
Octofitter.skytrackplot(model, chain; gaia_id=4318465066420528000, keplerian_mult=20)

# …or, with no network access, name the direction yourself:
Octofitter.skytrackplot(model, chain; ra=294.83, dec=14.93, keplerian_mult=20)
```

(When the model declares a full absolute frame, `ra`/`dec` are taken from it and neither
keyword is needed.)

And of course, you can make a corner plot as usual:
```@example 1
using PairPlots
octocorner(model, chain, small=true)
```

### Cross-validation
You can use the full suite of tools for construting models that subset different amounts of data in different ways. See "Cross Validataion".

The most powerful is exhaustive leave-one-out cross validataion plus calculation of expected log pointwise density (ELPD) to score the fit.

```julia
likelihood_mat, epochs = Octofitter.pointwise_like(model, chain)
# `likelihood_mat` is now a N_sample x N_data matrix. 
using ParetoSmooth
result = psis_loo(
    collect(likelihood_mat'),
    chain_index=ones(Int,size(chain,1))
)
```

### Simulate data from a posterior draw, and re-fit with or without noise
Optional consistency checks---could be used in a loop as part of e.g. simulation based calibration.
```julia
template = Octofitter.mcmcchain2result(model, chain, 1)
sim_system = Octofitter.generate_from_params(model.system, template; add_noise=true)
sim_model = Octofitter.LogDensityModel(sim_system)

# Optional initialization speed up hack:
sim_model.starting_points = model.starting_points;

# then re-fit...
sim_chain, pt = octofit_pigeons(sim_model, n_rounds=5)
octoplot(sim_model, sim_chain)
```

## Gaia BH 3

The following tutorial reproduces the fit to Gaia BH3. This one can take longer to run
since the orbit is ultra well constrained. For that reason, we don't run it automatically when building the docs. Please go ahead and run the code on your own computer; ETA=approx 20 minutes.

As a first step, we will load the astrometry data for Gaia-BH3 and plot it:
```@example 1
headers = [
    :transit_id
    :ccd_id
    :obs_time_tcb
    :centroid_pos_al
    :centroid_pos_error_al
    :parallax_factor_al
    :scan_pos_angle
    :outlier_flag
]
df = CSV.read(joinpath(@__DIR__, "astrom.dat"), DataFrame, skipto=7, header=headers, delim=' ', ignorerepeated=true)
df.epoch = jd2mjd.(df.obs_time_tcb)

# The table published with the BH3 paper gives scan_pos_angle in DEGREES, but
# GaiaDR4AstromObs takes sincos(scan_pos_angle) directly and so expects RADIANS.
df.scan_pos_angle = deg2rad.(df.scan_pos_angle)
nothing # hide
```

We can now construct a likelihood object for this data:
```@example 1
ref_epoch_mjd = 57936.375
gaia_bh3_astrom_obs = GaiaDR4AstromObs(
    df,
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
```

Next, the two bodies. Note that both masses are in solar masses:
```@example 1
orbit_ref_epoch = mean(gaia_bh3_astrom_obs.table.epoch)

A_bh3 = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(0.76, 0.05), lower=0.1)   # Msol
        flux = 1.0     # the luminous star sets the photocentre's flux scale
    end
)

BH = Body(
    name="BH",
    about=A_bh3,
    variables=@variables begin
        mass ~ LogUniform(1.0, 1000.0)   # Msol
        flux = 0.0     # a black hole contributes no light
        a ~ Uniform(0, 1000)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0,2pi)
        i ~ Sine()
        Ω ~ Uniform(0,2pi)
        θ ~ Uniform(0,2pi)
        epoch = $orbit_ref_epoch
    end
)
```

This object also has published RV data from Gaia, which we can load and use as normal.


```@example 1
headers_rv = [
    :transit_id,
    :obs_time_tcb,
    :radial_velocity_kms,
    :radial_velocity_err_kms,
]
dfrv = CSV.read(joinpath(@__DIR__, "epochrv.dat"), DataFrame, skipto=7, header=headers_rv, delim=' ', ignorerepeated=true)
dfrv.epoch = jd2mjd.(dfrv.obs_time_tcb)
dfrv.rv = dfrv.radial_velocity_kms * 1e3
dfrv.σ_rv = dfrv.radial_velocity_err_kms * 1e3

# Calculate mean RV for the prior
mean_rv = mean(dfrv.rv)

rvlike = RadialVelocityObs(
    dfrv,
    target = A_bh3,
    ref = Barycentre,
    name="GaiaRV",
    variables=@variables begin
        offset ~ Normal(mean_rv, 10_000)  # wide prior on RV offset centred on mean RV  
        jitter ~ LogUniform(0.01, 100_000)  # RV jitter parameter
    end
)
errorbars(
    dfrv.obs_time_tcb,
    dfrv.rv,
    dfrv.σ_rv
)
```

Now, we assemble the system. Both observations are listed on the system; neither is attached
to a body:
```@example 1
sys = System(
    name="gaiadr4test",
    bodies=[A_bh3, BH],
    observations=[gaia_bh3_astrom_obs, rvlike],
    variables=@variables begin
        # Note: keep these physically plausible to prevent numerical errors
        plx ~ Uniform(0.01,100) # mas
    end
)
model = Octofitter.LogDensityModel(sys, verbosity=4)
```

We will initialize the model starting positions:
```@example 1
# Note: you can see the required format for parameter initialization by running:
# nt = Octofitter.drawfrompriors(model.system);
# println(nt)
# then remove any derived parameters (parameters in your model that are on the right of an `=`)

init_chain = initialize!(model, (;
    plx = 1.6686144513164856,
    bodies = (;
        A  = (; mass = 0.7792923132247755),
        BH = (; mass = 36.032664849109906,
                a = 18.905647598089196,
                e = 0.7583328001601555,
                i = 1.9216027029499319),
    ),
    observations = (;
        GaiaDR4 = (; astrometric_jitter = 0.027554101045898238),
        GaiaRV  = (; offset = -359481.6706770764, jitter = 2143.4793485877644),
    ),
))
```

If you don't pick the starting point, you can also just run Pigeons for 8-10 rounds,
which is recommended anyway for convergence, and the sampler will find this result.

Now, we can perform the fit. It is a little slow since we have many hundreds of RV and astrometry data points.
```@example 1
chain, pt = octofit_pigeons(model, n_rounds=6) # might need more rounds to converge
```

Pigeons can resume an existing run, so you can add rounds until the reported
`log(Z₁/Z₀)` and `Λ` settle instead of starting over:
```@example 1
increment_n_rounds!(pt, 1)
chain, pt = octofit_pigeons(pt)
```

Finally, we can visualize the results. `octoplot` produces the along-scan panel from the
astrometry and a radial-velocity panel from the RVs, so one
call covers the whole fit:
```@example 1
octoplot(model, chain)
```

```@example 1
rvplot(model, chain)
```


[`rvpostplot`](@ref) is the RV-only slice of that same figure, and
[`gaiastarplot`](@ref) / [`skytrackplot`](@ref) add the sky-plane views of the astrometry:
```@example 1
rvpostplot(model, chain)
```

```@example 1
Octofitter.gaiastarplot(model, chain)
```

```@example 1
Octofitter.skytrackplot(model, chain; gaia_id=4318465066420528000)
```

```@example 1
octocorner(model, chain, small=true)
```
