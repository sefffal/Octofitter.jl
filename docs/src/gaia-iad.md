# Fitting Gaia DR4 IAD

This tutorial shows how you can use Octofitter to fit preliminary RV and absolute astrometry data from DR4, using various already published data, including:
* Gaia-BH3 black hole (the only real DR4 data released so far)
* Data from https://dace.unige.ch/openData/?record=10.82180/dace-gaia-ohp

The final format of the Gaia IAD data may change, in which case this tutorial will be updated.

```@example 1
using CairoMakie
using Octofitter
using Distributions
using CSV, DataFrames
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




We can now construct a likelihood object for this data. We must also supply the Gaia ID, which will be used to query the full Gaia solution for this object (for now, using DR3):
```@example 1
ref_epoch_mjd = 57936.375

gaia_dr4_obs = GaiaDR4AstromObs(
    df, 
    # For plotting reasons, you need to supply a Gaia ID to know e.g. the absolute Ra and Dec
    # For these ~simulated examples, you can pick one at random
    gaia_id=4373465352415301632,
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10) # mas
        ra_offset_mas  ~ Normal(0, 10000)
        dec_offset_mas ~ Normal(0, 10000)
        pmra ~ Uniform(-1000, 1000) # mas/yr
        pmdec ~  Uniform(-1000, 1000) # mas/yr
        # ra_offset_mas = -0.00154
        # dec_offset_mas = 0.0019354
        # pmra = 5.421
        # pmdec = -24.121
        plx = system.plx
        ref_epoch = $ref_epoch_mjd
    end
)
```


```@example 1
mjup2msol = Octofitter.mjup2msol
ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd
orbit_ref_epoch = mean(gaia_dr4_obs.table.epoch)

b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[],
    variables=@variables begin
        a ~ LogUniform(0.01, 100)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0,2pi)
        i ~ Sine()
        Ω ~ Uniform(0,2pi)
        θ ~ Uniform(0,2pi)
        tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
        mass ~ LogUniform(0.01, 1000)
    end
)
sys = System(
    name="target_1",
    companions=[b],
    observations=[gaia_dr4_obs,],
    variables=@variables begin
        M = 1.0
        # Note: keep these physically plausible to prevent numerical errors
        plx ~ Uniform(0.01,100) # mas
    end
)
model = Octofitter.LogDensityModel(sys, verbosity=4)
```

Find a starting point via global optimization and variational approximation, and plot the initial points against the data:
```@example 1
init_chain = initialize!(model)
octoplot(model, init_chain)
```


Sample using parallel tempering (could also use HMC for these unimodel distributions):
```@example 1
using Pigeons
chain, pt = octofit_pigeons(model, n_rounds=10)
```

Plot the results:
```@example 1
octoplot(model, chain)
```

If we pick an individual draw, we can also plot the orbit against the Gaia data more directly. Like RV, this only works with individual draws because the Gaia points are "detrended" in a sense from the values of a particular draw (they move around the plot depending on the draw).

```@example 1
# Picks MAP sample
# gaiastarplot(model, chain)           

# Pick specific arbitrary sample
idx = rand(1:size(chain,1)) # pick an integer randomly
Octofitter.gaiastarplot(model, chain, idx)    
```

A good practice is to plot a few different values from the posterior, or e.g. plot draws from 5th, 50th, and 95th percentile in a key orbit parameter like

*  semi-major axis
*  eccentricity 
*  inclination 

Here, we show semi-major axis/period
```@example 1
fig = Figure(size=(920,345))

percentile_positions = round.(Int, [0.05, 0.50, 0.95] .* size(chain,1))
indices = [partialsortperm(chain["b_a"][:], k) for k in percentile_positions]

# hint! Try `"b_e", "b_a", and "b_i"

ax1 = Octofitter.gaiastarplot!(fig[1,1], model, chain, indices[1]) 
ax2 = Octofitter.gaiastarplot!(fig[1,2], model, chain, indices[2])
ax3 = Octofitter.gaiastarplot!(fig[1,3], model, chain, indices[3]) 

ax1.title = "5th percentile period"
ax2.title = "50th percentile period"
ax3.title = "95th percentile period"
hideydecorations!(ax2)
hideydecorations!(ax3)
linkaxes!(ax1,ax2,ax3)

fig
```




We can also plot a sky track, reconstructed against the official gaia solution (queried automatically when constructing the observation object).
```@example 1
# specific sample
idx = 42
# optional parameter keplerian_mult: exagerates the keplerian component for visualization
Octofitter.skytrackplot(model, chain, idx, keplerian_mult=10) 
```



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
template = Octofitter.mcmcchain2result(model,chain,1)
sim_system = Octofitter.generate_from_params(model.system, template; add_noise=true)
sim_model = Octofitter.LogDensityModel(sim_system)

# Optional initialization speed up hack:
sim_model.starting_points = model.starting_points;

# then re-fit...
sim_chain, pt = octofit_pigeons(sim_model, n_rounds=5)
Octofitter.gaiastarplot(sim_model, sim_chain)
```


## Gaia BH 3

The following tutorial reproduces the fit to Gaia BH3. This one can take longer to run
since the orbit is ultra well constrained. For that reason, we don't run it automatically when building the docs. Please go ahead and run the code on your own computer; ETA=approx 20 minutes.

As a first step, we will load the astrometry data for Gaia-BH3 and plot it:
```julia
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

scatter(
    df.obs_time_tcb,
    df.centroid_pos_al,
)
```

We can now construct a likelihood object for this data. We must also supply the Gaia ID, which will be used to query the full Gaia solution for this object (for now, using DR3):
```julia
ref_epoch_mjd = 57936.375
gaia_bh3_astrom_obs = GaiaDR4Astrom(
    df, 
    gaia_id=4373465352415301632,
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10) # mas
        ra_offset_mas  ~ Normal(0, 10000)
        dec_offset_mas ~ Normal(0, 10000)
        pmra ~ Uniform(-1000, 1000) # mas/yr
        pmdec ~  Uniform(-1000, 1000) # mas/yr
        # ra_offset_mas = -0.00154
        # dec_offset_mas = 0.0019354
        # pmra = 5.421
        # pmdec = -24.121
        plx = system.plx
        ref_epoch = $ref_epoch_mjd
    end
)
```

This object also has published RV data from Gaia, which we can load and use as normal:
```julia
using CSV, DataFrames
using OctofitterRadialVelocity

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

rvlike = StarAbsoluteRVObs(
    dfrv,
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

Now, we define a model that incorporates this data:
```julia
mjup2msol = Octofitter.mjup2msol
ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd
orbit_ref_epoch = mean(gaia_dr4_obs.table.epoch)

b = Planet(
    name="BH",
    basis=Visual{KepOrbit},
    observations=[],
    variables=@variables begin
        a ~ Uniform(0, 1000)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0,2pi)
        i ~ Sine()
        Ω ~ Uniform(0,2pi)
        θ ~ Uniform(0,2pi)
        tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
        mass = system.M_sec / mjup2msol
    end
)
# DR3 catalog position and reference epoch
ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd
sys = System(
    name="gaiadr4test",
    companions=[b],
    observations=[gaia_bh3_astrom_obs, rvlike,],
    variables=@variables begin
        M_pri ~ truncated(Normal(0.76,0.05),lower=0.1) # M sun
        M_sec ~ LogUniform(1, 1000) # M sun
        M = M_pri + M_sec
        # Note: keep these physically plausible to prevent numerical errors
        plx ~ Uniform(0.01,100) # mas
    end
)
model = Octofitter.LogDensityModel(sys, verbosity=4)
```


We will initialize the model starting positions and visualize them:
```julia
# Note: you can see the required format for paramter initialization by running:
# nt = Octofitter.drawfrompriors(model.system);
# println(nt)
# then remove any derived parameters (parameters in your model that are on the right of an `=`)

init_chain = initialize!(model, (;
    M_pri = 0.7792923132247755,
    M_sec = 36.032664849109906,
    plx = 1.6686144513164856,
    observations = (GaiaDR4 = (astrometric_jitter = 0.027554101045898238,),
    GaiaRV = (offset = -359481.6706770764,jitter = 2143.4793485877644)),
    planets = (BH = (
        a = 18.905647598089196,
        e = 0.7583328001601555,
        i = 1.9216027029499319,
    ),
)))
octoplot(model, init_chain, show_rv=true)
```



!!! note
    If you don't pick the starting point, you cabn also just run Pigeons for 8-10 rounds, which is recommended anyways for convergence, and the sampler will find this result.

Now, we can perform the fit. It is a little slow since we have many hundreds of RV and astrometry data points.
```julia
using Pigeons
chain, pt = octofit_pigeons(model, n_rounds=6) # might need more rounds to converge
```

```julia
increment_n_rounds!(pt,1)
chain,pt = octofit_pigeons(pt)
```

Finally, we can visualize the results:
```julia
octoplot(model, chain, show_rv=true, mark_epochs_mjd=mjd.([
    "2017"
    "2022"
    "2027"
]))
```


```julia
# Picks MAP sample
Octofitter.gaiastarplot(model, chain)
```


```julia
octocorner(model, chain, small=true)
```

```julia
Octofitter.rvpostplot(model, chain)
```

