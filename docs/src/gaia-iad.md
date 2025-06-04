# Fitting Gaia DR4 IAD

This tutorial shows how you can use Octofitter to fit preliminary RV and absolute astrometry data from DR4, using the already published data for the Gaia-BH3 black hole.

The final format of the Gaia IAD data may change, in which case this tutorial will be updated.

```@example 1
using CairoMakie
using Octofitter
using Distributions
using CSV, DataFrames
```


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

scatter(
    df.obs_time_tcb,
    df.centroid_pos_al,
)
```

We can now construct a likelihood object for this data. We must also supply the Gaia ID, which will be used to query the full Gaia solution for this object (for now, using DR3):
```@example 1
gaiaIADLike = GaiaDR4Astrom(df, gaia_id=4318465066420528000)
```

This object also has published RV data from Gaia, which we can load and use as normal:
```@example 1
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

rvlike = StarAbsoluteRVLikelihood(
    dfrv,
    offset=:offset_gaiarv,
    jitter=:jitter_gaiarv,
)
errorbars(
    dfrv.obs_time_tcb,
    dfrv.rv,
    dfrv.σ_rv
)
```

Now, we define a model that incorporates this data:
```@example 1
mjup2msol = Octofitter.mjup2msol
ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd
orbit_ref_epoch = mean(gaiaIADlike.table.epoch)
mean_rv = mean(dfrv.rv)
@planet b Visual{KepOrbit} begin
    a ~ Uniform(0, 1000)
    e ~ Uniform(0, 0.99)
    ω ~ UniformCircular()
    i ~ Sine()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system, b, $orbit_ref_epoch)
    mass = system.M_sec / $mjup2msol
end
sys = @system gaiadr4test begin
    M_pri = 0.76
    M_sec ~ LogUniform(1, 1000) # Msol
    M = system.M_pri + system.M_sec# 
    plx ~ Uniform(0.01, 100)
    pmra ~ Uniform(-1000, 1000)
    pmdec ~ Uniform(-1000, 1000)

    # Put a prior of the catalog value +- 10,000 mas on position
    # We could just fit the deltas directly, but we can propagate 
    # handling all non-linear effects if we know the actual ra, 
    # dec, and rv.
    ra_offset_mas ~ Normal(0, 100)
    dec_offset_mas ~ Normal(0, 100)
    dec = $gaiaIADlike.gaia_sol.dec + system.ra_offset_mas / 60 / 60 / 1000
    ra = $gaiaIADlike.gaia_sol.ra + system.dec_offset_mas / 60 / 60 / 1000 / cosd(system.dec)

    # Absolutely no idea what range is reasonable here
    astrometric_jitter ~ LogUniform(0.00001, 10) # mas
    offset_gaiarv ~ Normal(mean_rv, 10_000) # wide prior on RV offset centred on mean RV
    jitter_gaiarv ~ LogUniform(0.01, 100_000)
end gaiaIADlike rvlike b

model = Octofitter.LogDensityModel(sys, verbosity=4)
```


We will initialize the model starting positions and visualize them:
```@example 1
init_chain = initialize!(model, (;
    plx=gaiaIADlike.gaia_sol.parallax,
    pmra=gaiaIADlike.gaia_sol.pmra,
    pmdec=gaiaIADlike.gaia_sol.pmdec,
    ra_offset_mas=0.0,
    dec_offset_mas=0.0,
    astrometric_jitter=0.1,
    jitter_gaiarv=0.11,
))
octoplot(model, init_chain, show_rv=true)
```

Now, we can perform the fit. It is a little slow since we have many hundreds of RV and astrometry data points.
```@example 1
using Pigeons
chain, pt = octofit_pigeons(model, n_rounds=10) # might need more rounds to converge
```

Finally, we can visualize the results:
```@example 1
octoplot(model, chain, show_rv=true, mark_epochs_mjd=mjd.([
    "2017"
    "2022"
    "2027"
]))
```

```@example 1
octocorner(model, chain, small=true)
```

```@example 1
Octofitter.rvpostplot(model, chain)
```

