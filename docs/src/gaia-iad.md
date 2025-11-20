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
gaiaIADlike = GaiaDR4Astrom(
    df, 
    gaia_id=4318465066420528000,
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10) # mas
    end
)
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
```@example 1
mjup2msol = Octofitter.mjup2msol
ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd
orbit_ref_epoch = mean(gaiaIADlike.table.epoch)

b = Planet(
    name="BH",
    basis=AbsoluteVisual{KepOrbit},
    observations=[],
    variables=@variables begin
        a ~ Uniform(0, 1000)
        e ~ Uniform(0, 0.99)
        ω ~ UniformCircular()
        i ~ Sine()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
        mass = system.M_sec / mjup2msol
    end
)
# DR3 catalog position and reference epoch
gaia_ra = 294.82786250815144
gaia_dec = 14.930979608612361
ref_epoch_mjd = Octofitter.meta_gaia_DR3.ref_epoch_mjd
sys = System(
    name="gaiadr4test",
    companions=[b],
    observations=[gaiaIADlike, rvlike,],
    variables=@variables begin
        M_pri ~ truncated(Normal(0.76,0.05),lower=0.1) # M sun
        M_sec ~ LogUniform(1, 1000) # M sun
        M = M_pri + M_sec
        # Note: keep these physically plausible to prevent numerical errors
        plx ~ Uniform(0.01,100) # mas
        pmra ~ Uniform(-1000, 1000) # mas/yr
        pmdec ~  Uniform(-1000, 1000) # mas/yr
        rv = −333.2e3 # m/s

        # Put a prior of the catalog value +- 10,000 mas on position
        ra_offset_mas ~ Normal(0, 10000)
        dec_offset_mas ~ Normal(0, 10000)
        dec = $gaia_dec + ra_offset_mas / 60 / 60 / 1000
        ra = $gaia_ra + dec_offset_mas / 60 / 60 / 1000 / cosd(dec)
        # Important! This is the reference epoch for the ra and dec provided above, *not* necessarily DR4.
        ref_epoch = $ref_epoch_mjd
    end
)
model = Octofitter.LogDensityModel(sys, verbosity=4)
```


We will initialize the model starting positions and visualize them:
```@example 1
# Note: you can see the required format for paramter initialization by running:
# nt = Octofitter.drawfrompriors(model.system);
# println(nt)
# then remove any derived parameters (parameters in your model that are on the right of an `=`)

init_chain = initialize!(model, (;
    M_pri = 0.7792923132247755,
    M_sec = 36.032664849109906,
    plx = 1.6686144513164856,
    pmra = -27.89740759925553,
    pmdec = -156.1023951519146,
    ra_offset_mas = 236.8072112885035,
    dec_offset_mas = 45.781075653307376,
    observations = (GaiaDR4 = (astrometric_jitter = 0.027554101045898238,),
    GaiaRV = (offset = -359481.6706770764,jitter = 2143.4793485877644)),
    planets = (BH = (
        a = 18.905647598089196,
        e = 0.7583328001601555,
        ωx = -0.19433584569119122,
        ωy = -0.9414842877197981,
        i = 1.9216027029499319,
        Ωx = -0.9890745570284801,
        Ωy = 0.9268637554445821,
        θx = 0.4250634152645573,
        θy = -0.19794636356747858,
    ),
)); verbosity=4)
octoplot(model, init_chain, show_rv=true)
```

!!! note
    If you don't pick the starting point, you cabn also just run Pigeons for 8-10 rounds, which is recommended anyways for convergence, and the sampler will find this result.

Now, we can perform the fit. It is a little slow since we have many hundreds of RV and astrometry data points.
```@example 1
using Pigeons
chain, pt = octofit_pigeons(model, n_rounds=6) # might need more rounds to converge
```

```julia
increment_n_rounds!(pt,1)
chain,pt = octofit_pigeons(pt)
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

