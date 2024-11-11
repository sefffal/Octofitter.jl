# [Basic RV Fit](@id fit-rv)

You can use Octofitter to fit radial velocity data, either alone or in combination with other kinds of data.
Multiple instruments (any number) are supported, as are arbitrary trends, and gaussian processes to model stellar activity.

!!! note
    Radial velocity modelling is supported in Octofitter via the extension package OctofitterRadialVelocity. To install it, run 
    `pkg> add OctofitterRadialVelocity`

For this example, we will fit the orbit of the planet K2-131, and reproduce this [RadVel tutorial](https://radvel.readthedocs.io/en/latest/tutorials/GaussianProcess-tutorial.html).


We will use the following packages:
```@example 1
using Octofitter
using OctofitterRadialVelocity
using PlanetOrbits
using CairoMakie
using PairPlots
using CSV
using DataFrames
using Distributions
```

We will start by downloading and preparing a table of radial velocity measurements, and create a [`StarAbsoluteRVLikelihood`](@ref) object to hold them.


The following functions allow you to directly load data from various public RV databases:
* `HARPS_DR1_rvs("star-name")`
* `HARPS_RVBank_observations("star-name")`
* `Lick_rvs("star-name")`
* `HIRES_rvs("star-name")`

## Basic Fit
Start by fitting a 1-planet Keplerian model.

Calling those functions with the name of a star will return a [`StarAbsoluteRVLikelihood`](@ref) table. 

If you would like to manually specify RV data, use the following format:
```julia
rv_like = StarAbsoluteRVLikelihood(
    # epoch is in units of MJD. `jd` is a helper function to convert.
    # you can also put `years2mjd(2016.1231)`.
    # rv and σ_rv are in units of meters/second
    (epoch=jd(2455110.97985),  rv=−6.54, σ_rv=1.30),
    (epoch=jd(2455171.90825),  rv=−3.33, σ_rv=1.09),
    # ... etc.
    instrument_name="insert name here"
)
```

For this example, to replicate the results of RadVel, we will download their example data for K2-131 and format it for Octofitter:
```@example 1
rv_file = download("https://raw.githubusercontent.com/California-Planet-Search/radvel/master/example_data/k2-131.txt")
rv_dat_raw = CSV.read(rv_file, DataFrame, delim=' ')
rv_dat = DataFrame();
rv_dat.epoch = jd2mjd.(rv_dat_raw.time)
rv_dat.rv = rv_dat_raw.mnvel
rv_dat.σ_rv = rv_dat_raw.errvel
tels = sort(unique(rv_dat_raw.tel))

# This table includes data from two insturments. We create a separate
# likelihood object for each:
rvlike_harps = StarAbsoluteRVLikelihood(
    rv_dat[rv_dat_raw.tel .== "harps-n",:],
    instrument_name="harps-n",
    offset=:rv0_harps,
    jitter=:jitter_harps,
)
rvlike_pfs = StarAbsoluteRVLikelihood(
    rv_dat[rv_dat_raw.tel .== "pfs",:],
    instrument_name="pfs",
    offset=:rv0_pfs,
    jitter=:jitter_pfs,
)
```


Now, create a planet. We can use the [`RadialVelocityOrbit`](https://sefffal.github.io/PlanetOrbits.jl/dev/api/#Required-Parameters) type from PlanetOrbits.jl that requires fewer parameters (eg no inclination or longitude of ascending node). We could instead use a `Visual{KepOrbit}` or similar
if we wanted to include these parameters and visualize the orbit in the plane of the sky.


```@example 1
@planet b RadialVelocityOrbit begin
    e = 0
    ω = 0.0

    # To match RadVel, we set a prior on Period and calculate semi-major axis from it
    P ~ truncated(
        Normal(0.3693038/365.256360417, 0.0000091/365.256360417),
        lower=0.0001
    )
    a = cbrt(system.M * b.P^2) # note the equals sign. 

    τ ~ UniformCircular(1.0)
    tp = b.τ*b.P*365.256360417 + 57782 # reference epoch for τ. Choose an MJD date near your data.
    
    # minimum planet mass [jupiter masses]. really m*sin(i)
    mass ~ LogUniform(0.001, 10)
end


@system k2_132 begin
    # total mass [solar masses]
    M ~ truncated(Normal(0.82, 0.02),lower=0.1) # (Baines & Armstrong 2011).

    rv0_harps ~ Normal(-6693,100) # m/s
    rv0_pfs ~ Normal(0,100) # m/s

    jitter_harps ~ LogUniform(0.1,100) # m/s
    jitter_pfs ~ LogUniform(0.1,100) # m/s

end rvlike_harps rvlike_pfs b

```

Note how the `rvlike` object was attached to the `k2_132` system instead of the planet. This is because
the observed radial velocity is of the star, and is caused by any/all orbiting planets.

The `rv0` and `jitter` parameters specify priors for the instrument-specific offset and white noise jitter standard deviation. The `_i` index matches the `inst_idx` used to create the observation table.

Note also here that the `mass` variable is really `msini`, or the minimum mass of the planet.

We can now prepare our model for sampling.
```@example 1
model = Octofitter.LogDensityModel(k2_132)
```

Sample:
```@example 1
using Random
rng = Random.Xoshiro(0)

chain = octofit(rng, model)
```

Excellent! Let's plot an orbit sampled from the posterior:
```@example 1
using CairoMakie: Makie
fig = Octofitter.rvpostplot(model, chain) # saved to "k2_132-rvpostplot.png"
```

We can also plot many number of samples from the posterior:
```@example 1
using CairoMakie: Makie
octoplot(model, chain)
```


And create a corner plot:
```@example 1
using PairPlots
using CairoMakie: Makie
octocorner(model, chain)
```

This example continues in [Fit Gaussian Process](@ref fit-rv-gp).
