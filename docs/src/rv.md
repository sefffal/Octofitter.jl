# [Fit RV and Astrometry](@id fit-rv-pma)

In this example, we will fit an orbit model to a combination of radial velocity and Hipparcos-GAIA proper motion anomaly for the star $\epsilon$ Eridani. We will use some of the radial velocity data collated in [Mawet et al 2019](https://iopscience.iop.org/article/10.3847/1538-3881/aaef8a).

!!! note
    Radial velocity modelling is supported in Octofitter via the extension package OctofitterRadialVelocity. To install it, run 
    `pkg> add OctofitterRadialVelocity`

Datasets from two different radial velocity insturments are included and modelled together with separate jitters and instrumental offsets.


```@example 1


using Octofitter, OctofitterRadialVelocity, Distributions, PlanetOrbits

gaia_id = 5164707970261890560 


@planet b Visual{KepOrbit} begin
    e = 0
    ω=0.0
    mass ~ Uniform(0, 3)
    a~Uniform(1, 5)
    i~Sine()
    Ω~UniformCircular()
    
    τ ~ UniformCircular(1.0)
    P = √(b.a^3/system.M)
    tp =  b.τ*b.P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
end # No planet astrometry is included since it has not yet been directly detected


# Convert from JD to MJD
# Data tabulated from Mawet et al
jd(mjd) = mjd - 2400000.5

# # You could also specify it manually like so:
# rvs = StarAbsoluteRVLikelihood(
#     (;inst_idx=1, epoch=jd(2455110.97985),  rv=−6.54, σ_rv=1.30),
#     (;inst_idx=1, epoch=jd(2455171.90825),  rv=−3.33, σ_rv=1.09),
#     ...
# )

# We will load in data from two RV
rvharps = OctofitterRadialVelocity.HARPS_RVBank_rvs("HD22049")
rvharps.table.inst_idx .= 1
rvhires = OctofitterRadialVelocity.HIRES_rvs("HD22049")
rvhires.table.inst_idx .= 2
rvs_merged = StarAbsoluteRVLikelihood(
    cat(rvhires.table, rvharps.table,dims=1),
    instrument_names=["HARPS", "HIRES"]
)
scatter(rvs_merged.table.epoch[:], rvs_merged.table.rv[:])


@system ϵEri begin
    M ~ truncated(Normal(0.78, 0.01),lower=0)
    plx ~ gaia_plx(;gaia_id)
    pmra ~ Normal(-975, 10)
    pmdec ~ Normal(20,  10)

    rv0_1 ~ Normal(0,10)
    rv0_2 ~ Normal(0,10)
    jitter_1 ~ truncated(Normal(0,10),lower=0)
    jitter_2 ~ truncated(Normal(0,10),lower=0)
end HGCALikelihood(;gaia_id) rvs_merged b

# Build model
model = Octofitter.LogDensityModel(ϵEri; autodiff=:ForwardDiff, verbosity=4) # defaults are ForwardDiff, and verbosity=0
```

Now sample:
```@example 1
using Random
rng = Random.Xoshiro(0)

results = octofit(rng, model)
```

We can now plot the results with a multi-panel plot:
```@example 1
## Save results plot
octoplot(model, results)
```


We can also plot just the RV fit:
```@example 1
fig = OctofitterRadialVelocity.rvpostplot(model, results)
```

And a corner plot:
```@example 1
using CairoMakie, PairPlots
octocorner(model, results, small=true)
```