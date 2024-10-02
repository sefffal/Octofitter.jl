# [Fit RV and Proper Motion Anomaly](@id fit-rv-pma)

In this example, we will fit an orbit model to a combination of radial velocity and Hipparcos-GAIA proper motion anomaly for the star $\epsilon$ Eridani. We will use some of the radial velocity data collated in [Mawet et al 2019](https://iopscience.iop.org/article/10.3847/1538-3881/aaef8a).

!!! note
    Radial velocity modelling is supported in Octofitter via the extension package OctofitterRadialVelocity. To install it, run 
    `pkg> add OctofitterRadialVelocity`

Datasets from two different radial velocity insturments are included and modelled together with separate jitters and instrumental offsets.


```@example 1


using Octofitter, OctofitterRadialVelocity, Distributions, PlanetOrbits, CairoMakie

gaia_id = 5164707970261890560 


@planet b Visual{KepOrbit} begin
    # For speed of example, we are fitting a circular orbit only.s
    e = 0
    ω=0.0
    mass ~ Uniform(0, 3)
    a ~ truncated(Normal(3.48,0.1),lower=0)
    i ~ Sine()
    Ω ~ UniformCircular()
    
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
#     (;inst_idx=1, epoch=jd(2455171.90825),  rv=−3.33, σ_rv=1.09), # units in meters/second
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
scatter(
    rvs_merged.table.epoch[:],
    rvs_merged.table.rv[:],
    axis=(
        xlabel="time [mjd]",
        ylabel="RV [m/s]",
    )
)
```

We load the HGCA data for this target:
```@example
hgca_like = HGCALikelihood(;gaia_id, N_ave=1)
```
In the interests of time, we set `N_ave=1` to speed up the computation. This parameter controls how the model smears out the simulated Gaia and Hipparcos measurements. For a real target, leave it at the default value once you have completed testing.


```@example 1
@system ϵEri begin
    M ~ truncated(Normal(0.78, 0.01),lower=0)
    plx ~ gaia_plx(;gaia_id)
    pmra ~ Normal(-975, 10)
    pmdec ~ Normal(20,  10)

    # Radial velocity zero point per instrument
    rv0 ~ Product([
        Normal(0,10), # m/s
        Normal(0,10),
    ])
    # Radial jitter per instrument. 
    jitter ~ Product([
        truncated(Normal(0,10),lower=0), # m/s
        truncated(Normal(0,10),lower=0),
    ])
    # You can also set both instruments to the same jitter, eg by putting instead (with = not ~):
    # jitter_2 = system.jitter_1
end hgca_like rvs_merged b

# Build model
model = Octofitter.LogDensityModel(ϵEri)
```


Now sample. You could use HMC via `octofit` or tempered sampling via `octofit_pigeons`. When using tempered sampling, make sure to start julia with `julia --thread=auto`. Each additional round doubles the number of posterior samples, so `n_rounds=10` gives 1024 samples. You should adjust `n_chains` to be roughly double the `Λ` value printed out during sample, and `n_chains_variational` to be roughly double the `Λ_var` column. 
```@example 1
using Pigeons
results, pt = octofit_pigeons(model, n_rounds=8, n_chains=10, n_chains_variational=10);
nothing # hide
```

We can now plot the results with a multi-panel plot:
```@example 1
octoplot(model, results[1:200,:,:], show_mass=true)
```


We can also plot just the RV curve from the maximum *a-posteriori* fit:
```@example 1
fig = OctofitterRadialVelocity.rvpostplot(model, results)
```

We can see what the visual orbit looks like for the maximum a-posteriori orbit:
```@example 1
i_max = argmax(results[:logpost][:])
fig = octoplot(
    model,
    results[i_max,:,:],
    # change the colour map a bit:
    colormap=Makie.cgrad([Makie.wong_colors()[1], "#FAFAFA"]),
    show_astrom=true,
    show_astrom_time=false,
    show_rv=false,
    show_hgca=false,
    mark_epochs_mjd=[
        mjd("2030")
    ]
)
Label(fig[0,1], "Maximum a-posteriori orbit sample")
fig
```


And a corner plot:
```@example 1
using CairoMakie, PairPlots
octocorner(model, results, small=true)
```