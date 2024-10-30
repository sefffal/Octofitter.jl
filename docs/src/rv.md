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
    a ~ Uniform(3, 10)
    i ~ Sine()
    Ω ~ Uniform(0,2pi)
    τ ~ Uniform(0,1.0)
    P = √(b.a^3/system.M)
    tp =  b.τ*b.P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
end # No planet astrometry is included since it has not yet been directly detected


# We will load in data from one RV instruments.
# We use `MarginalizedStarAbsoluteRVLikelihood` instead of 
# `StarAbsoluteRVLikelihood` to automatically marginalize out
# the radial velocity zero point of each instrument, saving one parameter.
hires_data = OctofitterRadialVelocity.HIRES_rvs("HD22049")
rvlike_hires = MarginalizedStarAbsoluteRVLikelihood(
    hires_data,
    instrument_name="HIRES",
    jitter=:jitter_hires,
)
```

We load the HGCA data for this target:
```@example 1
hgca_like = HGCALikelihood(;gaia_id, N_ave=1)
```
In the interests of time, we set `N_ave=1` to speed up the computation. This parameter controls how the model smears out the simulated Gaia and Hipparcos measurements. For a real target, leave it at the default value once you have completed testing.


```@example 1
@system ϵEri begin
    M ~ truncated(Normal(0.82, 0.02),lower=0.5, upper=1.5) # (Baines & Armstrong 2011).
    plx ~ gaia_plx(;gaia_id)
    pmra ~ Normal(-975, 10)
    pmdec ~ Normal(20,  10)

    # Jitter per instrument
    jitter_hires ~ LogUniform(0.1, 100) # m/s

end hgca_like rvlike_hires b
# Build model
model = Octofitter.LogDensityModel(ϵEri)
```


Now sample. You could use HMC via `octofit` or tempered sampling via `octofit_pigeons`. When using tempered sampling, make sure to start julia with `julia --thread=auto`. Each additional round doubles the number of posterior samples, so `n_rounds=10` gives 1024 samples. You should adjust `n_chains` to be roughly double the `Λ` value printed out during sample, and `n_chains_variational` to be roughly double the `Λ_var` column. 
```@example 1
using Pigeons
results, pt = octofit_pigeons(model, n_rounds=10, n_chains=10, n_chains_variational=0, explorer=SliceSampler());
nothing # hide
```

We can now plot the results with a multi-panel plot:
```@example 1
octoplot(model, results, show_mass=true)
```


We can also plot just the RV curve from the maximum *a-posteriori* fit:
```@example 1
fig = Octofitter.rvpostplot(model, results)
```

We can see what the visual orbit looks like for the maximum a-posteriori sample (note, we would need to run an optimizer to get the true MAP value; this is just the MCMC sample with higest posterior density):
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
        mjd("2037")
    ]
)
Label(fig[0,1], "Maximum a-posteriori orbit sample")
Makie.resize_to_layout!(fig)
fig
```


And a corner plot:
```@example 1
using CairoMakie, PairPlots
octocorner(model, results, small=true)
```