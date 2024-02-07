# [Fit Radial Velocity](@id fit-rv)

You can use Octofitter to fit radial velocity data, either alone or in combination with other kinds of data.
Multiple instruments (up to 10) are supported, as are gaussian processes to model stellar activity..

!!! note
    Radial velocity modelling is supported in Octofitter via the extension package OctofitterRadialVelocity. To install it, run 
    `pkg> add OctofitterRadialVelocity`

For this example, we will fit the orbit of the planet K2-131 to perform the same fit as in the RadVel [Gaussian Process Fitting](https://radvel.readthedocs.io/en/latest/tutorials/GaussianProcess-tutorial.html) tutorial.


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
    (;inst_idx=1, epoch=jd(2455110.97985),  rv=−6.54, σ_rv=1.30),
    (;inst_idx=1, epoch=jd(2455171.90825),  rv=−3.33, σ_rv=1.09),
    # ... etc.
)
```

For this example, to replicate the results of RadVel, we will download their example data for K2-131.
and format it for Octofitter:
```@example 1
rv_file = download("https://raw.githubusercontent.com/California-Planet-Search/radvel/master/example_data/k2-131.txt")
rv_dat_raw = CSV.read(rv_file, DataFrame, delim=' ')
rv_dat = DataFrame();
rv_dat.epoch = jd2mjd.(rv_dat_raw.time)
rv_dat.rv = rv_dat_raw.mnvel
rv_dat.σ_rv = rv_dat_raw.errvel
tels = sort(unique(rv_dat_raw.tel))
rv_dat.inst_idx = map(rv_dat_raw.tel) do tel
    findfirst(==(tel), tels)
end
rvlike = StarAbsoluteRVLikelihood(
    rv_dat,
    instrument_names=["harps-n", "psf"],
)
```


Now, create a planet. We can use the [`RadialVelocityOrbit`](https://sefffal.github.io/PlanetOrbits.jl/dev/api/#Required-Parameters) type from PlanetOrbits.jl that requires fewer parameters (eg no inclination or longitude of ascending node). We could instead use a `Visual{KepOrbit}` or similar
if we wanted to include these parameters and visualize the orbit in the plane of the sky.


```@example 1


@planet b RadialVelocityOrbit begin
    e = 0
    ω = 0.0

    # To match RadVel, we set a prior on Period and calculate semi-major axis from it
    P ~ truncated(Normal(0.3693038/365.25, 0.0000091/365.25),lower=0.00001)
    a = cbrt(system.M * b.P^2) # note the equals sign. 

    τ ~ UniformCircular(1.0)
    tp =  b.τ*b.P*365.25 + 57782 # reference epoch for τ. Choose an MJD date near your data.
    
    # planet mass [jupiter masses]
    mass ~ LogUniform(0.001, 10)
end


@system k2_132 begin
    # total mass [solar masses]
    M ~ truncated(Normal(0.82, 0.02),lower=0) # (Baines & Armstrong 2011).

    # HARPS-N
    rv0_1 ~ Normal(0,10000) # m/s
    jitter_1 ~ LogUniform(0.01,100) # m/s

    # FPS
    rv0_2 ~ Normal(0,10000) # m/s
    jitter_2 ~ LogUniform(0.01,100) # m/s
end rvlike b

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

Excellent! Let's plot the maximum likelihood orbit:
```@example 1
using CairoMakie: Makie
fig = OctofitterRadialVelocity.rvpostplot(model, chain) # saved to "k2_132-rvpostplot.png"
```


Create a corner plot:
```@example 1
using PairPlots
using CairoMakie: Makie
octocorner(model, chain, small=true)
```



## Gaussian Process Fit
Let us now add a Gaussian process to model stellar activity. This should improve the fit.

We start by writing a functino that creates a Gaussian process kernel from a set of system 
parameters. We will create a quasi-periodic kernel.
```@example 1
using AbstractGPs

function gauss_proc_quasi_periodic(θ_system)
    # θ_system is a named tuple of variable values for this posterior draw.
    η₁ = θ_system.gp_η₁ # h
    η₂ = θ_system.gp_η₂ # λ
    η₃ = θ_system.gp_η₃ # θ
    η₄ = θ_system.gp_η₄ # ω
    kernel = η₁^2 *  
        (SqExponentialKernel() ∘ ScaleTransform(1/(η₂))) *
        (PeriodicKernel(r=[η₄]) ∘ ScaleTransform(1/(η₃)))
    gp = GP(kernel)
    return gp
end
```

Now provide this function when we construct our likelihood object:
```@example 1
rvlike = StarAbsoluteRVLikelihood(rv_dat,
    instrument_names=["harps-n", "psf"],
    gaussian_process = gauss_proc_quasi_periodic
)
```

We create a new system model with added parameters `gp_η₁` etc for the Gaussian process fit.

```@example 1

gp_explength_mean = 9.5*sqrt(2.) # sqrt(2)*tau in Dai+ 2017 [days]
gp_explength_unc = 1.0*sqrt(2.)
gp_perlength_mean = sqrt(1. /(2. *3.32)) # sqrt(1/(2*gamma)) in Dai+ 2017
gp_perlength_unc = 0.019
gp_per_mean = 9.64 # T_bar in Dai+ 2017 [days]
gp_per_unc = 0.12

# No change to planet model
@planet b RadialVelocityOrbit begin
    e = 0
    ω = 0.0
    P ~ truncated(Normal(0.3693038/365.25, 0.0000091/365.25),lower=0.00001)
    a = cbrt(system.M * b.P^2)
    mass ~ LogUniform(0.001, 10)
    
    τ ~ UniformCircular(1.0)
    tp =  b.τ*b.P*365.25 + 57782 # reference epoch for τ. Choose an MJD date near your data.
end


@system k2_132 begin
    M ~ truncated(Normal(0.82, 0.02),lower=0) # (Baines & Armstrong 2011).

    # HARPS-N
    rv0_1 ~ Normal(0,10000)
    jitter_1 ~ LogUniform(0.01,10)
    # FPS
    rv0_2 ~ Normal(0,10000)
    jitter_2 ~ LogUniform(0.01,10)

    # Add priors on GP kernel hyper-parameters.
    gp_η₁ ~ truncated(Normal(25,10),lower=0)
    gp_η₂ ~ truncated(Normal(gp_explength_mean,gp_explength_unc),lower=0)
    gp_η₃ ~ truncated(Normal(gp_per_mean,1),lower=0)
    gp_η₄ ~ truncated(Normal(gp_perlength_mean,gp_perlength_unc),lower=0)

end rvlike b

model = Octofitter.LogDensityModel(k2_132; autodiff=:ForwardDiff)
```

Sample from the model as before:
```@example 1
using Random
rng = Random.Xoshiro(0)

chain = octofit(
    rng, model,
    adaptation = 100,
    iterations = 100,
)
```
For real data, we would want to increase the adaptation and iterations to about 1000 each.


Plot the results:
```@example 1
fig = OctofitterRadialVelocity.rvpostplot(model, chain) # saved to "k2_132-rvpostplot.png"
```


Corner plot:
```@example 1
octocorner(model, chain, small=true) # saved to "k2_132-pairplot-small.png"
```

It is expensive (per iteration) to fit Guassian processes. If you have many cores (or a cluster) available you might also try starting julia with multiple threads (`julia --threads=auto`) and using `octofit_pigeons`. Pigeons is less efficient with single-core, but has better parallelization scaling:
```julia
using Pigeons
chain, pt = octofit_pigeons(
    model,
    n_rounds=10
)
display(chain)
```

