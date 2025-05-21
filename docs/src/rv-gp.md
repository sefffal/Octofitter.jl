# [Fit Gaussian Process](@id fit-rv-gp)


This example shows how to fit a Gaussian process to model stellar activity in RV data. It continues from [Basic RV Fit](@ref fit-rv).


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


If you would like to manually specify RV data, use the following format:
```@example 1
rv_like = StarAbsoluteRVLikelihood(
    # epoch is in units of MJD. `jd` is a helper function to convert.
    # you can also put `years2mjd(2016.1231)`.
    # rv and σ_rv are in units of meters/second
    (epoch=jd2mjd(2455110.97985),  rv=−6.54, σ_rv=1.30),
    (epoch=jd2mjd(2455171.90825),  rv=−3.33, σ_rv=1.09),
    # ... etc.
    instrument_name="insert name here",
    offset=:rv0_harps,
    jitter=:jitter_harps,
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

## Gaussian Process Fit
Let us now add a Gaussian process to model stellar activity. This should improve the fit.

We start by writing a function that creates a Gaussian process kernel from a set of system parameters.
We will create a quasi-periodic kernel.
```@example 1
using AbstractGPs

function gauss_proc_quasi_periodic(θ_system)
    # θ_system is a named tuple of variable values for this posterior draw.
    # Note: for multiple instruments, you can supply different gaussian process
    # functions that e.g. use different hyper-parameter names.
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
rvlike_harps = StarAbsoluteRVLikelihood(
    rv_dat[rv_dat_raw.tel .== "harps-n",:],
    instrument_name="harps-n",
    offset=:rv0_harps,
    jitter=:jitter_harps,
    # NEW:
    gaussian_process = gauss_proc_quasi_periodic
)
rvlike_pfs = StarAbsoluteRVLikelihood(
    rv_dat[rv_dat_raw.tel .== "pfs",:],
    instrument_name="pfs",
    offset=:rv0_pfs,
    jitter=:jitter_pfs,
    # NEW:
    gaussian_process = gauss_proc_quasi_periodic
)
```

Note that the two instruments do not need to use the same Gaussian process kernels, nor the same hyper parameter names.

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


@system k2_132_gp begin
    M ~ truncated(Normal(0.82, 0.02),lower=0.1) # (Baines & Armstrong 2011).

    rv0_harps ~ Normal(-6693,100) # m/s
    rv0_pfs ~ Normal(0,100) # m/s

    jitter_harps ~ LogUniform(0.1,100) # m/s
    jitter_pfs ~ LogUniform(0.1,100) # m/s

    # Add priors on GP kernel hyper-parameters.
    gp_η₁ ~ truncated(Normal(25,10),lower=0.1,upper=100)
    # Important: ensure the period and exponential length scales
    # have physically plausible lower and upper limits to avoid poor numerical conditioning
    gp_η₂ ~ truncated(Normal(gp_explength_mean,gp_explength_unc),lower=5,upper=100)
    gp_η₃ ~ truncated(Normal(gp_per_mean,1),lower=2, upper=100)
    gp_η₄ ~ truncated(Normal(gp_perlength_mean,gp_perlength_unc),lower=0.2, upper=10)

end rvlike_harps rvlike_pfs b

model = Octofitter.LogDensityModel(k2_132_gp)
```

Initialize the starting points, and confirm the data are entered correcly:
```@example 1
init_chain = initialize!(model)
fig = Octofitter.rvpostplot(model, init_chain)
```


Sample from the model using MCMC (the no U-turn sampler)
```@example 1
# Seed the random number generator
using Random
rng = Random.Xoshiro(0)

chain = octofit(
    rng, model,
    adaptation = 100,
    iterations = 100,
)
```
For real data, we would want to increase the adaptation and iterations to about 1000 each.


Plot one sample from the results:
```@example 1
fig = Octofitter.rvpostplot(model, chain) # saved to "k2_132-rvpostplot.png"
```

Plot many samples from the results:
```@example 1
fig = octoplot(
    model,
    chain,
    # Some optional tweaks to the appearance:
    N=50, # only plot 50 samples
    figscale=1.5, # make it larger
    alpha=0.05, # make each sample more transparent
    colormap="#0072b2",
) # saved to "k2_132-plot-grid.png"
```


Corner plot:
```@example 1
octocorner(model, chain, small=true) # saved to "k2_132-pairplot-small.png"
```
