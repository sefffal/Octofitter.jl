# [Fit Gaussian Process](@id fit-rv-gp)


This example shows how to fit a Gaussian process to model stellar activity in RV data. It continues from [Basic RV Fit](@ref fit-rv).


!!! note
    The Celerite backend and the public RV archive loaders live in the extension package OctofitterRadialVelocity. To install it, run 
    `pkg> add OctofitterRadialVelocity`. The AbstractGPs backend needs only core Octofitter plus AbstractGPs.

There are two different GP packages supported by OctofitterRadialVelocity: AbstractGPs, and Celerite. Important note: Celerite.jl does not support Julia 1.0+, so we currently bundle a fork that has been patched to work. When / if Celerite.jl is updated we will switch back to the public package.

!!! note "Computational cost"
    GP models are significantly more computationally expensive than non-GP models. Plan for longer sampling times when using Gaussian processes.

!!! warning "GP models need strictly increasing epochs"
    A Gaussian process is fit to the residuals *at this observation's own epochs,
    in table order*, and two rows at the same epoch make the covariance singular
    for any stationary kernel. This is checked in the constructor and raises a
    clear error; merge duplicated rows, or nudge them apart by the exposure time.
    Observations without a `gaussian_process` still accept duplicated epochs.

!!! note "Correlated noise is available for relative RV too"
    There is one [`RadialVelocityObs`](@ref), so `gaussian_process=` and
    `trend_function=` work whether you point it at the star
    (`target=A, ref=Barycentre`) or at a companion (`target=b, ref=A`).


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
using Pigeons
```

We will pick up from our tutorial [Basic RV Fit](@ref fit-rv) with the data already downloaded and available as a table called `rv_dat`:
```@example 1
rv_file = download("https://raw.githubusercontent.com/California-Planet-Search/radvel/master/example_data/k2-131.txt")
rv_dat_raw = CSV.read(rv_file, DataFrame, delim=' ')
rv_dat = DataFrame();
rv_dat.epoch = jd2mjd.(rv_dat_raw.time)
rv_dat.rv = rv_dat_raw.mnvel
rv_dat.σ_rv = rv_dat_raw.errvel
tels = sort(unique(rv_dat_raw.tel))
nothing # hide
```

The bodies are the same as in the basic fit, so we define them once here and
reuse them for both backends:
```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(0.82, 0.02), lower=0.1) # M⊙ (Baines & Armstrong 2011)
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        # Radial-velocity-only fit: fix the inclination and the node.
        i = pi/2
        Ω = 0.0
        e = 0.0
        ω = 0.0
        # To match RadVel we put the prior on the period directly. `P` is an
        # orbital element, in days.
        P ~ truncated(Normal(0.3693038, 0.0000091), lower=0.0001)
        τ ~ UniformCircular(1.0)
        tp = τ * P + 57782 # reference epoch for τ. Choose an MJD date near your data.
        # Minimum planet mass (really m·sin i). Masses are solar masses, and
        # `mjup` is a plain multiplicative constant.
        mass ~ LogUniform(0.001mjup, 10mjup)
    end
)
nothing # hide
```


## Gaussian Process Fit with AbstractGPs
Let us now add a Gaussian process to model stellar activity. This should improve the fit.

We start by writing a function that creates a Gaussian process kernel from a set of observation parameters.
We will create a quasi-periodic kernel. We provide this function as an argument `gaussian_process` to the likelihood constructor:

```@example 1
using AbstractGPs

gp_explength_mean = 9.5*sqrt(2.) # sqrt(2)*tau in Dai+ 2017 [days]
gp_explength_unc = 1.0*sqrt(2.)
gp_perlength_mean = sqrt(1. /(2. *3.32)) # sqrt(1/(2*gamma)) in Dai+ 2017
gp_perlength_unc = 0.019
gp_per_mean = 9.64 # T_bar in Dai+ 2017 [days]
gp_per_unc = 0.12

quasiperiodic = θ_obs -> GP(
    θ_obs.η_1^2 *
    (SqExponentialKernel() ∘ ScaleTransform(1/(θ_obs.η_2))) *
    (PeriodicKernel(r=[θ_obs.η_4]) ∘ ScaleTransform(1/(θ_obs.η_3)))
)

rvlike_harps = RadialVelocityObs(
    rv_dat[rv_dat_raw.tel .== "harps-n", :];
    target=A, ref=Barycentre,
    name="harps-n",
    gaussian_process = quasiperiodic,
    variables=@variables begin
        offset ~ Normal(-6693,100) # m/s
        jitter ~ LogUniform(0.1,100) # m/s
        # Add priors on GP kernel hyper-parameters.
        η_1 ~ truncated(Normal(25,10),lower=0.1,upper=100)
        # Important: ensure the period and exponential length scales
        # have physically plausible lower and upper limits to avoid poor numerical conditioning
        η_2 ~ truncated(Normal(gp_explength_mean,gp_explength_unc),lower=5,upper=100)
        η_3 ~ truncated(Normal(gp_per_mean,1),lower=2, upper=100)
        η_4 ~ truncated(Normal(gp_perlength_mean,gp_perlength_unc),lower=0.2, upper=10)
    end
)
rvlike_pfs = RadialVelocityObs(
    rv_dat[rv_dat_raw.tel .== "pfs", :];
    target=A, ref=Barycentre,
    name="pfs",
    gaussian_process = quasiperiodic,
    variables=@variables begin
        offset ~ Normal(0,100) # m/s
        jitter ~ LogUniform(0.1,100) # m/s
        η_1 ~ truncated(Normal(25,10),lower=0.1,upper=100)
        η_2 ~ truncated(Normal(gp_explength_mean,gp_explength_unc),lower=5,upper=100)
        η_3 ~ truncated(Normal(gp_per_mean,1),lower=2, upper=100)
        η_4 ~ truncated(Normal(gp_perlength_mean,gp_perlength_unc),lower=0.2, upper=10)
    end
)

## No change to the rest of the model

sys = System(
    name = "k2_131",
    bodies=[A, b],
    observations=[rvlike_harps, rvlike_pfs],
)

model = Octofitter.LogDensityModel(sys)
```

Note that the two instruments do not need to use the same Gaussian process kernels, nor the same hyper parameter names. 

!!! note
    Tip: If you want the instruments to *share* the Gaussian process kernel hyper parameters, move the variables up to the system's `@variables` block, and forward them to the observation variables block e.g. `η_1 = system.η_1`, `η_2 = system.η_2`.

!!! note "Three hooks, if you want a different GP library"
    The AbstractGPs path above works by duck typing through three functions you
    can add methods to for any backend: `Octofitter.gp_condition(gp, epochs, σ²)`,
    `Octofitter.gp_ln_like(fx, residuals)`, and
    `Octofitter.gp_predict(fx, residuals, epochs) -> (mean, var)` (which must
    return the *latent* predictive variance — the held-out point's own white
    noise is added by the caller). `OctofitterRadialVelocity` plugs Celerite in
    exactly this way.


Initialize the starting points, and confirm the data are entered correcly:
```@example 1
init_chain = initialize!(model)
octoplot(model, init_chain)
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


Plot the fit: the RV time series, the residual strip, and a phase-folded panel
for the planet. [`octoplot`](@ref) draws everything the model has;
[`rvplot`](@ref) is the single-draw radial-velocity figure.

The fitted Gaussian process is part of the picture. Its conditioned prediction
is drawn as a band around the model curve, subtracted from the residuals before
they are plotted or folded, and its predictive variance is added into `σ_eff` —
so the residual strip shows what the fit is *left with* rather than the
correlated structure the GP was fitted to absorb, and the whitened z-scores are
standard normal for a model that is working. `gpband=false` turns the band off.

Conditioning a GP once per posterior draw is not free; `ndraws=` is the lever if
a figure is taking too long.

```@example 1
octoplot(model, chain)
```

Some optional tweaks to the appearance:
```@example 1
octoplot(
    model,
    chain,
    N=50,          # only plot 50 samples
    figscale=1.5,  # make it larger
)
```


Corner plot:
```@example 1
octocorner(model, chain, small=true)
```


## Gaussian Process Fit with Celerite

We now demonstrate an approximate quasi-static kernel implemented using Celerite. 
For the class of kernels supported by Celerite, the performance scales
much better with the number of data points. This makes it a good choice
for modelling large RV datasets.

!!! warning
    Make sure that you type `using OctofitterRadialVelocity.Celerite` and not `using Celerite`. 
    Celerite.jl does not support Julia 1.0+, so we currently bundle a fork that has been patched to work. When / if Celerite.jl is updated we will switch back to the public package.



```@example 1
using OctofitterRadialVelocity.Celerite

quasistatic = θ_obs -> Celerite.CeleriteGP(
    Celerite.RealTerm(
        #=log_a=# log(θ_obs.B*(1+θ_obs.C)/(2+θ_obs.C)),
        #=log_c=# log(1/θ_obs.L)
    ) + Celerite.ComplexTerm(
        #=log_a=#  log(θ_obs.B/(2+θ_obs.C)),
        #=log_b=#  -Inf,
        #=log_c=#  log(1/θ_obs.L),
        #=log_d=#  log(2pi/θ_obs.Prot)
    )
)

rvlike_harps_cel = RadialVelocityObs(
    rv_dat[rv_dat_raw.tel .== "harps-n", :];
    target=A, ref=Barycentre,
    name="harps-n",
    gaussian_process = quasistatic,
    variables=@variables begin
        offset ~ Normal(-6693,100) # m/s
        jitter ~ LogUniform(0.1,100) # m/s
        # Add priors on GP kernel hyper-parameters.
        B ~ Uniform(0.00001, 2000000)
        C ~ Uniform(0.00001, 200)
        L ~ Uniform(2, 200)
        Prot ~ Uniform(8.5, 20)
    end
)
rvlike_pfs_cel = RadialVelocityObs(
    rv_dat[rv_dat_raw.tel .== "pfs", :];
    target=A, ref=Barycentre,
    name="pfs",
    gaussian_process = quasistatic,
    variables=@variables begin
        offset ~ Normal(0,100) # m/s
        jitter ~ LogUniform(0.1,100) # m/s
        B ~ Uniform(0.00001, 2000000)
        C ~ Uniform(0.00001, 200)
        L ~ Uniform(2, 200)
        Prot ~ Uniform(8.5, 20)
    end
)

## No change to the rest of the model

sys_cel = System(
    name = "k2_131_celerite",
    bodies=[A, b],
    observations=[rvlike_harps_cel, rvlike_pfs_cel],
)

using DifferentiationInterface
using FiniteDiff
model_cel = Octofitter.LogDensityModel(sys_cel, autodiff=AutoFiniteDiff())
```


The Celerite implementation doesn't support our default autodiff backend (ForwardDiff.jl), so
either (A) switch the gradients over to finite differences, as above, and then sample with
the Pigeons slice sampler, which doesn't require gradients at all, or (B) use Enzyme
autodiff.


Initialize the starting points, and confirm the data are entered correcly:
```@example 1
init_chain_cel = initialize!(model_cel)
octoplot(model_cel, init_chain_cel)
```


```@example 1
chain_cel, pt = octofit_pigeons(model_cel, n_rounds=7)
octoplot(model_cel, chain_cel)
```

!!! tip "Why this fit is tempered and the AbstractGPs one above is not"
    [`octofit_pigeons`](@ref) explores with a slice sampler, which needs only log-density
    *evaluations* — no gradient. That is exactly what the Celerite backend can supply
    cheaply and ForwardDiff cannot supply at all, so tempering sidesteps the autodiff
    problem instead of working around it with finite differences. It also handles the
    multi-modality a quasi-periodic GP kernel introduces: the rotation-period
    hyperparameter has real posterior mass at harmonics and aliases of the true period,
    and those are separated modes.

    The AbstractGPs model earlier on the page differentiates cleanly, so [`octofit`](@ref)
    (HMC) is the cheaper choice there. If you prefer HMC here too, keep the
    `AutoFiniteDiff()` backend, seed it with [`initialize!`](@ref), and run several chains
    from different starting points so a missed period harmonic shows up as disagreement.

!!! note "Cross-validating a GP model"
    `Octofitter.likeobj_from_epoch_subset` works on a GP-backed
    `RadialVelocityObs`: the requested rows become the held-out set and the
    remaining rows are silently retained as the conditioning set, which is what a
    correlated model needs in order to score a point at all. Prediction is
    implemented for the Celerite backend only — the AbstractGPs case raises a
    clear error naming the missing `gp_predict` method.
