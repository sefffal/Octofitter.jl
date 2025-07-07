# [Generating and Fitting Simulated Data](@id data-simulation)

This tutorial demonstrates how to use Octofitter to simulate synthetic data, fit models to that data, and validate the results. This is a crucial workflow for understanding model performance and testing analysis pipelines.

We'll simulate data from a known set of parameters, fit a model to recover those parameters, and compare the posterior to the true values used for simulation.


## Setup

```@example 1
using Octofitter, OctofitterRadialVelocity, Distributions
using CairoMakie, PairPlots, Pigeons
using CSV, DataFrames
```

## Define the Model

In order to simulate data from Octofitter, you start by defining a real model with data. The simulate step will use the provided epochs, data types, and uncertainties when simulating data. If you don't have real data, you can enter in arbitrary values for e.g. delta R.A., but use expected epochs and uncertainties.

For this example, we'll use the combined astrometry, proper motion anomaly, and radial velocity model from the [Astrometry, PMA, and RV](@ref astrom-pma-rv) tutorial to define the epochs and uncertainties.

```@example 1
# HGCA likelihood for HD 91312
hgca_like = HGCALikelihood(gaia_id=756291174721509376)

# Simulated relative astrometry data (from discovery paper)
astrom_dat = Table(;
    epoch = [mjd("2016-12-15"), mjd("2017-03-12"), mjd("2017-03-13"), mjd("2018-02-08"), mjd("2018-11-28"), mjd("2018-12-15")],
    ra    = [133., 126., 127., 083., 058., 056.],
    dec   = [-174., -176., -172., -133., -122., -104.],
    σ_ra  = [07.0, 04.0, 04.0, 10.0, 10.0, 08.0],
    σ_dec = [07.0, 04.0, 04.0, 10.0, 20.0, 08.0],
    cor   = [0.2, 0.3, 0.1, 0.4, 0.3, 0.2]
)

astrom_like = PlanetRelAstromLikelihood(
    astrom_dat,
    name = "SCExAO",
    variables = @variables begin
        jitter = 0
        northangle = 0
        platescale = 1
    end
)

# Simulated RV data
rv_dat = Table(;
    epoch = [mjd("2008-05-01"), mjd("2010-02-15"), mjd("2016-03-01")],
    rv    = [1300, 700, -2700],
    σ_rv  = [150, 150, 150]
)

rvlike = StarAbsoluteRVLikelihood(
    rv_dat,
    name="SOPHIE",
    variables=@variables begin
        jitter ~ truncated(Normal(10, 5), lower=0)
        offset ~ Normal(0, 1000)
    end
)

# Planet model
planet_b = Planet(
    name="b",
    basis=AbsoluteVisual{KepOrbit},
    likelihoods=[ObsPriorAstromONeil2019(astrom_like)],
    variables=@variables begin
        a ~ LogUniform(0.1,400)
        e ~ Uniform(0,0.999)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        mass = system.M_sec
        θ ~ Uniform(0, 2pi)
        M = system.M
        tp = θ_at_epoch_to_tperi(θ, 57737.0; M, e, a, i, ω, Ω)
        F = 0.0
    end
)

# System model
ra = 158.30707896392835
dec = 40.42555422701387

sys = System(
    name="HD91312_simulation",
    companions=[planet_b],
    likelihoods=[hgca_like, rvlike],
    variables=@variables begin
        M_pri ~ truncated(Normal(1.61, 0.1), lower=0.1)
        M_sec ~ LogUniform(0.5, 1000)
        M = M_pri + M_sec*Octofitter.mjup2msol
        plx ~ gaia_plx(gaia_id=756291174721509376)
        pmra ~ Normal(-137, 10)
        pmdec ~ Normal(2, 10)
        ra = $ra
        dec = $dec
        rv = 0*1e3
        ref_epoch = Octofitter.meta_gaia_DR3.ref_epoch_mjd
    end
)

model = Octofitter.LogDensityModel(sys)
nothing # hide
```

## Generate Synthetic Data

We have three choices for generating simulated data:
1. Draw values from the priors
2. Use values from a fitted posterior 
3. Specifying values manually

We will look at each.

### 1. Draw values from the priors

We can draw a value from the priors like so:
```@example 1
params_to_simulate = Octofitter.drawfrompriors(model.system)
```


### 2. Use values from a fitted posterior 

We can select a particular draw from the posterior and use this to generate new data
```julia
# Perform MCMC fitting on real data
using Pigeons
chain_real, pt = octofit_pigeons(model, n_rounds=5)

# Use one particular draw as the basis of our simulation
draw_number = 1
params_to_simulate = Octofitter.mcmcchain2result(model, chain_real, draw_number)
```

### Specifying values manually
We can also specify all values for the simulation manually. This process is a bit more involved. 

Start by drawing parameters from the priors:
```@example 1
template = Octofitter.drawfrompriors(model.system)
```

Copy this output as a template, and replace values as needed. Note that if some parameters are calculated based on others in your model, you will have to repeat those calculations here.

!!! warning
    Note that the output below is just an example, you must generate your own template from your model and modify it as needed. The exact structure is not garuanteed to be stable between versions of Octofitter.

```@example 1
# Define our "true" parameter values for simulation
M_pri = 1.61
M_sec = 85.0
M = M_pri + M_sec*Octofitter.mjup2msol
params_to_simulate = (
    M_pri = M_pri, 
    M_sec = M_sec,  # Jupiter masses
    M = M,
    plx = 21.5,
    pmra = -137.0,
    pmdec = 2.0,
    ra = ra,
    dec = dec,
    rv = 0.0,
    ref_epoch = Octofitter.meta_gaia_DR3.ref_epoch_mjd,
    observations = (
        SOPHIE = (jitter = 15.0, offset = 0.0),
    ),
    planets = (
        b = (
            a = 45.0,
            e = 0.15,
            ω = 1.2,
            mass = M_sec,
            i = 0.8,
            Ω = 2.1,
            θ = 1.5,
            M = M,
            tp = θ_at_epoch_to_tperi(1.5, 57737.0; M=M, e=0.15, a=45.0, i=0.8, ω=1.2, Ω=2.1),
            F = 0.0,
            observations = NamedTuple()
        ),
    )
)
```


## Generate synthetic system with simulated data

```@example 1
sim_system = Octofitter.generate_from_params(model.system, params_to_simulate)
sim_model = Octofitter.LogDensityModel(sim_system)
```

Let's plot the simulated orbit and data to see what we generated:

```@example 1
# Convert true parameters to chain format for plotting
true_chain = Octofitter.result2mcmcchain([params_to_simulate])

# Plot the simulated system with true parameters
fig = octoplot(sim_model, true_chain, colormap=:red)
fig
```

## Fit the Simulated Data

Now we'll sample from the posterior using the simulated data:

```@example 1
# Sample from the simulated data
chain, pt = octofit_pigeons(sim_model, n_rounds=8, explorer=SliceSampler())
display(chain)
```

## Compare Results

Let's visualize how well our sampling recovered the true parameters:

```@example 1
# Plot the posterior from fitting simulated data
fig = octoplot(sim_model, chain)
fig
```

## Overlay True and Recovered Parameters

For a direct comparison, we can overlay the true orbit with the posterior samples:

```@example 1
# Create astrometry plot showing both posterior and true orbit
fig = Octofitter.astromplot(sim_model, chain, use_arcsec=false, ts=1:2)
ax = fig.content[1]

# Add true orbit in red
Octofitter.astromplot!(ax, sim_model, true_chain, use_arcsec=false, ts=1:2, colormap=Makie.cgrad([:red]))

# Make true orbit line more visible
ax.scene.plots[6].linewidth = 6

fig
```

## Corner Plot Comparison

Compare the parameters used to generate the simulated data, and the recovered posterior:

```@example 1
# Create corner plot showing both posterior and true values
octocorner(
    sim_model,
    chain,
    small=true,
    truth=(PairPlots.Truth((;
        M=collect(true_chain[:M][:]),
        b_a=collect(true_chain[:b_a][:]),
        b_e=collect(true_chain[:b_e][:]),
        b_i=collect(true_chain[:b_i][:]),
        b_mass=collect(true_chain[:b_mass][:]),
    ),label="Simulated Orbit", color=:red)=>(
        PairPlots.MarginLines(),
        PairPlots.BodyLines(),
    ),)
)
```
