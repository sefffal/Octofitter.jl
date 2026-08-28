# [Basic Astrometry Fit](@id fit-astrometry)

Here is a worked example of a one-planet model fit to relative astrometry (positions measured between the planet and the host star). 

Start by loading the Octofitter and Distributions packages:
```@example 1
using Octofitter, Distributions
```

### Creating the bodies

We start by specifying what bodies (stars, planets, etc) we want to model.
In Octofitter, we specify models using a "probabilistic programming language"
Quantities with a `~` are random variables. The distributions on the right hand sides are **priors**. You must specify a proper prior for any quantity which is allowed to vary.
Quantities with an `=` are *derived*: constants, or fixed mathematical functions of the other variables.
You can also apply priors to derived quantities.

The host star is a body like any other --- start by creating it and specifying its mass in solar masses:
```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.2, 0.1), lower=0.1) # [M⊙]
    end
)
```

Now we create another body to represent the planet. We indicate that the planet orbits the star with the `about=A` argument. Then, the orbit parameters in the `@variables` block are describe its orbit versus `A`.

```@example 1
planet_b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass = 0.0             # [M⊙] fix to zero for this example
        a ~ Uniform(0, 100)    # [AU]
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()             # [rad]
        ω ~ UniformCircular()  # [rad]
        Ω ~ UniformCircular()  # [rad]
        θ ~ UniformCircular()  # [rad] position angle at `epoch`
        epoch = 50420.0        # [MJD] reference epoch for θ. Choose a date near your data.
    end
)
nothing # hide
```

**`name`**: Try to give each body a short name consisting only of letters and/or trailing numbers. It is used to name that body's columns in the output chain (`b_a`, `b_e`, ...).

**`about`**: which body (or barycentre) this one orbits. Exactly one body, normally the host star, must be unplaced (no `about` argument). 
`about=(A, b)` would place a third body about the *barycentre* of A and b creating a **Jacobi** chain. This is what you usually want for outer planets. See 
[Jacobi vs. astrocentric](https://sefffal.github.io/PlanetOrbits.jl/dev/hierarchies/#Jacobi-vs.-astrocentric)
in the PlanetOrbits manual for which to pick.

**`variables`**: pass a block of variables: priors, derived quantities, and so on.
You need to provide the `mass`, and a sufficient set of orbital parameters.
Some models also want `flux` or `flux_<band>`.


For orbital elements, you must supply exactly one of these options:
| group | alternatives |
|---|---|
| size | `a` [AU] or `P` [**days**] |
| shape | (`e`, `ω`) or (`secosω`, `sesinω`) or (`ecosω`, `esinω`) |
| phase | `tp` or `M0` + `epoch` or `θ` + `epoch` |
| orientation | `i`, `Ω` |

Here we used `θ` (the planet's position angle on the sky at `epoch`) to fix the phase,
which is usually much better constrained by relative data than the epoch of periastron.

Priors can be any continuous univariate distribution from the Distributions.jl package.
Many are supported, including `Uniform`, `LogNormal`, `LogUniform`, `Sine`, and `Beta`.
See the section on [Priors](@ref priors) for more information.
The variables can be specified in any order.

You can also hardcode a particular value for any parameter if you don't want it to vary. Simply replace eg. `e ~ Uniform(0, 0.999)` with `e = 0.1`.
This `=` syntax works for arbitrary mathematical expressions and even functions.
The `=` syntax also works to access variables from the system level, e.g. `plx = system.plx`.

!!! warning
    You must specify a proper prior for any quantity which is allowed to vary. 
    "Uninformative" priors like `1/x` must be given bounds, and can be specified with `LogUniform(lower, upper)`.

!!! warning
    Make sure that variables like mass and eccentricity can't be negative. You can pass a distribution to `truncated` to prevent this, e.g. `mass ~ truncated(Normal(1, 0.1),lower=0)`.

### Creating the observations


We will create an observation object to contain our relative astrometry data. We can specify this data in several formats. It can be listed in the code or loaded from a file (eg. a CSV file, FITS table, or SQL database). You can use any Julia table object.

```@example 1
astrom_dat_1 = Table(;
    epoch= [50000,  50120, 50240, 50360,50480, 50600, 50720, 50840,], # MJD (days)
    ra   = [-505.764, -502.57, -498.209, -492.678,-485.977, -478.11, -469.08, -458.896,], # mas
    dec  = [-66.9298, -37.4722, -7.92755, 21.6356, 51.1472,  80.5359,  109.729,  138.651, ], # mas
    # Tip! Type this as \sigma + <TAB key>!
    σ_ra = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, ],  # mas
    σ_dec = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, ], # mas
    cor =  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ]
)
astrom_obs_1 = RelAstromObs(astrom_dat_1; target=planet_b, ref=A, name="relastrom")
nothing # hide
```
Each observation must be given a `name`, which is used to label it in plots and the chains.

The arguments `target=planet_b` and  `ref=A` are required. They explain to Octofitter that the relative astrometry positions are measured *between planet b and the star A*. You could just as easily supply observations measured between one planet and a second planet, or an entirely different star.

In Octofitter, `epoch` is always the modified Julian date (measured in days). If you're not sure what this is, you can get started by just putting in arbitrary time offsets measured in days.

In this case, we specified `ra` and `dec` offsets in milliarcseconds. We could instead specify `sep` (projected separation) in milliarcseconds and `pa` in radians. You cannot mix the two formats in a single [`RelAstromObs`](@ref) but you can create two different observation objects, one for each format, and add them both to your model:

```@example 1
astrom_dat_2 = Table(
    epoch = [42000, ], # MJD
    sep = [505.7637580573554, ], # mas
    pa = [deg2rad(24.1), ], # radians
    # Tip! Type this as \sigma + <TAB key>!
    σ_sep = [70, ],
    σ_pa = [deg2rad(10.2), ],
)
astrom_obs_2 = RelAstromObs(astrom_dat_2; target=planet_b, ref=A, name="relastrom2")
nothing # hide
```

!!! note
    Tip: You can load data from a CSV file:

```julia
    using CSV
    astrom_dat = CSV.read("mydata.csv", Table)
```

Every observation says what it observes (`target`) and what it is measured against
(`ref`). Relative astrometry of the planet with respect to the star is `target=planet_b,
ref=A`:


#### Advanced Options
You can group your data in different observation objects, each with their own instrument name. Each group can have its own `platescale`, `northangle`, and astrometric `jitter` variables for modelling instrument-specific systematics.

```@example 1
astrom_obs_1 = RelAstromObs(
    astrom_dat_1;
    target = planet_b,
    ref = A,
    name = "GPI astrom",
    variables = @variables begin
        jitter ~ Uniform(0, 10) # mas [optional]
        northangle ~ Normal(0, deg2rad(1)) # radians of offset [optional]
        platescale ~ truncated(Normal(1, 0.01), lower=0) # 1% relative platescale uncertainty 
    end
)

astrom_obs_2 = RelAstromObs(
    astrom_dat_2;
    target = planet_b,
    ref = A,
    name = "SPHERE astrom",
    variables = @variables begin
        jitter ~ Uniform(0, 10) # mas [optional]
        northangle ~ Normal(0, deg2rad(1)) # radians of offset [optional]
        platescale ~ truncated(Normal(1, 0.01), lower=0) # 1% relative platescale uncertainty 
    end
)
nothing # hide
```

### Creating a system

Now we assemble the bodies and observations into a "system". Properties of the whole
system or reference frame are specified here, like parallax distance.
For multi-planet systems, it makes sense to create shared variables here — for example a single inclination used by two planets.

```@example 1
sys = System(
    name = "Tutoria",
    bodies=[A, planet_b],
    observations=[astrom_obs_1, astrom_obs_2],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

nothing #hide
```

The name of your system will be used for output file names by default. We suggest naming it something like `"PDS70-astrom-model"`.


### Prepare model
We now convert our model into efficient, compiled code:
```@example 1
model = Octofitter.LogDensityModel(sys)
```

This object implements the julia LogDensityProblems.jl interface and can be passed to a wide variety of samplers.


### Initialize starting points for chains

Run the `initialize!` function to find good starting points for the chain. You can provide guesses for parameters if you want to. Body variables are nested under `bodies`:
```julia
init_chain = initialize!(model) # No guesses provided, slower global optimization will be used
```

```@example 1
init_chain = initialize!(model, (;
    plx = 50,
    bodies = (;
        A = (; mass = 1.21),
        b = (;
            a = 10.0,
            e = 0.01,
        )
    )
))
nothing # hide
```

!!! warning
    Never initialize a value on the bounds of the prior. For example, exactly 0.00000 eccentricity is disallowed by the `Uniform(0,1)` prior. 

### Visualize the starting points

Plot the inital values to make sure that they are reasonable, and match your data. This is a great time to confirm that your data were entered in correctly.

```@example 1
using CairoMakie
octoplot(model, init_chain)
```

The starting points for sampling look reasonable!

!!! note
    The return value from `initialize!` is a "variational approximation". You can pass that chain to any function expecting a `chain` argument, like `Octofitter.savechain` or `octocorner`. It gives a very rough approximation of the posterior we expect.

### Sampling
Now we are ready to draw samples from the posterior:
```@example 1
octofit(model, verbosity = 0,iterations=2,adaptation=2,); # hide
chain = octofit(model, iterations=1000)
```

You will get an output that looks something like this with a progress bar that updates every second or so. You can reduce or completely silence the output by reducing the `verbosity` value down to 0 from a default of 2 (or get more info with `verbosity=4`).

Once complete, the `chain` object will hold the posterior samples. Displaying it prints out a summary table like the one shown above.

For a basic model like this with few epochs and well-specified uncertainties, sampling should take less than a minute on a typical laptop.

Sampling can take much longer when you have measurements with very small uncertainties (e.g. VLTI-GRAVITY).

### Diagnostics
The first thing you should do with your results is check a few diagnostics to make sure the sampler converged as intended.

A few things to watch out for: check that you aren't getting many numerical errors (`ratio_divergent_transitions`). 
This likely indicates a problem with your model: either invalid values of one or more parameters are encountered (e.g. the prior on semi-major axis includes negative values) or that there is a region of very high curvature that is failing to sample properly. This latter issue can lead to a bias in your results.

One common mistake is to use a distribution like `Normal(10,3)` for semi-major axis. This left tail of this distribution includes negative values, and our orbit model is not defined for negative semi-major axes. A better choice is a `truncated(Normal(10,3), lower=0.1)` distribution (not including zero, since a=0 is not defined).

Next, you can make a trace plot of different variabes to visually inspect the chain:
```@example 1
using CairoMakie
lines(
    chain["b_a"][:],
    axis=(;
        xlabel="iteration",
        ylabel="semi-major axis (AU)"
    )
)
```

And an auto-correlation plot:
```@example 1
using StatsBase
using CairoMakie
lines(
    autocor(chain["b_e"][:], 1:500),
    axis=(;
        xlabel="lag",
        ylabel="autocorrelation",
    )
)
```
This plot shows that these samples are not correlated after only about 5 iterations. No thinning is necessary.

To confirm convergence, you may also examine the `rhat` column from chains. This diagnostic approaches 1 as the chains converge and should at the very least equal `1.0` to one significant digit (3 recommended).

Finaly, you might consider running multiple chains. Simply run `octofit` multiple times, and store the result in different variables. Then you can combine the chains using `chainscat` and run additional inter-chain convergence diagnostics:
```@example 1
using MCMCChains
chain1 = octofit(model)
chain2 = octofit(model)
chain3 = octofit(model)
merged_chains = chainscat(chain1, chain2, chain3)
gelmandiag(merged_chains)
```

This will check that the means and variances are similar between chains that were initialized at different starting points.

### Analysis
As a first pass, let's plot a sample of orbits drawn from the posterior.
The function `octoplot` is a conveninient way to generate a multi-panel plot of the orbits and every data channel in the model:
```@example 1
using CairoMakie
octoplot(model,merged_chains)
```
This function draws orbits from the posterior and displays them in a plot. Any astrometry points are overplotted. 

You can control how many orbits are drawn, the figure scale, and which panels appear. See [`octoplot`](@ref) for more details.

#### Predicting quantities you did not observe

`octoplot` is not restricted to the data in the model. Ask `channels=` for an
**observable** the fit has no data for at all, and it draws the posterior's
prediction for that quantity instead — the curves alone, with nothing to overlay:

```julia
# What radial velocity would a spectrograph see, given this astrometry?
octoplot(model, merged_chains; channels=radvel, show_sky=false)

# The star's reflex proper motion, likewise.
octoplot(model, merged_chains; channels=pmra, show_sky=false)
```

This is the figure to make when deciding whether a target is worth
spectroscopic or absolute-astrometry time. Note that it needs the companion to
have mass: with `mass = 0.0`, as in the model above, the host's reflex signal is
identically zero. Give `planet_b` a `mass` prior first.

The natural follow-up is to take the observation, keep it *out* of the fit, and
check the prediction against it — see
[Predicting data the model never saw](@ref heldout-prediction), which also shows
how to score the result numerically.

### Pair Plot
A very useful visualization of our results is a pair-plot, or corner plot. We can use the `octocorner` function and our PairPlots.jl package for this purpose:
```@example 1
using CairoMakie
using PairPlots
octocorner(model, merged_chains, small=true)
```
Remove `small=true` to display all variables.

In this case, the sampler was able to resolve the complicated degeneracies between eccentricity, the longitude of the ascending node, and argument of periapsis.


### Working with the fitted orbits

Chain columns are named `<body>_<variable>` for body variables, `<observation>_<variable>`
for observation variables, and are unprefixed for system variables:

```@example 1
sma_planet_b = chain["b_a"]      # a (samples × chains) matrix
nothing # hide
```

To evaluate the orbits themselves, rebuild the whole `PlanetOrbits.System` for a draw and
query it. Every observable takes `(solution, target, reference)`:

```@example 1
posys = construct_system(model, chain, 1)        # draw #1
traj = orbitsolve(posys, [mjd("2025-01-01"), mjd("2030-01-01")])
raoff(traj[1], :b, :A), decoff(traj[1], :b, :A)   # [mas]
```

### Saving your chain

You can save your chain in FITS table format by running:
```julia
Octofitter.savechain("mychain.fits", chain)
```

You can load it back via:
```julia
chain = Octofitter.loadchain("mychain.fits"; model)
```

Passing `model` is recommended: `loadchain` then checks the chain's column names against
the model, and errors with a rename hint instead of silently giving you `missing` values.

### Saving your model

You may choose to save your model so that you can reload it later to make plots, etc:
```@example 1
using Serialization
serialize("model1.jls", model)
```

Which can then be loaded at a later time using:
```julia
using Serialization
using Octofitter # must include all the same imports as your original script
model = deserialize("model1.jls")
```

!!! warning
    Serialized models are only loadable/restorable on the same computer, version of Octofitter, and version of Julia. They are not intended for long-term archiving. For reproducibility, make sure to keep your original model definition script.


### Comparing chains
We can compare two different chains by passing them both to `octocorner`. Let's compare the `init_chain` with the full results from `octofit`:
```@example 1
octocorner(model, chain, init_chain, small=true)
```
