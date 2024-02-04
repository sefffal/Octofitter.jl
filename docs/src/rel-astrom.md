# [Fit Relative Astrometry](@id fit-astrometry)

Here is a worked example of a one-planet model fit to relative astrometry (positions measured between the planet and the host star). 

Start by loading the Octofitter and Distributions packages:
```@example 1
using Octofitter, Distributions
```

## Specifying the data
We will create a likelihood object to contain our relative astrometry data. We can specify this data in several formats. It can be listed in the code or loaded from a file (eg. a CSV file, FITS table, or SQL database).

```@example 1
astrom_like = PlanetRelAstromLikelihood(
    (epoch = 50000, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50120, ra = -502.570356287689, dec = -37.47217527025044, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50240, ra = -498.2089148883798, dec = -7.927548139010479, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50360, ra = -492.67768482682357, dec = 21.63557115669823, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50480, ra = -485.9770335870402, dec = 51.147204404903704, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50600, ra = -478.1095526888573, dec = 80.53589069730698, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50720, ra = -469.0801731788123, dec = 109.72870493064629, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50840, ra = -458.89628893460525, dec = 138.65128697876773, σ_ra = 10, σ_dec = 10, cor=0),
)
```
In Octofitter, `epoch` is always the modified Julian date (measured in days). If you're not sure what this is, you can get started by just putting in arbitrary time offsets measured in days.

In this case, we specified `ra` and `dec` offsets in milliarcseconds. We could instead specify `sep` (projected separation) in milliarcseconds and `pa` in radians. You cannot mix the two formats in a single `PlanetRelAstromLikelihood` but you can create two different likelihood objects, one for each format.

Another way we could specify the data is by column:
```@example 1
astrom_like = PlanetRelAstromLikelihood(Table(;
    epoch= [
        50000,
        50120,
        50240,
        50360,
        50480,
        50600,
        50720,
        50840,
    ],
    ra = [
        -505.764,
        -502.57,
        -498.209,
        -492.678,
        -485.977,
        -478.11,
        -469.08,
        -458.896,
    ],
    dec = [
        -66.9298,
        -37.4722,
        -7.92755,
        21.6356, 
        51.1472, 
        80.5359, 
        109.729, 
        138.651, 
    ],
    σ_ra = fill(10.0, 8),
    σ_dec = fill(10.0, 8),
    cor = fill(0.0, 8)
))
```

Finally we could also load the data from a file somewhere. Here is an example 
of loading a CSV:
```julia
using CSV # must install CSV.jl first
astrom_like = CSV.read("mydata.csv", PlanetRelAstromLikelihood)
```

You can also pass a DataFrame or any other table format.

## Creating a planet

We now create our first planet model. Let's name it planet `b`. 
The name of the planet will be used in the output results.

In Octofitter, we specify planet and system models using a "probabilistic
programming language". Quantities with a `~` are random variables. The distributions on the right hand sides are **priors**. You must specify a 
proper prior for any quantity which is allowed to vary. 

We now create our planet `b` model using the `@planet` macro.
```@example 1
@planet b Visual{KepOrbit} begin
    a ~ truncated(Normal(10, 4), lower=0, upper=100)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()

    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50420)
end astrom_like
nothing # hide
```

There's a lot going on here, so let's break it down.

First, `Visual{KepOrbit}` is a kind of orbit parameterization from PlanetOrbits.jl that we will use for this. A `KepOrbit` uses the traditional Keplerian orbital elements like semi-major axis and inclination.
The `Visual{...}` part surrounding it adds support for parallax distance and allows a Keplerian orbit to be projected into the plane of the sky, which is where our relative astrometry data lives! A `KepOrbit` by itself can only be used to fit position measurements in astronomical units. The `Visual{...}` part makes it so we can calculate and specify observed angles between the star and planet.

You can read about orbit parameterizations in the [PlanetOrbits.jl documentation page](https://sefffal.github.io/PlanetOrbits.jl/dev/api/).

Other options include the similar `ThieleInnesOrbit` which uses a different parameterization, as well as `RadVelOrbit` which is useful for modelling radial velocity data. These can also be wrapped in `Visual{...}` if we want to e.g. project a `RadVelOrbit` onto the plane of the sky.

After the `begin` comes our variable specification. Quantities with a `~` are random variables. The distribution son the right hand sides are **priors**.You must specify a proper prior for any quantity which is allowed to vary. 
"Uninformative" priors like `1/x` must be given bounds, and can be specified with `LogUniform(lower, upper)`.

Priors can be any univariate distribution from the Distributions.jl package.

For a `KepOrbit` you must specify the following parameters:
* `a`: Semi-major axis, astronomical units (AU)
* `i`: Inclination, radius
* `e`: Eccentricity in the range [0, 1)
* `tp`: Epoch of periastron passage
* `ω`: Argument of periastron, radius
* `Ω`: Longitude of the ascending node, radians.


Many different distributions are supported as priors, including `Uniform`, `LogNormal`, `LogUniform`, `Sine`, and `Beta`. See the section on [Priors](@ref priors) for more information.
The parameters can be specified in any order.

You can also hardcode a particular value for any parameter if you don't want it to vary. Simply replace eg. `e ~ Uniform(0, 0.999)` with `e = 0.1`.
This `=` syntax works for arbitrary mathematical expressions and even functions. We use it here to reparameterize `tp`.

`tp` is a date which sets the location of the planet around its orbit. It repeats every orbital period and the orbital period depends on the scale of the orbit. This makes it quite hard to sample from. We therefore reparameterize using `θ` parameter, representing the position angle of the planet at a given reference epoch. This parameterization speeds up sampling quite a bit!

After the variables block are zero or more `Likelihood` objects. These are observations specific to a given planet that you would like to include in the model. If you would like to sample from the priors only, don't pass in any observations.

For this example, we specify `PlanetRelAstromLikelihood` block. This is where you can list the position of a planet relative to the star at different epochs.

When you have created your planet, you should see the following output. If you don't, you can run `display(b)` to force the text to be output:
```@example 1
b # hide
```


## Creating a system

A system represents a host star with one or more planets. Properties of the whole system are specified here, like parallax distance and mass of the star. This is also where you will supply data like images, astrometric acceleration, or stellar radial velocity since they don't belong to any planet in particular.

```@example 1
@system Tutoria begin
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end b
nothing #hide
```

`Tutoria` is the name we have given to the system. It could be eg `PDS70`, or anything that will help you keep track of the results.

The variables block works just like it does for planets. Here, the two parameters you must provide are:
* `M`: Gravitational parameter of the central body, expressed in units of Solar mass.
* `plx`: Distance to the system expressed in milliarcseconds of parallax.

`M` is always required for all choices of parameterizations supported by PlanetOrbits.jl. `plx` is needed due to our choice to use `Visual{...}` orbit and relative astrometry.
The prior on `plx` can be looked up from GAIA for many targets by using the function `gaia_plx`:
```julia
    plx ~ gaia_plx(;gaia_id=12345678910)
```

After that, just list any planets that you want orbiting the star. Here, we pass planet `b`.

This is also where we could pass likelihood objects for system-wide data like stellar radial velocity.

You can display your system object by running `display(Tutoria)` (or whatever you chose to name your system).


## Prepare model
We now convert our declarative model into efficient, compiled code.
The `autodiff` flag specifies what Julia automatic differentiation package we should use to calculate the gradients of our model.


```@example 1
model = Octofitter.LogDensityModel(Tutoria)
```

This type implements the julia LogDensityProblems.jl interface and can be passed to a wide variety of samplers.

### A note on automatic differentiation
Julia has many packages for calculating the gradients of arbitrary code, and several are supported with Octofitter. `autodiff=:ForwardDiff` is a very robust choice for forward-mode auto-diff, and it works well on most codes and models. `:Enzyme` is a state-of-the-art auto-diff forward and reverse mode package that works with many, but not all models. Give it a try for a free speed-boost, but fall back to `ForwardDiff` if you see an error message or crash.
The `:Zygote` reverse diff package is also partially supported, but usually only of interest when fitting gaussian-process models.
If you absolutely must, you can fallback to `:FiniteDiff` which uses classical finite differencing methods to approximate gradients for un-differentiable code. 

## Sampling
Great! Now we are ready to draw samples from the posterior.

Start sampling:
```@example 1
# Provide a seeded random number generator for reproducibility of this example.
# This is not necessary in general: you may simply omit the RNG parameter if you prefer.
using Random
rng = Random.Xoshiro(1234)
octofit(rng, model, verbosity = 2,iterations=2,adaptation=2,); # hide
chain = octofit(rng, model)
```

You will get an output that looks something like this with a progress bar that updates every second or so. You can reduce or completely silence the output by reducing the `verbosity` value down to 0.


The sampler will begin by drawing orbits randomly from the priors (50,000 by default). It will then pick the orbit with the highest posterior density as a starting point. These are then passed to AdvancedHMC to adapt following the Stan windowed adaption scheme.

Once complete, the `chain` object will hold the sampler results. Displaying it prints out a summary table like the one shown above.

For a basic model like this, sampl]ing should take less than a minute on a typical laptop.

## Diagnostics
The first thing you should do with your results is check a few diagnostics to make sure the sampler converged as intended.

A few things to watch out for: check that you aren't getting many (any, really) numerical errors (`ratio_divergent_transitions`). 
This likely indicates a problem with your model: either invalid values of one or more parameters are encountered (e.g. the prior on semi-major axis includes negative values) or that there is a region of very high curvature that is failing to sample properly. This latter issue can lead to a bias in your results.

One common mistake is to use a distribution like `Normal(10,3)` for semi-major axis. This left hand side of this distribution includes negative values which are not physically possible. A better choice is a `truncated(Normal(10,3), lower=0)` distribution.

You may see some warnings during initial step-size adaptation. These are probably nothing to worry about if sampling proceeds normally afterwards.

You should also check the acceptance rate (`mean_accept`) is reasonably high and the mean tree depth (`mean_tree_depth`) is reasonable (~4-8). 
Lower than this and the sampler is taking steps that are too large and encountering a U-turn very quicky. Much larger than this and it might be being too conservative. The default maximum tree depth is 10. These don't affect the accuracy of your results, but do affect the efficiency of the sampling.

Next, you can make a trace plot of different variabes to visually inspect the chain:
```@example 1
using Plots: Plots
Plots.plot(
    chain["b_a"],
    xlabel="iteration",
    ylabel="semi-major axis (AU)"
)
```

And an auto-correlation plot:
```@example 1
using StatsBase
Plots.plot(
    autocor(chain["b_e"], 1:500),
    xlabel="lag",
    ylabel="autocorrelation",
)
```
This plot shows that these samples are not correlated after only above 5 steps. No thinning is necessary.

To confirm convergence, you may also examine the `rhat` column from chains. This diagnostic approaches 1 as the chains converge and should at the very least equal `1.0` to one significant digit (3 recommended).

Finaly, if you ran multiple chains (see later tutorials to learn how), you can run 
```julia
using MCMCChains
gelmandiag(chain)
```
As an additional convergence test.

## Analysis
As a first pass, let's plot a sample of orbits drawn from the posterior.

```@example 1
using Plots: Plots
plotchains(chain, :b, kind=:astrometry, color="b_a")
```
This function draws orbits from the posterior and displays them in a plot. Any astrometry points are overplotted. 

We can overplot the astrometry data like so:
```@example 1
Plots.plot!(astrom_like, label="astrometry")
```

The function `octoplot` is a conveninient way to generate a 9-panel plot of velocities and position:
```@example 1
octoplot(model, chain)
```


## Pair Plot
A very useful visualization of our results is a pair-plot, or corner plot. We can use the `octocorner` function and our PairPlots.jl package for this purpose:
```@example 1
using CairoMakie: Makie
using PairPlots
octocorner(model, chain, small=true)
```
Remove `small=true` to display all variables, or run `pairplot(chain)` to include internal variables.

In this case, the sampler was able to resolve the complicated degeneracies between eccentricity, the longitude of the ascending node, and argument of periapsis.

## Notes on Hamiltonian Monte Carlo
Unlike most other astrometry modelling code, Octofitter uses Hamiltonian Monte Carlo instead of Affine Invariant MCMC (e.g. emcee in Python) This sampling method makes use of derivative information, and is much more efficient. This package by default uses the No U-Turn sampler, as implemented in AdvancedHMC.jl.

Derviatives for a complex model are usualy tedious to code, but Octofitter uses ForwardDiff.jl to generate them automatically.

When using HMC, only a few chains are necessary. This is in contrast to Affine Invariant MCMC based packages where hundreds or thousands of walkers are required.
One chain should be enough to cover the whole posterior, but you can run a few different chains to make sure each has converged to the same distribution.

Similarily, many fewer samples are required. This is because unlike Affine Invariant MCMC, HMC produces samples that are much less correlated after each step (i.e. the autocorrelation time is much shorter).
