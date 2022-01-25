# [Fitting Astrometry](@id fit-astrometry)

Here is a worked example of a basic model. It contains a star with a single planet, and several astrometry points.

The full code is available on [GitHub](https://github.com/sefffal/DirectDetections.jl/examples/basic-example.jl)

Start by loading the DirectDetections and Plots packages:
```julia
using DirectDetections, Distributions, Plots
```

## Creating a planet

Create our first planet. Let's name it planet X.
```julia
@named X = Planet(
    Priors(
        a = TruncatedNormal(1, 0.5, 0, Inf),
        e = TruncatedNormal(0.0, 0.2, 0, 1.0),
        τ = Normal(0.5, 1),
        ω = Normal(deg2rad(250.), deg2rad(80.)),
        i = Normal(deg2rad(20.), deg2rad(10.)),
        Ω = Normal(deg2rad(200.), deg2rad(30.)),
    ),
    Astrometry(
        (epoch=5000.,  ra=-364., dec=-1169., σ_ra=70., σ_dec=30.),
        (epoch=5014.,  ra=-493., dec=-1104., σ_ra=70., σ_dec=30.),
        (epoch=5072.,  ra=-899., dec=-629., σ_ra=10., σ_dec=50.),
    )
)
```

There's a lot going on here, so let's break it down.

The `Priors` block accepts the priors that you would like for the orbital parameters of this planet. Priors can be any univariate distribution from the Distributions.jl package.
You will want to always specify the following parameters:
* `a`: Semi-major axis, astronomical units (AU)
* `i`: Inclination, radius
* `e`: Eccentricity in the range [0, 1)
* `τ`: Epoch of periastron passage, in fraction of orbit \[0,1] (periodic outside these bounds)
* `ω`: Argument of periastron, radius
* `Ω`: Longitude of the ascending node, radians.

The parameter τ represents the epoch of periastron passage as a fraction of the planet's orbit between 0 and 1. This follows the same convention as Orbitize! and you can read more about their choice in ther FAQ.

The parameters can be specified in any order.

The `Astrometry` block is optional. This is where you can list the position of a planet at different epochs if it known. `epoch` is a modified Julian date that the observation was taken. the `ra`, `dec`, `σ_ra`, and `σ_dec` parameters are the position of the planet at that epoch, relative to the star. All values in milliarcseconds (mas).


Many different distributions are supported as priors, including `Uniform`, `LogNormal`, `LogUniform`, `Sine`, and `Beta`. See the section on [Priors](@ref priors) for more information.

## Creating a system

A system represents a host star with one or more planets. Properties of the whole system are specified here, like parallax distance and mass of the star. This is also where you will supply data like images and astrometric acceleration in later tutorials, since those don't belong to any planet in particular.

```julia
@named HD82134 = System(
    Priors(
        μ = Normal(1.0, 0.01),
        plx =Normal(1000.2, 0.02),
    ),  
    X,
)
```

The `Priors` block works just like it does for planets. Here, the two parameters you must provide are:
* `μ`: Gravitational parameter of the central body, expressed in units of Solar mass.
* `plx`: Distance to the system expressed in milliarcseconds of parallax.

After that, just list any planets that you want orbiting the star. Here, we pass planet X.

You can name the system and planets whatever you like.
NB: the `@named` convenience macro just passes in the name as a keyword argument, e.g. `name=:HD82134`. This makes sure that the variable name matches what gets displayed in the package output, and saved a few keystrokes. (taken from ModellingToolkit.jl)

## Sampling
Great! Now we are ready to draw samples from the posterior.

Start sampling:
```julia
chain = DirectDetections.hmc(
    HD82134,
    adaptation =   8_000,
    iterations = 100_000,
);
```

You will get an output that looks something like with a progress bar that updates every second or so:
```
┌ Info: Guessing a good starting location by sampling from priors
└   N = 500000
┌ Info: Found initial stepsize
└   initial_ϵ = 2.44140625e-5
[ Info: Will adapt step size and mass matrix
[ Info: progress logging is enabled globally
Sampling100%|███████████████████████████████| Time: 0:03:43
[ Info: Sampling compete. Building chains.
Sampling report for chain 1:
mean_accept =         0.926384435746952
num_err_frac =        0.0
mean_tree_depth =     7.072923809523809
max_tree_depth_frac = 0.0

Chains MCMC chain (105000×8×1 Array{Float64, 3}):

Iterations        = 1:1:105000
Number of chains  = 1
Samples per chain = 105000
Wall duration     = 215.28 seconds
Compute duration  = 215.28 seconds
parameters        = μ, plx, X[a], X[e], X[τ], X[ω], X[i], X[Ω]

Summary Statistics
  parameters        mean       std   naive_se      mcse          ess      rhat   ess_per_sec 
      Symbol     Float64   Float64    Float64   Float64      Float64   Float64       Float64 

           μ      0.9991    0.0100     0.0000    0.0000   73497.2356    1.0000      341.4109
         plx   1000.2001    0.0201     0.0001    0.0001   71398.1409    1.0000      331.6602
        X[a]      1.1400    0.0073     0.0000    0.0000   73252.1068    1.0000      340.2722
        X[e]      0.1300    0.0715     0.0002    0.0005   22658.1132    1.0005      105.2519
        X[τ]      0.5019    1.0032     0.0031    0.0039   72030.3424    1.0000      334.5969
        X[ω]      3.5913    0.7692     0.0024    0.0057   16317.8944    1.0000       75.8002
        X[i]      0.7151    0.0596     0.0002    0.0003   49773.0387    1.0000      231.2068
        X[Ω]      3.3891    0.2884     0.0009    0.0023   14232.4927    1.0001       66.1131

Quantiles
  parameters        2.5%       25.0%       50.0%       75.0%       97.5% 
      Symbol     Float64     Float64     Float64     Float64     Float64 

           μ      0.9795      0.9924      0.9991      1.0059      1.0188
         plx   1000.1607   1000.1865   1000.2000   1000.2136   1000.2395
        X[a]      1.1257      1.1350      1.1399      1.1449      1.1544
        X[e]      0.0419      0.0793      0.1090      0.1636      0.3154
        X[τ]     -1.4665     -0.1724      0.5005      1.1751      2.4657
        X[ω]      1.8839      3.0722      3.8094      4.2234      4.5469
        X[i]      0.5929      0.6765      0.7172      0.7559      0.8261
        X[Ω]      2.7051      3.2433      3.4295      3.5671      3.8877
```

The sampler will begin by drawing orbits randomly from the priors (50,000 by default). It will then pick the orbit with the highest posterior density as a starting point for HMC adaptation. This recipe is a good way to find a point somewhat close to the typical set. Starting at the global maximum on the other hand, has at times not led to good sampling.

Once complete, the `chain` object will hold the sampler results. Displaying it prints out a summary table like the one shown above.

For a basic model like this, sampling should take less than a minute on a typical laptop.

## Diagnostics
The first thing you should do with your results is check a few diagnostics to make sure the sampler converged as intended.

A few things to watch out for: check that you aren't getting many (any, really) numerical errors (`num_err_frac`). This likely means your priors include values that are not possible/supported like an eccentricity outside [0,1), or a negative semi-major axis. 

One common mistake is to use a distribution like `Normal(10,3)` for semi-major axis. This left hand side of this distribution includes negative values which are not physically possible. A better choice is a `TruncatedNormal(10,3,0,Inf)` distribution.

You may see some warnings during initial step-size adaptation. These are probably nothing to worry about if sampling proceeds normally afterwards.

You should also check the acceptance rate (`mean_accept`) is reasonably high and the mean tree depth (`mean_tree_depth`) is reasonable (~4-8). 
Lower than this and the sampler is taking steps that are too large and encountering a U-turn very quicky. Much larger than this and it might be being too conservative. The default maximum tree depth is 10. It should not average anything close to this value, but occasional high values are okay.

The rate at which the maximum tree depth is hit is `max_tree_dept_frac` and should be quite low. Reaching the maximum tree depth does not invalidate the results, but does lower the sampling efficiency.

Next, you can make a trace plot of different variabes to visually inspect the chain:
```julia
plot(
    chain["X[a]"],
    xlabel="iteration",
    ylabel="semi-major axis (aU)"
)
```
![trace plot](assets/astrometry-trace-plot.png)

And an auto-correlation plot:
```julia
using StatsBase
plot(
    autocor(chain["X[e]", 1:500),
    xlabel="lag",
    ylabel="autocorrelation",
)
```
This plot shows that these samples are not correlated after only above 5 steps. No thinning is necessary.
![autocorrelation plot](assets/astrometry-autocor-plot.png)

To confirm convergence, you may also examine the `rhat` column from chains. This diagnostic approaches 1 as the chains converge and should at the very least equal `1.0` to one significant digit (3 recommended).

Finnaly, if you ran multiple chains (see later tutorials to learn how), you can run 
```julia
using MCMCChains
gelmandiag(chain)
```
As an additional convergence test.

## Analysis
As a first pass, let's plot a sample of orbits drawn from the posterior.

```julia 
using Plots
plotmodel(HD82134)
```
This function draws orbits from the posterior and displays them in a plot. Any astrometry points are overplotted. If other data like astrometric acceleration is provided, additional panels will appear.
![model plot](assets/astrometry-model-plot.png)

## Pair Plot
A very useful visualization of our results is a pair-plot, or corner plot. We can use our PairPlots.jl package for this purpose:
```julia
using Plots, PairPlots
table = (;
    a=         chain["X[a]"],
    e=         chain["X[e]"],
    i=rad2deg.(chain["X[i]"]),
    Ω=rad2deg.(chain["X[Ω]"]),
    ω=rad2deg.(chain["X[ω]"]),
    τ=         chain["X[τ]"],
)
labels=["a", "e", "i", "\\Omega", "\\omega", "\\tau"]
units = ["(au)", "", "(\\degree)", "(\\degree)", "(\\degree)", ""]
PairPlots.corner(table, labels, units)
```
You can read more about the syntax for creating pair plots in the PairPlots.jl documentation page.
[![corner plot](assets/astrometry-corner-plot.png)](assets/astrometry-corner-plot.svg)
In this case, the sampler was able to resolve the complicated degeneracies between eccentricity, the longitude of the ascending node, and argument of periapsis.

## Notes on Hamiltonian Monte Carlo
Traditional Affine Invariant MCMC is supported (similar to the python `emcee` package), but it is recommended that you use Hamiltonian Monte Carlo. This sampling method makes use of derivative information, and is much more efficient. This package by default uses the No U-Turn sampler, as implemented in AdvancedHMC.jl.

Derviatives for a complex model are usualy tedious to code, but DirectDetections uses ForwardDiff.jl to generate them automatically.

When using HMC, only a few chains are necessary. This is in contrast to Affine Invariant MCMC based packages where hundreds or thousands of walkers are required.
One chain should be enough to cover the whole posterior, but you can run a few different chains to make sure each has converged to the same distribution.

Similarily, many fewer samples are required. This is because unlike Affine Invariant MCMC, HMC produces samples that are much less correlated after each step (i.e. the autocorrelation time is much shorter).
