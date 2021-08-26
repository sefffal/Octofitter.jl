# [Basic Model](@id fit-astrometry)

Here is a worked example of a basic model. It contains a star with a single planet, and several astrometry points.

Start by loading the DirectDetections and Plots packages:
```julia
using DirectDetections, Plots
```

## Creating a planet

Create our first planet. Let's name it planet X.
```julia
@named X = DirectDetections.Planet(
    Priors(
        a = TruncatedNormal(15, 2, 0, Inf),
        e = TruncatedNormal(0.1, 0.1, 0, 0.5),
        τ = Normal(0.5,0.5),
        ω = Normal(1π, 1π),
        i = Normal(1π, 1π),
        Ω = Normal(1π, 1π),
    ),
    Astrometry(
        (epoch=5123.,  ra=531.0, dec=542, σ_ra=2., σ_dec=2.),
        (epoch=5832.,  ra=489.1, dec=752., σ_ra=2., σ_dec=2.),
    )
)
```

There's a lot going on here, so let's break it down.

The `Priors` block accepts the priors that you would like for the orbital parameters of this planet. Priors can be any univariate distribution from the Distributions.jl package.
You will want to always specify the following parameters:
* `a`: Semi-major axis in astronomical units (AU)
* `i`: Inclination in radians
* `e`: Eccentricity in the range [0, 1)
* `τ`: Epoch of periastron passage, in fraction of orbit \[0,1]
* `ω`: Argument of periastron
* `Ω`: Longitude of the ascending node, radians.

The parameter τ represents the epoch of periastron passage as a fraction of the planet's orbit between 0 and 1. This follows the same convention as Orbitize! and you can read more about their choice in ther FAQ.

The parameters can be specified in any order.

The `Astrometry` block is optional. This is where you can list the position of a planet at different epochs if it known. `epoch` is a modified Julian date that the observation was taken. the `ra`, `dec`, `σ_ra`, and `σ_dec` parameters are the position of the planet at that epoch, relative to the star. All values in milliarcseconds (mas).

## Creating a system

A system represents a host star with one or more planets. Properties of the whole system are specified here, like parallax distance and mass of the star. This is also where you will supply data like images and astrometric acceleration in later tutorials, since those don't belong to any planet in particular.

```julia
@named HD82134 = System(
    Priors(
        μ = Normal(1.53, 0.01),
        plx =Normal(24.2, 0.02),
    ),  
    X,
)
```

The `Priors` block works just like it does for planets. Here, the two parameters you must provide are:
* `μ`: Graviataion parameter of the central body, expressed in units of Solar mass.
* `plx`: Distance to the system expressed in milliarcseconds of parallax.

After that, just list any planets that you want orbiting the star. Here, we pass planet X.

You can name the system and planets whatever you like.
NB: the `@named` convenience macro just passes in the name as a keyword argument, e.g. `name=:HD82134`. This makes sure that the variable name matches what gets displayed in the package output, and saved a few keystrokes. (taken from ModellingToolkit.jl)

## Sampling
Great! Now we are ready to draw samples from the posterior.

Start sampling:
```julia
chains, stats = DirectDetections.hmctf(
    HD82134;
    burnin=2_000,
    numwalkers=1,
    numsamples_perwalker=10_000
);
```

You will get an output that looks something like with a progress bar that updates every second or so:
```
┌ Info: Guessing a good starting location by sampling from priors
└   N = 100000
┌ Info: Found good location
│   mapv = 298.85965496281216
│   a =
│    1-element Vector{Float64}:
└     17.526543775440597
Sampling 28%|█████████                      |  ETA: 0:11:32
  iterations:                    422
  n_steps:                       8191
  is_accept:                     true
  acceptance_rate:               0.8361569510892266
  log_density:                   335.34337289077916
  hamiltonian_energy:            -333.56163860438176
  hamiltonian_energy_error:      0.6794974484806744
  max_hamiltonian_energy_error:  1.0270049717873349
  tree_depth:                    13
  numerical_error:               false
  step_size:                     0.00024442197084181946
  nom_step_size:                 0.00024442197084181946
  is_adapt:                      true
  mass_matrix:                   DenseEuclideanMetric(diag=[7.836792441654134e-5, 5.04...])
```

A few things to watch out for: check that you aren't getting many (any, really) `numerical_error=true`. This likely indicates that the priors are too restrictive, and the sampler keeps taking steps outside of their valid range. It could also indicate a problem with DirectDetections, e.g. if the sampler is picking negative eccentricities.

## Analysis
As a first pass, let's plot a sample of orbits drawn from the posterior.

```julia 
using Plots
plotmodel(chains[1], HD82134)
```

## Notes on Hamiltonian Monte Carlo
Traditional Affine Invariant MCMC is supported (similar to the python `emcee` package), but it is recommended that you use Hamiltonian Monte Carlo. This sampling method makes use of derivative information, and is much more efficient. 

Derviatives for a complex model are usualy tedious to code, but DirectDetections uses ForwardDiff.jl to generate them automatically.

When using HMC, only a few chains are necessary. This is in contrast to Affine Invariant MCMC based packages where hundreds or thousands of walkers are required.
One chain should be enough to cover the whole posterior, but you can run a few different chains to make sure each has converged to the same distribution.

Similarily, many fewer samples are required. This is because unlike Affine Invariant MCMC, HMC produces samples that are much less correlated after each step (i.e. the autocorrelation time is much shorter).