# [Quick Start](@id quick-start)

This guide introduces the key concepts in Octofitter:
* Observation objects to hold your data
* `Body` and `System` models to specify variables, priors, and system architecture
* Sampling from the posterior using MCMC
* Plotting the results
* Saving the chain

For installation instructions, see [Installation](@ref install).

If you are porting a script written for Octofitter v1 (v8 or earlier), read
[Migrating to Octofitter v2](@ref v2-migration) first — the model syntax changed.

## Example: Fit a Single Planet Orbit to Relative Astrometry

Load the required packages:
```@example 1
using Octofitter, Distributions, CairoMakie, PairPlots
```

Create a [`RelAstromObs`](@ref) object containing your observational data. In this case
it is the position of the planet relative to the star, but many other kinds of data are
supported:
```@example 1
astrom_dat = Table(
    epoch = [50000.0, 50120.0, 50240.0],  # Dates in MJD
    ra    = [-505.7, -502.5, -498.2],     # [mas] East positive
    dec   = [-66.9, -37.4, -7.9],         # [mas] North positive
    σ_ra  = [10.0, 10.0, 10.0],           # [mas] Uncertainties
    σ_dec = [10.0, 10.0, 10.0],           # [mas] Uncertainties
    cor   = [0.0, 0.0, 0.0]               # RA/Dec correlations
)
nothing # hide
```

A model is a flat list of **bodies** and a flat list of **observations**. The host star is
a body like any other, so define it first:
```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.2, 0.1), lower=0.1)  # [M⊙]
    end
)
nothing # hide
```

Now the planet. `about=A` says the planet orbits the star; its orbital elements and their
[prior distributions](@ref priors) go in its own variables block:
```@example 1
b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass = 0.0                 # [M⊙] relative astrometry alone can't measure it
        a ~ Uniform(0, 100)        # Semi-major axis [AU]
        e ~ Uniform(0.0, 0.5)      # Eccentricity
        i ~ Sine()                 # Inclination [rad]
        ω ~ UniformCircular()      # Argument of periastron [rad]
        Ω ~ UniformCircular()      # Longitude of ascending node [rad]
        θ ~ UniformCircular()      # Position angle at the reference epoch [rad]
        epoch = 50000.0            # The reference epoch for θ [MJD]
    end
)
nothing # hide
```

!!! note
    Make sure to adjust the epoch `50000.0` above to match your most constraining data
    epoch. `θ` and `epoch` together fix where the planet is on its orbit; you could
    equally supply `tp` (epoch of periastron passage) or `M0` and `epoch`.

!!! note "Where did `M` go?"
    There is no `M` variable any more. The gravitating mass of an orbit is the total mass
    of the bodies it binds, computed from the model rather than declared. Here that is
    `A.mass + b.mass`, so fixing `b.mass = 0.0` makes `A.mass` play exactly the role
    v1's `M` did.

Now build the observation, saying what it observes (`target`) and what it is measured
against (`ref`), and assemble the system:
```@example 1
astrom = RelAstromObs(astrom_dat; target=b, ref=A, name="GPI astrom")

sys = System(
    name="HD1234",
    bodies=[A, b],
    observations=[astrom],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)  # Parallax [mas]
    end
)
```

Which frame variables the system block defines chooses the frame: `plx` alone gives
angular observables in milliarcseconds, which is what relative astrometry needs. See
[System Construction](@ref derived) for more options.

Compile the model into efficient sampling code:
```@example 1
model = Octofitter.LogDensityModel(sys)
```

Initialize the starting points for the chains. You can optionally provide starting values
for particular variables (`UniformCircular` priors are a special case — see
[Priors](@ref priors)). Body variables are nested under `bodies`:
```@example 1
init_chain = initialize!(model, (;
    plx = 50.001,
    bodies = (;
        A = (; mass = 1.18),
        b = (;
            a = 10.0,
            e = 0.01,
        )
    )
))
nothing # hide
```

Visualize the starting point. You can use this plot to make absolutely sure your data was
entered correctly:
```@example 1
octoplot(model, init_chain)
```

Sample from the posterior using Hamiltonian Monte Carlo (see [Samplers](@ref samplers) for
other options):
```@example 1
octofit(model, verbosity=0, iterations=2, adaptation=2); # hide
chain = octofit(model, iterations=1000)
```

Visualize the results with orbit plots and a corner plot:
```@example 1
octoplot(model, chain)     # Plot orbits and data
```

```@example 1
octocorner(model, chain, small=true)   # Corner plot of posterior
```

Save the results to a FITS file (see [Loading and Saving Data](@ref loading-saving) for
other formats):
```julia
Octofitter.savechain("output.fits", chain)
chain = Octofitter.loadchain("output.fits"; model)
```

Passing `model` to `loadchain` is recommended: it checks that the chain's columns match
the model you are about to use it with, instead of silently returning `missing` for any
that don't.

## Working with Dates

These helper functions convert dates to and from Modified Julian Days:
```@example 1
mjd("2020-01-01")     # Date string to MJD
years2mjd(2020.0)     # Decimal year to MJD
mjd2date(50000)       # MJD to date
```

## Next Steps
See the Tutorials section for complete examples, starting with the
[Basic Astrometry Fit](@ref fit-astrometry).
