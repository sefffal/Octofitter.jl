# Quick Start (@id quick-start)

This guide introduces the key concepts in Octofitter: 
* Likelihood objects to hold your data
* `Planet` and `System` models to specify variables, priors, and system architecture
* Sampling from the posterior using MCMC
* Plotting the results
* Saving the chain

For installation instructions, see [Installation](@ref install).


## Example: Fit a Single Planet Orbit to Relative Astrometry

Load the required packages:
```julia
using Octofitter, Distributions, CairoMakie, PairPlots
```

Create a [`PlanetRelAstromLikelihood`](@ref) object containing your observational data. In this case its the position of the planet relative to the star, but many other kinds of data are supported:
```julia
astrom_dat = Table(
    epoch = [50000, 50120, 50240],      # Dates in MJD
    ra = [-505.7, -502.5, -498.2],      # [mas] East positive
    dec = [-66.9, -37.4, -7.9],         # [mas] North positive
    σ_ra = [10.0, 10.0, 10.0],          # [mas] Uncertainties
    σ_dec = [10.0, 10.0, 10.0],         # [mas] Uncertainties
    cor = [0.0, 0.0, 0.0]               # RA/Dec correlations
)
astrom = PlanetRelAstromLikelihood(astrom_dat)
```

Define a planet model with orbital elements and their [prior distributions](@ref priors):
```julia
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[astrom],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)  # Total mass (solar masses) for this orbit
        a ~ Uniform(0, 100)        # Semi-major axis [AU]
        e ~ Uniform(0.0, 0.5)      # Eccentricity  
        i ~ Sine()                 # Inclination [rad]
        ω ~ UniformCircular()      # Argument of periastron [rad]
        Ω ~ UniformCircular()      # Longitude of ascending node [rad]
        θ ~ UniformCircular()      # Position angle at reference epoch [rad]
        # Epoch of periastron passage
        # We calculate it from the position angle above
        tp = θ_at_epoch_to_tperi(θ, 50000; M, e, a, i, ω, Ω)  
    end
)
```

!!! note
    Make sure to adjust the epoch `50000` above to match your most constraining data epoch.

Define the system with its mass and distance - see [System Construction](@ref derived) for more options:
```julia
sys = System(
    name="HD1234",
    companions=[planet_b],
    likelihoods=[],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)  # Parallax (mas)
    end
)
```

That there are many different orbit parameterizations, each requiring different of parameters names. The `KepOrbit` is a full 3D keplerian orbit with Campbell parameters. `Visual` means that we have a defined parallax distance `plx` that can map separations in AU to arcseconds.

Compile the model into efficient sampling code:
```julia
model = Octofitter.LogDensityModel(sys)
```

Initialize the starting points for the chains. You can optionally provide starting values directly for certain variables.
```julia
init_chain = initialize!(model, (;
    plx = 50.001,
    planets = (;
        b=(;
            M = 1.18,
            a = 10.0,
            e = 0.01,
        )
    )
)) 
```

Visualize the starting point. You can use this plot to make absolutely sure your data is entered in correctly:
```julia
octoplot(model, init_chain)
```

Sample from the posterior using Hamiltonian Monte Carlo (see [Samplers](@ref samplers) for other options):
```julia
chain = octofit(model, iterations=1000)
```

Visualize the results with orbit plots and a corner plot:
```julia
octoplot(model, chain)     # Plot orbits and data
```

```julia
octocorner(model, chain)   # Corner plot of posterior
```

Save the results to a FITS file (see [Loading and Saving Data](@ref loading-saving) for other formats):
```julia
Octofitter.savechain("output.fits", chain)
```

## Working with Dates

Convert dates to and from Modified Julian Days using these helper functions:
```julia
mjd("2020-01-01")     # Date string to MJD
years2mjd(2020.0)     # Decimal year to MJD
mjd2date(50000)       # MJD to date
```

## Next Steps
See the Tutorials section for complete examples.