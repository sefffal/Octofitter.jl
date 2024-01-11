# *Octofitter*

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/sefffal/Octofitter.jl)

Welcome to the documentation page for Octofitter.jl. 
This page includes tutorial and an API reference for using this package.

Octofitter is a Julia package for performing Bayesian inference 
against a wide variety of exoplanet / binary star data.
You can also use Octofitter from Python using the [Python guide](@ref python).

The package provides a simple but powerful modelling language which is used to generate
efficient, differentiable code. You can then plug it into a variety of samplers.
The package also contains analysis and visualization tools for understanding your results.


!!! note
    Octofitter is under active development and is only tested against Julia 1.9+

**Supported data:**
* Fit exoplanet orbits to relative astrometry
* Fit radial velocity data
* Model stellar activity with Gaussian processes
* Model stellar astrometric accerlation (Gaia-Hipparcos proper motion anomaly)
* "De-orbiting": combine a sequence of images with orbital motion to detect planets
* Sample directly from images and interferometric visibilities
* experimental support for transit data based on Transits.jl

You can freely combine any of the above data types. 
Any and all combinations work together.

**Modelling features:**
* multiple planets (one or more)
* hyperbolic orbits
* co-planar, and non-coplanar systems
* arbitrary priors and parameterizations
* link mass to photometry via atmosphere models
* hierarchical models (with a bit of work from the user)

**Speed:**

Fit astrometry on your laptop in seconds!

* Highly optimized code and derivatives are generated from your model
* Higher order sampler (No U-Turn sampler) which explores the parameter space very efficiently 
* The sampler is automatically warmed up using a variational approximation from the Pathfinder algorithm (Pathfinder.jl) 

Multi-body physics is not currently supported. A Pull-request to PlanetOrbits.jl implementing this functionality would be welcome.

See also: the python libraries [Orbitize!](https://orbitize.readthedocs.io/en/latest/), [orvara](https://github.com/t-brandt/orvara), and [exoplanet](https://docs.exoplanet.codes/en/latest/).

### Read the paper
In addition to these documentation and tutorial pages, you can read the paper published in the [Astronomical Journal](https://dx.doi.org/10.3847/1538-3881/acf5cc) (open-access).

Please cite this paper if you use Octofitter in your work.


## Rady?
Ready to get started? Follow our [installation guide](@ref getting-started) and then follow our [first tutorial](@fit-astrometry).