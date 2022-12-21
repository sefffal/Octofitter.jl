
# DirectDetections.jl
[GitHub](https://github.com/sefffal/DirectDetections.jl)

Welcome to the documentation page for DirectDetections.jl. 
This page includes tutorial and an API reference for using this package.

DirectDetections is a Julia package for performing Bayesian inference
against direct images of exoplanets, relative astrometry, astrometric acceleration
of the host star, and radial velocity (future).

The package provides a simple but powerful modelling language which is used to generate
efficient, differentiable code for your system.
The package also contains analysis and visualization tools for understanding your results.

**Supported data:**
* sample directly from images
* exoplanet astrometry 
* stellar astrometric acceleration
* radial velocity
* experimental support for transit data based on Transits.jl

Any and all combinations also work together.

**Modelling features:**
* multiple planets (one or more)
* co-planar, and non-coplanar systems
* arbitrary priors and parameterizations
* link mass to photometry via atmosphere models

**Speed:**

Fit astrometry on your laptop in minutes!

* Highly optimized code and derivatives are generated from your model
* Higher order sampler (No U-Turn sampler) which explores the parameter space very efficiently 
* Run on a single core, multiple threads, or hundreds of nodes by changing just a single line of code

The package supports only bound, 2-body Keplerian orbits. Support for hyperbolic orbits and multi-body physics are not currently planned. Pull-requests to PlanetOrbits implementing this functionality would be welcome.

See also: the python libraries [Orbitize!](https://orbitize.readthedocs.io/en/latest/), [orvara](https://github.com/t-brandt/orvara), and [exoplanet](https://docs.exoplanet.codes/en/latest/).


### Getting Started
```@contents
Pages = ["getting-started.md"]
Depth = 5
```

### Tutorials
```@contents
Pages = ["modelling.md", "pma.md",  "images.md", "derived.md", "mass-photometry.md"]
Depth = 5
```

### Documentation
```@contents
Pages = ["samplers.md", "chains.md", "kepler.md", "api.md"]
Depth = 5
```

