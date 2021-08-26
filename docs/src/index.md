
# DirectDetections.jl

Welcome to the documentation page for DirectDetections.jl. 
This page includes tutorial and an API reference for using this package.

DirectDetections is a Julia package for performing Bayesian inference
against direct images of exoplanets, exoplanet astrometry, astrometric acceleration
of the host star, and radial velocity (future).

You build a model of the system using the functions described below, list any
data you might have, and start the sampler. The package also contains analysis
and visualization tools for understanding your results.

**Supported data:**
* sample directly from images
* exoplanet astrometry 
* stellar astrometric acceleration
Any and all combinations also work together.

**Modelling features:**
* multiple planets (one or more)
* co-planar, and non-coplanar systems
* hierarchical models
* link mass to photometry via an atmosphere model

## Table of Contents
```@contents
```