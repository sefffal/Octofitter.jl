# Octofitter.jl

Octofitter is a Julia package for performing Bayesian inference
against direct images of exoplanets, relative astrometry, astrometric acceleration
of the host star, and radial velocity (future).

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sefffal.github.io/Octofitter.jl/dev)

The package provides a simple but powerful modelling language which is used to generate
efficient, differentiable code for your system.
The package also contains analysis and visualization tools for understanding your results.

**Supported data:**
* sample directly from images
* exoplanet astrometry 
* stellar astrometric acceleration

Any combination of the above.

**Modelling features:**
* multiple planets (one or more)
* co-planar, and non-coplanar systems
* arbitrary priors and parameterizations
* link mass to photometry via atmosphere models

**Speed:**
<p>Fit astrometry on your laptop in seconds!</p>

* Highly optimized code and derivatives are generated from your model
* Higher order sampler (No U-Turn sampler) which explores the parameter space very efficiently 

The package supports only bound, 2-body Keplerian orbits. Support for hyperbolic orbits and multi-body physics are not currently planned.

See also: [Orbitize!](https://orbitize.readthedocs.io/en/latest/), [orvara](https://github.com/t-brandt/orvara), and [exoplanet](https://docs.exoplanet.codes/en/latest/).

For instructions, see the documentation page:

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sefffal.github.io/Octofitter.jl/dev)

A particularily unique feature of this package is that you can model the flux and orbital motion of planets using a sequence of image(s) without any obvious detections. Using this tool, the SNR of a planet can grow with roughly the square root of the number of images. You can spot planets even if there is orbtial motion between the images, or constrain orbits using images with no detection.

![](images/readme-example.png)
