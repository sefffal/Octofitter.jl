# Fit Relative RV Data

Octofitter includes support for fitting relative radial velocity data. Currently this is only tested with a single companion. Please open an issue if you would like to fit multiple companions simultaneously.

The convention we adopt is that positive relative radial velocity is the velocity of the companion (exoplanets) minus the velocity of the host (star).

To fit relative RV data, start by creating a likelihood object:
```@example 1
using Octofitter
using OctofitterRadialVelocity
using CairoMakie

rel_rv_like = PlanetRelativeRVLikelihood(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11),
    )
```
The columns in the table should be . See the standard radial velocity tutorial for examples on how this data can be loaded from a CSV file.

The relative RV likelihood does not incorporate an instrument-specific RV offset. A jitter parameter called `jitter` can still be specified in the `@planet` block, as can parameters for a gaussian process model of stellar noise. Unlike the `StarAbsoluteRVLikelihood`, only a single instrument jitter parameter is supported. If you need to model relative radial velocities from multiple instruments with different jitters, please open an issue on GitHub.

Next, create a planet model. If you are modelling **only** relative RVs, you can use a `RadialVelocityOrbit` parameterization. If, however, you also are modelling absolute radial velocities, then it is possible for the model to constrain inclination and you should use a full `Visual{KepOrbit}` parameterization.
