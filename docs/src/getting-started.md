# Getting Started

The first step to using Octofitter.jl is to install Julia. If you're used to Python, don't worry --- Julia is easy to install, and you won't need to code anything other than changing your input data.


## Installing Julia
Visit the [julialang.org](https://julialang.org/downloads/) Downloads page, and select the latest stable version for your operating system. Currently, this is 1.8.3. Click the `[help]` links next to your operating system if you require more detailed instructions.

## Installing Octofitter
Normally, Julia packages are installed from the General registry. Since Octofitter isn't quite ready for prime time, you can instead install it directly from GitHub.

1. Start julia in a terminal by running `julia`
2. Type `]` to enter package-mode (see Julia documentation for more details)
3. Type `add https://github.com/sefffal/Octofitter.jl Distributions`

You will need the Distributions package so that you can specify priors for different parameters in your models.

If you would like to visualize your results, you can also install the Plots package:
4. Type `add Plots`

This will take a little while to download all the required packages and precompile for your system.

## Extension Packages
Some Octofitter functionality exists in extension packages, including radial velocity fitting.
If you need one of these packages you can install them like so:
```
pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterImages
pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterRadialVelocity
```
That is, specify the extension package you want to install after the colon (`:`).

These aren't included by default since they may include a number of heavier dependencies that aren't needed by all users.

## Fitting your first model
Start with the [Fit AstrometryLikelihood](@ref fit-astrometry) tutorial which shows how to model of one planet with some astrometry points.
