# [Installation](@id install)

The first step to using Octofitter.jl is to install Julia. If you're used to Python, don't worry --- Julia is easy to install, and you won't need to code anything other than changing your input data.


## Installing Julia
Visit the [julialang.org](https://julialang.org/downloads/) Downloads page, and select the latest stable version for your operating system. This is 1.10.0 at the time of writing. Click the `[help]` links next to your operating system if you require more detailed instructions.

## Installing Octofitter

1. Start julia in a terminal by running `julia`
2. Type `]` to enter package-mode (see Julia documentation for more details)
3. Type `add Octofitter Distributions CairoMakie PairPlots`

You will need the Distributions,jl package so that you can specify priors for different parameters in your models.
[CairoMakie.jl](http://makie.juliaplots.org/) is used for generating plots and isn't needed if you only want text-based summary outputs. [PairPlots.jl](https://sefffal.github.io/PairPlots.jl/dev/) (in combination with CairoMakie) is used for generating corner plots and can also be skipped if these aren't of interest.

## Extension Packages
Some Octofitter functionality exists in extension packages, including radial velocity fitting.
If you need one of these packages you can install them like so:
```
pkg> add OctofitterRadialVelocity
pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterImages
pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterInterferometry
```

These aren't included by default since they may include a number of heavier dependencies that aren't needed by all users.
They are descibed further in relevant sections of the documentation.

## Fitting your first model
Start with the [Quick Start](@ref fit-astrometry) tutorial. It shows how one can model the orbit of one planet based on relative astrometry points.
