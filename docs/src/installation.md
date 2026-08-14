# [Installation](@id install)

The first step to using Octofitter.jl is to install Julia. If you're used to Python, don't worry --- Julia is easy to install, and you won't need to code anything other than changing your input data.


## Installing Julia
Visit the [julialang.org](https://julialang.org/downloads/) Downloads page, and select the latest stable version for your operating system. Click the `[help]` links next to your operating system if you require more detailed instructions.

Octofitter requires Julia **1.11 or newer**, and **1.12 is recommended**.

## Installing Octofitter

1. Start julia in a terminal by running `julia`
2. Type `]` to enter package-mode (see Julia documentation for more details)
3. Type `add Octofitter Distributions CairoMakie PairPlots`

You will need the Distributions.jl package so that you can specify priors for different parameters in your models.

## Extension Packages
Some Octofitter functionality exists in extension packages.
If you need one of these packages you can install them like so:
```
pkg> add OctofitterRadialVelocity
pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterImages
pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterInterferometry
```

These aren't included by default since they may include a number of heavier dependencies that aren't needed by all users.
They are described further in relevant sections of the documentation.


## Companion Packages by Use Case

Depending on your analysis, you may need additional packages. Here's a guide organized by task:

**For plotting and corner plots:**
```julia
pkg> add CairoMakie PairPlots
```

**For radial velocities:**
```julia
pkg> add OctofitterRadialVelocity
```

**For direct imaging:**
```julia
pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterImages
```

**For interferometry:**
```julia
pkg> add http://github.com/sefffal/Octofitter.jl:OctofitterInterferometry
```


!!! note "Parallel tempering (Pigeons)"
    [`octofit_pigeons`](@ref) — for multimodal posteriors, discrete variables, and
    Bayesian model comparison via `stepping_stone` Add
    Pigeons and run `using Pigeons`. See [Samplers](@ref samplers).

## Fitting your first model
Start with the [Quick Start](@ref quick-start) guide, then the [Basic Astrometry Fit](@ref fit-astrometry) tutorial. It shows how one can model the orbit of one planet based on relative astrometry points.
