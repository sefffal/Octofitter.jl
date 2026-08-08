# [Installation](@id install)

The first step to using Octofitter.jl is to install Julia. If you're used to Python, don't worry --- Julia is easy to install, and you won't need to code anything other than changing your input data.


## Installing Julia
Visit the [julialang.org](https://julialang.org/downloads/) Downloads page, and select the latest stable version for your operating system. Click the `[help]` links next to your operating system if you require more detailed instructions.

Octofitter requires Julia **1.11 or newer**, and **1.12 is recommended** — it is what the test suite and the documentation build track, and several of the allocation-free code paths depend on escape analysis that older versions do not perform.

## Installing Octofitter

1. Start julia in a terminal by running `julia`
2. Type `]` to enter package-mode (see Julia documentation for more details)
3. Type `add Octofitter Distributions CairoMakie PairPlots`

You will need the Distributions.jl package so that you can specify priors for different parameters in your models.
[CairoMakie.jl](http://makie.juliaplots.org/) is used for generating plots and isn't needed if you only want text-based summary outputs. [PairPlots.jl](https://sefffal.github.io/PairPlots.jl/dev/) (in combination with CairoMakie) is used for generating corner plots and can also be skipped if these aren't of interest.

Octofitter re-exports [PlanetOrbits.jl](https://sefffal.github.io/PlanetOrbits.jl/dev/), so
`orbitsolve`, `raoff`, `radvel`, `Barycentre`, `Photocentre`, the `mjup`/`mearth`/`msun`
constants and the rest of the orbit machinery are available as soon as you have loaded
Octofitter. You only need to `add PlanetOrbits` separately if you want to use it on its own.

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

!!! note "Basic radial velocity does not need a separate package"
    `RadialVelocityObs` — which covers both stellar reflex RV (`target=A, ref=Barycentre`)
    and relative RV of a companion (`target=b, ref=A`) — lives in core Octofitter.
    You only need `OctofitterRadialVelocity` for the Celerite Gaussian-process backend,
    [`MarginalizedRVObs`](@ref), and the archive loaders (HARPS, HIRES, Lick, CES).

## Companion Packages by Use Case

Depending on your analysis, you may need additional packages. Here's a guide organized by task:

**For plotting and corner plots:**
```julia
pkg> add CairoMakie PairPlots
```

**For radial velocity archive loaders, Celerite GPs, or the marginalized RV likelihood:**
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

!!! tip "Quick install for common workflows"
    For a typical multi-planet RV analysis:
    ```julia
    pkg> add Octofitter OctofitterRadialVelocity Distributions CairoMakie PairPlots
    ```

!!! note "Parallel tempering (Pigeons)"
    [`octofit_pigeons`](@ref) — for multimodal posteriors, discrete variables, and
    Bayesian model comparison via `stepping_stone` — lives in a package extension. Add
    Pigeons and `using Pigeons` before calling it, or the function will exist with no
    methods. See [Samplers](@ref samplers).

## Fitting your first model
Start with the [Quick Start](@ref quick-start) guide, then the [Basic Astrometry Fit](@ref fit-astrometry) tutorial. It shows how one can model the orbit of one planet based on relative astrometry points.
