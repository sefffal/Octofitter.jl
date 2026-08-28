# [Installation](@id install)

The first step to using Octofitter.jl is to install Julia. If you're used to Python, don't worry --- Julia is easy to install, and you won't need to code anything other than changing your input data.


## Installing Julia
Visit the [julialang.org](https://julialang.org/downloads/) Downloads page, and select the latest stable version for your operating system. Click the `[help]` links next to your operating system if you require more detailed instructions.

Octofitter requires Julia **1.11 or newer**, and **1.12 is recommended**.

## Installing Octofitter

These docs describe **Octofitter v9**. If you are coming from v8, your scripts
need updating --- see the [migration guide](@ref v9-migration).

Start Julia in a terminal by running `julia`, then copy this in. It creates a
project environment dedicated to your analysis and installs everything into it:

```julia
using Pkg
Pkg.activate("my-orbit-fit")     # a fresh, dedicated environment
Pkg.add(["Octofitter", "Distributions", "CairoMakie", "PairPlots"])
```

You will need the Distributions.jl package so that you can specify priors for different parameters in your models.
CairoMakie and PairPlots are used for all the plots in these tutorials.

Every following session, start Julia and run `using Pkg; Pkg.activate("my-orbit-fit")`
before `using Octofitter` --- otherwise you are back in your default environment.

If you have an **existing** environment from the v8 era (or from the v9
prerelease branches), the cleanest upgrade is a fresh environment as above.
`Pkg.status()` should then show `Octofitter v9.x` and `PlanetOrbits v1.x`; if
it shows `Octofitter v8.x` or `PlanetOrbits v0.11.x`, something in that
environment is holding the old line back --- start fresh.

### Extension packages

Some Octofitter functionality lives in extension packages, which are not included
by default since they pull in heavier dependencies. `OctofitterRadialVelocity`
is a registered package --- add it by name:

```julia
pkg> add OctofitterRadialVelocity
```

`OctofitterImages` and `OctofitterInterferometry` are unregistered and live in
subdirectories of the Octofitter.jl repository, so they need a `PackageSpec`:

```julia
using Pkg
Pkg.activate("my-orbit-fit")
Pkg.add([
    # direct imaging:
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl", subdir="OctofitterImages"),
    # interferometry:
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl", subdir="OctofitterInterferometry"),
])
```

Pick only the ones you need; they are described further in the relevant sections
of the documentation.

## [Installing a development version](@id install-prerelease)

To try unreleased work from a branch (for example, the tip of `main`), install
the whole stack from branches **in a fresh environment and a single `Pkg.add`
call** --- one resolve sees the whole set at once, so Pkg cannot pair a
development Octofitter with a released PlanetOrbits or vice versa:

```julia
using Pkg
Pkg.activate("my-orbit-fit-dev")     # fresh, dedicated environment
Pkg.add([
    PackageSpec(url="https://github.com/sefffal/PlanetOrbits.jl", rev="master"),
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl",   rev="main"),
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl",   rev="main", subdir="OctofitterRadialVelocity"),
])
```

`Pkg.update()` moves branch checkouts to the current branch tip. If precompilation
fails with an `UndefVarError` mentioning a symbol from another Octofitter-family
package, the environment has mixed a development package with a released
dependency --- start a fresh environment and repeat the single `Pkg.add`.
`Pkg.status()` is the quickest diagnosis: every Octofitter-family entry should
show a git revision (`#main`), or none of them should.

## Companion Packages by Use Case

Depending on your analysis, you may need additional packages. These are ordinary
registered packages, so you can `pkg> add` them by name at any time:

**For plotting and corner plots:**
```julia
pkg> add CairoMakie PairPlots
```

**For parallel tempering, Bayesian evidence, and multimodal posteriors:**
```julia
pkg> add Pigeons
```

!!! note "Parallel tempering (Pigeons)"
    [`octofit_pigeons`](@ref) — for multimodal posteriors, discrete variables, and
    Bayesian model comparison via `stepping_stone` — lives in a package extension.
    Add Pigeons and run `using Pigeons`. See [Samplers](@ref samplers).

## Fitting your first model
Start with the [Quick Start](@ref quick-start) guide, then the [Basic Astrometry Fit](@ref fit-astrometry) tutorial. It shows how one can model the orbit of one planet based on relative astrometry points.
