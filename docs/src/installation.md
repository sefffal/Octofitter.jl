# [Installation](@id install)

The first step to using Octofitter.jl is to install Julia. If you're used to Python, don't worry --- Julia is easy to install, and you won't need to code anything other than changing your input data.


## Installing Julia
Visit the [julialang.org](https://julialang.org/downloads/) Downloads page, and select the latest stable version for your operating system. Click the `[help]` links next to your operating system if you require more detailed instructions.

Octofitter requires Julia **1.11 or newer**, and **1.12 is recommended**.

## Installing Octofitter

!!! warning "Octofitter v9 is an unreleased prerelease"
    These docs describe **Octofitter v9**, which is built on **PlanetOrbits v2**.
    Neither is registered yet, so `pkg> add Octofitter` gives you the *released*
    v8 line and the code on these pages will not run. Until v9 is registered,
    install from the `v2` branches with the single command below --- and read
    [Installing the v9 prerelease](@ref install-prerelease) if anything goes wrong.

Start Julia in a terminal by running `julia`, then copy this in. It creates a
project environment dedicated to your analysis and installs everything into it
in **one** resolve:

```julia
using Pkg
Pkg.activate("my-orbit-fit")     # a fresh, dedicated environment
Pkg.add([
    PackageSpec(url="https://github.com/sefffal/PlanetOrbits.jl", rev="v2"),
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl",   rev="v2"),
    PackageSpec(name="Distributions"),
    PackageSpec(name="CairoMakie"),
    PackageSpec(name="PairPlots"),
])
```

You will need the Distributions.jl package so that you can specify priors for different parameters in your models.
CairoMakie and PairPlots are used for all the plots in these tutorials.

Every following session, start Julia and run `using Pkg; Pkg.activate("my-orbit-fit")`
before `using Octofitter` --- otherwise you are back in your default environment,
which does not have the prerelease installed.

## [Installing the v9 prerelease](@id install-prerelease)

Two rules keep this painless. They matter because PlanetOrbits `v2` carries the
version number `2.0.0-DEV` and Octofitter `v2` carries `9.0.0-DEV`; neither is in
the General registry, so Pkg can only find them if you name their branches.

1. **Install into a fresh environment.** Adding a `v2` branch on top of an
   environment that already has a released Octofitter is the main way people get
   stuck (see below). `Pkg.activate("some-new-directory")` first.
2. **Add every Octofitter package in the same `Pkg.add` call**, together with
   PlanetOrbits. One `Pkg.add` with several `PackageSpec`s is a single resolve,
   so Pkg sees the whole prerelease stack at once. Adding them one at a time also
   works *provided* PlanetOrbits goes first, but there is no reason to risk it.

### Extension packages

Some Octofitter functionality lives in extension packages, which are not included
by default since they pull in heavier dependencies. They live in subdirectories of
the Octofitter.jl repository, so they need `rev="v2"` **and** `subdir=` --- add
them in the same call as everything else:

```julia
using Pkg
Pkg.activate("my-orbit-fit")
Pkg.add([
    PackageSpec(url="https://github.com/sefffal/PlanetOrbits.jl", rev="v2"),
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl",   rev="v2"),
    # radial velocity:
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl", rev="v2", subdir="OctofitterRadialVelocity"),
    # direct imaging:
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl", rev="v2", subdir="OctofitterImages"),
    # interferometry:
    PackageSpec(url="https://github.com/sefffal/Octofitter.jl", rev="v2", subdir="OctofitterInterferometry"),
])
```

Pick only the ones you need; they are described further in the relevant sections
of the documentation.

!!! warning "Don't `add OctofitterRadialVelocity` by name"
    The registered `OctofitterRadialVelocity` belongs to the v8 line. Adding it by
    name into a v9 environment either fails to resolve or drags the released
    Octofitter back in.

### Updating

`Pkg.update()` moves branch checkouts to the current tip of `v2`. Update the whole
stack at once so PlanetOrbits and Octofitter stay in step:

```julia
using Pkg
Pkg.activate("my-orbit-fit")
Pkg.update()
```

### If something went wrong

Two failures account for nearly all reports.

**`Unsatisfiable requirements detected for package Octofitter`**, mentioning
`PlanetOrbits [fd6f9641] is fixed to version 2.0.0-DEV`:

```
Unsatisfiable requirements detected for package Octofitter [daf3887e]:
 ├─possible versions are: 1.0.0 - 8.3.0 or uninstalled
 └─restricted by compatibility requirements with PlanetOrbits [fd6f9641] to versions: uninstalled — no versions left
```

You have a **released** Octofitter in the environment and are adding PlanetOrbits
`v2` on top of it. Released Octofitter cannot use PlanetOrbits 2. Start a fresh
environment and use the single `Pkg.add` above.

**`UndefVarError: 'mjup' not defined in 'PlanetOrbits'`** (or
`UndefVarError: 'gp_condition' not defined in 'Octofitter'`) while precompiling:

You have a v9 package resolved against a **v1/v8** dependency --- the symptom of
adding one prerelease package without the others. Again: fresh environment, one
`Pkg.add`.

In both cases `Pkg.status()` is the quickest diagnosis. A correct v9 environment
shows a git revision after every Octofitter-family entry:

```
  [daf3887e] Octofitter v9.0.0-DEV `https://github.com/sefffal/Octofitter.jl#v2`
  [fd6f9641] PlanetOrbits v2.0.0-DEV `https://github.com/sefffal/PlanetOrbits.jl#v2`
```

If you see `Octofitter v8.x` or `PlanetOrbits v0.11.x` with no `#v2` after it,
that environment is on the released line.

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
