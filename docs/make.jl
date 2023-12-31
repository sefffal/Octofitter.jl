using Pkg
cd(@__DIR__)
Pkg.activate(".")
pkg"dev .. ../OctofitterRadialVelocity ../OctofitterImages ../OctofitterWhereistheplanet ../OctofitterInterferometry"
Pkg.instantiate()
Pkg.precompile()

using Documenter, Octofitter, OctofitterRadialVelocity

# Increase resolution of figures
using CairoMakie
CairoMakie.activate!(px_per_unit=4)


makedocs(
    sitename="Octofitter.jl",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Tutorials" => [
            "Fit Relative Astrometry" => "rel-astrom.md",
            "Fit Proper Motion Anomaly" => "pma.md",
            "Fit Images" => "images.md",
            "Fit Radial Velocity" => "rv-1.md",
            "Fit RV and Astrometry" => "rv.md",
            "Connecting Mass with Photometry" => "mass-photometry.md",
            "Loading and Saving Data" => "loading-saving.md",
            "Custom Likelihoods" => "custom-likelihood.md",
            "Calling from Python" => "python.md",
            # "Multiple Planets" => "multi-planets.md",
            # "" => "multi-planets.md",
        ],
        "Documentation" => [
            "Priors" => "priors.md",
            "Derived Variables" => "derived.md",
            "Sampler" => "samplers.md",
            "Parallel Sampling" => "parallel-sampling.md",
            "Chains" => "chains.md",
            "Kepler Solver" => "kepler.md",
            "API" => "api.md"
        ]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)


deploydocs(
    repo = "github.com/sefffal/Octofitter.jl.git",
    devbranch = "main"
)
