using Pkg
cd(@__DIR__)
Pkg.activate(".")


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
            "Fit with Observable Priors" => "rel-astrom-obs.md",
            "Fit Radial Velocity" => "rv-1.md",
            "Fit Proper Motion Anomaly" => "pma.md",
            "Fit Images" => "images.md",
            "Fit RV and Proper Motion Anomaly" => "rv.md",
            "Fit RV and Relative Astrometry" => "fit-rv-astrom.md",
            "Fit Resonant Co-Planar Model" => "fit-coplanar.md",
            "Loading and Saving Data" => "loading-saving.md",
            "Connecting Mass with Photometry" => "mass-photometry.md",
            "Custom Likelihoods" => "custom-likelihood.md",
            "Python Guide" => "python.md",
            "Simulation Based Calibration" => "sbc.md",
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
    ),
    pagesonly=true
)


deploydocs(
    repo = "github.com/sefffal/Octofitter.jl.git",
    devbranch = "main"
)
