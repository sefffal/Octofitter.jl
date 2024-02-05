using Pkg
cd(@__DIR__)
Pkg.activate(".")


using Documenter, Octofitter, OctofitterRadialVelocity

# Increase resolution of figures
using CairoMakie
CairoMakie.activate!(px_per_unit=2)


makedocs(
    sitename="Octofitter.jl",
    pages = [
        "Home" => "index.md",
        "Installation" => "getting-started.md",
        "Tutorials" => [
            "Fit Relative Astrometry" => "rel-astrom.md",
            "Fit with Observable Priors" => "rel-astrom-obs.md",
            "Fit Radial Velocity" => "rv-1.md",
            "Fit Proper Motion Anomaly" => "pma.md",
            "Fit Orbits to Images" => "images.md",
            "Fit RV and Proper Motion Anomaly" => "rv.md",
            "Fit RV and Relative Astrometry" => "fit-rv-astrom.md",
            "Fit Resonant Co-Planar Model" => "fit-coplanar.md",
            "Fit Interferometer Data" => "fit-interfere.md",
            "Fit Likelihood Map" => "fit-likemap.md",
            "Fit Thiele-Innes Paramerers" => "thiele-innes.md"
        ],
        "Statistical Methods" => [
            "Prior Predictive Checks" => "prior-pred.md",
            "Posterior Predictive Checks" => "post-pred.md",
            "Simulation Based Calibration" => "sbc.md",
            # "Leave-One-Out Cross Validation" => "psis-loo.md",
        ],
        "Documentation" => [
            "Using Python" => "python.md",
            "Priors" => "priors.md",
            "Derived Variables" => "derived.md",
            "Loading and Saving Data" => "loading-saving.md",
            "Connecting Mass with Photometry" => "mass-photometry.md",
            "Custom Likelihoods" => "custom-likelihood.md",
            "Sampler" => "samplers.md",
            "Parallel Sampling" => "parallel-sampling.md",
            "Chains" => "chains.md",
            "Kepler Solver" => "kepler.md",
            "Full API Documentation" => "api.md"
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
