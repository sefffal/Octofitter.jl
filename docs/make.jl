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
            "Fit Thiele-Innes Paramerers" => "thiele-innes.md",
            "Fit Relative RV Data" => "fit-rv-rel.md",
            "Fit GRAVITY Wide Data" => "fit-grav-wide.md",
            # TODO:
            "Detection Limits" => "limits.md",
        ],
        "Statistical Methods" => [
            "Prior Predictive Checks" => "prior-pred.md",
            "Posterior Predictive Checks" => "post-pred.md",
            # TODO:
            "Cross Validataion" => "cross-validation.md",
            "Simulation Based Calibration" => "sbc.md",
        ],
        "Compatibility" => [
            # "RadVel" => "compat-radvel.md",
            "Orbitize!" => "compat-orbitize.md",
        ],
        "Documentation" => [
            "Using Python" => "python.md",
            "Chains" => "chains.md",
            "Loading and Saving Data" => "loading-saving.md",
            "Sampler" => "samplers.md",
            "Parallel Sampling" => "parallel-sampling.md",
            "Priors" => "priors.md",
            "Derived Variables" => "derived.md",
            "Connecting Mass with Photometry" => "mass-photometry.md",
            "Custom Likelihoods" => "custom-likelihood.md",
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
