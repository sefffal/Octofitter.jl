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
        "Getting Started" => [
            "Installation" => "installation.md",
            "Quick Start" => "quick-start.md",
            "FAQ" => "faq.md",
            "Migration Guide" => "migration.md",
        ],
        "Tutorials" => [
            "Relative Astrometry" => [
                "Basic Astrom Fit" => "rel-astrom.md",
                "Observable Priors" => "rel-astrom-obs.md",
                "Resonant Co-Planar Model" => "fit-coplanar.md",
                "Thiele-Innes Parameters" => "thiele-innes.md",    
            ],
            "Radial Velocity" => [
                "Basic RV Fit" => "rv-1.md",
                "Gaussian Process" => "rv-gp.md",
                "Multiple Planets" => "rv-multi-planet.md",
                "Relative RV Data" => "fit-rv-rel.md",
            ],
            "Absolute Astrometry" => [
                "Proper Motion Anomaly" => "pma.md",
                "Hipparcos IAD" => "hipparcos.md",
                "Gaia DR4 IAD" => "gaia-iad.md"
            ],
            "Images and More" => [
                "Image Data (de-orbiting)" => "images.md",
                "Extract Astrom. and Photometry" => "extract-phot-astrom.md",
                "Connect Mass and Photometry" => "mass-photometry.md",
                "Interferometer Data" => "fit-interfere.md",
                "Likelihood Map" => "fit-likemap.md",
                "GRAVITY Wide Data" => "fit-grav-wide.md",
            ],
            "Joint Models" => [
                "Astrometry, PMA, and RV" => "astrom-pma-rv.md",
                "RV and Relative Astrometry" => "fit-rv-astrom.md",
                "RV and Proper Motion Anomaly" => "rv.md",
                "Calculate Detection Limits" => "limits.md",
            ],
            "Bayesian Workflows" => [
                "Circular or Eccentric? Model Comparison" => "eccentric-or-circular.md",
                "Generating and Fitting Simulated Data" => "data-simulation.md",
                "Prior Predictive Checks" => "prior-pred.md",
                "Posterior Predictive Checks" => "post-pred.md",
                "Cross Validation" => "cross-validation.md",
                "Simulation Based Calibration" => "sbc.md",
            ],
        ],
        "Documentation" => [
            "Using Python" => "python.md",
            "Chains" => "chains.md",
            "Orbit plots with `octoplot`"=>"octoplot.md",
            "RV plots with `rvpostplot`"=>"rvpostplot.md",
            "Loading and Saving Data" => "loading-saving.md",
            "Sampler" => "samplers.md",
            "Distributed Sampling" => "parallel-sampling.md",
            "Priors" => "priors.md",
            "Derived Variables" => "derived.md",
            "Custom Likelihoods" => "custom-likelihood.md",
            "Kepler Solver" => "kepler.md",
            "Orbitize! Compatibility" => "compat-orbitize.md",
            "Full API Documentation" => "api.md"
        ],
        "Developer Documentation" => [
            "Architecture Overview" => "dev/architecture.md",
            "Epoch Tables and Kepler Pre-Solving" => "dev/epoch-tables-kepler.md"
        ]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold=nothing
    ),
    pagesonly=true,
    warnonly=:example_block
)


deploydocs(
    repo = "github.com/sefffal/Octofitter.jl.git",
    devbranch = "main",
    push_preview = true
)
