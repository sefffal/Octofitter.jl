using Pkg
cd(@__DIR__)
Pkg.activate(".")


using Documenter, Octofitter, OctofitterRadialVelocity

# Increase resolution of figures
using CairoMakie
CairoMakie.activate!(px_per_unit=2)

# =============================================================================
# Page Definitions - Split into partitions for parallel CI builds
# =============================================================================

# Partition 1: Getting Started + Documentation (fast, mostly text)
const PAGES_PARTITION_1 = [
    "Home" => "index.md",
    "Getting Started" => [
        "Installation" => "installation.md",
        "Quick Start" => "quick-start.md",
        "FAQ" => "faq.md",
        "Migration Guide" => "migration.md",
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
        "Architecture Overview" => "dev/architecture.md"
    ]
]

# Partition 2: Relative Astrometry + Images (medium complexity)
const PAGES_PARTITION_2 = [
    "Tutorials" => [
        "Relative Astrometry" => [
            "Basic Astrom Fit" => "rel-astrom.md",
            "Observable Priors" => "rel-astrom-obs.md",
            "Resonant Co-Planar Model" => "fit-coplanar.md",
            "Thiele-Innes Parameters" => "thiele-innes.md",
        ],
        "Images and More" => [
            "Image Data (de-orbiting)" => "images.md",
            "Extract Astrom. and Photometry" => "extract-phot-astrom.md",
            "Connect Mass and Photometry" => "mass-photometry.md",
            "Interferometer Data" => "fit-interfere.md",
            "Likelihood Map" => "fit-likemap.md",
            "GRAVITY Wide Data" => "fit-grav-wide.md",
        ],
    ],
]

# Partition 3: Absolute Astrometry + Limits (heavy - pma, g23h, gaia-iad, limits)
const PAGES_PARTITION_3 = [
    "Tutorials" => [
        "Absolute Astrometry" => [
            "Proper Motion Anomaly" => "pma.md",
            "Hipparcos IAD" => "hipparcos.md",
            "Joint Gaia-Hipparcos (G23H)" => "g23h.md",
            "G23H Full Example" => "g23h-example.md",
            "Gaia DR4 Epoch Astrometry" => "gaia-iad.md",
            "Gaia DR4 Simulation" => "gaia-dr4-simulation.md",
        ],
        "Joint Models" => [
            "Calculate Detection Limits" => "limits.md",
        ],
    ],
]

# Partition 4: RV + Joint Models + Bayesian (heavy - astrom-pma-rv, rv tutorials)
const PAGES_PARTITION_4 = [
    "Tutorials" => [
        "Radial Velocity" => [
            "Basic RV Fit" => "rv-1.md",
            "Gaussian Process" => "rv-gp.md",
            "Multiple Planets" => "rv-multi-planet.md",
            "Relative RV Data" => "fit-rv-rel.md",
        ],
        "Joint Models" => [
            "Astrometry, PMA, and RV" => "astrom-pma-rv.md",
            "RV and Relative Astrometry" => "fit-rv-astrom.md",
            "RV and Proper Motion Anomaly" => "rv.md",
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
]

# Full pages structure for local builds and final merge
const PAGES_ALL = [
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
            "Joint Gaia-Hipparcos (G23H)" => "g23h.md",
            "G23H Full Example" => "g23h-example.md",
            "Gaia DR4 Epoch Astrometry" => "gaia-iad.md",
            "Gaia DR4 Simulation" => "gaia-dr4-simulation.md",
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
        "Architecture Overview" => "dev/architecture.md"
    ]
]

# =============================================================================
# Select pages based on DOCS_PARTITION environment variable
# =============================================================================

partition = get(ENV, "DOCS_PARTITION", "all")

pages_to_build = if partition == "1"
    @info "Building partition 1: Getting Started + Documentation"
    PAGES_PARTITION_1
elseif partition == "2"
    @info "Building partition 2: Relative Astrometry + Images"
    PAGES_PARTITION_2
elseif partition == "3"
    @info "Building partition 3: Absolute Astrometry + Limits"
    PAGES_PARTITION_3
elseif partition == "4"
    @info "Building partition 4: RV + Joint Models + Bayesian"
    PAGES_PARTITION_4
else
    @info "Building all pages"
    PAGES_ALL
end

# For partition builds, we skip deploydocs (handled by merge job)
skip_deploy = partition != "all"

makedocs(
    sitename="Octofitter.jl",
    pages = pages_to_build,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold=nothing
    ),
    pagesonly=true,
    warnonly=:example_block
)

# Only deploy when building all pages (final merge job)
if !skip_deploy
    deploydocs(
        repo = "github.com/sefffal/Octofitter.jl.git",
        devbranch = "main",
        push_preview = true
    )
end
