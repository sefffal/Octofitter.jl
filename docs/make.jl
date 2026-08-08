using Pkg
cd(@__DIR__)
Pkg.activate(".")

# Draft mode (`julia --project=docs docs/make.jl draft`, or OCTOFITTER_DOCS_DRAFT=true)
# skips `@example` execution so the build finishes in minutes instead of hours,
# while still failing on unresolved cross-references. Used by the pre-commit
# hook in .githooks/ to catch bad `@ref`s before they reach CI.
const DRAFT = "draft" in ARGS || get(ENV, "OCTOFITTER_DOCS_DRAFT", "false") == "true"

# The subpackages are loaded so that `api.md` can `@docs` their observation
# types (`ImageObs`, `LogLikelihoodMapObs`, `InterferometryObs`,
# `GRAVITYWideKPObs`) — Documenter can only resolve a docstring for a module
# that is in scope here.
using Documenter, Octofitter, OctofitterRadialVelocity,
      OctofitterImages, OctofitterInterferometry

# CairoMakie is loaded even in draft mode: the Makie package extension
# provides some docstrings (e.g. `rvplot_animated`), and without it their
# `@docs` entries and `@ref`s fail to resolve.
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
            "Migrating to v9" => "v9-migration.md",
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
                "Gaia DR4 Pre-Release Data" => "gaia-dr4-prerelease.md",
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
                "Detection Completeness Mapping" => "completeness.md",
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
            "Radial velocity figures"=>"rvplot.md",
            "Loading and Saving Data" => "loading-saving.md",
            "Sampler" => "samplers.md",
            "Distributed Sampling" => "parallel-sampling.md",
            "Priors" => "priors.md",
            "Corrections & Data Provenance" => "corrections.md",
            "Derived Variables" => "derived.md",
            "Custom Likelihoods" => "custom-likelihood.md",
            "Kepler Solver" => "kepler.md",
            "Orbitize! Compatibility" => "compat-orbitize.md",
            "Full API Documentation" => "api.md"
        ],
        "Developer Documentation" => [
            "Architecture Overview" => "dev/architecture.md"
        ]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold=nothing
    ),
    pagesonly=true,
    warnonly=:example_block,
    draft=DRAFT,
)


if !DRAFT
    deploydocs(
        repo = "github.com/sefffal/Octofitter.jl.git",
        devbranch = "main",
        push_preview = true
    )
end
