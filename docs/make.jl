using Documenter, DirectDetections


makedocs(
    sitename="DirectDetections.jl",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Tutorials" => [
            "Fit Astrometry" => "modelling.md",
            "Fit Images" => "images.md",
            "Fit Astrometric Acceleration" => "pma.md",
            "Derived Variables" => "derived.md",
            "Connecting Mass with Photometry" => "mass-photometry.md",
            # "Multiple Planets" => "multi-planets.md",
            # "" => "multi-planets.md",
        ],
        "Documentation" => [
            "Samplers" => "samplers.md",
            "Chains" => "chains.md"
        ]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)


deploydocs(
    repo = "github.com/sefffal/DirectDetections.jl.git",
    devbranch = "main"
)