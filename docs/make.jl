using Documenter, DirectDetections


makedocs(
    sitename="DirectDetections.jl",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Tutorials" => [
            "Fit Astrometry" => "modelling.md",
            # "Multiple Planets" => "multi-planets.md",
            # "" => "multi-planets.md",
        ]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)


deploydocs(
    repo = "github.com/sefffal/DirectDetections.jl.git",
)