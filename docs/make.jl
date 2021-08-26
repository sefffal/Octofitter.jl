using Documenter, DirectDetections


makedocs(
    sitename="Direct Detections",
    pages = [
        "Getting Started" => "getting-started.md",
        "Tutorials" => [
            "Basic Model" => "modelling.md",
        ]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)


