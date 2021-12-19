using DirectDetections, Distributions, Plots

@named X = DirectDetections.Planet(
    Priors(
        a = Normal(1, 0.5),
        e = TruncatedNormal(0.0, 0.2, 0, 1.0),
        τ = Normal(0.5, 1),
        ω = Normal(deg2rad(250.), deg2rad(80.)),
        i = Normal(deg2rad(20.), deg2rad(10.)),
        Ω = Normal(deg2rad(200.), deg2rad(30.)),
    ),
    Astrometry(
        (epoch=5000.,  ra=-364., dec=-1169., σ_ra=70., σ_dec=30.),
        (epoch=5014.,  ra=-493., dec=-1104., σ_ra=70., σ_dec=30.),
        (epoch=5072.,  ra=-899., dec=-629., σ_ra=20., σ_dec=50.),
    )
)

@named HD82134 = System(
    Priors(
        μ = Normal(1.0, 0.01),
        plx =Normal(1000.2, 0.02),
    ),  
    X,
)

##

chain = DirectDetections.hmc(
    HD82134,
    adaptation =   5_000,
    iterations = 100_000,
);

##
plotmodel(chain, color=:e)

##
using PairPlots
corner(chain)

## Or for more customized output:
using PairPlots
table = (;
    a=         chain["X[a]"],
    H=         chain["X[GPI_H]"],
    e=         chain["X[e]"],
    i=rad2deg.(chain["X[i]"]),
    Ω=rad2deg.(chain["X[Ω]"]),
    ω=rad2deg.(chain["X[ω]"]),
    τ=         chain["X[τ]"],
)
labels=[
    "a",
    "H",
    "e",
    "i",
    "\\Omega",
    "\\omega",
    "\\tau",
]
units = [
    "(au)",
    "(arb.)",
    "",
    "(\\degree)",
    "(\\degree)",
    "(\\degree)",
    "",
]
corner(table, labels, units)

##
cd(@__DIR__)
savefig("../docs/src/assets/astrometry-corner-plot.svg")