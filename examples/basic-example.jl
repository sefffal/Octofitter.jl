using DirectDetections, Distributions, Plots

@named X = DirectDetections.Planet(
    Priors(
        a = TruncatedNormal(1, 0.5, 0, Inf),
        e = TruncatedNormal(0.0, 0.2, 0, 1.0),
        τ = Uniform(0, 1),
        ω = TruncatedNormal(deg2rad(250.), deg2rad(80.), 0, 2pi),
        i = TruncatedNormal(deg2rad(20.), deg2rad(10.), 0, pi),
        Ω = TruncatedNormal(deg2rad(200.), deg2rad(30.), 0, 2pi),
    ),
    Astrometry(
        (epoch=5000.,  ra=-364., dec=-1169., σ_ra=70., σ_dec=30.),
        (epoch=5014.,  ra=-493., dec=-1104., σ_ra=70., σ_dec=30.),
        (epoch=5072.,  ra=-899., dec=-629., σ_ra=20., σ_dec=50.),
    )
)

@named HD82134 = System(
    Priors(
        μ =  TruncatedNormal(1.0, 0.01, 0, Inf),
        plx = TruncatedNormal(1000.2, 0.02, 0, Inf),
    ),  
    X,
)

##
scatter(X.astrometry, label="X Astrometry", aspectratio=1)
xlims!(-1500,1500)
ylims!(-1500,1500)
##

chain = DirectDetections.hmc(
    HD82134,
    adaptation =   5_000,
    iterations =   5_000,
)

##
plotmodel(chain, color=:i)
xlims!(-1500,1500)
ylims!(-1500,1500)

## Or for more customized output:
using PairPlots
table = (;
    a=         chain["X[a]"],
    e=         chain["X[e]"],
    i=rad2deg.(chain["X[i]"]),
    Ω=rad2deg.(chain["X[Ω]"]),
    ω=rad2deg.(chain["X[ω]"]),
    τ=         chain["X[τ]"],
)
labels=[
    "a",
    "e",
    "i",
    "\\Omega",
    "\\omega",
    "\\tau",
]
units = [
    "(au)",
    "",
    "(\\degree)",
    "(\\degree)",
    "(\\degree)",
    "",
]
PairPlots.corner(table, labels, units)

##
cd(@__DIR__)
savefig("../docs/src/assets/astrometry-corner-plot.svg")