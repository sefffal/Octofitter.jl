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
        (epoch=5072.,  ra=-899., dec=-629., σ_ra=10., σ_dec=50.),
    )
)

@named HD82134 = System(
    Priors(
        μ = Normal(1.0, 0.01),
        plx =Normal(1000.2, 0.02),
    ),  
    X,
)

chains, stats = DirectDetections.hmc(
    HD82134;
    burnin=8_000,
    numwalkers=1,
    numsamples_perwalker=100_000,
);

##
plotmodel(chains[1], HD82134)

##
using PairPlots
table = (;
    a=chains[1].planets[1].a,
    e=chains[1].planets[1].e,
    i=rad2deg.(chains[1].planets[1].i),
    Ω=rad2deg.(chains[1].planets[1].Ω),
    ω=rad2deg.(chains[1].planets[1].ω),
    t=(chains[1].planets[1].τ),
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
corner(table, labels, units, plotscatter=false)#, hist2d_kwargs=(;nbins=15))

##
corner(table, labels, units)#, hist2d_kwargs=(;nbins=15))
cd(@__DIR__)
savefig("../docs/src/assets/astrometry-corner-plot.png")