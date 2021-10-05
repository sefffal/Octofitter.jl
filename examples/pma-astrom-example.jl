using DirectDetections, Distributions, Plots
using DirectOrbits: mjd

@named B = DirectDetections.Planet(
    Deterministic(
        e = (sys, pl) -> 10^pl.loge,
        a = (sys, pl) -> 10^pl.loga,
        mass = (sys, pl) -> 10^pl.logm,
    ),
    Priors(
        # a = Uniform(0, 25),
        # e = TruncatedNormal(0.0, 0.2, 0, 1.0),
        # Note: priors with sharp edges (e.g. Uniform priors) are challenging for HMC samplers.
        # An alternative could be wide Gaussians, for example.
        τ = Normal(0.5, 1),
        ω = Normal(pi, 2pi),
        i = Normal(pi, 2pi),
        Ω = Normal(pi, 2pi),
        # mass = Uniform(0, 1),

        # Reparameterize a few properties for better sampling
        loge = TruncatedNormal(-2, 1.5, -Inf, 0),
        loga = Uniform(-1, 1.5),
        logm = Uniform(-2, 0),

    ),
    Astrometry(
        (epoch=mjd("2016-12-15"), ra=0.133*1e3, dec=-0.174*1e3, σ_ra=0.007*1e3, σ_dec=0.007*1e3),
        (epoch=mjd("2017-03-12"), ra=0.126*1e3, dec=-0.176*1e3, σ_ra=0.004*1e3, σ_dec=0.004*1e3),
        (epoch=mjd("2017-03-13"), ra=0.127*1e3, dec=-0.172*1e3, σ_ra=0.004*1e3, σ_dec=0.004*1e3),
        (epoch=mjd("2018-02-08"), ra=0.083*1e3, dec=-0.133*1e3, σ_ra=0.010*1e3, σ_dec=0.010*1e3),
        (epoch=mjd("2018-11-28"), ra=0.058*1e3, dec=-0.122*1e3, σ_ra=0.010*1e3, σ_dec=0.020*1e3),
        (epoch=mjd("2018-12-15"), ra=0.056*1e3, dec=-0.104*1e3, σ_ra=0.008*1e3, σ_dec=0.008*1e3),
    )
)

@named HD82134 = System(
    # Deterministic(
    #     j = Returns(1),
    #     k = Returns(1)
    # ),
    Priors(
        μ = Normal(1.61, 0.05),
        plx = gaia_plx(gaia_id=756291174721509376),
    ),  
    ProperMotionAnomHGCA(gaia_id=756291174721509376),
    B,
)

## Sampling from the posterior

chain, stats = DirectDetections.hmc(
    HD82134,
    adaptation =  1_000,
    iterations = 50_000,
    tree_depth =     12,
);

## Plot the input data and samples from the posterior
plotmodel(chain, HD82134, clims=(5, 20))

## Create a corner plot / pair plot.
# We can access any property from the chain specified in Priors or in Deterministic.
using PairPlots
table = (;
    a=         chain["B[a]"],
    μ=         chain["μ"],
    m=         chain["B[mass]"],
    e=         chain["B[e]"],
    i=rad2deg.(chain["B[i]"]),
    Ω=rad2deg.(chain["B[Ω]"]),
    ω=rad2deg.(chain["B[ω]"]),
    τ=         chain["B[τ]"],
)
labels=[
    "a",
    "\\mu",
    "m",
    "e",
    "i",
    "\\Omega",
    "\\omega",
    "\\tau",
]
units = [
    "(au)",
    "(_\\odot)",
    "(_\\odot)",
    "",
    "(\\degree)",
    "(\\degree)",
    "(\\degree)",
    "",
]
corner(table, labels, units);
nothing

