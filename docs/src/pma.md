# [Fit Astrometric Acceleration](@id fit-pma)

One of the features of DirectDetections.jl is support for proper motion anomaly / astrometric acceleration.
These data points are typically calculated by finding the difference between a long term proper motion of a star between the Hipparcos and GAIA catalogs, and their proper motion calculated within the windows of each catalog.

For Hipparcos/GAIA this gives four data points that can constrain the dynamical mass & orbits of planetary companions (assuming we subtract out the net trend).

You can specify these quantities manually, but the easiest way is to use the Hipparcos-GAIA Catalog of Accelerations (HGCA, [https://arxiv.org/abs/2105.11662](https://arxiv.org/abs/2105.11662)). Support for loading this catalog is built into DirectDetections.jl.

Let's look at the star and companion [HD 91312 A & B](https://arxiv.org/abs/2109.12124), discovered by SCExAO. We will use their published astrometry and proper motion anomaly extracted from the HGCA.

The first step is to find the GAIA source ID for your object. For HD 91312, SIMBAD tells us the GAIA DR2 ID is `756291174721509376` (which we will assume is the same in eDR3).

## Planet Model
For this model, we will have to add the variable `mass` as a prior.
The units used on this variable are Jupiter masses, in contrast to `μ`, the primary's mass, in solar masses.  A reasonable uninformative prior for `mass` is `Uniform(0,1000)` or `LogUniform(1,1000)` depending on the situation.

```@setup
using DirectDetections, Distributions, Plots
```


```@example
@named B = DirectDetections.Planet(
    Priors(
        a = Uniform(1, 25),
        e = Beta(2,15),
        # Note: priors with sharp edges (e.g. Uniform priors) are challenging for HMC samplers.
        # An alternative could be wide Gaussians, for example.
        τ = Normal(0.5, 1),
        ω = Normal(pi, 2pi),
        i = Normal(pi, 2pi),
        Ω = Normal(pi, 2pi),
        # mass = LogUniform(0.5, 2000),
        mass = Uniform(0.5, 2000),
    ),
    Astrometry(
        (epoch=mjd("2016-12-15"), ra=133., dec=-174., σ_ra=07.0, σ_dec=07.),
        (epoch=mjd("2017-03-12"), ra=126., dec=-176., σ_ra=04.0, σ_dec=04.),
        (epoch=mjd("2017-03-13"), ra=127., dec=-172., σ_ra=04.0, σ_dec=04.),
        (epoch=mjd("2018-02-08"), ra=083., dec=-133., σ_ra=10.0, σ_dec=10.),
        (epoch=mjd("2018-11-28"), ra=058., dec=-122., σ_ra=10.0, σ_dec=20.),
        (epoch=mjd("2018-12-15"), ra=056., dec=-104., σ_ra=08.0, σ_dec=08.),
    )
);
```


## System Model & Specifying Proper Motion Anomaly
Now that we have our planet model, we create a system model to contain it.

```@example
@named HD91312 = System(
    Priors(
        μ = Normal(1.61, 0.05),
        plx = gaia_plx(gaia_id=756291174721509376),
    ),  
    ProperMotionAnomHGCA(gaia_id=756291174721509376),
    B,
)
```

We specify priors on `μ` and `plx` as usual, but here we use the `gaia_plx` helper function to read the parallax and uncertainty directly from the HGCA using its source ID.

After the priors, we add the proper motion anomaly measurements from the HGCA. If this is your first time running this code, you will be prompted to automatically download and cache the catalog which may take around 30 seconds.


## Sampling from the posterior
Ssample from our model as usual:

```@example
chain, stats = DirectDetections.hmc(
    HD91312, 0.65,
    adaptation =  1_000,
    iterations = 10_000,
);
display(chain)
```

This is quite quick even on an older laptop, single core.
## Analysis

We can use the `plotmodel` helper to summarize the results.


```@example
plotmodel(chain, HD91312, color=:mass, pma_scatter=:mass)
```
<!-- [![mass histogram](assets/pma-astrometry-posterior.png)](assets/pma-astrometry-posterior.svg) -->


### Pair Plot
Visualize all the parameters as a pair-plot:

```julia
##Create a corner plot / pair plot.
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
corner(table, labels, units)
```
[![pair plot](assets/pma-astrometry-mass-corner.png)](assets/pma-astrometry-mass-corner.svg)
