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

Initial setup:
```julia
using DirectDetections, Distributions, Plots
```


```julia
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

```julia
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

```julia
chain, stats = DirectDetections.hmc(
    HD91312, 0.65,
    adaptation =  1_000,
    iterations = 10_000,
);
display(chain)
```

Output:
```
Sampling report:
mean_accept = 0.8538973684799083
num_err_frac = 0.0102
mean_tree_depth = 3.2545
max_tree_depth_frac = 0.0
Chains MCMC chain (10000×9×1 Array{Float64, 3}):

Iterations        = 1:1:10000
Number of chains  = 1
Samples per chain = 10000
Wall duration     = 11.43 seconds
Compute duration  = 11.43 seconds
parameters        = μ, plx, B[a], B[e], B[τ], B[ω], B[i], B[Ω], B[mass]

Summary Statistics
  parameters       mean       std   naive_se      mcse         ess      rhat   ess_per_sec ⋯
      Symbol    Float64   Float64    Float64   Float64     Float64   Float64       Float64 ⋯

           μ     1.6123    0.0496     0.0005    0.0006   6567.9890    1.0001      574.5769 ⋯
         plx    29.1438    0.1411     0.0014    0.0015   8176.9434    0.9999      715.3305 ⋯
        B[a]     6.7216    0.0928     0.0009    0.0017   3061.4284    1.0002      267.8181 ⋯
        B[e]     0.2056    0.0440     0.0004    0.0008   3378.2624    1.0003      295.5352 ⋯
        B[τ]    -0.8507    0.0265     0.0003    0.0005   3789.5866    0.9999      331.5184 ⋯
        B[ω]    -6.6247    0.1028     0.0010    0.0019   3071.5674    1.0009      268.7050 ⋯
        B[i]     1.4768    0.0292     0.0003    0.0005   3296.3292    1.0001      288.3675 ⋯
        B[Ω]    11.9009    0.0147     0.0001    0.0002   5376.1724    1.0000      470.3151 ⋯
     B[mass]   271.7459   92.4623     0.9246    2.4161   1195.0382    1.0022      104.5436 ⋯

Quantiles
  parameters       2.5%      25.0%      50.0%      75.0%      97.5% 
      Symbol    Float64    Float64    Float64    Float64    Float64

           μ     1.5160     1.5786     1.6120     1.6462     1.7085
         plx    28.8621    29.0483    29.1438    29.2381    29.4209
        B[a]     6.5367     6.6589     6.7218     6.7851     6.9002
        B[e]     0.1241     0.1753     0.2034     0.2352     0.2973
        B[τ]    -0.8941    -0.8695    -0.8538    -0.8351    -0.7918
        B[ω]    -6.8125    -6.6942    -6.6294    -6.5626    -6.4023
        B[i]     1.4156     1.4579     1.4775     1.4967     1.5310
        B[Ω]    11.8709    11.8911    11.9012    11.9110    11.9284
     B[mass]   165.5598   209.0807   246.8737   306.0517   536.6443
```

This is quite quick even on an older laptop, single core.

## Analysis

We can use the `plotmodel` helper to summarize the results.


```julia
plotmodel(chain, HD91312, color=:mass, pma_scatter=:mass)
```
[![mass histogram](assets/pma-astrometry-posterior.png)](assets/pma-astrometry-posterior.svg)


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
