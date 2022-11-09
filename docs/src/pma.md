# [Fit Astrometric Acceleration](@id fit-pma)

One of the features of DirectDetections.jl is support for proper motion anomaly / astrometric acceleration.
These data points are typically calculated by finding the difference between a long term proper motion of a star between the Hipparcos and GAIA catalogs, and their proper motion calculated within the windows of each catalog.

For Hipparcos/GAIA this gives four data points that can constrain the dynamical mass & orbits of planetary companions (assuming we subtract out the net trend).

You can specify these quantities manually, but the easiest way is to use the Hipparcos-GAIA Catalog of Accelerations (HGCA, [https://arxiv.org/abs/2105.11662](https://arxiv.org/abs/2105.11662)). Support for loading this catalog is built into DirectDetections.jl.

Let's look at the star and companion [HD 91312 A & B](https://arxiv.org/abs/2109.12124), discovered by SCExAO. We will use their published astrometry and proper motion anomaly extracted from the HGCA.

The first step is to find the GAIA source ID for your object. For HD 91312, SIMBAD tells us the GAIA DR2 ID is `756291174721509376` (which we will assume is the same in eDR3).

## Planet Model
For this model, we will have to add the variable `mass` as a prior.
The units used on this variable are Jupiter masses, in contrast to `M`, the primary's mass, in solar masses.  A reasonable uninformative prior for `mass` is `Uniform(0,1000)` or `LogUniform(1,1000)` depending on the situation.

Initial setup:
```julia
using DirectDetections, Distributions, Plots
```


```julia
@named B = Planet{VisualOrbit}(
    Variables(
        a = Uniform(1, 25),
        e = Beta(2,15),
        τ = Uniform(0, 1),
        ω = Uniform(0, 2pi),
        i = Sine(), # The Sine() distribution is defined by DirectDetections
        Ω = Uniform(0, pi),
        mass = Uniform(0.5, 2000),
        # Anoter option would be:
        # mass = LogUniform(0.5, 2000),
    ),
    Astrometry(
        (epoch=mjd("2016-12-15"), ra=133., dec=-174., σ_ra=07.0, σ_dec=07.),
        (epoch=mjd("2017-03-12"), ra=126., dec=-176., σ_ra=04.0, σ_dec=04.),
        (epoch=mjd("2017-03-13"), ra=127., dec=-172., σ_ra=04.0, σ_dec=04.),
        (epoch=mjd("2018-02-08"), ra=083., dec=-133., σ_ra=10.0, σ_dec=10.),
        (epoch=mjd("2018-11-28"), ra=058., dec=-122., σ_ra=10.0, σ_dec=20.),
        (epoch=mjd("2018-12-15"), ra=056., dec=-104., σ_ra=08.0, σ_dec=08.),
    )
)
```

The `@named` macro simply passes the variable name you give as an arugment to the function i.e. `Planet(..., name=:B)`. This ensures the parameter names in the output are consistent with the code.


## System Model & Specifying Proper Motion Anomaly
Now that we have our planet model, we create a system model to contain it.

```julia
@named HD91312 = System(
    Variables(
        M = LogNormal(1.61, 1),
        plx = gaia_plx(gaia_id=756291174721509376),
    ),  
    ProperMotionAnomHGCA(gaia_id=756291174721509376),
    B,
)
```

We specify priors on `M` and `plx` as usual, but here we use the `gaia_plx` helper function to read the parallax and uncertainty directly from the HGCA using its source ID.

After the priors, we add the proper motion anomaly measurements from the HGCA. If this is your first time running this code, you will be prompted to automatically download and cache the catalog which may take around 30 seconds.


## Sampling from the posterior
Ssample from our model as usual:

```julia
chain = DirectDetections.hmc(
    HD91312, 0.85,
    adaptation =  1_000,
    iterations = 10_000,
)
```

Output:
```        
┌ Info: Guessing a good starting location by sampling from priors
└   N = 50000
┌ Info: Found initial stepsize
└   initial_ϵ = 0.0125
[ Info: Will adapt step size and mass matrix
[ Info: progress logging is enabled globally
Sampling100%|███████████████████████████████| Time: 0:00:24
  iterations:                    30000
  n_steps:                       7
  is_accept:                     true
  acceptance_rate:               0.863239872503996
  log_density:                   -51.71436291701389
  hamiltonian_energy:            54.325211045874674
  hamiltonian_energy_error:      0.19872154840668088
  max_hamiltonian_energy_error:  0.3162241419617615
  tree_depth:                    3
  numerical_error:               false
  step_size:                     0.5200876630871248
  nom_step_size:                 0.5200876630871248
  is_adapt:                      false
  mass_matrix:                   DenseEuclideanMetric(diag=[0.002497517511145372, 0.02 ...])
[ Info: Resolving derived variables
[ Info: Constructing chains
Sampling report:
mean_accept = 0.7834391801902499
num_err_frac = 0.017433333333333332
mean_tree_depth = 2.8884666666666665
max_tree_depth_frac = 0.0
Chains MCMC chain (30000×9×1 Array{Float64, 3}):

Iterations        = 1:1:30000
Number of chains  = 1
Samples per chain = 30000
Wall duration     = 24.15 seconds
Compute duration  = 24.15 seconds
parameters        = M, plx, B[a], B[e], B[τ], B[ω], B[i], B[Ω], B[mass]

Summary Statistics
  parameters       mean       std   naive_se      mcse          ess      rhat   ess_per_sec 
      Symbol    Float64   Float64    Float64   Float64      Float64   Float64       Float64

           M     1.6107    0.0499     0.0003    0.0003   21449.1180    1.0000      888.1255
         plx    29.1444    0.1396     0.0008    0.0010   23713.4280    1.0000      981.8818
        B[a]     6.7402    0.0890     0.0005    0.0010    7643.8317    1.0000      316.5017
        B[e]     0.2043    0.0449     0.0003    0.0005    5958.7172    1.0000      246.7276
        B[τ]     1.1459    0.0254     0.0001    0.0003    9732.8386    1.0000      402.9994
        B[ω]    -0.3675    0.0984     0.0006    0.0011    8379.4363    1.0000      346.9602
        B[i]     1.4706    0.0280     0.0002    0.0003   10928.8907    1.0000      452.5233
        B[Ω]    -0.6672    0.0147     0.0001    0.0001   19642.2568    1.0001      813.3103
     B[mass]   245.0888   69.1839     0.3994    1.2972    2404.3514    1.0000       99.5549

Quantiles
  parameters       2.5%      25.0%      50.0%      75.0%      97.5% 
      Symbol    Float64    Float64    Float64    Float64    Float64

           M     1.5131     1.5777     1.6109     1.6441     1.7083
         plx    28.8718    29.0491    29.1440    29.2396    29.4166
        B[a]     6.5706     6.6794     6.7398     6.7994     6.9165
        B[e]     0.1229     0.1730     0.2023     0.2332     0.2989
        B[τ]     1.1036     1.1281     1.1433     1.1608     1.2032
        B[ω]    -0.5503    -0.4335    -0.3705    -0.3075    -0.1589
        B[i]     1.4133     1.4525     1.4715     1.4902     1.5217
        B[Ω]    -0.6967    -0.6769    -0.6670    -0.6573    -0.6390
     B[mass]   159.7716   196.9441   228.3291   274.9650   431.3844
```

This takes about a minute on the first run due to JIT startup latency; subsequent runs are very quick even on e.g. an older laptop.

## Analysis

The first step is to look at the table output above generated by MCMCChains.jl.
The `rhat` column gives a convergence measure. Each parameter should have an `rhat` very close to 1.000.
If not, you may need to run the model for more iterations or tweak the parameterization of the model to improve sampling.
The `ess` column gives an estimate of the effective sample size.
The `mean` and `std` columns give the mean and standard deviation of each parameter.

The second table summarizes the 2.5, 25, 50, 75, and 97.5 percentiles of each parameter in the model.

Since this chain is well converged, we can begin examining the results.
Use the `plotmodel` function to display orbits from the posterior against the input data:


```julia
plotmodel(chain, color=:mass, pma_scatter=:mass)
```
[![mass histogram](assets/pma-astrometry-posterior.png)](assets/pma-astrometry-posterior.svg)


### Pair Plot
If we wish to examine the covariance between parameters in more detail, we can construct a pair-plot (aka. corner plot).

For a quick look, you can just run `corner(chain)`, but for more professional output you may wish to customize the labels, units, unit labels, etc:


```julia
##Create a corner plot / pair plot.
# We can access any property from the chain specified in Variables
using PairPlots
table = (;
    a=         chain["B[a]"],
    M=         chain["M"],
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
    "(M_\\odot)",
    "(M_{jup})",
    "",
    "(\\degree)",
    "(\\degree)",
    "(\\degree)",
    "",
]
corner(table, labels, units)
```
[![pair plot](assets/pma-astrometry-mass-corner.png)](assets/pma-astrometry-mass-corner.svg)


## Fitting Astrometric Acceleration Only
If you wish to look at the possible locations/masses of a planet around a star using onnly GAIA/HIPPARCOS,
you can follow a simplified approach.

As a start, you can restrict the orbital parameters to just semi-major axis, epoch of periastron passage, and mass.

```julia
@named b = Planet{VisualOrbit}(
    Variables(
        a = LogUniform(0.1, 100),
        τ = Uniform(0, 1),
        mass = LogUniform(1, 2000),
        i = Returns(0),
        e = Returns(0),
        Ω = Returns(0),
        ω = Returns(0),
    )
)
@named HD91312 = System(
    Variables(
        M = TruncatedNormal(1.61, 0.2, 0, Inf),
        plx = gaia_plx(gaia_id=756291174721509376),
    ),  
    ProperMotionAnomHGCA(gaia_id=756291174721509376),
    b,
)
```

This models assumes a circular, face-on orbit.

The Julia helper function `Returns(x)` is just a function that always returns the same value no matter the arguments. It  is equivalent to `(args...; kwargs...) -> x`.

```julia
chains = DirectDetections.hmc(
    HD91312, 0.85,
    MCMCThreads(),
    num_chains=4,
    adaptation =  1_000,
    iterations = 50_000,
)
```
```
┌ Info: Guessing a good starting location by sampling from priors
└   N = 50000
┌ Warning: The current proposal will be rejected due to numerical error(s).
│   isfinite.((θ, r, ℓπ, ℓκ)) = (true, false, false, false)
└ @ AdvancedHMC ~/.julia/packages/AdvancedHMC/HQHnm/src/hamiltonian.jl:47
┌ Warning: The current proposal will be rejected due to numerical error(s).
│   isfinite.((θ, r, ℓπ, ℓκ)) = (true, false, false, false)
└ @ AdvancedHMC ~/.julia/packages/AdvancedHMC/HQHnm/src/hamiltonian.jl:47
┌ Info: Found initial stepsize
└   initial_ϵ = 0.0001953125
[ Info: Will adapt step size and mass matrix
[ Info: progress logging is enabled globally
[ Info: Sampling compete. Building chains.
Sampling report for chain 1:
mean_accept =         0.9386094717407892
num_err_frac =        0.0
mean_tree_depth =     5.732058823529412
max_tree_depth_frac = 0.0

Sampling report for chain 2:
mean_accept =         0.9562432650363307
num_err_frac =        0.0
mean_tree_depth =     5.930529411764706
max_tree_depth_frac = 0.0

Sampling report for chain 3:
mean_accept =         0.941867350117908
num_err_frac =        0.0
mean_tree_depth =     5.783862745098039
max_tree_depth_frac = 0.0

Sampling report for chain 4:
mean_accept =         0.9504366454930078
num_err_frac =        0.0
mean_tree_depth =     5.86243137254902
max_tree_depth_frac = 0.0

Chains MCMC chain (51000×9×4 Array{Float64, 3}):

Iterations        = 1:1:51000
Number of chains  = 4
Samples per chain = 51000
Wall duration     = 64.1 seconds
Compute duration  = 64.1 seconds
parameters        = M, plx, b[a], b[τ], b[mass], b[i], b[e], b[Ω], b[ω]

Summary Statistics
  parameters      mean       std   naive_se      mcse           ess      rhat   ess_per_sec 
      Symbol   Float64   Float64    Float64   Float64       Float64   Float64       Float64 

           M    1.6085    0.2006     0.0004    0.0005   132909.5151    1.0000     2073.3743
         plx   29.1452    0.1408     0.0003    0.0004   123574.6795    1.0000     1927.7519
        b[a]    1.1321    0.0477     0.0001    0.0001   128054.3813    1.0000     1997.6348
        b[τ]    0.8804    0.0027     0.0000    0.0000   129890.5148    1.0000     2026.2783
     b[mass]   74.0519    7.6350     0.0169    0.0208   124168.8358    1.0000     1937.0207
        b[i]    0.0000    0.0000     0.0000    0.0000           NaN       NaN           NaN
        b[e]    0.0000    0.0000     0.0000    0.0000           NaN       NaN           NaN
        b[Ω]    0.0000    0.0000     0.0000    0.0000           NaN       NaN           NaN
        b[ω]    0.0000    0.0000     0.0000    0.0000           NaN       NaN           NaN

Quantiles
  parameters      2.5%     25.0%     50.0%     75.0%     97.5% 
      Symbol   Float64   Float64   Float64   Float64   Float64 

           M    1.2152    1.4730    1.6086    1.7439    2.0012
         plx   28.8705   29.0501   29.1454   29.2403   29.4210
        b[a]    1.0329    1.1013    1.1342    1.1651    1.2198
        b[τ]    0.8751    0.8786    0.8804    0.8822    0.8856
     b[mass]   59.3819   68.8615   73.9459   79.1446   89.3102
        b[i]    0.0000    0.0000    0.0000    0.0000    0.0000
        b[e]    0.0000    0.0000    0.0000    0.0000    0.0000
        b[Ω]    0.0000    0.0000    0.0000    0.0000    0.0000
        b[ω]    0.0000    0.0000    0.0000    0.0000    0.0000

```

With such simple models, the mean tree depth is often very low and sampling proceeds very quickly.

A good place to start is a histogram of planet mass vs. semi-major axis:
```julia
histogram2d(chains["b[a]"], chains["b[mass]"], color=:plasma, xguide="sma (au)", yguide="mass (Mjup)")
```
[![2d histogram](assets/pma-a-vs-mass.svg)](assets/pma-a-vs-mass.svg)


You can also visualize the orbits and proper motion:
```julia
plotmodel(chains)
```
[![2d histogram](assets/pma-model.png)](assets/pma-model.svg)
