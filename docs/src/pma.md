# [Fit Proper Motion Anomaly](@id fit-pma)

Octofitter.jl supports fitting orbit models to astrometric motion in the form of GAIA-Hipparcos proper motion anomaly (HGCA; [https://arxiv.org/abs/2105.11662](https://arxiv.org/abs/2105.11662)).
These data points are calculated by finding the difference between a long term proper motion of a star between the Hipparcos and GAIA catalogs, and their proper motion calculated within the windows of each catalog. This gives four data points that can constrain the dynamical mass & orbits of planetary companions (assuming we subtract out the net trend).

If your star of interest is in the HGCA, all you need is it's GAIA DR3 ID number. You can find this number by searching for your target on [SIMBAD](http://simbad.cds.unistra.fr).

For this tutorial, we will examine the star and companion [HD 91312 A & B](https://arxiv.org/abs/2109.12124) discovered by SCExAO. We will use their published astrometry and proper motion anomaly extracted from the HGCA.

We will also perform a model comparison: we will fit the same model to four different subsets of data to see how each dataset are impacting the final constraints. This is an important consistency check, especially with proper motion / absolute astrometry data which be susceptible to systematic errors.

The first step is to find the GAIA source ID for your object. For HD 91312, SIMBAD tells us the GAIA DR3 ID is `756291174721509376`.

## Fitting Astrometric Motion Only


Initial setup:
```@example 1
using Octofitter, Distributions, Random
```

We begin by finding orbits that are consistent with the astrometric motion. Later, we will add in relative astrometry to the fit from direct imaging to further constrain the planet's orbit and mass.

Compared to previous tutorials, we will now have to add a few additional variables to our model. The first is a prior on the mass of the companion, called `mass`.
The units used on this variable are Jupiter masses, in contrast to `M`, the primary's mass, in solar masses.  A reasonable uninformative prior for `mass` is `Uniform(0,1000)` or `LogUniform(1,1000)` depending on the situation.

For this model, we also want to place a prior on the host star mass rather than system total mass. For exoplanets there is litte difference between these two values, but in this example we have a reasonably informative prior on the host mass, and know from the paper that the companion is has a non-neglible effect on the total system mass.

To make this parameterization change, we specify priors on both masses in the `@system` block, and connect it to the planet.

### Retrieving the HGCA
To start, we retrieve the HGCA data for this object.
```julia
hgca_like = HGCALikelihood(
    gaia_id=3937211745905473024,
    variables=@variables begin
        # Optional: flux ratio for luminous companions, one entry per companion
        # fluxratio ~ Product([Uniform(0, 1), Uniform(0, 1), ])  # uncomment if needed for unresolved companions
    end
)
```

You can optionally provide flux ratio priors in the variables block to represent the flux ratio of the companions to the host, if you don't want to approximate it as zero. This is to handle luminous companions that are unresolved by gaia.


If you're in a hurry, and you're study orbits with periods much longer than the mission durations of Gaia or Hipparcos (>> 4 years) then you might consider using a faster approximation that the Gaia and Hipparcos measurements were instantaneous. You can do so as follows:

```@example 1
hgca_like = HGCAInstantaneousLikelihood(gaia_id=756291174721509376, N_ave=1) 
```
`N_ave` is an optional argument to control over how many epochs the measurements are approximated, e.g. N_ave=10 implies that the position and proper motion was measured instantaneously 10 times over each mission and averaged.

### Planet Model

```@example 1
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    variables=@variables begin
        a ~ LogUniform(0.1,20)
        e ~ Uniform(0,0.999)
        ω ~ UniformCircular()
        i ~ Sine() # The Sine() distribution is defined by Octofitter
        Ω ~ UniformCircular()

        mass = super.M_sec

        M = super.M
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ, 57423.0; M, e, a, i, ω, Ω) # epoch of GAIA measurement
    end
)
```


### System Model & Specifying Proper Motion Anomaly
Now that we have our planet model, we create a system model to contain it.

We specify priors on `plx` as usual, but here we use the `gaia_plx` helper function to read the parallax and uncertainty directly from the HGCA catalog using its source ID.

We also add parameters for the star's long term proper motion. This is usually close to the long term trend between the Hipparcos and GAIA measurements. If you're not sure what to use here, try `Normal(0, 1000)`; that is, assume a long-term proper motion of 0 +- 1000 milliarcseconds / year.


```@example 1
sys = System(
    name="HD91312_pma",
    companions=[planet_b],
    likelihoods=[hgca_like],
    variables=@variables begin
        M_pri ~ truncated(Normal(1.61, 0.1), lower=0.1) # Msol
        M_sec ~ LogUniform(0.5, 1000) # MJup
        M = M_pri + M_sec*Octofitter.mjup2msol # Msol

        plx ~ gaia_plx(gaia_id=756291174721509376)
                
        # Priors on the center of mass proper motion
        pmra ~ Normal(-137, 10)
        pmdec ~ Normal(2,  10)
    end
)

model_pma = Octofitter.LogDensityModel(sys)
```


After the priors, we add the proper motion anomaly measurements from the HGCA. If this is your first time running this code, you will be prompted to automatically download and cache the catalog which may take around 30 seconds.


### Sampling from the posterior (PMA only)

Because proper motion anomaly data is quite sparse, it can often produce multi-modal posteriors. If your orbit already has several relative astrometry or RV data points, this is less of an issue. But in many cases it is recommended to use the `Pigeons.jl` sampler instead of Octofitter's default. This sampler is less efficient for unimodal distributions, but is more robust at exploring posteriors with distinct, widely separated peaks. 

To install and use `Pigeons.jl` with Octofitter, type `using Pigeons` at in the terminal and accept the prompt to install the package. You may have to restart Julia.

!!! note
    `octofit_pigeons` scales very well across multiple cores. Start julia with `julia --threads=auto` to make sure you have multiple threads available for sampling.

We now sample from our model using Pigeons:
```@example 1
using Pigeons
chain_pma, pt = octofit_pigeons(model_pma, n_rounds=13, explorer=SliceSampler()) 
display(chain_pma)
```

Note that `octofit_pigeons` took somewhat longer to run than `octofit` typically does; however, as we will see, it sampled successfully from severally completely disconnected modes in the posterior. That makes it a good fit for sampling from proper motion anomaly and relative astrometry with limited orbital coverage.

### Analysis

The first step is to look at the table output above generated by MCMCChains.jl.
The `rhat` column gives a convergence measure. Each parameter should have an `rhat` very close to 1.000.
If not, you may need to run the model for more iterations or tweak the parameterization of the model to improve sampling.
The `ess` column gives an estimate of the effective sample size.
The `mean` and `std` columns give the mean and standard deviation of each parameter.

The second table summarizes the 2.5, 25, 50, 75, and 97.5 percentiles of each parameter in the model.


### Pair Plot
If we wish to examine the covariance between parameters in more detail, we can construct a pair-plot (aka. corner plot).

```@example 1
# Create a corner plot / pair plot.
# We can access any property from the chain specified in Variables
using CairoMakie: Makie
using PairPlots
octocorner(model_pma, chain_pma, small=true)
```

Notice how there are completely separated peaks? The default Octofitter sample (Hamiltonian Monte Carlo) is capabale of jumping 2-3σ gaps between modes, but such widely separated peaks can cause issues (hence why we used Pigeons in this example).

### Posterior Mass vs. Semi-Major Axis

Given that this posterior is quite unconstrained, it is useful to make a simplified plot marginalizing over all orbital 
parameters besides semi-major axis. We can do this using PairPlots.jl:
```@example 1
using CairoMakie, PairPlots
pairplot(
    (; a=chain_pma["b_a"][:], mass=chain_pma["b_mass"][:]) =>
        (
            PairPlots.Scatter(color=:red),
            PairPlots.MarginHist(),
            PairPlots.MarginQuantileLines(),
            PairPlots.MarginQuantileText(),
        ),
    labels=Dict(:mass=>"mass [Mⱼᵤₚ]", :a=>"sma. [au]"),
    axis = (;
        a = (;
            scale=Makie.pseudolog10,
            ticks=2 .^ (0:1:6)
        )
    )
)
```
