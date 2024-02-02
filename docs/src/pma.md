# [Fit Proper Motion Anomaly](@id fit-pma)

Octofitter.jl supports fitting orbit models to astrometric motion in the form of GAIA-Hipparcos proper motion anomaly (HGCA; [https://arxiv.org/abs/2105.11662](https://arxiv.org/abs/2105.11662)).
These data points are calculated by finding the difference between a long term proper motion of a star between the Hipparcos and GAIA catalogs, and their proper motion calculated within the windows of each catalog. This gives four data points that can constrain the dynamical mass & orbits of planetary companions (assuming we subtract out the net trend).

If your star of interest is in the HGCA, all you need is it's GAIA DR3 ID number. You can find this number by searching for your target on [SIMBAD](http://simbad.cds.unistra.fr).

For this tutorial, we will examine the star and companion [HD 91312 A & B](https://arxiv.org/abs/2109.12124) discovered by SCExAO. We will use their published astrometry and proper motion anomaly extracted from the HGCA.

The first step is to find the GAIA source ID for your object. For HD 91312, SIMBAD tells us the GAIA DR3 ID is `756291174721509376`.

## Fitting Astrometric Motion Only


Initial setup:
```@example 1
using Octofitter, Distributions, Plots, Random
```

We begin by finding orbits that are consistent with the astrometric motion. Later, we will add in relative astrometry to the fit from direct imaging to further constrain the planet's orbit and mass.

Compared to previous tutorials, we will now have to add a few additional variables to our model. The first is a prior on the mass of the companion, called `mass`.
The units used on this variable are Jupiter masses, in contrast to `M`, the primary's mass, in solar masses.  A reasonable uninformative prior for `mass` is `Uniform(0,1000)` or `LogUniform(1,1000)` depending on the situation.

For this model, we also want to place a prior on the host star mass rather than system total mass. For exoplanets there is litte difference between these two values, but in this example we have a reasonably informative prior on the host mass, and know from the paper that the companion is has a non-neglible effect on the total system mass.

To make this parameterization change, we specify priors on both masses in the `@system` block, and connect it to the planet.


### Planet Model

```@example 1
@planet B Visual{KepOrbit} begin
    a ~ LogUniform(0.1,20)
    e ~ Uniform(0,0.999)
    ω ~ UniformCircular()
    i ~ Sine() # The Sine() distribution is defined by Octofitter
    Ω ~ UniformCircular()

    mass = system.M_sec

    τ ~ UniformCircular(1.0)
    P = √(B.a^3/system.M)
    tp =  B.τ*B.P*365.25 + 57737 # reference epoch for τ. Choose an MJD date near your data.
end
```


### System Model & Specifying Proper Motion Anomaly
Now that we have our planet model, we create a system model to contain it.

We specify priors on `plx` as usual, but here we use the `gaia_plx` helper function to read the parallax and uncertainty directly from the HGCA catalog using its source ID.

We also add parameters for the star's long term proper motion. This is usually close to the long term trend between the Hipparcos and GAIA measurements. If you're not sure what to use here, try `Normal(0, 1000)`; that is, assume a long-term proper motion of 0 +- 1000 milliarcseconds / year.

```@example 1
@system HD91312 begin
    
    M_pri ~ truncated(Normal(1.61, 0.1), lower=0) # Msol
    M_sec ~ LogUniform(0.5, 1000) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol

    plx ~ gaia_plx(gaia_id=756291174721509376)
            
    # Priors on the center of mass proper motion
    pmra ~ Normal(-137, 10)
    pmdec ~ Normal(2,  10)
end HGCALikelihood(gaia_id=756291174721509376) B

model = Octofitter.LogDensityModel(HD91312)
```


After the priors, we add the proper motion anomaly measurements from the HGCA. If this is your first time running this code, you will be prompted to automatically download and cache the catalog which may take around 30 seconds.


### Sampling from the posterior

Because proper motion anomaly data is quite sparse, it can often produce multi-modal posteriors. If your orbit already has several relative astrometry or RV data points, this is less of an issue. But in many cases it is recommended to use the `Pigeons.jl` sampler instead of Octofitter's default. This sampler is less efficient for unimodal distributions, but is more robust at exploring posteriors with distinct, widely separated peaks. 

To install and use `Pigeons.jl` with Octofitter, type `using Pigeons` at in the terminal and accept the prompt to install the package. You may have to restart Julia.

We now sample from our model using Pigeons:
```@example 1
using Pigeons
model = Octofitter.LogDensityModel(HD91312)
Random.seed!(1)
chain, pt = octofit_pigeons(model, n_rounds=12) 
display(chain)
```

Note that `octofit_pigeons` took considerably longer to run than `octofit` typically does; however, as we will see, it sampled successfully from severally completely disconnected modes in the posterior. That makes it a good fit for sampling from proper motion anomaly and relative astrometry with limited orbital coverage.

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
octocorner(model, chain, small=true)
```

Notice how there are completely separated peaks? The default Octofitter sample (Hamiltonian Monte Carlo) is capabale of jumping 2-3σ gaps between modes, but such widely separated peaks can cause issues (hence why we used Pigeons in this example).

### Posterior Mass vs. Semi-Major Axis

Given that this posterior is quite unconstrained, it is useful to make a simplified plot marginalizing over all orbital 
parameters besides semi-major axis. We can do this using PairPlots.jl:
```@example 1
using CairoMakie, PairPlots
pairplot(
    (; a=chain["B_a"][:], mass=chain["B_mass"][:]) =>
        (
            PairPlots.Scatter(color=:red),
            PairPlots.MarginHist(),
            PairPlots.MarginConfidenceLimits()
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

## Adding Relative Astrometry

The first orbit fit to only Hipparcos/GAIA data was very unconstrained. We will now add six epochs of
relative astrometry (measured from direct images) gathered from the [discovery paper](https://arxiv.org/abs/2109.12124).


```@example 1
astrom_like = PlanetRelAstromLikelihood(
    (epoch=mjd("2016-12-15"), ra=133., dec=-174., σ_ra=07.0, σ_dec=07., cor=0.2),
    (epoch=mjd("2017-03-12"), ra=126., dec=-176., σ_ra=04.0, σ_dec=04., cor=0.3),
    (epoch=mjd("2017-03-13"), ra=127., dec=-172., σ_ra=04.0, σ_dec=04., cor=0.1),
    (epoch=mjd("2018-02-08"), ra=083., dec=-133., σ_ra=10.0, σ_dec=10., cor=0.4),
    (epoch=mjd("2018-11-28"), ra=058., dec=-122., σ_ra=10.0, σ_dec=20., cor=0.3),
    (epoch=mjd("2018-12-15"), ra=056., dec=-104., σ_ra=08.0, σ_dec=08., cor=0.2),
)
Plots.plot(astrom_like, lims=:symmetric)
```

We use the same model as before, but now condition the planet model `B` on the astrometry data by
adding `astrom_like` to the end of the `@planet` defintion.
```@example 1
@planet B Visual{KepOrbit} begin
    a ~ LogUniform(0.1,200)
    e ~ Uniform(0,0.999)
    ω ~ UniformCircular()
    i ~ Sine()
    Ω ~ UniformCircular()

    mass = system.M_sec

    τ ~ UniformCircular(1.0)
    P = √(B.a^3/system.M)
    tp =  B.τ*B.P*365.25 + 57737 # reference epoch for τ. Choose an MJD date near your data.
end astrom_like # Note the relative astrometry added here!

@system HD91312 begin
    
    M_pri ~ truncated(Normal(1.61, 0.1), lower=0)
    M_sec ~ LogUniform(0.5, 1000) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol

    plx ~ gaia_plx(gaia_id=756291174721509376)
            
    # Priors on the centre of mass proper motion
    pmra ~ Normal(-137, 10)
    pmdec ~ Normal(2,  10)
end HGCALikelihood(gaia_id=756291174721509376) B

model = Octofitter.LogDensityModel(HD91312)
```


Sample as before:
```@example 1
using Pigeons
model = Octofitter.LogDensityModel(HD91312)
Random.seed!(1)
chain, pt = octofit_pigeons(model, n_rounds=12) 
display(chain)
```


The mass vs. semi-major axis posterior is now much more constrained:
```@example 1
using CairoMakie, PairPlots
pairplot(
    (; a=chain["B_a"][:], mass=chain["B_mass"][:]) =>
        (
            PairPlots.Scatter(color=:red),
            PairPlots.MarginHist(),
            PairPlots.MarginConfidenceLimits()
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


It is now useful to display the orbits projected onto the plane of the sky using `octoplot`. This function produces a nine-panel figure
showing posterior predictive distributions for velocity (in three dimensions), projected positions vs. time in the plane of the sky, 
and various other two and three-dimensional views.
```@example 1
octoplot(model, chain)
```


Finally we display the posterior as a corner plot:
```@example 1
# Create a corner plot / pair plot.
# We can access any property from the chain specified in Variables
using CairoMakie: Makie
using PairPlots
octocorner(model, chain, small=true)
```



Notice that our data is consitent with orbits with three completely different orbital periods. It might be worth separating the posterior and quoting uncertainties on the three modes separately:
```@example 1
mode1 = chain[       chain["B_a"][:] .< 5.5]
mode2 = chain[5.5 .< chain["B_a"][:] .< 7.5]
mode3 = chain[7.5 .< chain["B_a"][:]       ]
octocorner(model,mode1, mode2, mode3, small=true)
```

Here we plot just the semi-major axis and eccentricity:
```@example 1
viz_layers = (
    PairPlots.Scatter(),
    PairPlots.MarginDensity(),
    PairPlots.MarginConfidenceLimits()
)
pairplot(
    (;a=mode1["B_a"][:], e=mode1["B_e"][:])=>viz_layers,
    (;a=mode2["B_a"][:], e=mode2["B_e"][:])=>viz_layers,
    (;a=mode3["B_a"][:], e=mode3["B_e"][:])=>viz_layers,
)
```

