# Hierarchical Co-Planar, Resonant Model

This example shows how you can fit a two planet model to relative astrometry data. This functionality would work equally well with RV, images, etc.

We will restrict the two planets to being exactly co-planar, and have a near 2:1 period resonance, and zero-eccentricity.

For this example, we will use astrometry from the HR8799 system collated by Jason Wang and retrieved from the website [Whereistheplanet.com](http://whereistheplanet.com).


```@example 1
using Octofitter
using CairoMakie
using PairPlots
using Distributions
using PlanetOrbits
```

```@example 1
astrom_b = PlanetRelAstromLikelihood(
    (epoch=53200.0, object=1, ra=1471.0, σ_ra=6.0, dec=887.0, σ_dec=6.0, cor=0, ),
    (epoch=54314.0, object=1, ra=1504.0, σ_ra=3.0, dec=837.0, σ_dec=3.0, cor=0, ),
    (epoch=54398.0, object=1, ra=1500.0, σ_ra=7.0, dec=836.0, σ_dec=7.0, cor=0, ),
    (epoch=54727.0, object=1, ra=1516.0, σ_ra=4.0, dec=818.0, σ_dec=4.0, cor=0, ),
    (epoch=55042.0, object=1, ra=1526.0, σ_ra=4.0, dec=797.0, σ_dec=4.0, cor=0, ),
    (epoch=55044.0, object=1, ra=1531.0, σ_ra=7.0, dec=794.0, σ_dec=7.0, cor=0, ),
    (epoch=55136.0, object=1, ra=1524.0, σ_ra=10.0, dec=795.0, σ_dec=10.0, cor=0, ),
    (epoch=55390.0, object=1, ra=1532.0, σ_ra=5.0, dec=783.0, σ_dec=5.0, cor=0, ),
    (epoch=55499.0, object=1, ra=1535.0, σ_ra=15.0, dec=766.0, σ_dec=15.0, cor=0, ),
    (epoch=55763.0, object=1, ra=1541.0, σ_ra=5.0, dec=762.0, σ_dec=5.0, cor=0, ),
    (epoch=56130.0, object=1, ra=1545.0, σ_ra=5.0, dec=747.0, σ_dec=5.0, cor=0, ),
    (epoch=56226.0, object=1, ra=1549.0, σ_ra=4.0, dec=743.0, σ_dec=4.0, cor=0, ),
    (epoch=56581.0, object=1, ra=1545.0, σ_ra=22.0, dec=724.0, σ_dec=22.0, cor=0, ),
    (epoch=56855.0, object=1, ra=1560.0, σ_ra=13.0, dec=725.0, σ_dec=13.0, cor=0, ),
    (epoch=58798.03906, object=1, ra=1611.002, σ_ra=0.133, dec=604.893, σ_dec=0.199, cor=-0.406, ),
    (epoch=59453.245, object=1, ra=1622.924, σ_ra=0.32, dec=570.534, σ_dec=0.296, cor=-0.905, ),
    (epoch=59454.231, object=1, ra=1622.872, σ_ra=0.204, dec=571.296, σ_dec=0.446, cor=-0.79, ),
)

astrom_c = PlanetRelAstromLikelihood(
    (epoch=53200.0, object=2, ra=-739.0, σ_ra=6.0, dec=612.0, σ_dec=6.0, cor=0, ),
    (epoch=54314.0, object=2, ra=-683.0, σ_ra=4.0, dec=671.0, σ_dec=4.0, cor=0, ),
    (epoch=54398.0, object=2, ra=-678.0, σ_ra=7.0, dec=678.0, σ_dec=7.0, cor=0, ),
    (epoch=54727.0, object=2, ra=-663.0, σ_ra=3.0, dec=693.0, σ_dec=3.0, cor=0, ),
    (epoch=55042.0, object=2, ra=-639.0, σ_ra=4.0, dec=712.0, σ_dec=4.0, cor=0, ),
    (epoch=55136.0, object=2, ra=-636.0, σ_ra=9.0, dec=720.0, σ_dec=9.0, cor=0, ),
    (epoch=55390.0, object=2, ra=-619.0, σ_ra=4.0, dec=728.0, σ_dec=4.0, cor=0, ),
    (epoch=55499.0, object=2, ra=-607.0, σ_ra=12.0, dec=744.0, σ_dec=12.0, cor=0, ),
    (epoch=55763.0, object=2, ra=-595.0, σ_ra=4.0, dec=747.0, σ_dec=4.0, cor=0, ),
    (epoch=56130.0, object=2, ra=-578.0, σ_ra=5.0, dec=761.0, σ_dec=5.0, cor=0, ),
    (epoch=56226.0, object=2, ra=-572.0, σ_ra=3.0, dec=768.0, σ_dec=3.0, cor=0, ),
    (epoch=56581.0, object=2, ra=-542.0, σ_ra=22.0, dec=784.0, σ_dec=22.0, cor=0, ),
    (epoch=56855.0, object=2, ra=-540.0, σ_ra=12.0, dec=799.0, σ_dec=12.0, cor=0, ),
)

fig = Makie.scatter(astrom_b.table.ra, astrom_b.table.dec, axis=(;autolimitaspect=1))
Makie.scatter!(astrom_c.table.ra, astrom_c.table.dec)
Makie.scatter!([0], [0], marker='⋆', markersize=50, color=:black)
fig
```


We now specify our two planet model for planets b & c.

```@example 1
@planet b Visual{KepOrbit} begin
    e = 0.0
    ω = 0.0

    # Use the system inclination and longitude of ascending node
    # variables
    i = system.i
    Ω = system.Ω

    # Specify the period as ~ 1% around 2X the P_nominal variable
    P_mul ~ Normal(1, 0.1)
    P = 2*system.P_nominal * b.P_mul

    a = cbrt(system.M * b.P^2)
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,59454.231)  # reference epoch for θ. Choose an MJD date near your data.
end astrom_b
@planet c Visual{KepOrbit} begin
    e = 0.0
    ω = 0.0

    # Use the system inclination and longitude of ascending node
    # variables
    i = system.i
    Ω = system.Ω

    # Specify the period as ~ 1% the P_nominal variable
    P_mul ~ truncated(Normal(1, 0.1), lower=0.1)
    P = system.P_nominal * c.P_mul

    a = cbrt(system.M * c.P^2)

    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,c,59454.231)  # reference epoch for θ. Choose an MJD date near your data.
end astrom_c
@system HR8799_res_co begin
    plx ~ gaia_plx(;gaia_id=2832463659640297472)
    M ~ truncated(Normal(1.5, 0.02), lower=0.1)
    # We create inclination and longitude of ascending node variables at the
    # system level.
    i ~ Sine()
    Ω ~ UniformCircular()
    # We create a nominal period of planet c variable. 
    P_nominal ~ Uniform(50, 300) # years
end b c


model = Octofitter.LogDensityModel(HR8799_res_co)
```


Let's now sample from the model:
```@example 1
using Pigeons
results,pt = octofit_pigeons(model, n_rounds=10);
nothing # hide
```

Plots the orbits:
```@example 1
octoplot(model, results)
```

Corner plot:
```@example 1
octocorner(model, results, small=true)
```

Now examine the period ratio:
```@example 1
hist(
    results[:b_P][:] ./ results[:c_P][:],
    axis=(;
        xlabel="period ratio",
        ylabel="counts",
    )
)
```
