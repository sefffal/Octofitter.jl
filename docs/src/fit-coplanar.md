# Hierarchical Co-Planar, Near-Resonant Model

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
astrom_dat_b = Table(;
    epoch = [53200.0, 54314.0, 54398.0, 54727.0, 55042.0, 55044.0, 55136.0, 55390.0, 55499.0, 55763.0, 56130.0, 56226.0, 56581.0, 56855.0, 58798.03906, 59453.245, 59454.231],
    ra    = [1471.0, 1504.0, 1500.0, 1516.0, 1526.0, 1531.0, 1524.0, 1532.0, 1535.0, 1541.0, 1545.0, 1549.0, 1545.0, 1560.0, 1611.002, 1622.924, 1622.872],
    dec   = [887.0, 837.0, 836.0, 818.0, 797.0, 794.0, 795.0, 783.0, 766.0, 762.0, 747.0, 743.0, 724.0, 725.0, 604.893, 570.534, 571.296],
    σ_ra  = [6.0, 3.0, 7.0, 4.0, 4.0, 7.0, 10.0, 5.0, 15.0, 5.0, 5.0, 4.0, 22.0, 13.0, 0.133, 0.32, 0.204],
    σ_dec = [6.0, 3.0, 7.0, 4.0, 4.0, 7.0, 10.0, 5.0, 15.0, 5.0, 5.0, 4.0, 22.0, 13.0, 0.199, 0.296, 0.446],
    cor   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.406, -0.905, -0.79]
)

astrom_dat_c = Table(;
    epoch = [53200.0, 54314.0, 54398.0, 54727.0, 55042.0, 55136.0, 55390.0, 55499.0, 55763.0, 56130.0, 56226.0, 56581.0, 56855.0],
    ra    = [-739.0, -683.0, -678.0, -663.0, -639.0, -636.0, -619.0, -607.0, -595.0, -578.0, -572.0, -542.0, -540.0],
    dec   = [612.0, 671.0, 678.0, 693.0, 712.0, 720.0, 728.0, 744.0, 747.0, 761.0, 768.0, 784.0, 799.0],
    σ_ra  = [6.0, 4.0, 7.0, 3.0, 4.0, 9.0, 4.0, 12.0, 4.0, 5.0, 3.0, 22.0, 12.0],
    σ_dec = [6.0, 4.0, 7.0, 3.0, 4.0, 9.0, 4.0, 12.0, 4.0, 5.0, 3.0, 22.0, 12.0],
    cor   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)

astrom_b = PlanetRelAstromLikelihood(
    astrom_dat_b,
    instrument_name = "GPI",
    variables = @variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)

astrom_c = PlanetRelAstromLikelihood(
    astrom_dat_c,
    instrument_name = "GPI",
    variables = @variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)

fig = Makie.scatter(astrom_b.table.ra, astrom_b.table.dec, axis=(;autolimitaspect=1))
Makie.scatter!(astrom_c.table.ra, astrom_c.table.dec)
Makie.scatter!([0], [0], marker='⋆', markersize=50, color=:black)
fig
```


We now specify our two planet model for planets b & c.

```@example 1
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[astrom_b],
    variables=@variables begin
        e = 0.0
        ω = 0.0
        M_pri = super.M_pri
        M_b = super.M_b
        M_c = super.M_c
        M = M_pri + M_b*Octofitter.mjup2msol + M_c*Octofitter.mjup2msol
        mass = M_b

        # Use the system inclination and longitude of ascending node
        # variables
        i = super.i
        Ω = super.Ω

        # Specify the period as ~ 1% around 2X the P_nominal variable
        P_mul ~ Normal(1, 0.1)
        P_nominal = super.P_nominal
        P = 2*P_nominal * P_mul

        a = cbrt(M * P^2)
        θ_x ~ Normal()
        θ_y ~ Normal()
        θ = atan(θ_y, θ_x)
        tp = θ_at_epoch_to_tperi(θ, 59454.231; M, e, a, i, ω, Ω)  # reference epoch for θ. Choose an MJD date near your data.
    end
)

planet_c = Planet(
    name="c",
    basis=Visual{KepOrbit},
    likelihoods=[astrom_c],
    variables=@variables begin
        e = 0.0
        ω = 0.0
        M_pri = super.M_pri
        M_b = super.M_b
        M_c = super.M_c
        M = M_pri + M_b*Octofitter.mjup2msol
        mass = M_c

        # Use the system inclination and longitude of ascending node
        # variables
        i = super.i
        Ω = super.Ω

        # Specify the period as ~ 1% the P_nominal variable
        P_mul ~ truncated(Normal(1, 0.1), lower=0.1)
        P_nominal = super.P_nominal
        P = P_nominal * P_mul

        a = cbrt(M * P^2)

        θ_x ~ Normal()
        θ_y ~ Normal()
        θ = atan(θ_y, θ_x)
        tp = θ_at_epoch_to_tperi(θ, 59454.231; M, e, a, i, ω, Ω)  # reference epoch for θ. Choose an MJD date near your data.
    end
)

sys = System(
    name="HR8799_res_co",
    companions=[planet_b, planet_c],
    likelihoods=[],
    variables=@variables begin
        plx ~ gaia_plx(;gaia_id=2832463659640297472)
        M_pri ~ truncated(Normal(1.5, 0.02), lower=0.1)
        M_b ~ Uniform(0, 12)
        M_c ~ Uniform(0, 12)
        # We create inclination and longitude of ascending node variables at the
        # system level.
        i ~ Sine()
        Ω_x ~ Normal()
        Ω_y ~ Normal()
        Ω = atan(Ω_y, Ω_x)
        # We create a nominal period of planet c variable. 
        P_nominal ~ Uniform(50, 300) # years
    end
)

model = Octofitter.LogDensityModel(sys)
```


Initialize the starting points, and confirm the data are entered correcly:
```@example 1
init_chain = initialize!(model, (;
    plx =24.4549,
  M_pri = 1.48,
    M_b = 5.73,
    M_c = 5.14,
    P_nominal = 230,
))
octoplot(model, init_chain)
```


Now sample from the model using Pigeons parallel tempering:
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
