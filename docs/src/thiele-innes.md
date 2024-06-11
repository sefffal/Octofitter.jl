# Fit with a Thiele-Innes Basis

This example shows how to fit relative astrometry using a Thiele-Innes orbital basis instead of the traditional Campbell basis used in other tutorials. The Thiele-Innes basis is more suitable then Campbell for fitting low-eccentricity orbits, because it does not have the issues where `ω`, `Ω`, and `tp` become degenerate as eccentricity and/or inclination fall to zero.

At the end, we will convert our results back into the Campbell basis to compare.

```@example 1
using Octofitter
using CairoMakie
using PairPlots
using Distributions

astrom_like = PlanetRelAstromLikelihood(
    (epoch = 50000, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50120, ra = -502.570356287689, dec = -37.47217527025044, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50240, ra = -498.2089148883798, dec = -7.927548139010479, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50360, ra = -492.67768482682357, dec = 21.63557115669823, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50480, ra = -485.9770335870402, dec = 51.147204404903704, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50600, ra = -478.1095526888573, dec = 80.53589069730698, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50720, ra = -469.0801731788123, dec = 109.72870493064629, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 50840, ra = -458.89628893460525, dec = 138.65128697876773, σ_ra = 10, σ_dec = 10, cor=0),
)
@planet b ThieleInnesOrbit begin
    e ~ Uniform(0.0, 0.5)
    A ~ Normal(0, 10000) # milliarcseconds
    B ~ Normal(0, 10000) # milliarcseconds
    F ~ Normal(0, 10000) # milliarcseconds
    G ~ Normal(0, 10000) # milliarcseconds
    
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50000.0)  # reference epoch for θ. Choose an MJD date near your data.
end astrom_like

@system TutoriaPrime begin
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end b

model = Octofitter.LogDensityModel(TutoriaPrime)
```

We now sample from the model as usual:
```@example 1
using Random
Random.seed!(0)

results = octofit(model)
```
Notice that the fit was very very fast! The Thiele-Innes orbital paramterization is easier to explore than the default Campbell in 
many cases.

We now display the results:
```@example 1
octoplot(model,results)
```

```@example 1
octocorner(model, results, small=false)
```

## Conversion back to Campbell Elements
To convert our chain into the more familiar Campbell parameterization, we have to do a few steps. We start by turning the chain table into a an array of orbit objects, and then convert their type:

```@example 1
orbits_ti = Octofitter.construct_elements(results, :b, :) # colon means all rows
```

Here is one of those entries:
```@example 1
display(orbits_ti[1])
```

We can now make a table of results (and visualize them in a corner plot) by querying properties of these objects:
```@example 1
table = (;
    B_a = semimajoraxis.(orbits_ti),
    B_e = eccentricity.(orbits_ti),
    B_i = rad2deg.(inclination.(orbits_ti)),
)
pairplot(table)
```

We can also convert the orbit objects into Campbell parameters:
```@example 1
orbits_campbell = Visual{KepOrbit}.(orbits_ti)
orbits_campbell[1]
```
