# Fit with a Thiele-Innes Basis

This example shows how to fit relative astrometry using a Thiele-Innes orbital basis instead of the traditional Campbell basis used in other tutorials. The Thiele-Innes basis is more suitable then Campbell for fitting low-eccentricity orbits, because it does not have the issues where `ω`, `Ω`, and `tp` become degenerate as eccentricity and/or inclination fall to zero.

At the end, we will convert our results back into the Campbell basis to compare.

```@example 1
using Octofitter
using CairoMakie
using PairPlots
using Distributions

astrom_dat = Table(;
    epoch = [50000, 50120, 50240, 50360, 50480, 50600, 50720, 50840],
    ra    = [-505.7637580573554, -502.570356287689, -498.2089148883798, -492.67768482682357, -485.9770335870402, -478.1095526888573, -469.0801731788123, -458.89628893460525],
    dec   = [-66.92982418533026, -37.47217527025044, -7.927548139010479, 21.63557115669823, 51.147204404903704, 80.53589069730698, 109.72870493064629, 138.65128697876773],
    σ_ra  = [10, 10, 10, 10, 10, 10, 10, 10],
    σ_dec = [10, 10, 10, 10, 10, 10, 10, 10],
    cor   = [0, 0, 0, 0, 0, 0, 0, 0]
)

astrom_like = PlanetRelAstromLikelihood(
    astrom_dat,
    name = "GPI",
    variables = @variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)

planet_b = Planet(
    name="b",
    basis=ThieleInnesOrbit,
    likelihoods=[astrom_like],
    variables=@variables begin
        e ~ Uniform(0.0, 0.5)
        A ~ Normal(0, 1000) # milliarcseconds
        B ~ Normal(0, 1000) # milliarcseconds
        F ~ Normal(0, 1000) # milliarcseconds
        G ~ Normal(0, 1000) # milliarcseconds
        
        M = super.M
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ, 50000.0; super.plx, M, e, A, B, F, G)  # reference epoch for θ. Choose an MJD date near your data.
    end
)

sys = System(
    name="TutoriaPrime",
    companions=[planet_b],
    likelihoods=[],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
```


Initialize the starting points, and confirm the data are entered correcly:
```@example 1
init_chain = initialize!(model)
octoplot(model, init_chain)
```

We now sample from the model as usual:
```@example 1
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
orbits_ti = Octofitter.construct_elements(model, results, :b, :) # colon means all rows
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
