# Observale-Based Priors

This tutorial shows how to fit an orbit to relative astrometry using the observable-based priors of[O'Neil et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....158....4O). Please cite that paper if you use this functionality.

We will fit the same astrometry as in the [previous tutorial](@ref fit-astrometry), and just change our priors.


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
@planet b Visual{KepOrbit} begin
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    # Results will be sensitive to the prior on period
    P ~  LogUniform(0.1, 50)
    a = ∛(system.M * b.P^2)
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50420)
end astrom_like ObsPriorAstromONeil2019(astrom_like)

@system TutoriaPrime begin
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end b

model = Octofitter.LogDensityModel(TutoriaPrime)
```


Now run the fit:
```@example 1
using Random
Random.seed!(0)

results_obspri = octofit(model,iterations=5000,)
```


```@example 1
using Plots: Plots
plotchains(results_obspri, :b, kind=:astrometry, color="b_e")
Plots.plot!(astrom_like, label="astrometry")
Plots.xlims!(-1000,1000)
Plots.ylims!(-1000,1000)
```

Compare this with the previous fit using uniform priors:
```@example 1
@planet b Visual{KepOrbit} begin # hide
    a ~ Uniform(0, 100) # hide
    e ~ Uniform(0.0, 0.5) # hide
    i ~ Sine() # hide
    ω ~ UniformCircular() # hide
    Ω ~ UniformCircular() # hide
    P = √(b.a^3/system.M) # hide
    θ ~ UniformCircular() # hide
    tp = θ_at_epoch_to_tperi(system,b,50420) # hide
end astrom_like # hide
@system Tutoria begin # hide
    M ~ truncated(Normal(1.2, 0.1), lower=0) # hide
    plx ~ truncated(Normal(50.0, 0.02), lower=0) # hide
end b # hide
model = Octofitter.LogDensityModel(Tutoria) # hide
Random.seed!(0) # hide
results = octofit(model,iterations=5000,verbosity=0) # hide
plotchains(results, :b, kind=:astrometry, color="b_e")
Plots.plot!(astrom_like, label="astrometry")
Plots.xlims!(-1000,1000)
Plots.ylims!(-1000,1000)
```

We can compare the results in a corner plot:
```@example 1
octocorner(model,results,results_obspri,small=true)
```