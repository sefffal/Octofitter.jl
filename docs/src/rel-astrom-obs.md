# Observable-Based Priors

This tutorial shows how to fit an orbit to relative astrometry using the observable-based priors of[O'Neil et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....158....4O). Please cite that paper if you use this functionality.

We will fit the same astrometry as in the [previous tutorial](@ref fit-astrometry), and just change our priors.


```@example 1
using Octofitter
using CairoMakie
using PairPlots
using Distributions

astrom_like = PlanetRelAstromLikelihood(
    (epoch = 50000, ra = -494.4, dec = -76.7, σ_ra =  12.6, σ_dec =  12.6, cor=  0.2),
    (epoch = 50120, ra = -495.0, dec = -44.9, σ_ra =  10.4, σ_dec =  10.4, cor=  0.5),
    (epoch = 50240, ra = -493.7, dec = -12.9, σ_ra =   9.9, σ_dec =   9.9, cor=  0.1),
    (epoch = 50360, ra = -490.4, dec =  19.1, σ_ra =   8.7, σ_dec =   8.7, cor= -0.8),
    (epoch = 50480, ra = -485.2, dec =  51.0, σ_ra =   8.0, σ_dec =   8.0, cor=  0.3),
    (epoch = 50600, ra = -478.1, dec =  82.8, σ_ra =   6.9, σ_dec =   6.9, cor= -0.0),
    (epoch = 50720, ra = -469.1, dec = 114.3, σ_ra =   5.8, σ_dec =   5.8, cor=  0.1),
    (epoch = 50840, ra = -458.3, dec = 145.3, σ_ra =   4.2, σ_dec =   4.2, cor= -0.2),
)
obs_pri = ObsPriorAstromONeil2019(astrom_like)
@planet b Visual{KepOrbit} begin
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    # Results will be sensitive to the prior on period
    P ~  LogUniform(0.1, 150) # Period, years
    a = ∛(system.M * b.P^2)
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50420)
end astrom_like obs_pri

@system TutoriaPrime begin
    M ~ truncated(Normal(1.2, 0.1), lower=0.1)
    plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
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
octoplot(model, results_obspri)
```

Compare this with the previous fit using uniform priors:
```@example 1
@planet b Visual{KepOrbit} begin # hide
    a ~ Uniform(0, 100) # hide
    e ~ Uniform(0.0, 0.5) # hide
    i ~ Sine() # hide
    ω ~ UniformCircular() # hide
    # Ω ~ UniformCircular() # hide
    # P = √(b.a^3/system.M) # hide
    # θ ~ UniformCircular() # hide
    # ω ~ Uniform(0,2pi) # hide
    Ω ~ Uniform(0,2pi) # hide
    P = √(b.a^3/system.M) # hide
    θ ~ Uniform(0,2pi) # hide
    tp = θ_at_epoch_to_tperi(system,b,50420) # hide
end astrom_like # hide
@system Tutoria begin # hide
    M ~ truncated(Normal(1.2, 0.1), lower=0.1) # hide
    plx ~ truncated(Normal(50.0, 0.02), lower=0.1) # hide
end b # hide
model_with_uniform_priors = Octofitter.LogDensityModel(Tutoria) # hide
Random.seed!(0) # hide
results_unif_pri = octofit(model_with_uniform_priors,iterations=5000,verbosity=0) # hide
octoplot(model_with_uniform_priors, results_unif_pri)
```

We can compare the results in a corner plot:
```@example 1
octocorner(model,results_unif_pri,results_obspri,small=true)
```
