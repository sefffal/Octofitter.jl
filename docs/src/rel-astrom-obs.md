# Observable-Based Priors

This tutorial shows how to fit an orbit to relative astrometry using the observable-based priors of[O'Neil et al. 2019](https://ui.adsabs.harvard.edu/abs/2019AJ....158....4O). Please cite that paper if you use this functionality.

We will fit the same astrometry as in the [previous tutorial](@ref fit-astrometry), and just change our priors.


```@example 1
using Octofitter
using CairoMakie
using PairPlots
using Distributions

astrom_dat = Table(;
    epoch = [50000, 50120, 50240, 50360, 50480, 50600, 50720, 50840],
    ra    = [-494.4, -495.0, -493.7, -490.4, -485.2, -478.1, -469.1, -458.3],
    dec   = [-76.7, -44.9, -12.9, 19.1, 51.0, 82.8, 114.3, 145.3],
    σ_ra  = [12.6, 10.4, 9.9, 8.7, 8.0, 6.9, 5.8, 4.2],
    σ_dec = [12.6, 10.4, 9.9, 8.7, 8.0, 6.9, 5.8, 4.2],
    cor   = [0.2, 0.5, 0.1, -0.8, 0.3, -0.0, 0.1, -0.2]
)

astrom_obs = PlanetRelAstromObs(
    astrom_dat,
    name = "obs_prior_example",
    variables = @variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)
# We wrap the likelihood in this prior
obs_pri_astrom_obs = ObsPriorAstromONeil2019(astrom_obs)

planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    # NOTE! We only provide the wrapped obs_pri_astrom_obs
    likelihoods=[obs_pri_astrom_obs],
    variables=@variables begin
        M = system.M
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        # Results will be sensitive to the prior on period
        P ~  LogUniform(0.1, 150) # Period, years
        a = ∛(M * P^2)
        θ_x ~ Normal()
        θ_y ~ Normal()
        θ = atan(θ_y, θ_x)
        tp = θ_at_epoch_to_tperi(θ, 50420; M, e, a, i, ω, Ω)
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

Initialize the model starting points and confirm the data are entered correctly:
```@example 1
init_chain = initialize!(model,)
octoplot(model, init_chain)
```

Now run the fit:
```@example 1
results_obspri = octofit(model,iterations=5000,)
```

Plot the MCMC results:
```@example 1
octoplot(model, results_obspri)
```

Compare this with the previous fit using uniform priors:
```@example 1
astrom_obs_uniform = PlanetRelAstromObs( # hide
    astrom_dat, # hide
    name = "uniform_prior_example", # hide
    variables = @variables begin # hide
        jitter = 0        # mas # hide
        northangle = 0    # radians # hide
        platescale = 1    # relative # hide
    end # hide
) # hide
planet_b_uniform = Planet( # hide
    name="b", # hide
    basis=Visual{KepOrbit}, # hide
    likelihoods=[astrom_obs_uniform], # hide
    variables=@variables begin # hide
        M = system.M # hide
        a ~ Uniform(0, 100) # hide
        e ~ Uniform(0.0, 0.5) # hide
        i ~ Sine() # hide
        ω_x ~ Normal() # hide
        ω_y ~ Normal() # hide
        ω = atan(ω_y, ω_x) # hide
        Ω ~ Uniform(0,2pi) # hide
        P = √(a^3/M) # hide
        θ ~ Uniform(0,2pi) # hide
        tp = θ_at_epoch_to_tperi(θ, 50420; M, e, a, i, ω, Ω) # hide
    end # hide
) # hide
sys_uniform = System( # hide
    name="Tutoria", # hide
    companions=[planet_b_uniform], # hide
    likelihoods=[], # hide
    variables=@variables begin # hide
        M ~ truncated(Normal(1.2, 0.1), lower=0.1) # hide
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1) # hide
    end # hide
) # hide
model_with_uniform_priors = Octofitter.LogDensityModel(sys_uniform) # hide
results_unif_pri = octofit(model_with_uniform_priors,iterations=5000,verbosity=0) # hide
octoplot(model_with_uniform_priors, results_unif_pri)
```

We can compare the results in a corner plot:
```@example 1
octocorner(model,results_unif_pri,results_obspri,small=true)
```
