# [Multi-Planet RV Fits](@id fit-rv-multi)


This tutorial shows how to perform a multi-planet RV fit, and compare the bayesian evidence between the two models.


```@example 1
using Octofitter
using OctofitterRadialVelocity
using CairoMakie
using PairPlots
using Distributions
using PlanetOrbits
```

To begin, we create simulated data. We imagine that we have two different instruments.
```@example 1
using Random
Random.seed!(1)

orb_template_1 = orbit(a = 1.0,e = 0.05,ω = 1π/4,M = 1.0,tp =58800)
mass_1 = 0.25*1e-3
orb_template_2 = orbit(a = 5.0,e = 0.4,ω = 1π/4,M = 1.0,tp =59800)
mass_2 = 1.0*1e-3

epochs = (58400:150:69400) .+ 10 .* randn.()
rv = radvel.(orb_template_1, epochs, mass_1) .+ radvel.(orb_template_2, epochs, mass_2)
rvlike1 = MarginalizedStarAbsoluteRVLikelihood(
    Table(epoch=epochs, rv=rv .+ 4 .* randn.(), σ_rv=[4 .* abs.(randn.()) .+ 1 for _ in 1:length(epochs)]),
    name="DATA 1",
    variables=@variables begin
        jitter ~ LogUniform(0.1, 100) # m/s
    end
)

epochs = (65400:100:71400) .+ 10 .* randn.()
rv = radvel.(orb_template_1, epochs, mass_1) .+ radvel.(orb_template_2, epochs, mass_2)
rvlike2 = MarginalizedStarAbsoluteRVLikelihood(
    Table(epoch=epochs, rv=rv .+ 2 .* randn.() .+ 7, σ_rv=[2 .* abs.(randn.()) .+ 1 for _ in 1:length(epochs)]),
    name="DATA 2",
    variables=@variables begin
        jitter ~ LogUniform(0.1, 100) # m/s
    end
)

fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel="epoch [mjd]",
    ylabel="rv [m/s]"
)
Makie.scatter!(ax, rvlike1.table.epoch, rvlike1.table.rv)
Makie.errorbars!(ax, rvlike1.table.epoch, rvlike1.table.rv, rvlike1.table.σ_rv)
Makie.scatter!(ax, rvlike2.table.epoch, rvlike2.table.rv)
Makie.errorbars!(ax, rvlike2.table.epoch, rvlike2.table.rv, rvlike2.table.σ_rv)
fig
```

## Two Planet Model

```@example 1
planet_b = Planet(
    name="b",
    basis=RadialVelocityOrbit,
    likelihoods=[],
    variables=@variables begin
        M_pri = system.M_pri
        M_b = system.M_b
        M_c = system.M_c
        M = M_pri + (M_b + M_c) * Octofitter.mjup2msol
        e ~ Uniform(0,0.999999)
        mass = M_b
        ω ~ Uniform(0,2pi)
        τ ~ Uniform(0,1.0)

        P_kep_yrs ~ Uniform(0, 100)
        a = ∛(M * P_kep_yrs^2)
        tp = τ*P_kep_yrs*365.25 + 58400
    end
)

planet_c = Planet(
    name="c",
    basis=RadialVelocityOrbit,
    likelihoods=[],
    variables=@variables begin
        M_pri = system.M_pri
        M_c = system.M_c
        M = M_pri + M_c * Octofitter.mjup2msol
        e ~ Uniform(0,0.999999)
        mass = M_c
        ω ~ Uniform(0,2pi)
        τ ~ Uniform(0,1.0)

        P_kep_yrs ~ Uniform(0, 100)
        a = ∛(M * P_kep_yrs^2)
        tp = τ*P_kep_yrs*365.25 + 58400
    end
)

sim_2p = System(
    name="sim_2p",
    companions=[planet_b, planet_c],
    likelihoods=[rvlike1, rvlike2],
    variables=@variables begin
        M_pri = 1.0
        M_b ~ Uniform(0, 10)
        M_c ~ Uniform(0, 10)
    end
)

model_2p = Octofitter.LogDensityModel(sim_2p)
```

Sample from the posterior
```@example 1
using Pigeons
results_2p, pt_2p = octofit_pigeons(model_2p, n_rounds=10)
```

Plot RV curve, phase folded curve, and binned residuals:
```@example 1
Octofitter.rvpostplot(model_2p, results_2p)
```


## One Planet Model

We now create a new system object that only includes one planet (we dropped c, in this case).
```@example 1
sim_1p = System(
    name="sim_1p",
    companions=[planet_b],
    likelihoods=[rvlike1, rvlike2],
    variables=@variables begin
        M_pri = 1.0
        M_b ~ Uniform(0, 10)
        M_c = 0.0
    end
)

model_1p = Octofitter.LogDensityModel(sim_1p)
```

Sample from the posterior
```@example 1
using Pigeons
results_1p, pt_1p = octofit_pigeons(model_1p, n_rounds=10)
```

Plot RV curve, phase folded curve, and binned residuals:
```@example 1
Octofitter.rvpostplot(model_1p, results_1p)
```

## Model Comparison: Bayesian Evidence

Octofitter with Pigeons directly calculates the log Bayesian evidence using the "stepping stone" method. This should be more reliable than even nested sampling, and certainly more reliable than approximate methods like the BIC etc.

```@example 1
Z1 = stepping_stone(pt_1p)
Z2 = stepping_stone(pt_2p)

log_BF₁₀ = Z2-Z1
```

Here is a standard guideline you can use to interpret the evidence:

| Log Bayes Factor log(BF₁₀) | Interpretation                 |
|----------------------------|--------------------------------|
| > 4.61                     | Extreme evidence for $H_A$     |
| 3.40 - 4.61                | Very strong evidence for $H_A$ |
| 2.30 - 3.40                | Strong evidence for $H_A$      |
| 1.10 - 2.30                | Moderate evidence for $H_A$    |
| 0 - 1.10                   | Anecdotal evidence for $H_A$   |
| 0                          | No evidence                    |
| -1.10 - 0                  | Anecdotal evidence for $H_B$   |
| -2.30 - -1.10              | Moderate evidence for $H_B$    |
| -3.40 - -2.30              | Strong evidence for $H_B$      |
| -4.61 - -3.40              | Very strong evidence for $H_B$ |
| < -4.61                    | Extreme evidence for $H_B$     |

As you can see, the evidence for there being two planets is "extreme" in this case.
Try adjusting the masses of the two planets and see how this changes!

## Parameterizations

When using the evidence for model comparisons, a model with more specific priors will have more evidence than an quivalent model with broad priors.

In our two planet model above, we made two exactly equivalent planets. If you inspect the chains, you may notice that the two planets often flip back and forth -- sometimes `b` has the longer period, and sometimes `c` does. 

For example, here is a histogram of the period of planet b:
```@example 1
hist(vec(results_2p[:b_P_kep_yrs]), bins=100)
```

We can refine the two planet model a bit by adjusting the priors such that planet `c` always has a longer period than planet `b`.

This will make analysis a little more straightforward, but crucially it will also increase the evidence of this model, by approximately halving the prior volume---thus making a more specific prediction.

There are several ways we could do this. Here, we add a "nominal period" variable and reparameterize the two planets as ratios of this nominal period.


```@example 1
planet_b_v2 = Planet(
    name="b",
    basis=RadialVelocityOrbit,
    likelihoods=[],
    variables=@variables begin
        M_pri = system.M_pri
        M_b = system.M_b
        M_c = system.M_c
        M = M_pri + (M_b + M_c) * Octofitter.mjup2msol
        e ~ Uniform(0,0.999999)
        mass = M_b
        ω ~ Uniform(0,2pi)
        τ ~ Uniform(0,1.0)

        P_yrs_nom = system.P_yrs_nom
        P_ratio_b = system.P_ratio_b
        P_kep_yrs = P_yrs_nom * P_ratio_b
        a = ∛(M * P_kep_yrs^2)
        tp = τ*P_kep_yrs*365.25 + 58400
    end
)

planet_c_v2 = Planet(
    name="c",
    basis=RadialVelocityOrbit,
    likelihoods=[],
    variables=@variables begin
        M_pri = system.M_pri
        M_c = system.M_c
        M = M_pri + M_c * Octofitter.mjup2msol
        e ~ Uniform(0,0.999999)
        mass = M_c
        ω ~ Uniform(0,2pi)
        τ ~ Uniform(0,1.0)

        P_yrs_nom = system.P_yrs_nom
        P_ratio_c = system.P_ratio_c
        P_kep_yrs = P_yrs_nom * P_ratio_c
        a = ∛(M * P_kep_yrs^2)
        tp = τ*P_kep_yrs*365.25 + 58400
    end
)


sim_2p_v2 = System(
    name="sim_2p_v2",
    companions=[planet_b_v2, planet_c_v2],
    likelihoods=[rvlike1, rvlike2],
    variables=@variables begin
        M_pri = 1.0
        M_b ~ Uniform(0, 10)
        M_c ~ Uniform(0, 10)
        
        P_yrs_nom ~ Uniform(0, 100)
        P_ratio_b ~ Uniform(0, 0.5)
        P_ratio_c ~ Uniform(0.5, 1)
    end
)

model_2p_v2 = Octofitter.LogDensityModel(sim_2p_v2)
```

Sample from the posterior
```@example 1
using Pigeons
results_2p_v2, pt_2p_v2 = octofit_pigeons(model_2p_v2, n_rounds=10)
```

The planet with the wider orbit is now consistently plotted in the bottom panel (meaning that planet b and c are no longer trading back and forth):
```@example 1
Octofitter.rvpostplot(model_2p_v2, results_2p_v2)
```

If we look again at the log-evidence, we see that this parameterization (Z3) is even more favoured. This is because this small change in parameterization makes considerably 
```@example 1
Z1 = stepping_stone(pt_1p)
Z2 = stepping_stone(pt_2p)
Z3 = stepping_stone(pt_2p_v2)

Z1, Z2, Z3
```


As a final treat, let's animate the orbit plots. All the previous images were visualizing a single posterior draw. In this animation, we'll loop over many different samples:

```@example 1
Octofitter.rvpostplot_animated(model_2p_v2, results_2p_v2)
```

```@raw html
<video src="rv-posterior.mp4" autoplay loop width=300 height=300>
```



## Note about the evidence ratio
The pigeons method returns the log evidence ratio. If the priors are properly normalized, this is equal to the log evidence.

In other cases (e.g. if using `ObsPriorAstromONeil2019` or `UniformCircular`) you may need to calculate the log_Z0 term yourself. This can be done as follows:
```@example 1
prior_model = Octofitter.LogDensityModel(Octofitter.prior_only_model(model_1p.system, exclude_all=true))
_, pt_prior = octofit_pigeons(prior_model, n_rounds=10) # should be very quick!
log_Z0 = stepping_stone(pt_prior)
```

Subtract this from the stepping stone value to get the true evidence:
```@example 1
log_Z1_over_Z0 = stepping_stone(pt_1p)
log_Z1 = log_Z1_over_Z0 - log_Z0
```
