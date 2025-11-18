# Marginalizing over Circular or Eccentric Hypotheses

This tutorial demonstrates how to simultaneously model both circular and eccentric orbit possibilities in a single fit, allowing you to rigorously evaluate whether eccentricity is necessary to explain your data.

## Motivation

When fitting orbits, a common question arises: **Is the orbit circular, or is eccentricity required by the data?**

A naive approach might be to simply fit with a prior like `e ~ Uniform(0, 0.99)` and check if the posterior excludes zero. However, this approach has a subtle problem: it doesn't properly account for the different model complexities. An eccentric orbit has more free parameters (eccentricity `e` and argument of periastron `ω`) compared to a circular orbit where these parameters are undefined or fixed.

A more rigorous approach is to **marginalize over both hypotheses simultaneously** using a discrete indicator variable. This allows us to:

1. Let the data determine which hypothesis is preferred
2. Calculate a Bayes factor comparing the two models
3. Avoid artificially penalizing the circular model by forcing it to explore unnecessary parameter space

## The Spike-and-Slab Prior

The technique we'll use is called a **spike-and-slab prior**. The idea is simple but powerful:

- The **"spike"** represents the circular orbit hypothesis (e = 0)
- The **"slab"** represents the eccentric orbit hypothesis (e > 0)
- We use a discrete indicator variable to switch between them

Mathematically, we can write:
```
eccentric ~ Bernoulli(0.5)        # 50% prior probability for each hypothesis
e′ ~ Uniform(0.0, 0.99)           # If eccentric, what is e?
e = eccentric × e′                 # Spike at zero when eccentric=0, slab when eccentric=1
```

When `eccentric = 0`, we get `e = 0` (the spike). When `eccentric = 1`, we get `e = e′` (the slab, uniformly distributed).

This prior structure naturally encodes both hypotheses and allows the posterior probability `P(eccentric = 1 | data)` to tell us which model the data prefer.

## Why This Differs from Just Using a Uniform Prior

You might wonder: why not just use `e ~ Uniform(0, 0.99)` and see if the posterior includes zero?

The key difference is in **model comparison**. With a simple uniform prior:
- You're always fitting an eccentric model, even when e ≈ 0
- The parameter `ω` is being fit even when it's physically meaningless (for circular orbits)
- You can't as easily calculate a Bayes factor between circular and eccentric hypotheses

With the spike-and-slab approach:
- When `eccentric = 0`, the model *actually becomes circular* with fewer degrees of freedom
- The posterior mean of the `eccentric` parameter directly gives you the probability of the eccentric hypothesis
- You can calculate a proper Bayes factor to quantify the evidence

## Setting Up the Model

Let's start by loading the necessary packages:

```@example 1
using Octofitter, Distributions
using CairoMakie, PairPlots
```

We'll use the same synthetic astrometry data from the basic tutorial:

```@example 1
astrom_dat = Table(;
    epoch = [50000, 50120, 50240, 50360, 50480, 50600, 50720, 50840,], # MJD (days)
    ra    = [-505.764, -502.57, -498.209, -492.678, -485.977, -478.11, -469.08, -458.896,], # mas
    dec   = [-66.9298, -37.4722, -7.92755, 21.6356, 51.1472, 80.5359, 109.729, 138.651,], # mas
    σ_ra  = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,], # mas
    σ_dec = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,], # mas
    cor   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,]
)
astrom_obs = PlanetRelAstromObs(astrom_dat, name="relastrom")
```

Now we define our planet model with the spike-and-slab prior on eccentricity:

```@example 1
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[astrom_obs],
    variables=@variables begin
        M = system.M
        plx = system.plx
        a ~ Uniform(0, 100)

        # Spike-and-slab prior for eccentricity
        eccentric ~ Bernoulli(0.5)      # 50% prior odds for eccentric vs circular
        e′ ~ Uniform(0.0, 0.99)         # If eccentric, what eccentricity?
        e = eccentric * e′              # Spike at 0 when eccentric=0, slab when eccentric=1


        # ω is also multiplied by eccentric because it's undefined for circular orbits
        ω′ ~ UniformCircular()
        ω = eccentric * ω′              # When e=0, ω has no physical meaning

        i ~ Sine()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ, 50420; M, e, a, i, ω, Ω)
    end
)
nothing # hide
```

!!! note "Why multiply ω by the indicator variable?"
    When an orbit is circular (`e = 0`), the argument of periastron `ω` becomes undefined—there is no periastron! By multiplying `ω` by the `eccentric` indicator, we ensure that when the model is circular, `ω` is automatically set to zero and doesn't waste computational effort exploring meaningless values. 

Now we complete the system definition:

```@example 1
sys = System(
    name="CircularOrEccentric",
    companions=[planet_b],
    likelihoods=[],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)

model = Octofitter.LogDensityModel(sys)
```

## Sampling with Pigeons

!!! warning "Important: Use Pigeons for Discrete Variables"
    The default HMC sampler (`octofit`) is **not compatible** with discrete variables like our `eccentric` indicator. You **must** use the Pigeons sampler via `octofit_pigeons` for models with discrete parameters.

    Make sure you have Pigeons installed:
    ```julia
    using Pkg
    Pkg.add("Pigeons")
    ```

Let's sample from our model using Pigeons:

```@example 1
using Pigeons
chain, pt = octofit_pigeons(model, n_rounds=10)
```

The Pigeons sampler with the default SliceSampler can explore both discrete and continuous parameter spaces. It automatically handles the switching between circular and eccentric hypotheses.

## Interpreting the Results

Let's look at the posterior distribution of our indicator variable:

```@example 1
using Statistics
mean_eccentric = mean(chain[:b_eccentric][:])
println("Posterior probability of eccentric orbit: ", round(mean_eccentric, digits=3))
```

The mean of the `eccentric` indicator variable directly gives us the posterior probability that the orbit is eccentric (given the data and our priors).

## Calculating the Bayes Factor

One of the most powerful aspects of this approach is that we can calculate a **Bayes factor** comparing the two hypotheses. The Bayes factor is defined as:

```math
\text{BF} = \frac{P(\text{data} | \text{eccentric model})}{P(\text{data} | \text{circular model})}
```

It turns out that for a spike-and-slab model with equal prior odds (Bernoulli(0.5)), the Bayes factor can be calculated directly from the posterior mean of the indicator variable:

```math
\text{BF} = \frac{\bar{I}}{1 - \bar{I}}
```

where $\bar{I}$ is the posterior mean of the `eccentric` indicator.

!!! note "Mathematical Intuition"
    Why does this work? By Bayes' theorem, the posterior odds equal the prior odds times the Bayes factor:

    ```math
    \frac{P(\text{eccentric} | \text{data})}{P(\text{circular} | \text{data})} = \frac{P(\text{eccentric})}{P(\text{circular})} \times \text{BF}
    ```

    Since we used `Bernoulli(0.5)`, our prior odds are 1:1. Therefore, the posterior odds *are* the Bayes factor. The posterior mean $\bar{I}$ gives us $P(\text{eccentric} | \text{data})$, so $1 - \bar{I}$ gives us $P(\text{circular} | \text{data})$, and their ratio is the Bayes factor.

Let's calculate it:

```@example 1
bayes_factor = mean_eccentric / (1 - mean_eccentric)
println("Bayes factor (eccentric vs circular): ", round(bayes_factor, digits=2))
println("Bayes factor (circular vs eccentric): ", round(1/bayes_factor, digits=2))
```

**Interpreting Bayes Factors:**
- BF > 10: Strong evidence for eccentric orbit
- BF = 3-10: Moderate evidence for eccentric orbit
- BF = 1/3-3: Data are ambiguous
- BF < 1/10: Strong evidence for circular orbit

In this case, you can interpret the Bayes factor to determine whether eccentricity is justified by your data.

## Visualizing Results

### All Samples

Let's first plot all posterior samples together:

```@example 1
octoplot(model, chain)
```

This plot includes both circular (when `eccentric = 0`) and eccentric (when `eccentric = 1`) samples, which is why you might see orbits with a range of eccentricities.

### Eccentric Samples Only

We can subset the chain to show only the eccentric orbit samples:

```@example 1
chain_eccentric = chain[chain[:b_eccentric][:] .> 0]
println("Number of eccentric samples: ", length(chain_eccentric))
octoplot(model, chain_eccentric)
```

### Circular Samples Only

Similarly, we can plot only the circular orbit samples:

```@example 1
chain_circular = chain[chain[:b_eccentric][:] .== 0]
println("Number of circular samples: ", length(chain_circular))
octoplot(model, chain_circular)
```

These subsetted plots help you visualize what each hypothesis predicts for your data.

!!! tip "Understanding the Split"
    The ratio of eccentric to circular samples in your chain reflects the posterior probability of each hypothesis. If you have roughly equal numbers, the data don't strongly prefer one model over the other. If one dominates, that tells you the data have a clear preference.

## Corner Plots

Let's examine the parameter correlations with a corner plot:

```@example 1
octocorner(model, chain, small=true)
```

Notice how the eccentricity parameter `b_e` has a spike at zero (circular orbits) and a continuous distribution above zero (eccentric orbits). This is the "spike and slab" structure!

You can also compare the eccentric and circular subsets:

```@example 1
octocorner(model, chain_eccentric, chain_circular, small=true)
```

## Generallization

5. **Generalization**: This technique can be extended to other discrete model choices, such as:
   - Coplanar vs non-coplanar multi-planet systems
   - Including vs excluding a Gaussian process for stellar activity
   - Different numbers of planets
