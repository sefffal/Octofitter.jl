# Octofitter.jl v7 Migration Guide

This guide helps you migrate your code from Octofitter.jl v6 to v7, which introduced a significant API redesign.

## Overview of Changes

The v7 API redesign eliminates the `@planet` and `@system` macros in favor of explicit `Planet()` and `System()` constructors. This change provides better error handling, clearer variable scoping, and more flexible model composition.

## Key Migration Steps

### 1. Model Definition Syntax

#### Old Syntax (v6):
```julia
@planet b Visual{KepOrbit} begin
    a ~ Uniform(0, 100)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system, b, 50000)
end astrom

@system HD1234 begin
    M ~ truncated(Normal(1.2, 0.1), lower=0.1)
    plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
end b
```

#### New Syntax (v7):
```julia
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[astrom],
    variables=@variables begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)  # Total mass for this orbit
        a ~ Uniform(0, 100)
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(θ, 50000; M, e, a, i, ω, Ω)
    end
)

sys = System(
    name="HD1234",
    companions=[planet_b],
    likelihoods=[],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)
```

### 2. Likelihood Construction

#### Old Syntax (v6):
```julia
astrom = PlanetRelAstromLikelihood(Table(
    epoch = [50000, 50120, 50240],
    ra = [-505.7, -502.5, -498.2],
    dec = [-66.9, -37.4, -7.9],
    σ_ra = [10.0, 10.0, 10.0],
    σ_dec = [10.0, 10.0, 10.0],
    cor = [0.0, 0.0, 0.0]
))
```

#### New Syntax (v7):

Your likelihood objects must be given a name in most cases:
```julia
astrom_dat = Table(
    epoch = [50000, 50120, 50240],
    ra = [-505.7, -502.5, -498.2],
    dec = [-66.9, -37.4, -7.9],
    σ_ra = [10.0, 10.0, 10.0],
    σ_dec = [10.0, 10.0, 10.0],
    cor = [0.0, 0.0, 0.0]
)
astrom = PlanetRelAstromLikelihood(astrom_dat, name="GPI astrom")
```

### 3. Radial Velocity Models

#### Old Syntax (v6):
```julia
rvlike_hires = MarginalizedStarAbsoluteRVLikelihood(
    hires_data,
    instrument_name="HIRES",
    jitter=:jitter_hires,
)
```

#### New Syntax (v7):
Instrument-specific variables are handled directly in the likelihood definition.
```julia
rvlike_hires = MarginalizedStarAbsoluteRVLikelihood(
    hires_data,
    name="HIRES",
    variables=@variables begin
        jitter ~ LogUniform(0.1, 100) # m/s
    end
)
```

### 4. Variable Access in Derived Parameters

#### Old Syntax (v6):
```julia
P = √(b.a^3/system.M)

tp = θ_at_epoch_to_tperi(system, b, 50000)
```

#### New Syntax (v7):

You no longer have to prefix with the planet name or `system`. Just use variabels directly.

The `θ_at_epoch_to_tperi` function syntax has changed. You now provided the necessary parameters for the calculation directly.

```julia
P = √(a^3/M)

# In planet variables block:
tp = θ_at_epoch_to_tperi(θ, 50000; M, e, a, i, ω, Ω)

# Access system variables with system. prefix when needed:
M = system.M  # if planet needs system mass
```


## Key Differences Summary

1. **Model Construction**: Replace `@planet` and `@system` macros with explicit `Planet()` and `System()` constructors
2. **Likelihood Names**: Most likelihoods now require a `name` parameter
3. **Variable Scoping**: Use direct variable names in derived expressions instead of qualified access
4. **Observation Variables**: Likelihood-specific variables are now defined in the likelihood's `@variables` block

