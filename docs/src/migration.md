# [Migration Guide] (@id migration)

## Table of Contents
- [Migrating to v8](#migrating-to-v8)
- [Migrating to v7](#migrating-to-v7)

---

## Migrating to v8

Octofitter.jl v8 introduces a cleaner naming convention and several new features. This section covers the changes you need to make to update your code.

### Overview of Changes

1. **Observation type renaming**: All `*Likelihood` types renamed to `*Obs` (e.g., `PlanetRelAstromLikelihood` → `PlanetRelAstromObs`)
2. **Constructor keyword change**: `likelihoods=[...]` → `observations=[...]` in `Planet()` and `System()`
3. **New Gaia DR4 support**: `GaiaDR4AstromObs` for fitting DR4 individual astrometric data
4. **New stability priors**: `NonCrossingPrior()` and `HillStabilityPrior()` for multi-planet systems

### Global Find & Replace

Apply these replacements to update your code:

**Constructor keyword:**

| Old (v7) | New (v8) |
|----------|----------|
| `likelihoods=[...]` | `observations=[...]` |

**Core Octofitter types:**

| Old Name (v7) | New Name (v8) |
|---------------|---------------|
| `PhotometryLikelihood` | `PhotometryObs` |
| `PlanetRelAstromLikelihood` | `PlanetRelAstromObs` |
| `HGCALikelihood` | `HGCAObs` |
| `HGCAInstantaneousLikelihood` | `HGCAInstantaneousObs` |
| `HipparcosIADLikelihood` | `HipparcosIADObs` |
| `GaiaHipparcosUEVAJointLikelihood` | `GaiaHipparcosUEVAJointObs` |
| `GaiaDifferenceLikelihood` | `GaiaDifferenceObs` |
| `GaiaCatalogFitLikelihood` | `GaiaCatalogFitObs` |
| `GaiaUEVALikelihood` | `GaiaUEVAObs` |
| `LogLikelihoodMap` | `LogLikelihoodMapObs` |

**OctofitterRadialVelocity types:**

| Old Name (v7) | New Name (v8) |
|---------------|---------------|
| `StarAbsoluteRVLikelihood` | `StarAbsoluteRVObs` |
| `MarginalizedStarAbsoluteRVLikelihood` | `MarginalizedStarAbsoluteRVObs` |
| `StarAbsoluteRVMarginLikelihood` | `MarginalizedStarAbsoluteRVObs` |
| `PlanetRelativeRVLikelihood` | `PlanetRelativeRVObs` |

**OctofitterImages types:**
| Old Name (v7) | New Name (v8) |
|---------------|---------------|
| `ImageLikelihood` | `ImageObs` |

**OctofitterInterferometry types:**
| Old Name (v7) | New Name (v8) |
|---------------|---------------|
| `InterferometryLikelihood` | `InterferometryObs` |
| `AbstractInterferometryLikelihood` | `AbstractInterferometryObs` |
| `GRAVITYWideKPLikelihood` | `GRAVITYWideKPObs` |

**OctofitterTransits types:**
| Old Name (v7) | New Name (v8) |
|---------------|---------------|
| `LightCurveLikelihood` | `LightCurveObs` |

**Note:** Types ending in `Prior` (e.g., `NonCrossingPrior`, `HillStabilityPrior`) and `UserLikelihood` keep their original names.

### Example Migration

**Old code (v7):**
```julia
astrom = PlanetRelAstromLikelihood(astrom_data, name="GPI")

planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    likelihoods=[astrom],
    variables=@variables begin
        # ...
    end
)

sys = System(
    name="HD1234",
    companions=[planet_b],
    likelihoods=[hgca],
    variables=@variables begin
        # ...
    end
)
```

**New code (v8):**
```julia
astrom = PlanetRelAstromObs(astrom_data, name="GPI")

planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[astrom],
    variables=@variables begin
        # ...
    end
)

sys = System(
    name="HD1234",
    companions=[planet_b],
    observations=[hgca],
    variables=@variables begin
        # ...
    end
)
```

### New Features in v8

#### Gaia DR4 Support

Fit Gaia DR4 individual astrometric data with the new `GaiaDR4AstromObs` type:

```julia
gaia_obs = GaiaDR4AstromObs(
    observations_table,
    gaia_id=5064625130502952704,
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)
    end
)
```

New plotting functions `gaiatimeplot` and `gaiastarplot` are available for visualizing DR4 fits.

#### Stability Priors

For multi-planet systems, enforce dynamical stability constraints:

```julia
sys = System(
    name="MultiPlanet",
    companions=[planet_b, planet_c],
    observations=[NonCrossingPrior(), astrometry_obs],
    variables=@variables begin
        M ~ truncated(Normal(1.0, 0.1), lower=0.1)
        plx ~ truncated(Normal(50, 1), lower=0.1)
    end
)
```

Available priors:
- `NonCrossingPrior()` - Prevents orbit crossing
- `LimitClosestApproachAUPrior(soft_au)` - Soft/hard closest approach constraints
- `HillStabilityPrior()` - Enforces Hill stability criterion

---

## Migrating to v7

This guide helps you migrate your code from Octofitter.jl v6 to v7, which introduced a significant API redesign.

### Overview of Changes

The v7 API redesign eliminates the `@planet` and `@system` macros in favor of explicit `Planet()` and `System()` constructors. This change provides better error handling, clearer variable scoping, and more flexible model composition.

This upgrade is particularly useful for large batch processing systems, and for models with large numbers of instruments.

!!! note
    Pro-tip: paste your old Octofitter scripts and this migration guide into an LLM to help automate the conversion.

### Key Migration Steps

#### 1. Model Definition Syntax

##### Old Syntax (v6):
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

##### New Syntax (v7+):
```julia
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[astrom],
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
    observations=[],
    variables=@variables begin
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end
)
```

#### 2. Observation Construction

##### Old Syntax (v6):
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

##### New Syntax (v7+):

Your observation objects must be given a name in most cases:
```julia
astrom_dat = Table(
    epoch = [50000, 50120, 50240],
    ra = [-505.7, -502.5, -498.2],
    dec = [-66.9, -37.4, -7.9],
    σ_ra = [10.0, 10.0, 10.0],
    σ_dec = [10.0, 10.0, 10.0],
    cor = [0.0, 0.0, 0.0]
)
astrom = PlanetRelAstromObs(astrom_dat, name="GPI astrom")
```

**Alternative: Vector-of-NamedTuples (still supported):**
You can still use the vector-of-namedtuples syntax from v6, but now it must be wrapped in a Table first:
```julia
astrom = PlanetRelAstromObs(
    Table([
        (epoch=50000, ra=-505.7, dec=-66.9, σ_ra=10.0, σ_dec=10.0, cor=0.0),
        (epoch=50120, ra=-502.5, dec=-37.4, σ_ra=10.0, σ_dec=10.0, cor=0.0),
        (epoch=50240, ra=-498.2, dec=-7.9,  σ_ra=10.0, σ_dec=10.0, cor=0.0),
    ]),
    name="GPI astrom"
)
```

#### 3. Radial Velocity Models

##### Old Syntax (v6):
```julia
rvlike_hires = MarginalizedStarAbsoluteRVLikelihood(
    hires_data,
    instrument_name="HIRES",
    jitter=:jitter_hires,
)
```

##### New Syntax (v7+):
Instrument-specific variables are handled directly in the observation definition.
```julia
rvlike_hires = MarginalizedStarAbsoluteRVObs(
    hires_data,
    name="HIRES",
    variables=@variables begin
        jitter ~ LogUniform(0.1, 100) # m/s
    end
)
```

#### 4. Variable Access in Derived Parameters

You no longer have to prefix with the planet name or `system`. Just use variabels directly.

The `θ_at_epoch_to_tperi` function syntax has changed. You now provided the necessary parameters for the calculation directly.

##### Old Syntax (v6):
```julia
P = √(b.a^3/system.M)

tp = θ_at_epoch_to_tperi(system, b, 50000)
```

##### New Syntax (v7):


```julia
P = √(a^3/M)

# In planet variables block:
tp = θ_at_epoch_to_tperi(θ, 50000; M, e, a, i, ω, Ω)

# Access system variables with system. prefix when needed:
M = system.M  # if planet needs system mass
```


### Key Differences Summary

1. **Model Construction**: Replace `@planet` and `@system` macros with explicit `Planet()` and `System()` constructors
2. **Observation Names**: Most observation types now require a `name` parameter
3. **Variable Scoping**: Use direct variable names in derived expressions instead of qualified access
4. **Observation Variables**: Observation-specific variables are now defined in the observation's `@variables` block

