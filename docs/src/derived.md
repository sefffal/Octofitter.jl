#  [Derived Variables] (@id derived)

Octofitter has a concept called "derived variables" that are inspired by PyMC3.
Derived variables are quantities that either have a fixed value, or a fixed mathematical relationship with the main variables in a model.

This concept is extremely powerful, as it lets you quickly create very sophisticated models.

Derived variables allow you to mark certain properties as constants, reparameterize models, link properties between planets in multi-planet systems, plug in physical models, and more.

## System Variables
Derived variables for the system as a whole can be created as follows:

```julia
sys = System(
    name="HD12345",
    companions=[],
    observations=[],
    variables=@variables begin
        M = 1.0
        plx ~ Normal(45., 0.02)
    end
)
```
In this case, instead of including `M` as a variable in the model, we define it as a function that always returns `1.0`. This is equivalent to passing `M=1.0`.

In the following case, let's define `M` as being calculated based on another variable in the model. This is how you can do reparameterizations in Octofitter.jl
```julia
sys = System(
    name="HD12345",
    companions=[],
    observations=[],
    variables=@variables begin
        plx ~ Normal(45., 0.02)
        logM ~ Normal(45., 0.02)
        M = 10^logM
    end
)
```
We defined a new variable `logM` as a prior, and then calculate M from it.

In general, you can write any function you want to map from any of combination of constants and variables in the model to new variables. The only constraints are that your functions always return the same outputs for the same inputs, and are differentiable. These functions will be called in a tight loop, so you should try to make sure they are as efficient as possible.


## Planet Variables
Derived variables for an individual planet are similar, but have access to both the planet's variables and the system as a whole.

Here is an example of reparameterizing `e` and `a` on a planet to be logarithmic quantities:
```julia
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[],
    variables=@variables begin
        ω ~ Normal(0.1, deg2rad(30.))
        i ~ Normal(0.1, deg2rad(30.))
        Ω ~ Normal(0.0, deg2rad(30.))
        loge ~ Uniform(-4, 1)
        loga ~ Normal(1, 1)
        e = 10^loge
        a = 10^loga

        M = system.M
        τ ~ UniformCircular(1.0)
        P = √(a^3/M)
        tp = τ*P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
    end
)
```
Here `e` is defined as log-uniform, and `a` as log-normal.


## Linking Planets
Planets can have Derived variables that are calculated from variables defined on the system as a whole.
This makes it easy to, for example, create a system of two planets that are co-planar.

```julia
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[],
    variables=@variables begin
        a ~ Uniform(0, 15)
        e ~ Uniform(0,0.99)
        ω ~ Normal(0.1, deg2rad(30.))
        i = system.i
        Ω = system.Ω

        M = system.M
        τ ~ UniformCircular(1.0)
        P = √(a^3/M)
        tp = τ*P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
    end
)

planet_c = Planet(
    name="c",
    basis=Visual{KepOrbit},
    observations=[],
    variables=@variables begin
        a ~ Uniform(15, 45)
        e ~ Uniform(0,0.99)
        ω ~ Normal(0.1, deg2rad(30.))
        i = system.i
        Ω = system.Ω

        M = system.M
        τ ~ UniformCircular(1.0)
        P = √(a^3/M)
        tp = τ*P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
    end
)

sys = System(
    name="HD12345",
    companions=[planet_b, planet_c],
    observations=[],
    variables=@variables begin
        plx ~ Normal(45., 0.02)
        M ~ Normal(45., 0.02)
        i ~ Normal(0.1, deg2rad(30.))
        Ω ~ Normal(0.0, deg2rad(30.))
    end
)
```
Notice how `i` and `Ω` are defined as variables on the System. The two planets B & C instead just take their values from the System. This way we can enforce co-planarity between planets without e.g. rejection sampling.

## Resolution Order
The order that variables are resolved is as follows:
* Variables defined as priors for the system and planets
* Derived variables on the system
* Derived variables on each planet

You can use one derived variable from another based on their order in the `@variables` block within `System()` or `Planet()` constructors.
You cannot access variables from a different planet inside a `Planet` variables block. If you need to do this, move the variable up to the `System` variables block.

!!! warning "Common Pitfall: Accessing Variables"
    Inside a `Planet` variables block, you **cannot** use `planet.variable_name` to access the current planet's variables. Instead:

    - Use the variable name directly if it's defined earlier in the same `@variables` block
    - Use `system.variable_name` to access system-level variables

    ```julia
    # ❌ WRONG - cannot use planet.P_days inside the same planet block
    variables=@variables begin
        P_days ~ Uniform(1, 100)
        tp = planet.P_days * τ  # This will NOT work!
    end

    # ✅ CORRECT - use the variable name directly
    variables=@variables begin
        P_days ~ Uniform(1, 100)
        tp = P_days * τ  # Use P_days directly
    end

    # ✅ CORRECT - access system variables via system.
    variables=@variables begin
        P_inner = system.P_inner  # Access from system
        tp = P_inner * τ
    end
    ```

## Interpolation Syntax with `$`

The `$` symbol in `@variables` blocks allows you to inject external values and functions into derived variable expressions. This is useful for referencing constants, interpolation functions, or any value defined outside the model.

```julia
# Define constants and functions outside the model
HOST_L_MAG = 4.5
mass_to_contrast(m) = m^0.5  # Simple example; use real atmosphere models in practice

planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[...],
    variables=@variables begin
        mass ~ Uniform(0, 10)
        flux = $mass_to_contrast(mass)  # Interpolate the function
    end
)
```

!!! warning "Interpolation Limitations"
    The `$` interpolation only works for **simple references**. Complex expressions inside `$()` will cause a syntax error:

    ```julia
    # ❌ This will FAIL with "syntax: $ expression outside quote"
    flux = $mass_to_L_contrast(planet.mass, system.age, $HOST_L_MAG)
    ```

    **Workaround**: Create a wrapper function that captures any constants:

    ```julia
    # Define constants
    HOST_L_MAG = 4.5
    TRUE_AGE = 10.0  # Myr

    # Create wrapper function that bakes in constants
    function mass_to_L_contrast_wrapper(mass)
        return mass_to_L_contrast(mass, TRUE_AGE, HOST_L_MAG)
    end

    # ✅ Now use the simple wrapper
    variables=@variables begin
        mass ~ Uniform(0, 10)
        flux = $mass_to_L_contrast_wrapper(mass)
    end
    ```

    Alternatively, if you need to use model variables (not constants), include them as arguments without `$`:

    ```julia
    # ✅ Model variables don't need $ - they're resolved at runtime
    variables=@variables begin
        mass ~ Uniform(0, 10)
        flux = $mass_to_L_contrast(mass, system.age, planet.temp)
    end
    ```

## Variable Scoping in Observations

When you define derived variables in an observation's `@variables` block, the available variable namespaces depend on **where the observation is attached**:

### Observations Attached to a Planet

When an observation (like `PhotometryObs` or `PlanetRelAstromObs`) is attached to a `Planet`, you can access:
- `planet.X` - variables defined on the planet
- `system.X` - variables defined on the system

```julia
H_band_data = PhotometryObs(
    data_table,
    name="H_band",
    variables=@variables begin
        # planet.mass comes from the Planet's variables block
        # system.age comes from the System's variables block
        flux = $H_band_contrast_interp(planet.mass, system.age)
    end
)

planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[H_band_data],  # Observation attached to planet
    variables=@variables begin
        mass ~ Uniform(0, 10)  # Accessed as planet.mass in observations
        # ...
    end
)

sys = System(
    name="HD12345",
    companions=[planet_b],
    observations=[],
    variables=@variables begin
        age ~ Normal(15, 1)  # Accessed as system.age in observations
        # ...
    end
)
```

### Observations Attached to a System

When an observation (like `StarAbsoluteRVObs`) is attached directly to the `System`, you can access:
- `system.X` - variables defined on the system

```julia
rv_obs = StarAbsoluteRVObs(
    rv_table,
    name="HARPS",
    variables=@variables begin
        offset ~ Normal(0, 10)
        jitter ~ LogUniform(0.1, 10)
        trend_slope = system.trend_slope  # From system variables
    end
)

sys = System(
    name="HD12345",
    companions=[...],
    observations=[rv_obs],  # Observation attached to system
    variables=@variables begin
        trend_slope ~ Uniform(-1, 1)  # Accessed as system.trend_slope
        # ...
    end
)
```

