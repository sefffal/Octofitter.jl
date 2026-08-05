#  [Derived Variables](@id derived)

Octofitter has a concept called "derived variables" that are inspired by PyMC3.
Derived variables are quantities that either have a fixed value, or a fixed mathematical relationship with the main variables in a model.

This concept is extremely powerful, as it lets you quickly create very sophisticated models.

Derived variables allow you to mark certain properties as constants, reparameterize models, link properties between bodies in multi-planet systems, plug in physical models, and more.

Inside a `@variables` block, `~` introduces a **prior** and `=` introduces a **derived
variable**.

## System Variables
Derived variables for the system as a whole can be created as follows:

```julia
sys = System(
    name="HD12345",
    bodies=[A],
    observations=[],
    variables=@variables begin
        age = 15.0
        plx ~ Normal(45., 0.02)
    end
)
```
In this case, instead of including `age` as a variable in the model, we define it as a constant that always returns `15.0`.

In the following case, let's define `plx` as being calculated from another variable in the model. This is how you can do reparameterizations in Octofitter.jl
```julia
sys = System(
    name="HD12345",
    bodies=[A],
    observations=[],
    variables=@variables begin
        logplx ~ Normal(1.65, 0.02)
        plx = 10^logplx
    end
)
```
We defined a new variable `logplx` as a prior, and then calculate `plx` from it.

In general, you can write any function you want to map from any of combination of constants and variables in the model to new variables. The only constraints are that your functions always return the same outputs for the same inputs, and are differentiable. These functions will be called in a tight loop, so you should try to make sure they are as efficient as possible.


## Body Variables
Derived variables for an individual body are similar, but have access to both the body's own variables and the system as a whole (through the `system.` prefix).

Here is an example of reparameterizing `e` and `a` on a planet to be logarithmic quantities:
```julia
planet_b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass = 0.0
        ω ~ Normal(0.1, deg2rad(30.))
        i ~ Normal(0.1, deg2rad(30.))
        Ω ~ Normal(0.0, deg2rad(30.))
        loge ~ Uniform(-4, 1)
        loga ~ Normal(1, 1)
        e = 10^loge
        a = 10^loga

        θ ~ UniformCircular()
        epoch = 58849.0   # reference epoch for θ. Choose an MJD date near your data.
    end
)
```
Here `e` is defined as log-uniform, and `a` as log-normal.

Three names in a body block are read specially: `mass` [M⊙], `flux` / `flux_<band>`, and
the orbital element keywords. Every other name is an ordinary local, free for later lines
in the same block to use — which is exactly what makes reparameterizations like the above
work.

## Linking Bodies
Bodies can have derived variables that are calculated from variables defined on the
system as a whole. This makes it easy to, for example, create a system of two planets
that are exactly co-planar:

```julia
planet_b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass = 0.0
        a ~ Uniform(0, 15)
        e ~ Uniform(0,0.99)
        ω ~ Normal(0.1, deg2rad(30.))
        i = system.i
        Ω = system.Ω

        θ ~ UniformCircular()
        epoch = 58849.0
    end
)

planet_c = Body(
    name="c",
    about=A,
    variables=@variables begin
        mass = 0.0
        a ~ Uniform(15, 45)
        e ~ Uniform(0,0.99)
        ω ~ Normal(0.1, deg2rad(30.))
        i = system.i
        Ω = system.Ω

        θ ~ UniformCircular()
        epoch = 58849.0
    end
)

sys = System(
    name="HD12345",
    bodies=[A, planet_b, planet_c],
    observations=[],
    variables=@variables begin
        plx ~ Normal(45., 0.02)
        i ~ Sine()
        Ω ~ Normal(0.0, deg2rad(30.))
    end
)
```
Notice how `i` and `Ω` are defined as variables on the System. The two planets b & c instead just take their values from the System. This way we can enforce co-planarity between planets without e.g. rejection sampling.

Hoisting a shared parameter to the system block like this is the *only* way to couple two
bodies exactly, and it is the intended one — a body block cannot see its siblings.

## Resolution Order

> Bodies look up and outward; the system block looks down and inward, after the fact;
> siblings never see each other.

1. System block, the lines that **do not** mention a body by name → visible everywhere as `system.*`
2. Every body's block, in declaration order
3. System block, the lines that **do** mention a body by name — these are **deferred automatically**
4. Observation blocks

You can use one derived variable to define another based on their order within a single
`@variables` block.

### Deferred system lines

A system line that mentions a body by name (`b.i`) is evaluated after every body block,
so couplings, period ratios, and constraints on derived quantities need no special
syntax:

```julia
sys = System(name="coplanar", bodies=[A, b, c], variables=@variables begin
    plx ~ Uniform(1, 100)
    i_shared ~ Sine()          # ordinary: hoisted, and read by both bodies
    mut_inc = b.Ω - c.Ω        # deferred: mentions a body, so it runs last
    mut_inc ~ Normal(0, 0.1)   # …and can then be constrained
end)
```

Deferral is transitive, and is detected by a static walk over the stored expressions. A
*body* that tries to read a deferred system variable is a genuine cycle, and is reported
as an error naming both sides.

!!! warning "Common Pitfall: Accessing Variables"
    Inside a `Body` variables block, you **cannot** use `planet.variable_name` or the
    body's own name to access its variables, and you cannot reach into a sibling.
    Instead:

    - Use the variable name directly if it's defined earlier in the same `@variables` block
    - Use `system.variable_name` to access system-level variables
    - Move anything shared between two bodies up to the `System` block

    ```julia
    # ❌ WRONG - cannot use b.P_days inside b's own block
    variables=@variables begin
        P_days ~ Uniform(1, 100)
        tp = b.P_days * 0.5  # This will NOT work!
    end

    # ✅ CORRECT - use the variable name directly
    variables=@variables begin
        P_days ~ Uniform(1, 100)
        tp = P_days * 0.5
    end

    # ✅ CORRECT - access system variables via system.
    variables=@variables begin
        P_inner = system.P_inner  # Access from system
        tp = P_inner * 0.5
    end
    ```

## Interpolation Syntax with `$`

The `$` symbol in `@variables` blocks allows you to inject external values and functions into derived variable expressions. This is useful for referencing constants, interpolation functions, or any value defined outside the model.

```julia
# Define constants and functions outside the model
HOST_L_MAG = 4.5
mass_to_contrast(m) = m^0.5  # Simple example; use real atmosphere models in practice

planet_b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass ~ LogUniform(1mjup, 50mjup)   # [M⊙]
        flux_L = $mass_to_contrast(mass)   # Interpolate the function
        a ~ Uniform(1, 100)
        e ~ Uniform(0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        tp ~ Uniform(50000, 60000)
    end
)
```

!!! warning "Interpolation Limitations"
    The `$` interpolation only works for **simple references**. Complex expressions inside `$()` will cause a syntax error:

    ```julia
    # ❌ This will FAIL with "syntax: $ expression outside quote"
    flux_L = $mass_to_L_contrast(mass, system.age, $HOST_L_MAG)
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
        mass ~ LogUniform(1mjup, 50mjup)
        flux_L = $mass_to_L_contrast_wrapper(mass)
    end
    ```

    Alternatively, if you need to use model variables (not constants), include them as arguments without `$`:

    ```julia
    # ✅ Model variables don't need $ - they're resolved at runtime
    variables=@variables begin
        mass ~ LogUniform(1mjup, 50mjup)
        flux_L = $mass_to_L_contrast(mass, system.age)
    end
    ```

## Variable Scoping in Observations

Observations are no longer attached to a planet; every observation lives on the `System`
and names its own `target` and `ref`. Consequently an observation's `@variables` block
sees exactly one namespace: `system.X`, including any deferred system variables.

```julia
rv_obs = RadialVelocityObs(
    rv_table;
    target = A,
    ref = Barycentre,
    name = "HARPS",
    variables = @variables begin
        offset ~ Normal(0, 100)
        jitter ~ LogUniform(0.1, 100)
        trend_slope = system.trend_slope  # From system variables
    end
)

sys = System(
    name="HD12345",
    bodies=[A, planet_b],
    observations=[rv_obs],
    variables=@variables begin
        plx ~ truncated(Normal(50, 1), lower=1)
        trend_slope ~ Uniform(-1, 1)      # Accessed as system.trend_slope
    end
)
```

If an observation needs a *body's* variable, expose it through the system block first —
which is a deferred line, and therefore already resolved by the time observation blocks
are evaluated:

```julia
variables=@variables begin
    plx ~ truncated(Normal(50, 1), lower=1)
    b_mass = b.mass          # deferred: mentions a body
end
```

Then read `system.b_mass` in the observation's block.

!!! note "Fluxes are body variables now"
    A per-band flux or contrast ratio is declared on the *body* as `flux_<band>`, not on
    the observation. Setting the host's `flux_H = 1.0` makes every other body's `flux_H` a
    contrast ratio against it. Observation types that used to own a `flux`, `fluxratio` or
    `flux` vector variable now read the bodies', and several of them raise an explicit
    error if you declare one on the observation. See
    [Migrating to Octofitter v2](@ref v2-migration).
