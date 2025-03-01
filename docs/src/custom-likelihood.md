# Adding a Custom Likelihood Function

It's fairly straightforward to add supprot for a new kind of observation to Octofitter.jl
You can also follow the same workflow if you want to handle an existing kind of observation in a new way—say, tweaking a calculation, or using Gaussian processes to better model noise in radial velocity data.

All the existing observation types are listed in  [`src/observations`](https://github.com/sefffal/Octofitter.jl/tree/main/src/observations)
and can be used as examples.

Note that these examples won't run if you copy and paste them, you'll need to modify them to suite your purposes.

## Creating a Likelihood type

The first step is to create a new data type to hold the observations. 

```julia
struct MyLikelihood{TTable<:Table} <: AbstractLikelihood
    table::TTable
    function MyLikelihood(observations...)
        table = Table(observations...)
        return new{typeof(table)}(table)
    end
end
```

Here we create a struct `MyLikelihood` that is a subtype of `AbstractLikelihood`. It's entirely up to you how you store the data in the struct. Here, we accept one or more arguments and forward them to the `TypedTables.Table` constructor so we can grab them out efficiently during sampling.

Try to follow the advice in the Julia Manual's performance tips section to ensure you've created a fully "concrete" type. This won't affect corectness, but will be important for performance down the road.

## Create likelihood functions

Now, create a method that extends `Octofitter.ln_like` for your custom observation type. 
If the likelihood function is specific to a planet (like astrometry, where the data is attached to a planet instead of the system) then the method signature should look like:

```julia
function Octofitter.ln_like(like::MyLikelihood, θ_planet, orbit,)

    # Access your data
    # like.table.col1[1]

    # Access planet variables
    # θ_planet.e

    # Access planet variables
    # θ_planet.e

    # Solve for position
    # x = raoff(orbit, like.table.date[1])

    return 0.0
end
```

If the data is attached to the system as a whole, like radial velocity, the method signature should look like:
```julia
function Octofitter.ln_like(like::MyLikelihood, θ_system, orbits)

    # Access your data
    # like.table.col1[1]

    # Access system variables
    # θ_system.M

    # Access planet variables
    # θ_system.planets.b.e
    
    return 0.0
end
```

Inside your method, you should calculate the log-likelihood of the data stored in `obs` given the parameters.

The `θ_planet` or `θ_system` variables hold the current parameters in a NamedTuple. For instance, `θ_planet` might look like `(a=12.4, e=0.1, ...)` and `θ_system` might look like `(M=1.52, plx=123.0, planets=(b=(a=12.4, e=0.1,), c=(a=5.3, e=0.4, ...)))`.

`orbit` or `orbits` will be set to a PlanetOrbits.jl orbit object or list of orbit objects. These can be used to efficiently solve orbits at a given time. They are created out of the `θ` parameters for you so that calculations of common factors can be shared across observation types and dates.

If your data requires a new parameter or parameters in the model, simply pass them in the `Variables` block of the planet or system and access them by name in `θ_planet` or `θ_system`. Everything will be handled for you.
If the parameter has a restricted domain where it is valid, ensure the prior passed in is truncated using `Distributions.truncate`. Then, the code will remap the variable automatically using Bijectors.jl to prevent invalid values.

## Bonus: Generative model
The above is sufficient to start sampling from the posterior. Ideally, you will also add a function that does the reverse: generate observations from a set of parameters. This is useful for a variety of statistical tests.

Simpling extend the `Octofitter.generate_from_params` function for your data type in much the same way:

```julia
# Generate new astrometry observations
function generate_from_params(like::MyLikelihood, elem::PlanetOrbits.AbstractOrbit, θ_planet)

    # Get epochs from observations (assuming you have an epoch column)
    epochs = like.table.e

    # Generate new data at that epoch based on current parameters. i.e.
    # "what would we see at date X if the true parameters were θ_planet"
    ras = raoff.(elem, epochs)
    decs = decoff.(elem, epochs)

    # Then, generate a new MyLikelihood with the simulated values
    astrometry_table = Table()

    return MyLikelihood(epoch=epochs, ra=ras, dec=decs)
end
```
