#  [Derived Variables] (@id derived)

DirectDetections has a concept called "derived variables" that are inspired by PyMC3.
Derived variables are quantities that either have a fixed value, or a fixed mathematical relationship with the main variables in a model.

This concept is extremely powerful, as it lets you quickly create very sophisticated models.

Derived variables allow you to mark certain properties as constants, reparameterize models, link properties between planets in multi-planet systems, plug in physical models, and more.

## System Variables
Derived variables for the system as a whole can be created as follows:

```julia
@named HD12345 = System(
    Variables(
        M = sys -> 1.0
        plx = Normal(45., 0.02),
    )
)
```
In this case, instead of including `M` as a variable in the model, we define it as a function that always returns `1.0`. This is equivalent to passing `M=1.0`.

In the following case, let's define `M` as being calculated based on another variable in the model. This is how you can do reparameterizations in DirectDetections.jl
```julia
@named HD12345 = System(
    Variables(
        plx =Normal(45., 0.02),
        logM =Normal(45., 0.02),
        M = sys -> 10^sys.logM,
    )
)
```
We defined a new variable `logM` under priors, and then used the `Derived` block to calculate `M` from the other variables in the model.

In general, you can write any function you want to map from any of combination of constants and variables in the model to new variables. The only constraints are that your functions always return the same outputs for the same inputs, and are differentiable. These functions will be called in a tight loop, so you should try to make sure they are as efficient as possible.


## Planet Variables
Derived variables for an individual planet are similar, but have access to both the planet's variables and the system as a whole.

Here is an example of reparameterizing `e` and `a` on a planet to be logarithmic quantities:
```julia
@named b = Planet{KeplerianElements}(
    Variables(
        τ = Normal(0.5, 1),
        ω = Normal(0.1, deg2rad(30.)),
        i = Normal(0.1, deg2rad(30.)),
        Ω = Normal(0.0, deg2rad(30.)),
        loge = Uniform(-4, 1),
        loga = Normal(1, 1)
        e = (sys, pl) -> 10^pl.loge,
        a = (sys, pl) -> 10^pl.loga,
    ),
)
```
Here `e` is defined as log-uniform, and `a` as log-normal.


## Linking Planets
Planets can have Derived variables that are calculated from variables defined on the system as a whole.
This makes it easy to, for example, create a system of two planets that are co-planar.

```julia
@named b = Planet{KeplerianElements}(
    Variables(
        a = Uniform(0, 15),
        e = TruncatedNormal(0, 0.1, 0, 1),
        ω = Normal(0.1, deg2rad(30.)),
        τ = Normal(0.5, 1),
        i = (sys, pl) -> sys.i,
        Ω = (sys, pl) -> sys.Ω,
    ),
)
@named c = Planet{KeplerianElements}(
    Variables(
        a = Uniform(15, 45),
        e = TruncatedNormal(0, 0.1, 0, 1),
        ω = Normal(0.1, deg2rad(30.)),
        τ = Normal(0.5, 1),
        i = (sys, pl) -> sys.i,
        Ω = (sys, pl) -> sys.Ω,
    ),
)
@named HD12345 = System(
    Variables(
        plx = Normal(45., 0.02),
        M = Normal(45., 0.02),
        i = Normal(0.1, deg2rad(30.)),
        Ω = Normal(0.0, deg2rad(30.)),
    ),
    b, c
)
```
Notice how `i` and `Ω` are defined as variables on the System. The two planets B & C omit those two variables from their priors, and instead just take their values from the System. This way we can enforce co-planarity between planets without e.g. rejection sampling.

## Resolution Order
The order that variables are resolved is as follows:
* Variables defined as priors for the system and planets
* Derived variables on the system
* Derived variables on each planet

It's not possible to observe one derived variable from another derived variable in a System or one Derived variable from another on a Planet but it is possible to see the System's Derived variables in the Planets'.
