# Fit Relative RV Data

Octofitter includes support for fitting relative radial velocity data — the
velocity of a companion measured against its host, rather than the reflex
velocity of the host itself.

The convention we adopt is that positive relative radial velocity is the velocity of the companion (exoplanets) minus the velocity of the host (star).

!!! note "Relative RV is the same likelihood as reflex RV"
    In v2 there is one radial velocity observation, [`RadialVelocityObs`](@ref).
    What distinguishes relative RV from the stellar reflex signal is *which
    references you name*:

    ```julia
    RadialVelocityObs(tab; target=A, ref=Barycentre)   # the star's reflex motion
    RadialVelocityObs(tab; target=b, ref=A)            # the companion, relative to the star
    ```

    That is the whole difference — v1's separate `PlanetRelativeRVObs` type is
    gone, and with it its own copy of the reflex-superposition bookkeeping.
    `ref` is a real reference, so `target=c, ref=b` (one companion measured
    against another) is expressible too.

To fit relative RV data, start by building the bodies and then the observation:
```@example 1
using Octofitter
using CairoMakie
using Distributions

rv_dat_1 = Table(
    epoch=55000:100:57400,
    rv = [
         -24022.74
        -18571.33
        14221.56
        26076.89
        -459.26
        -26319.26
        -13430.96
        19230.96
        23580.26
        -6786.28
        -27161.78
        -7548.58
        23177.95
        19780.94
        -12738.39
        -26503.74
        -1249.19
        25844.47
        14888.83
        -17986.76
        -24381.49
        5119.22
        27083.2
        9174.18
        -22241.45
    ],
    # Hint! Type as \sigma + <TAB>
    σ_rv= fill(15000.0, 25),
)
nothing # hide
```
See the [Basic RV Fit](@ref fit-rv) tutorial for examples on how this data can be loaded from a CSV file.

```@example 1
A = Body(
    name="A",
    variables=@variables begin
        mass = system.M_tot   # the companion is massless here, so the star carries it all
    end
)

b = Body(
    name="b",
    about=A,
    variables=@variables begin
        mass = 0.0            # M⊙; a test particle
        a ~ Uniform(0, 10)    # AU
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ Uniform(0, 2pi)
        Ω ~ Uniform(0, 2pi)
        # Kepler's third law, in years, from the total mass of this orbit.
        P_yrs = √(a^3 / system.M_tot)
        τ ~ Uniform(0, 1.0)
        tp = τ * P_yrs * 365.25 + 60000 # reference epoch for τ. Choose an MJD date near your data.
    end
)

rel_rv_obs = RadialVelocityObs(
    rv_dat_1;
    target=b, ref=A,           # <-- this is what makes it *relative* RV
    name="simulated data",
    variables=@variables begin
        jitter ~ LogUniform(0.1, 1000) # m/s
    end
)
nothing # hide
```

!!! note "Why `M_tot` lives in the system block"
    A body's `@variables` block can read `system.*` but never its siblings, so a
    quantity two bodies share is *hoisted* to the system. Here `M_tot` gives the
    star its mass and gives the orbit its period; had we instead sampled
    `A.mass` directly, `b`'s block could not have seen it.

The relative RV likelihood does not need an instrument-specific zero point —
the two stellar spectra are differenced against each other, so there is nothing
to offset. You may still declare an `offset` variable if your reduction has one;
as with reflex RV, nothing is added for you. A `jitter` parameter can be
specified in the observation's `@variables` block, as can parameters for a
Gaussian process model of correlated noise (new for relative RV in v2 — see
[Fit Gaussian Process](@ref fit-rv-gp)). Create one `RadialVelocityObs` per
instrument if you have several, each with its own jitter.

Next, assemble the system. Observations are a flat list on the system; they are
never attached to a body, because each one names its own `target` and `ref`.

```@example 1
sys = System(
    name="Example_System",
    bodies=[A, b],
    observations=[rel_rv_obs],
    variables=@variables begin
        M_tot ~ truncated(Normal(1.2, 0.1), lower=0.1) # total mass in solar masses
    end
)

model = Octofitter.LogDensityModel(sys)
```


### Initialize the model and verify starting point

```@example 1
init_chain = initialize!(model)

octoplot(model, init_chain)
```


```@example 1
using Random
rng = Random.Xoshiro(123)
chain = octofit(rng, model)
```


```@example 1
octoplot(model, chain)
```
