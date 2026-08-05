# Frequently Asked Questions



## How do you define position angle of the ascending node / longitude of the ascending node?

Ω is the **angle** at which the planet (or star) crosses the plane of the sky moving away from the observer, which is to say, moving from negative to positive on the z-axis (consistent with positive RV = redshift). See the more detailed entry "What is the definition of $\Omega$?" below for a discussion of why this is consistent for both the planet and the star.

## How do we calculate the position of a planet at a future epoch?

After fitting an orbit, rebuild the whole system for one or more draws from the chain and
query it at whatever epochs you like. `construct_system` replaces v1's
`construct_elements`: a draw is now one `PlanetOrbits.System` containing every body and
every orbit, rather than a per-planet orbit object.

```julia
posys = construct_system(model, chain, 1)   # draw #1; pass `:` for all draws
sols = orbitsolve(posys, [mjd("2028-01-01"), mjd("2029-01-01")])
```

Then read off whatever you need. Every observable takes `(solution, target, reference)`:

```julia
sol = sols[1]

# projected position in mas
X = raoff(sol, :b, :A)
Y = decoff(sol, :b, :A)

Proj_sep_mas = projectedseparation(sol, :b, :A)
PA_rad = posangle(sol, :b, :A)

# 3D position in AU
Z_au = posz(sol, :b, :A)
X_au = posx(sol, :b, :A)

# relative RV between star and companion
Rel_rv = radvel(sol, :b, :A)

# the star's reflex RV, relative to the system barycentre.
# Note the lowercase `barycentre`: it resolves the reference point against THIS
# sample's masses. The capitalized `Barycentre` is the declaration spec you write
# in an observation's `ref=`, and is not accepted by the observable functions.
Star_rv = radvel(sol, :A, PlanetOrbits.barycentre(posys))
```

Because the reference is explicit, the same call works for a companion measured against
another companion (`raoff(sol, :c, :b)`), against the barycentre, or against a
photocentre.

To build a predicted distribution rather than a single number, loop over draws:

```julia
ep = mjd("2028-01-01")
X = map(1:size(chain,1)) do i
    sol = orbitsolve(construct_system(model, chain, i), [ep])[1]
    raoff(sol, :b, :A)
end
scatter(X, axis=(;xlabel="draw", ylabel="ΔRA [mas]"))
```

!!! note "`mark_epochs_mjd` is gone"
    v1's `octoplot(model, chain, mark_epochs_mjd=[...])` has no v2 equivalent yet. Compute
    the predicted positions as above and add them to the figure yourself — `octoplot`
    returns named axes (`res.axes`) precisely so that you can.


## Can slope/GP parameters be shared between RV instruments?

You can share instrument parameters such as linear or quadratic terms between instruments by defining the variables at the system level, and forwarding them to each instrument's `@variables` block (see `<---` arrows):

```julia
# Instrument 1
rvlike_apf = RadialVelocityObs(
    rv_dat_apf,
    target = A,
    ref = Barycentre,
    name = "APF",

    # Linear trend
    trend_function = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000),

    variables=@variables begin
        offset ~ Normal(0, 100)     # m/s
        jitter ~ LogUniform(0.1,30) # m/s
        trend_slope = system.trend_slope  # <-----
    end
)
# Instrument 2
rvlike_hires = RadialVelocityObs(
    rv_dat_hires,
    target = A,
    ref = Barycentre,
    name = "HIRES",

    # Linear trend:
    trend_function = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000),

    variables=@variables begin
        offset ~ Normal(0, 100)     # m/s
        jitter ~ LogUniform(0.1,30) # m/s
        trend_slope = system.trend_slope  # <-----
    end
)
sys = System(
    name = "Star1",
    bodies=[A, planet_1],
    observations=[rvlike_apf, rvlike_hires],
    variables=@variables begin
        plx ~ truncated(Normal(50, 1), lower=1)

        trend_slope ~ Uniform(-1,1)  # <-----
    end
)
```

!!! warning "`offset` and `jitter` are no longer added for you"
    v1's `StarAbsoluteRVObs` silently injected `offset ~ Uniform(-1000, 1000)` and
    `jitter ~ LogUniform(0.001, 100)` when you gave it no `variables=` block. v2 never
    invents a prior. A v1 RV model copied across without those two lines fits with **no
    zero point and no jitter**, which will change your answers rather than error.

## What conventions does Octofitter use for orbital elements?

Octofitter solves orbits using [PlanetOrbits.jl](https://github.com/sefffal/PlanetOrbits.jl), which adopts the same conventions as Orbitize!. The full reference (including a derivation of projected position, velocity, and acceleration) lives in the [PlanetOrbits.jl coordinate conventions docs](https://sefffal.github.io/PlanetOrbits.jl/dev/conventions/), but the key points are summarized here for convenience.

**Coordinate system** (right-handed, observer-centric):

| Axis | Direction |
|------|-----------|
| $+x$ | Toward the East — increasing Right Ascension (note: this points to the **left** in the sky and in most plots) |
| $+y$ | Toward the North — increasing Declination |
| $+z$ | **Away** from the observer, so that $\partial z / \partial t > 0$ corresponds to a positive radial velocity (redshift) |

**Orbital elements:**

| Element | Meaning |
|---------|---------|
| $a$ | Semi-major axis (au) |
| $e$ | Eccentricity, range $[0, 1)$ |
| $i$ | Inclination (rad). $i = 0$ is **face-on**; $i = 90°$ is **edge-on** |
| $\omega$ | Argument of periastron of the **planet**, measured from the ascending node |
| $\Omega$ | Position angle / longitude of the ascending node (rad), measured counter-clockwise in the sky plane from North ($+y$) |
| $t_p$ | Epoch of periastron passage (MJD). Equivalent values differ by integer multiples of the period: $t_p' = t_p + iP$ |
| $P$ | Orbital period |
| $M$ | Total/central mass (solar masses) |
| `plx` | Parallax (mas), setting the system distance |

The individual elements are discussed in more detail in the entries below.

## What Coordinate System does Octofitter use?
Octofitter uses a coordinate system where
* $+x$ increases to the East (ie, x increases with increasing Right Ascension)
* $+y$ increases to the North
* $+z$ increases away from the observer.

This coordinate system has several nice properties: $+x$ increases with Right Ascension, $+y$ increases with Declination, and $\frac{\partial z}{\partial t}$ is positive for positive radial velocity / positive redshift.
**Note!** Since x increases with Right Ascension, that means +X is towards the **left** in the sky, and in most plots. If this bothers you, one thing to note is that it's also in the direction of increasing "time" as the sky rotates overhead AND this convention predates the modern concept of a cartesian graph :-) 

You can read more in the PlanetOrbits.jl docs.

## What is the definition of inclination ($i$)?
Octofitter (following PlanetOrbits.jl and Orbitize!) uses the convention where:

- $i = 0°$ is a **face-on** orbit — the orbital plane lies in the plane of the sky.
- $i = 90°$ is an **edge-on** orbit — the line of sight lies in the orbital plane.

Be aware that both conventions for $i=0$ appear in the literature, so it is worth double-checking when comparing across packages.

## What is the definition of $\Omega$?
$\Omega$ is the position angle of the ascending node, also known as the longitude of the ascending node. It is the **angle** at which the planet (or equivalently, the star) crosses the plane of the sky moving from a negative $z$ coordinate to a positive $z$ coordinate, i.e. moving **away from the observer**. It is measured counter-clockwise in the sky plane starting from North ($+y$).

Why "away" from the observer? That is because Octofitter uses a coordinate system where $+z$ increases away from the observer, such that radial velocity measured as a positive redshift corresponds to a positive velocity.

#### But isn't that contradictory for the star and the planet?

This is a common and very reasonable point of confusion. If $\Omega$ is defined where a body crosses the sky plane moving *away* from us (positive redshift), but the star and planet always move in opposite directions, then shouldn't $\Omega$ be different depending on whether you track the star or the planet?

There is actually no contradiction. The key is that $\Omega$ is the **angle** of the ascending node, and this angle is identical whether you track the star or the planet — even though:

- the two bodies pass through their respective ascending nodes at **different times** (180° out of phase), and
- for eccentric orbits, or when the masses are unequal, they cross the sky plane at **different points in space**.

Concretely: when the planet crosses the sky plane heading toward $+z$ (away from us, redshifted), that defines $\Omega$. At that *same instant* the star is crossing the sky plane heading toward $-z$ (toward us, blueshifted) — the star is at its *descending* node. Half an orbit later, the star reaches its own *ascending* node, crossing toward $+z$, and it does so at the same position angle $\Omega$. So both bodies share the same $\Omega$; they simply arrive there at opposite phases (and, in general, at different points along the line of nodes).

## What is the definition of $\omega$?
$\omega$ is the argument of periastron, which is the location where the **planet** makes its closest approach to the star, measured from the ascending node. This is consistent with most direct imaging conventions. 

**Note:** this is 180° offset from the typical definition used by codes that only fit radial velocity and/or transit, where the convention is to report the argument of periastron for the star. This is a significant potential source of confusion when comparing results between codes.

Note the contrast with $\Omega$ above. For $\Omega$, the star and planet share the *same angle* but reach the ascending node at *different times* (180° out of phase). For $\omega$, it is inverted: the planet and star reach periastron at the *same time* (they are always collinear with the barycenter), but the planet-referenced and star-referenced values of $\omega$ differ by *180° in angle*.

## I get a syntax error with `$` interpolation in `@variables`

If you see an error like `syntax: "$" expression outside quote` when using `$` interpolation in derived variables, it's likely because you have a complex expression inside `$()`.

The `$` interpolation only works for **simple references** to external variables or functions. For example:

```julia
# ❌ This FAILS - nested $ or complex expressions don't work
flux_L = $mass_to_L_contrast(mass, system.age, $HOST_L_MAG)

# ✅ This WORKS - simple function reference, model variables without $
flux_L = $mass_to_L_contrast(mass, system.age, tempK)
```

**Solution**: Create a wrapper function that captures your constants:

```julia
HOST_L_MAG = 4.5
TRUE_AGE = 10.0

function mass_to_L_contrast_wrapper(mass)
    return mass_to_L_contrast(mass, TRUE_AGE, HOST_L_MAG)
end

# Then use the simple wrapper
flux_L = $mass_to_L_contrast_wrapper(mass)
```

See [Derived Variables - Interpolation Syntax](@ref derived) for more details.

## Should I use `planet.X` or `system.X` in observation variables?

Always `system.X`. Observations are no longer attached to a planet — every observation
lives on the `System` and names its own `target` and `ref` — so an observation's
`@variables` block sees exactly one namespace, `system.`.

If an observation needs a *body's* variable, expose it through the system block first. A
system line that mentions a body by name is deferred until after every body block has
run, so it is already available by the time observation blocks are evaluated:

```julia
sys = System(
    name="HD12345",
    bodies=[A, planet_b],
    observations=[my_obs],
    variables=@variables begin
        age ~ Normal(15, 1)
        b_mass = b.mass        # deferred: mentions a body
    end
)
```

Then write `system.age` and `system.b_mass` inside `my_obs`'s variables block.

Note that per-band fluxes and contrast ratios are no longer observation variables at all:
they are declared on the body as `flux` or `flux_<band>`. Several observation types raise
an explicit error if you declare a `flux` variable on the observation. See
[Migrating to Octofitter v2](@ref v2-migration).

See [Derived Variables](@ref derived) for more details.

## What does the warning "Too many steps without any function evaluations" mean?

During model initialization with `initialize!()`, you may see a warning like:
```
Warning: Unrecognized stop reason: Too many steps (101) without any function evaluations
```

!!! info "This warning is safe to ignore"
    This warning comes from the underlying optimization library (Optim.jl via Pathfinder.jl) and indicates that the optimizer's line search has converged or reached a point where it cannot make further progress. **This warning is completely benign and does not affect your results.**

The warning typically appears when:
- The optimizer has already found a good starting point
- The line search algorithm determines no further step would improve the solution
- The optimization has effectively converged

**Your fit will proceed correctly.** Check the output of `initialize!()` - if it returns a valid chain with reasonable parameter values, the initialization was successful regardless of this warning.

## How to load past rounds saved with the Pigeons MCMC?
When using the Pigeons MCMC with Octofitter there is a checkpoint feature which is automatically set as false. When set as true, each round completed using the Pigeons MCMC will be saved to a results folder that will be automatically created. Inside this folder there will be two additional folders, all and latest. The former holders all past Pigeons runs while the latest folder contains the results from the most recent run. From these folders, you can reload old Pigeons MCMC rounds and start a run again from the saved checkpoint location (the last round found in a given folder). This can save time when running particularily long models.

When starting a new model with the Pigeons MCMC, compute the chain in the regular method but set the checkpoint variable to true,

```
chain, pt = octofit_pigeons(model, n_rounds=12, checkpoint=true)
```

Then, in future runs you can load past Pigeons MCMC rounds using the following:

```
pt = PT("path/to/results/folder")
pt = increment_n_rounds!(pt, number_of_rounds_to_add)
chain, pt = octofit_pigeons(pt)
```

For example, if I wanted to add 1 additional round to my model starting from the last round I ran, my code would look as follows:

```
#chain, pt = octofit_pigeons(model, n_rounds=12, checkpoint=true)
pt = PT("results/latest")
pt = increment_n_rounds!(pt, 1)
chain, pt = octofit_pigeons(pt)
```

Always be sure to comment out the initial creation of the chain and pt when you want to load and use past rounds. This will allow you to reload the last round and run additional rounds without having to restart the entire model.
