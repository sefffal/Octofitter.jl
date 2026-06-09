# Frequently Asked Questions



## How do you define position angle of the ascending node / longitude of the ascending node?

Ω is the **angle** at which the planet (or star) crosses the plane of the sky moving away from the observer, which is to say, moving from negative to positive on the z-axis (consistent with positive RV = redshift). See the more detailed entry "What is the definition of $\Omega$?" below for a discussion of why this is consistent for both the planet and the star.

## How do we calculate the position of a planet at a future epoch?

After fitting an orbit, there are two ways to calculate its projected position in the sky at a future epoch.

**Option 1**

The easiest way to get a report on the planet's future position is via `octoplot`
```
octoplot(model, chain, mark_epochs_mjd=[
    mjd("2028-01-01"),
    mjd("2029-01-01"),
    # etc
])
```
This will add a scatter point at the dates you requested and output a summary to the terminal. You can also control the alpha and number of orbits plotted. Something like this:
`octoplot(model, chain, mark_epochs_mjd=[mjd("2025-01-01"), mjd("2027-01-01")], N=100, alpha=1.0)`

**Option 2**

The other is to calculate the positions yourself and plot them however you like.
```
els = Octofitter.construct_elements(model, chain, :b, :)
```
Where `:b` is the name of the planet you chose, and `:` means all draws from the chain (you can also pass a particular iteration number of range of numbers).
Then, you can solve keplers equation and plot whatever you want:
```
sols = orbitsolve.(els, mjd("2025-01-01"))

# projected position in mas
X = raoff.(sols)
Y = decoff.(sols)

Proj_sep_mas = projectedseparation.(sols)
PA_rad = posangle.(sols)

# 3D position in AU
Z_au = posz.(sols)
X_au = posx.(sols)

# relative RV between star and companion
Rel_rv = radvel.(sols)
```

And so on. You can also query these values for the star.
Then, you can plot them however you like,
Eg.
```
scatter(X, Y, axis=(;xlabel="delta ra mas", autolimitaspect=1.0, xreversed = true))
```


## Can slope/GP parameters be shared between RV instruments?

You can share instrument parameters such as linear our quadratic terms between instruments by defining the variables at the system level, and forwarding them to each instrument's `@variables` block (see `<---` arrows):

```
# Instrument 1 likelihood
rvlike_apf = StarAbsoluteRVObs(
    rv_dat_apf,
    name="APF",

    # Linear trend
    trend_function = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000),
    
    variables=@variables begin
        offset = 0
        jitter ~ LogUniform(0.1,30) # m/s
        trend_slope = system.trend_slope  # <-----  
    end
)
# Instrument 2 likelihood
rvlike_hires = StarAbsoluteRVObs(
    rv_dat_hires,
    name="HIRES",
    
    # Linear trend:
    trend_function = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000), 
    
    variables=@variables begin
        offset = 0
        jitter ~ LogUniform(0.1,30) # m/s
        trend_slope = system.trend_slope  # <-----  
    end
)
sys = System(
    name = "Star1",
    companions=[planet_1],
    observations=[rvlike_apf, rvlike_hires],
    variables=@variables begin
        M ~ truncated(Normal(1.5, 0.06),lower=0.1, upper=10)

        trend_slope ~ Uniform(-1,1)  # <-----  
    end
)
```


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
flux = $mass_to_L_contrast(planet.mass, system.age, $HOST_L_MAG)

# ✅ This WORKS - simple function reference, model variables without $
flux = $mass_to_L_contrast(planet.mass, system.age, planet.temp)
```

**Solution**: Create a wrapper function that captures your constants:

```julia
HOST_L_MAG = 4.5
TRUE_AGE = 10.0

function mass_to_L_contrast_wrapper(mass)
    return mass_to_L_contrast(mass, TRUE_AGE, HOST_L_MAG)
end

# Then use the simple wrapper
flux = $mass_to_L_contrast_wrapper(mass)
```

See [Derived Variables - Interpolation Syntax](@ref derived) for more details.

## Should I use `planet.X` or `system.X` in observation variables?

It depends on **where the observation is attached** and **where the variable is defined**:

- **Observation attached to a Planet** (e.g., `PhotometryObs`, `PlanetRelAstromObs`):
  - Use `planet.X` for variables defined on the Planet
  - Use `system.X` for variables defined on the System

- **Observation attached to a System** (e.g., `StarAbsoluteRVObs`):
  - Use `system.X` for variables defined on the System

```julia
# PhotometryObs attached to a planet:
H_band_data = PhotometryObs(
    data_table,
    name="H_band",
    variables=@variables begin
        # mass is on planet, age is on system
        flux = $H_band_contrast_interp(planet.mass, system.age)
    end
)

planet_b = Planet(
    name="b",
    observations=[H_band_data],  # Attached to planet
    variables=@variables begin
        mass ~ Uniform(0, 10)  # Access via planet.mass
        # ...
    end
)
```

See [Derived Variables - Variable Scoping in Observations](@ref derived) for more details.

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
