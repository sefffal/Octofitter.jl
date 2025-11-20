# Frequently Asked Questions



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
rvlike_apf = StarAbsoluteRVLikelihood(
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
rvlike_hires = StarAbsoluteRVLikelihood(
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


## What Coordinate System does Octofitter use?
Octofitter uses a coordinate system where
* $+x$ increases to the East (ie, x increases with increasing Right Ascension)
* $+y$ increases to the North
* $+z$ increases away from the observer.

This coordinate system has several nice properties: $+x$ increases with Right Ascension, $+y$ increases with Declination, and $\frac{\partial z}{\partial t}$ is positive for positive radial velocity / positive redshift.
**Note!** Since x increases with Right Ascension, that means +X is towards the **left** in the sky, and in most plots. If this bothers you, one thing to note is that it's also in the direction of increasing "time" as the sky rotates overhead AND this convention predates the modern concept of a cartesian graph :-) 

You can read more in the PlanetOrbits.jl docs.

## What is the definition of $\omega$?
$\omega$ is the argument of periastron, which is the location where the **planet** makes its closest approach to the star, measured from the ascending node. This is consistent with most direct imaging conventions. 

**Note:** this is 180° offset from the typical definition used by codes that only fit radial velocity and/or transit, where the convention it to report the argument of periastron for the star. This is a significant potential source of confusion when comparing results between codes.

## What is the definition of $\Omega$?
$\Omega$ is the position angle of ascending node, also known as the longitude of ascending node. It is the point in an orbit where the planet (or equivalently, the star) moves from having a negative $z$ coordinate to having a positive $z$ coordinate. This happens where the planet (or star) moves cross the plane of the sky going **away from the observer**. 
Why "away" from the observer? That is because Octofitter uses a coordinate system where $+z$ increases away from the observer, such that radial velocity measured as a positive redshift corresponds to a positive velocity.
