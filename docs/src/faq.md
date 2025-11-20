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
