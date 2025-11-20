# Adding a Custom Likelihood Function

It's fairly straightforward to add support for a new kind of observation to Octofitter.jl
You can also follow the same workflow if you want to handle an existing kind of observation in a new way—say, tweaking a calculation, or using Gaussian processes to better model noise in radial velocity data.

All the existing observation types are listed in  [`src/likelihoods`](https://github.com/sefffal/Octofitter.jl/tree/main/src/likelihoods)
and can be used as examples.

Note that these examples won't run if you copy and paste them, you'll need to modify them to suite your purposes.

## Creating a Likelihood type

The first step is to create a new data type to hold the observations. 

```julia
"""
    data = Table(
        (epoch=50000, my_measurement=1.2, σ_measurement=0.1),
        (epoch=50100, my_measurement=1.5, σ_measurement=0.1),
    )
    MyLikelihood(
        data,
        name="MY_INSTRUMENT",
        variables=@variables begin
            my_parameter ~ Normal(0, 1)
        end
    )

A custom likelihood for my specific type of observations.
"""
struct MyLikelihood{TTable<:Table} <: AbstractLikelihood
    table::TTable
    name::String
    priors::Priors
    derived::Derived
    function MyLikelihood(
            observations;
            name="MY_LIKELIHOOD",
            variables::Tuple{Priors,Derived}=(@variables begin;end)
        )
        (priors,derived)=variables
        table = Table(observations)
        if !equal_length_cols(table)
            error("The columns in the input data do not all have the same length")
        end
        # Add any column validation here, e.g.:
        # expected_cols = (:epoch, :my_measurement, :σ_measurement)
        # if !issubset(expected_cols, Tables.columnnames(table))
        #     error("Expected columns $expected_cols")
        # end
        return new{typeof(table)}(table, name, priors, derived)
    end
end
export MyLikelihood
```

Here we create a struct `MyLikelihood` that is a subtype of `AbstractLikelihood`. The new API includes:

- **`table`**: The observational data as a TypedTables.Table
- **`name`**: Used for variable naming in MCMC chains
- **`priors`** and **`derived`**: Observation-specific variables from the `@variables` block

Try to follow the advice in the Julia Manual's performance tips section to ensure you've created a fully "concrete" type. This won't affect correctness, but will be important for performance down the road.

## Create likelihood functions

Now, create a method that extends `Octofitter.ln_like` for your custom observation type. 

If the likelihood function is specific to a planet (like astrometry, where the data is attached to a planet instead of the system) then the method signature should look like:

```julia
# MyLikelihood: attached to a planet
function Octofitter.ln_like(like::MyLikelihood, θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    # Access your data from the table
    for i_epoch in eachindex(like.table.epoch)
        epoch = like.table.epoch[i_epoch]
        measurement = like.table.my_measurement[i_epoch]
        σ_measurement = like.table.σ_measurement[i_epoch]
        
        # Access planet variables
        # θ_planet.e, θ_planet.a, etc.
        
        # Access observation-specific variables  
        my_param = θ_obs.my_parameter
        
        # Method 1: Use pre-solved orbit solutions (efficient!)
        # Get the pre-solved orbit solution for this planet at this epoch
        sol = orbit_solutions[i_planet][i_epoch + orbit_solutions_i_epoch_start]
        
        # Extract position from pre-solved solution
        predicted = raoff(sol) + my_param  # example calculation using pre-solved position
        
        # Method 2: Alternative - solve orbit on-the-fly (less efficient)
        # this_orbit = orbits[i_planet]
        # predicted = raoff(this_orbit, epoch) + my_param
        
        # Calculate likelihood contribution
        resid = predicted - measurement
        σ² = σ_measurement^2
        χ² = -(1/2) * resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    
    return ll
end
```

If the data is attached to the system as a whole, like radial velocity, the method signature should look like:
```julia
# MyLikelihood: attached to a system  
function Octofitter.ln_like(like::MyLikelihood, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    # Access your data from the table
    for i_epoch in eachindex(like.table.epoch)
        epoch = like.table.epoch[i_epoch]
        measurement = like.table.my_measurement[i_epoch]
        σ_measurement = like.table.σ_measurement[i_epoch]
        
        # Access system variables
        # θ_system.M, θ_system.plx, etc.
        
        # Access observation-specific variables
        my_param = θ_obs.my_parameter
        
        # Method 1: Use pre-solved orbit solutions for all planets (efficient!)
        predicted = zero(T)
        for planet_i in eachindex(orbits)
            # Get pre-solved solution for this planet at this epoch
            sol = orbit_solutions[planet_i][i_epoch + orbit_solutions_i_epoch_start[planet_i]]
            
            # Access planet-specific variables from θ_system.planets
            planet_keys = keys(θ_system.planets)
            planet_key = planet_keys[planet_i]
            θ_planet = θ_system.planets[planet_key]
            
            # Example: sum radial velocity contributions from all planets
            predicted += radvel(sol) * θ_planet.mass  # example calculation
        end
        predicted += my_param  # Add observation-specific offset
        
        # Method 2: Alternative - solve orbits on-the-fly (less efficient)
        # for planet_i in eachindex(orbits)
        #     orbit = orbits[planet_i]
        #     predicted += radvel(orbit, epoch) * θ_system.planets[planet_i].mass
        # end
        
        # Calculate likelihood contribution
        resid = predicted - measurement
        σ² = σ_measurement^2
        χ² = -(1/2) * resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end
    
    return ll
end
```

Inside your method, you should calculate the log-likelihood of the data stored in your likelihood object given the parameters.

The new API separates parameters into three categories:
- **`θ_system`**: System-level parameters like `M` (total mass), `plx` (parallax)
- **`θ_planet`**: Planet-specific parameters like `a` (semi-major axis), `e` (eccentricity) 
- **`θ_obs`**: Observation-specific parameters defined in the likelihood's `@variables` block

## Pre-solved Orbit Solutions (Performance Optimization)

For performance, Octofitter pre-solves orbits at all observation epochs and passes these solutions to your likelihood function:

- **`orbits`**: PlanetOrbits.jl orbit objects, one per planet
- **`orbit_solutions`**: Pre-solved orbit positions at each epoch
- **`orbit_solutions_i_epoch_start`**: Starting indices for epoch arrays (planet case) or per-planet (system case)

**Key advantages of using pre-solved solutions:**
- **Much faster**: Orbit solving is expensive, pre-solving avoids repeated calculations
- **Shared across observations**: Multiple observation objects can use the same pre-solved positions
- **Consistent**: All observations see the same orbital positions for the same parameters

**Usage patterns:**
```julia
# Planet case: Get solution for this planet at epoch i
sol = orbit_solutions[i_planet][i_epoch + orbit_solutions_i_epoch_start]

# System case: Get solution for planet j at epoch i  
sol = orbit_solutions[j][i_epoch + orbit_solutions_i_epoch_start[j]]

# Extract positions from solution
ra_offset = raoff(sol)      # RA offset in mas
dec_offset = decoff(sol)    # Dec offset in mas
radial_vel = radvel(sol)    # Radial velocity in m/s
# ... and many other functions available
```

If your likelihood requires new parameters, define them in the `variables` block of the likelihood constructor using the `@variables begin; end` syntax. These will be accessible in the `θ_obs` parameter.

If any parameter has a restricted domain where it is valid, ensure the prior is truncated using `Distributions.truncated()`. The code will automatically remap the variable using Bijectors.jl to prevent invalid values.

## Using your custom likelihood

Once you've defined your likelihood type and methods, you can use it in a model like any other likelihood:

```julia
# Create your data table
data = Table(
    (epoch=50000, my_measurement=1.2, σ_measurement=0.1),
    (epoch=50100, my_measurement=1.5, σ_measurement=0.1),
    (epoch=50200, my_measurement=1.1, σ_measurement=0.1),
)

# Create the likelihood with observation-specific variables
my_like = MyLikelihood(
    data,
    name="MY_INSTRUMENT",
    variables=@variables begin
        my_parameter ~ Normal(0, 0.5)  # Some calibration parameter
        jitter ~ LogUniform(0.01, 1.0)  # Additional uncertainty
    end
)

# Use it in a planet or system
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[my_like],  # Include your custom likelihood
    variables=@variables begin
        # Planet orbital parameters...
        a ~ Uniform(0, 100)
        e ~ Uniform(0.0, 0.5)
        # etc.
    end
)
```

## Bonus: Generative model
The above is sufficient to start sampling from the posterior. Ideally, you will also add a function that does the reverse: generate observations from a set of parameters. This is useful for a variety of statistical tests.

Simply extend the `Octofitter.generate_from_params` function for your data type:

```julia
# Generate new observations from model parameters
function Octofitter.generate_from_params(like::MyLikelihood, orbit::PlanetOrbits.AbstractOrbit, θ_planet, θ_obs)
    
    # Get epochs from original observations
    epochs = like.table.epoch
    
    # Generate new data at those epochs based on current parameters
    # i.e. "what would we observe at epoch X if the true parameters were θ_planet and θ_obs"
    simulated_measurements = []
    for epoch in epochs
        # Calculate predicted measurement using orbit and parameters
        predicted = raoff(orbit, epoch) + θ_obs.my_parameter  # example
        
        # Add noise based on original uncertainties
        σ = like.table.σ_measurement[findfirst(==(epoch), like.table.epoch)]
        noisy_measurement = predicted + σ * randn()
        push!(simulated_measurements, noisy_measurement)
    end
    
    # Create new table with simulated data
    simulated_table = Table(
        epoch=epochs,
        my_measurement=simulated_measurements,
        σ_measurement=like.table.σ_measurement  # Keep original uncertainties
    )
    
    # Return new likelihood object with simulated data
    return MyLikelihood(
        simulated_table,
        name=like.name,
        variables=(like.priors, like.derived)
    )
end
```

## Additional methods

You may also need to implement:

```julia
# For epoch subsetting (used in cross-validation)
function Octofitter.likeobj_from_epoch_subset(obs::MyLikelihood, obs_inds)
    return MyLikelihood(
        obs.table[obs_inds,:,1]; 
        name=obs.name, 
        variables=(obs.priors, obs.derived)
    )
end
```

This allows Octofitter to create subsets of your data for validation and testing purposes.
