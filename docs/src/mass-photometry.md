#  [Connecting Mass with Photometry] (@id mass-photometry)

You can make connections between model variables using physical models with the help of [Derived variables](@ref derived).

For this example, we will assume you are sampling from photometry data and have a model that takes a mass as input and returns the expected flux. The flux should be in the same units as the photometry measurements you included in the model.

For the sake of this example, we will assume you have photometry taken in two bands `H` and `J`. For demonstration purposes, we will just pretend that the flux can be calculated as the square root of the mass. In a real project, you would likely use an interpolation over some model grid. A good way to do this is Interpolations.jl.
```julia
H_band_contrast_interp(mass) = sqrt(mass) # your model or function here
J_band_contrast_interp(mass) = sqrt(mass) # your model or function here
```

First, create your photometry observations and likelihoods. Each PhotometryObs handles a single band with a `flux` variable defined in its variables block rather than in the planet definition. The `name` parameter is used for variable naming in the MCMC chain output (e.g., "b_H_band_flux"). This way, the flux variables will be calculated off of your model's mass parameter before getting compared to the photometry:
```julia
# Create separate photometry likelihoods for each band
H_band_table = Table(
    phot=[15.2], σ_phot=[0.5],
)
H_band_data = PhotometryObs(
    H_band_table,
    name="H_band",
    variables=@variables begin
        flux = $H_band_contrast_interp(system.mass)
    end
)

J_band_table = Table(
    phot=[14.8], σ_phot=[0.3],
)
J_band_data = PhotometryObs(
    J_band_table,
    name="J_band", 
    variables=@variables begin
        flux = $J_band_contrast_interp(system.mass)
    end
)

# Define the planet with orbital parameters and mass
planet_b = Planet(
    name="b",
    basis=Visual{KepOrbit},
    observations=[H_band_data, J_band_data],
    variables=@variables begin
        a ~ Normal(16, 3)
        e ~ truncated(Normal(0.2, 0.2), lower=0, upper=0.99)
        ω ~ Normal(0.6, 0.2)
        i ~ Normal(0.5, 0.2)
        Ω ~ Normal(0.0, 0.2)
        mass ~ Uniform(0, 1)
        
        M = system.M
        τ ~ UniformCircular(1.0)
        P = √(a^3/M)
        tp = τ*P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
    end
)
```

If your model grids contain more independent variables, like age, surface gravity, etc. you can create a multi-dimensional interpolator. I recommend `ThinPlate()` from Interpolations.jl as a starting point.

This might look a little like this:
```julia
# More complex interpolation function using multiple variables
K_band_contrast_interp(mass, age, temp) = sqrt(mass) * (age/10) * sqrt(temp/1000)

# Photometry likelihood with derived variable using system and planet properties
K_band_table = Table(
    phot=[13.5], σ_phot=[0.4],
)
K_band_data = PhotometryObs(
    K_band_table,
    name="K_band",
    variables=@variables begin
        flux = $K_band_contrast_interp(system.mass, system.age, system.temp)
    end
)

# Planet definition
planet_b = Planet(
    name="b", 
    basis=Visual{KepOrbit},
    observations=[K_band_data],
    variables=@variables begin
        a ~ Normal(16, 3)
        e ~ truncated(Normal(0.2, 0.2), lower=0, upper=0.99)
        ω ~ Normal(0.6, 0.2)
        i ~ Normal(0.5, 0.2)
        Ω ~ Normal(0.0, 0.2)
        mass ~ Uniform(0, 1)
        temp ~ Normal(1200, 500)
        age = system.age
        
        M = system.M
        τ ~ UniformCircular(1.0)
        P = √(a^3/M)
        tp = τ*P*365.25 + 58849 # reference epoch for τ. Choose an MJD date near your data.
    end
)

# System definition
sys = System(
    name="HD12345",
    companions=[planet_b],
    observations=[],
    variables=@variables begin
        M ~ Normal(1.0, 0.1)
        plx ~ Normal(12, 0.01)
        age ~ Normal(15, 1)
    end
)
```
Here the `K_band_contrast_interp` you supply accepts the mass of the planet, age of the system, and temperature of the planet as input, and returns the flux at K band in the same units as your photometry measurements.
