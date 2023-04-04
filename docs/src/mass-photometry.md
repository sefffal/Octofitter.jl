#  [Connecting Mass with Photometry] (@id mass-photometry)

You can make connections between model variables using physical models with the help of [Derived variables](@ref derived).

For this example, we will assume you are sampling from image data and have a model that takes a mass as input and returns the expected flux. The flux should be in the same units as the images you included in the model.

For the sake of this example, we will assume you have images taken in two bands `H` and `J`. For demonstration purposes, we will just that the flux can be calculated as the square root of the mass. In a real project, you would likely use an interpolation over some model grid. A good way to do this is Interpolations.jl.
```julia
H_band_contrast_interp(mass) = sqrt(mass)
J_band_contrast_interp(mass) = sqrt(mass)
```

Then, list your physical variable `mass` under priors. List the model functions under Derived. This way, the H and J band flux variables will be calculated off of your models mass parameter before getting compared to the images.
```julia
@planet b Visual{KepOrbit} begin
    a ~ Normal(16, 3)
    e ~ TruncatedNormal(0.2, 0.2, 0, 1.0)
    τ ~ Normal(0.5, 1)
    ω ~ Normal(0.6, 0.2)
    i ~ Normal(0.5, 0.2)
    Ω ~ Normal(0.0, 0.2)
    mass ~ Uniform(0, 1)
    H = H_band_contrast_interp(b.mass)
    J = J_band_contrast_interp(b.mass)
end
```

If your model grids contain more independent variables, like age, surface gravity, etc. you can create a multi-dimensional interpolator. I recomment `ThinPlate()` from Interpolations.jl as a starting point.

This might look a little like this:
models mass parameter before getting compared to the images.
```julia
@planet b Visual{KepOrbit} begin
    a ~ Normal(16, 3)
    e ~ TruncatedNormal(0.2, 0.2, 0, 1.0)
    τ ~ Normal(0.5, 1)
    ω ~ Normal(0.6, 0.2)
    i ~ Normal(0.5, 0.2)
    Ω ~ Normal(0.0, 0.2)
    mass ~ Uniform(0, 1)
    temp ~ Normal(1200, 500)
    K = K_band_contrast_interp(b.mass, system.age, b.temp)
end
@system HD12345 begin
    M ~ Normal(1.0, 0.1)
    plx ~ Normal(12, 0.01)
    age ~ Normal(15, 1)
end b
```
Here the `K_band_contrast_interp` you supply accepts the mass of the planet, age of the system, and temperature of the planet as input, and returns the flux at K band in the same units as your images.
