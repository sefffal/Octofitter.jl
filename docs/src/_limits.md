## Detection Limits

This guide shows how to calculate detection limits, in mass, or in photometry, as a function of orbital parameters.

There are a few use cases for this:

* Mass limit vs semi-major axis given one or more images and/or contrast curves
* Mass limit vs semi-major axis given an RV non-detection
* Mass limit vs semi-major axis given proper motion anomaly from the GAIA-Hipparcos Catalog of Accelerations
* Any combination of the above


We will once more use some sample data from the system [HD 91312 A & B](https://arxiv.org/abs/2109.12124) discovered by SCExAO. 


```@example 1
using Octofitter
using OctofitterImages
using OctofitterRadialVelocity
using Distributions
using Pigeons
using CairoMakie
```


```@example 1
const cooling_tracks = sonora_cooling_interpolator()
const sonora_temp_mass_L = sonora_photometry_interpolator(:Keck_L′)
```

```@example 1
@planet B Visual{KepOrbit} begin
    a ~ LogUniform(1, 65)
    e ~ Uniform(0,0.9)
    # e = 0
    ω ~ Uniform(0,2pi)#~ UniformCircular()
    i ~ Sine() # The Sine() distribution is defined by Octofitter
    Ω ~ Uniform(0,pi)# ~ UniformCircular()
    mass = system.M_sec

    # Calculate planet temperature from cooling track and planet mass variable
    tempK = cooling_tracks(system.age, B.mass)
    # Calculate absolute magnitude
    abs_mag_L = sonora_temp_mass_L(B.tempK, B.mass)
    # Deal with out-of-grid values by clamping to grid max and min
    abs_mal_L′ = if isfinite(B.abs_mag_L)
        B.abs_mag_L
    elseif B.mass > 10 
        8.2 # jump to absurdly bright
    else
        16.7 # jump to absordly dim
    end
    # Calculate relative magnitude
    rel_mag_L = B.abs_mal_L′ - system.rel_mag + 5log10(1000/system.plx)
    # Convert to contrast (same units as image)
    L = 10.0^(B.rel_mag_L/-2.5)

    θ ~ Uniform(0,2pi)#UniformCircular()
    tp = θ_at_epoch_to_tperi(system,B,57423.0) # epoch of GAIA measurement
end
@system HD91312_pma begin
    age = 10
    rel_mag = 5.65
    M_pri ~ truncated(Normal(0.95, 0.05), lower=0) # Msol
    M_sec ~ LogUniform(0.2, 65) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol

    plx ~ gaia_plx(gaia_id=6166183842771027328)
            
    # Priors on the center of mass proper motion
    pmra ~ Normal(0, 1000)
    pmdec ~ Normal(0,  1000)
end HGCALikelihood(gaia_id=6166183842771027328) B
model_pma = Octofitter.LogDensityModel(HD91312_pma)
```


```@example 1
chain_pma, pt = octofit_pigeons(model_pma, n_chains=8, n_rounds=14);
```


```@example 1
fig = Figure()
ax = Axis(
    fig[1,1],
    xscale=log10,
    yscale=log10,
    xticks = 2.0 .^ (0:2:6),
    yticks = 2.0 .^ (-1:2:6)
)
x = chain_pma[:B_a][:]
y = chain_pma[:B_mass][:]
scatter!(ax,x,y, markersize=5, color=Makie.wong_colors()[1])
fig
```



```@example 1
using AstroImages
download(
    "https://github.com/sefffal/Octofitter.jl/raw/main/docs/image-examples-1.fits",
    "image-examples-1.fits"
)

# Or multi-extension FITS (this example)
image = AstroImages.load("image-examples-1.fits").*5e-7 # units of contrast

image_data = ImageLikelihood(
    (band=:L, image=AstroImages.recenter(image), platescale=10.0, epoch=57423.6),
)
```



```@example 1
@planet B Visual{KepOrbit} begin
    a ~ LogUniform(1, 65)
    e ~ Uniform(0,0.9)
    # e = 0
    ω ~ Uniform(0,2pi)#~ UniformCircular()
    i ~ Sine() # The Sine() distribution is defined by Octofitter
    Ω ~ Uniform(0,pi)# ~ UniformCircular()
    mass = system.M_sec

    # Calculate planet temperature from cooling track and planet mass variable
    tempK = cooling_tracks(system.age, B.mass)
    # Calculate absolute magnitude
    abs_mag_L = sonora_temp_mass_L(B.tempK, B.mass)
    # Deal with out-of-grid values by clamping to grid max and min
    abs_mal_L′ = if isfinite(B.abs_mag_L)
        B.abs_mag_L
    elseif B.mass > 10 
        8.2 # jump to absurdly bright
    else
        16.7 # jump to absordly dim
    end
    # Calculate relative magnitude
    rel_mag_L = B.abs_mal_L′ - system.rel_mag + 5log10(1000/system.plx)
    # Convert to contrast (same units as image)
    L = 10.0^(B.rel_mag_L/-2.5)

    θ ~ Uniform(0,2pi)#UniformCircular()
    tp = θ_at_epoch_to_tperi(system,B,57423.0+1000) # epoch of GAIA measurement
end image_data 

@system HD91312_img begin
    # age ~ truncated(Normal(40, 15),lower=0, upper=200)
    age = 10
    M_pri ~ truncated(Normal(0.95, 0.05), lower=0) # Msol
    # Mass of secondary
    # Make sure to pick only a mass range that is covered by your models
    M_sec ~ LogUniform(0.55, 65) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol
    plx ~ gaia_plx(gaia_id=6166183842771027328)
    # Priors on the center of mass proper motion
    # pmra ~ Normal(0, 1000)
    # pmdec ~ Normal(0,  1000)
    rel_mag = 5.65
end  B
model_img = Octofitter.LogDensityModel(HD91312_img)
```

```@example 1
using Random; Random.seed!(1)

chain_img, pt = octofit_pigeons(model_img, n_rounds=14) ;
```




```julia

@planet B Visual{KepOrbit} begin
    a ~ LogUniform(1, 65)
    e ~ Uniform(0,0.9)
    # e = 0
    ω ~ Uniform(0,2pi)#~ UniformCircular()
    i ~ Sine() # The Sine() distribution is defined by Octofitter
    Ω ~ Uniform(0,pi)# ~ UniformCircular()
    mass = system.M_sec

    # Calculate planet temperature from cooling track and planet mass variable
    tempK = cooling_tracks(system.age, B.mass)
    # Calculate absolute magnitude
    abs_mag_L = sonora_temp_mass_L(B.tempK, B.mass)
    # Deal with out-of-grid values by clamping to grid max and min
    abs_mal_L′ = if isfinite(B.abs_mag_L)
        B.abs_mag_L
    elseif B.mass > 10 
        8.2 # jump to absurdly bright
    else
        16.7 # jump to absordly dim
    end
    # Calculate relative magnitude
    rel_mag_L = B.abs_mal_L′ - system.rel_mag + 5log10(1000/system.plx)
    # Convert to contrast (same units as image)
    L = 10.0^(B.rel_mag_L/-2.5)

    θ ~ Uniform(0,2pi)#UniformCircular()
    tp = θ_at_epoch_to_tperi(system,B,57423.0+1000) # epoch of GAIA measurement
end image_data 

@system HD91312_both begin
    # age ~ truncated(Normal(40, 15),lower=0, upper=200)
    age = 10
    M_pri ~ truncated(Normal(0.95, 0.05), lower=0) # Msol
    # Mass of secondary
    # Make sure to pick only a mass range that is covered by your models
    M_sec ~ LogUniform(0.55, 65) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol
    plx ~ gaia_plx(gaia_id=6166183842771027328)
    # Priors on the center of mass proper motion
    pmra ~ Normal(0, 1000)
    pmdec ~ Normal(0,  1000)
    rel_mag = 5.65
end  HGCALikelihood(gaia_id=6166183842771027328) B
model_both = Octofitter.LogDensityModel(HD91312_both)
using Random; Random.seed!(1)

chain_both, pt = octofit_pigeons(model_both, n_rounds=14) ;

```


```@example
fig = Figure()
ax = Axis(
    fig[1,1],
    xscale=Makie.pseudolog10,
    yscale=Makie.pseudolog10,
    xticks = 2.0 .^ (0:2:6),
    yticks = 2.0 .^ (-1:2:6),
    xlabel = "semi-major axis [au]",
    ylabel = Makie.rich("companion mass [m", Makie.subscript("jup"), "]"),
)



x = chain_img[:B_a][:]
y = chain_img[:B_mass][:]
scatter!(ax,x,y, markersize=3, color=Makie.wong_colors()[3], label="image", alpha=0.5)

x = chain_pma[:B_a][:]
y = chain_pma[:B_mass][:]
scatter!(ax,x,y, markersize=3, color=Makie.wong_colors()[1], label="PMA", alpha=0.5)


x = chain_both[:B_a][:]
y = chain_both[:B_mass][:]
scatter!(ax,x,y, markersize=3, color=Makie.wong_colors()[2], label="image+PMA", alpha=0.5)

Legend(fig[1,2], ax)

fig
```

We should also incorporate some kind of observability cut. Only consider planets that actually detectable ever.

```@example
ts = range(
    start=(57423+1000),
    stop =(57423+1000 + 15*365),
    length=100,
)

fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel="years after 1ˢᵗ image",
    ylabel="detection probability",
)
# Look at percentage of orbit separation + photometries that fall above contrast curve.
for (chain, name) in ((chain_pma, "PMA only"), (chain_img, "image only"), (chain_both, "PMA + image"))
    els  = Octofitter.construct_elements(chain, :B, :)
    sep = projectedseparation.(els, ts')

    platescale = 10.0
    cont = OctofitterImages.contrast_interp(AstroImages.recenter(image))
    mask = chain[:B_L] .> 3 .* cont.(sep./platescale)

    lines!(
        ax,
        (ts .- image_data.table.epoch[1])/365.25,
        mean(mask, dims=1)[:],
        label=name
    )
end
fig
# c=scatterlines(
#     ts,
#     mean(sep,dims=1)[:]
# )
# scatterlines!(
#     ts,
#     mean(sep,dims=1)[:]-std(sep,dims=1)[:]
# )
# scatterlines!(
#     ts,
#     mean(sep,dims=1)[:]+std(sep,dims=1)[:]
# )
```
