# Hipparcos Modelling

This tutorial explains how to model Hipparcos IAD data. The first example reproduces the catalog values of position, parallax, and proper motion. The second uses Hipparcos to constrain the mass of a directly imaged planet.

## Reproduce Catalog Values
This is the so-called "Nielsen" test from Nielsen et al (2020) and available in Orbitize!.

We start by using a system with a planet with zero mass to fit the straight line motion.

```@example 1
using Octofitter
using Distributions
using CairoMakie

hip_like = Octofitter.HipparcosIADLikelihood(;
    hip_id=21547,
    renormalize=true, # default: true
)

@planet b AbsoluteVisual{KepOrbit} begin
    mass = 0.
    e = 0. 
    ω = 0. 
    a = 1.
    i = 0
    Ω = 0.
    tp = 0.
end
@system c_Eri_straight_line begin
    M = 1.0 # Host mass not important for this example
    rv = 0.0 # system RV not significant for this example
    plx ~ Uniform(0,100)
    pmra ~ Uniform(-100, 100)
    pmdec ~  Uniform(-100, 100)


    # It is convenient to put a prior of the catalog value +- 10,000 mas on position
    ra_hip_offset_mas ~  Normal(0, 10000)
    dec_hip_offset_mas ~ Normal(0, 10000)
    dec = hip_like.hip_sol.dedeg + system.ra_hip_offset_mas/60/60/1000
    ra = hip_like.hip_sol.radeg + system.dec_hip_offset_mas/60/60/1000/cos(system.dec)

    ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd

end hip_like b

model = Octofitter.LogDensityModel(c_Eri_straight_line)
```

We can now sample from the model using Hamiltonian Monte Carlo. This should only take about 15 seconds.
```@example 1
# Typical depth is ~4, so limiting the max_depth down from default 12 speeds warm-up
chain = octofit(model, iterations=4_000, max_depth=6)
nothing # hide
```

Plot the posterior values:
```@example 1
octoplot(model,chain,show_astrom=false,show_astrom_time=false)
```


We now visualize the model fit compared to the Hipparcos catalog values:
```@example 1
using LinearAlgebra, StatsBase
fig = Figure(size=(1080,720))
j = i = 1
for prop in (
    (;chain=:ra, hip=:radeg, hip_err=:e_ra), 
    (;chain=:dec, hip=:dedeg, hip_err=:e_de),
    (;chain=:plx, hip=:plx, hip_err=:e_plx), 
    (;chain=:pmra, hip=:pm_ra, hip_err=:e_pmra), 
    (;chain=:pmdec, hip=:pm_de, hip_err=:e_pmde)
)
    global i, j, ax
    ax = Axis(
        fig[j,i],
        xlabel=string(prop.chain),
    )
    i+=1
    if i > 3
        j+=1
        i = 1
    end
    unc = hip_like.hip_sol[prop.hip_err]
    if prop.chain == :ra
        unc /= 60*60*1000 * cos(hip_like.hip_sol.dedeg)
    end
    if prop.chain == :dec
        unc /= 60*60*1000
    end
    if prop.hip == :zero
        n = Normal(0, unc)
    else
        mu = hip_like.hip_sol[prop.hip]
        n = Normal(mu, unc)
    end
    n0,n1=quantile.(n,(1e-4, 1-1e-4))
    nxs = range(n0,n1,length=200)
    h = fit(Histogram, chain[prop.chain][:], nbins=55)
    h = normalize(h, mode=:pdf)
    barplot!(ax, (h.edges[1][1:end-1] .+ h.edges[1][2:end])./2, h.weights, gap=0, color=:red, label="posterior")
    lines!(ax, nxs, pdf.(n,nxs), label="Hipparcos Catalog", color=:black, linewidth=2)
end
Legend(fig[i-1,j+1],ax,tellwidth=false)
fig
```


## Constrain Planet Mass

We now allow the planet to have a non zero mass and free orbit. We start by retrieving relative astrometry data on the planet, collated by Jason Wang and co on [whereistheplanet.com](whereistheplanet.com).

```@example 1
astrom_like1,astrom_like2 = Octofitter.Whereistheplanet_astrom("51erib",object=1)
nothing # hide
```

We specify our full model:
```@example 1
@planet b AbsoluteVisual{KepOrbit} begin
    a ~ truncated(Normal(10,1),lower=0)
    e ~ Uniform(0,0.99)
    ω ~ UniformCircular()
    i ~ Sine()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,57463) 
    mass = system.M_sec
end astrom_like1 

@system cEri begin
    M_pri ~ truncated(Normal(1.75,0.05), lower=0.03) # Msol
    M_sec ~ Uniform(0, 100) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol

    rv =  12.60e3 # m/s
    plx ~ Uniform(0,100)
    pmra ~ Uniform(-100, 100)
    pmdec ~  Uniform(-100, 100)

    # It is convenient to put a prior of the catalog value +- 1000 mas on position
    ra_hip_offset_mas ~  Normal(0, 10000)
    dec_hip_offset_mas ~ Normal(0, 10000)
    dec = hip_like.hip_sol.dedeg + system.ra_hip_offset_mas/60/60/1000
    ra = hip_like.hip_sol.radeg + system.dec_hip_offset_mas/60/60/1000/cos(system.dec)

    ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd
end hip_like b

model = Octofitter.LogDensityModel(cEri)
```


Now we sample:
```@example 1
chain = octofit(model)
nothing # hide
```

```@example 1
octoplot(model, chain, show_mass=true)
```

We see that we constrained both the orbit and the parallax. The mass is not strongly constrained by Hipparcos.