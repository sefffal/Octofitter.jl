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
    plx ~ Uniform(10,100)
    pmra ~ Uniform(-100, 100)
    pmdec ~  Uniform(-100, 100)


    # It is convenient to put a prior of the catalog value +- 10,000 mas on position
    ra_hip_offset_mas ~  Normal(0, 10000)
    dec_hip_offset_mas ~ Normal(0, 10000)
    dec = $hip_like.hip_sol.dedeg + system.ra_hip_offset_mas/60/60/1000
    ra = $hip_like.hip_sol.radeg + system.dec_hip_offset_mas/60/60/1000/cosd(system.dec)

    ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd

end hip_like b

model = Octofitter.LogDensityModel(c_Eri_straight_line)
```

Let's initialize the starting point for the chains to reasonable values
```@example 1
initialize!(model, (;
    plx=34.,
    pmra=44.25,
    pmdec=-64.5,
    ra_hip_offset_mas=0.,
    dec_hip_offset_mas=0.,
))
```

We can now sample from the model using Hamiltonian Monte Carlo. This should only take about 15 seconds.
```@example 1
using Pigeons
chain,pt = octofit_pigeons(model, n_rounds=6)
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
        unc /= 60*60*1000 * cosd(hip_like.hip_sol.dedeg)
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

We now allow the planet to have a non zero mass and have free orbit. We start by specifying relative astrometry data on the planet, collated by Jason Wang and co. on [whereistheplanet.com](http://whereistheplanet.com).

```@example 1
astrom_like1 = PlanetRelAstromLikelihood(
    (;epoch=57009.1, sep=454.24,  σ_sep=1.88, pa=2.98835, σ_pa=0.00401426),
    (;epoch=57052.1, sep=451.81,  σ_sep=2.06, pa=2.96723, σ_pa=0.00453786),
    (;epoch=57053.1, sep=456.8 ,  σ_sep=2.57, pa=2.97038, σ_pa=0.00523599),
    (;epoch=57054.3, sep=461.5 ,  σ_sep=23.9 ,pa=2.97404, σ_pa=0.0523599 ,),
    (;epoch=57266.4, sep=455.1 ,  σ_sep=2.23, pa=2.91994, σ_pa=0.00453786),
    (;epoch=57332.2, sep=452.88,  σ_sep=5.41, pa=2.89934, σ_pa=0.00994838),
    (;epoch=57374.2, sep=455.91,  σ_sep=6.23, pa=2.89131, σ_pa=0.00994838),
    (;epoch=57376.2, sep=455.01,  σ_sep=3.03, pa=2.89184, σ_pa=0.00750492),
    (;epoch=57415.0, sep=454.46,  σ_sep=6.03, pa=2.8962 , σ_pa=0.00890118),
    (;epoch=57649.4, sep=454.81,  σ_sep=2.02, pa=2.82394, σ_pa=0.00453786),
    (;epoch=57652.4, sep=451.43,  σ_sep=2.67, pa=2.82272, σ_pa=0.00541052),
    (;epoch=57739.1, sep=449.39,  σ_sep=2.15, pa=2.79357, σ_pa=0.00471239),
    (;epoch=58068.3, sep=447.54,  σ_sep=3.02, pa=2.70927, σ_pa=0.00680678),
    (;epoch=58442.2, sep=434.22,  σ_sep=2.01, pa=2.61171, σ_pa=0.00401426),
)
```

We specify our full model:
```@example 1
@planet b AbsoluteVisual{KepOrbit} begin
    a ~ truncated(Normal(10,1),lower=0.1)
    e ~ Uniform(0,0.99)
    ω ~ Uniform(0, 2pi)
    i ~ Sine()
    Ω ~ Uniform(0, 2pi)
    θ ~ Uniform(0, 2pi)
    tp = θ_at_epoch_to_tperi(system,b,58442.2) 
    mass = system.M_sec
end astrom_like1

@system cEri begin
    M_pri ~ truncated(Normal(1.75,0.05), lower=0.03) # Msol
    M_sec ~ LogUniform(0.1, 100) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol # Msol

    rv =  12.60e3 # m/s
    plx ~ Uniform(20,40)
    pmra ~ Uniform(-100, 100)
    pmdec ~  Uniform(-100, 100)

    # It is convenient to put a prior of the catalog value +- 1000 mas on position
    ra_hip_offset_mas ~  Normal(0, 1000)
    dec_hip_offset_mas ~ Normal(0, 1000)
    dec = $hip_like.hip_sol.dedeg + system.ra_hip_offset_mas/60/60/1000
    ra = $hip_like.hip_sol.radeg + system.dec_hip_offset_mas/60/60/1000/cos(system.dec)

    ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd
end hip_like b

model = Octofitter.LogDensityModel(cEri)
```


Now we sample:
```@example 1
using Pigeons
chain,pt = octofit_pigeons(model, n_rounds=8, explorer=SliceSampler())
chain
```

```@example 1
octoplot(model, chain, show_mass=true)
```

We see that we constrained both the orbit and the parallax. The mass is not strongly constrained by Hipparcos.