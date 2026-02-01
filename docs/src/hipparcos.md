# [Hipparcos IAD](@id fit-hipparcos)

This tutorial explains how to model Hipparcos IAD data. The first example reproduces the catalog values of position, parallax, and proper motion. The second uses Hipparcos to constrain the mass of a directly imaged planet.

## Reproduce Catalog Values
This is the so-called "Nielsen" test from Nielsen et al (2020) and available in Orbitize!.

We start by using a system with a planet with zero mass to fit the straight line motion.

```@example 1
using Octofitter
using Distributions
using CairoMakie

hip_obs = Octofitter.HipparcosIADObs(
    hip_id=21547,
    renormalize=true, # default: true
    variables=@variables begin
        # Optional: flux ratio for luminous companions, one entry per companion
        # fluxratio ~ Product([Uniform(0, 1), Uniform(0, 1)])  # uncomment if needed for unresolved companions
    end
)

planet_b = Planet(
    name="b",
    basis=AbsoluteVisual{KepOrbit},
    variables=@variables begin
        mass = 0.
        e = 0. 
        ω = 0. 
        a = 1.
        i = 0
        Ω = 0.
        tp = 0.
    end
)

sys = System(
    name="c_Eri_straight_line",
    companions=[planet_b],
    observations=[hip_obs],
    variables=@variables begin
        M = 1.0 # Host mass not important for this example
        rv = 0.0 # system RV not significant for this example
        plx ~ Uniform(10,100)
        pmra ~ Uniform(-100, 100)
        pmdec ~  Uniform(-100, 100)

        # It is convenient to put a prior of the catalog value +- 10,000 mas on position
        ra_hip_offset_mas ~  Normal(0, 10000)
        dec_hip_offset_mas ~ Normal(0, 10000)
        dec = $hip_obs.hip_sol.dedeg + ra_hip_offset_mas/60/60/1000
        ra = $hip_obs.hip_sol.radeg + dec_hip_offset_mas/60/60/1000/cosd(dec)

        ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd
    end
)

model = Octofitter.LogDensityModel(sys)
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
    unc = hip_obs.hip_sol[prop.hip_err]
    if prop.chain == :ra
        unc /= 60*60*1000 * cosd(hip_obs.hip_sol.dedeg)
    end
    if prop.chain == :dec
        unc /= 60*60*1000
    end
    if prop.hip == :zero
        n = Normal(0, unc)
    else
        mu = hip_obs.hip_sol[prop.hip]
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
astrom_dat = Table(;
    epoch = [57009.1, 57052.1, 57053.1, 57054.3, 57266.4, 57332.2, 57374.2, 57376.2, 57415.0, 57649.4, 57652.4, 57739.1, 58068.3, 58442.2],
    sep   = [454.24, 451.81, 456.8, 461.5, 455.1, 452.88, 455.91, 455.01, 454.46, 454.81, 451.43, 449.39, 447.54, 434.22],
    σ_sep = [1.88, 2.06, 2.57, 23.9, 2.23, 5.41, 6.23, 3.03, 6.03, 2.02, 2.67, 2.15, 3.02, 2.01],
    pa    = [2.98835, 2.96723, 2.97038, 2.97404, 2.91994, 2.89934, 2.89131, 2.89184, 2.8962, 2.82394, 2.82272, 2.79357, 2.70927, 2.61171],
    σ_pa  = [0.00401426, 0.00453786, 0.00523599, 0.0523599, 0.00453786, 0.00994838, 0.00994838, 0.00750492, 0.00890118, 0.00453786, 0.00541052, 0.00471239, 0.00680678, 0.00401426]
)

astrom_obs1 = PlanetRelAstromObs(
    astrom_dat,
    name="VLT/SPHERE",
    variables=@variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)
```

We specify our full model:
```@example 1
planet_b_mass = Planet(
    name="b",
    basis=AbsoluteVisual{KepOrbit},
    observations=[astrom_obs1],
    variables=@variables begin
        a ~ truncated(Normal(10,1),lower=0.1)
        e ~ Uniform(0,0.99)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        θ ~ Uniform(0, 2pi)
        M = system.M
        tp = θ_at_epoch_to_tperi(θ, 58442.2; M, e, a, i, ω, Ω) 
        mass = system.M_sec
    end
)

sys_mass = System(
    name="cEri",
    companions=[planet_b_mass],
    observations=[hip_obs],
    variables=@variables begin
        M_pri ~ truncated(Normal(1.75,0.05), lower=0.03) # Msol
        M_sec ~ LogUniform(0.1, 100) # MJup
        M = M_pri + M_sec*Octofitter.mjup2msol # Msol

        rv =  12.60e3 # m/s
        plx ~ Uniform(20,40)
        pmra ~ Uniform(-100, 100)
        pmdec ~  Uniform(-100, 100)

        # It is convenient to put a prior of the catalog value +- 1000 mas on position
        ra_hip_offset_mas ~  Normal(0, 1000)
        dec_hip_offset_mas ~ Normal(0, 1000)
        dec = $hip_obs.hip_sol.dedeg + ra_hip_offset_mas/60/60/1000
        ra = $hip_obs.hip_sol.radeg + dec_hip_offset_mas/60/60/1000/cos(dec)

        ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd
    end
)

model = Octofitter.LogDensityModel(sys_mass)
```

Initialize the starting points, and confirm the data are entered correcly:
```@example 1
init_chain = initialize!(model, (;
    plx=34.,
    pmra=44.25,
    pmdec=-64.5,
    ra_hip_offset_mas=0.,
    dec_hip_offset_mas=0.,
))
octoplot(model, init_chain)
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