# [Hipparcos IAD](@id fit-hipparcos)

This tutorial explains how to model Hipparcos intermediate astrometric data (IAD) with
[`HipparcosIADObs`](@ref). The first example reproduces the catalog values of position,
parallax, and proper motion. The second uses Hipparcos to constrain the mass of a directly
imaged planet.

You can also model Hipparcos IAD jointly with Gaia and HGCA proper motion anomaly with 
[`G23HObs`](@ref).

## The frame offset block

Hipparcos IAD fits allow you to specify variables for the 5 parameter solutions. This lets you marginalize over Hipparcos IAD reference frame systematics if you're combining this data with other absolute astrometry. If you only have Hipparcos IAD absolute astrometry, you can just set these to zero or forward the system reference frame ones (like `plx`):

| Variable | Meaning | Default |
|---|---|---|
| `iad_Δra`, `iad_Δdec` | offset from the catalog position, mas | `~ Uniform(-1000, 1000)` |
| `iad_Δpmra`, `iad_Δpmdec` | offset from the catalog proper motion, mas/yr | `~ Uniform(-1000, 1000)` |
| `iad_Δplx` | offset from the system's `plx`, mas | `= 0.0` |
| `iad_pmra`, `iad_pmdec` | catalog proper motion **plus** the offset | derived |
| `hip_iad_jitter` | excess per-transit dispersion, added in quadrature | `~ LogUniform(0.001, 100)` |


## Reproduce Catalog Values

This is the so-called "Nielsen" test from Nielsen et al (2020) and shown in the Orbitize! docs. It re-fits the IAD to get back the original Hipparcos 5-parameter solution as a test.

We start by using a system with a zero-mass companion, so that the only thing the model can
do is fit the straight-line motion of the star.

```@example 1
using Octofitter
using Distributions
using CairoMakie
using Pigeons

A = Body(
    name="A",
    variables=@variables begin
        mass = 1.0   # host mass is not important for this example
    end
)


hip_obs = Octofitter.HipparcosIADObs(
    hip_id=21547,
    host=A,
    companions=(),
    ref=Barycentre,
    renormalize=true, # default: true
)

sys = System(
    name="c_Eri_straight_line",
    bodies=[A,],
    observations=[hip_obs],
    variables=@variables begin
        plx ~ Uniform(10, 100)
    end
)

model = Octofitter.LogDensityModel(sys)
```




We can now sample from the model using parallel tempering. This should only take about 15 seconds.
```@example 1
init_chain = initialize!(model)
chain, pt = octofit_pigeons(model, n_rounds=7)
chain
```

Plot the a sample from the posterior:
```@example 1
Octofitter.hipparcosplot(model, chain)
```

We now visualize the model fit compared to the Hipparcos catalog values. The parallax is
compared directly; the position and proper motion are compared through the frame offsets,
which should recover zero offset and the catalog proper motion respectively:

```@example 1
using LinearAlgebra, StatsBase

# The catalog five-parameter solution this fit should reproduce.
hip_sol = hip_obs.hip_sol

comparisons = (
    (; chain=:plx,                      μ=hip_sol.plx,    σ=hip_sol.e_plx,  label="plx [mas]"),
    (; chain=:Hipparcos_IAD_iad_Δra,    μ=0.0,            σ=hip_sol.e_ra,   label="Δα⋆ [mas]"),
    (; chain=:Hipparcos_IAD_iad_Δdec,   μ=0.0,            σ=hip_sol.e_de,   label="Δδ [mas]"),
    (; chain=:Hipparcos_IAD_iad_pmra,   μ=hip_sol.pm_ra,  σ=hip_sol.e_pmra, label="μα⋆ [mas/yr]"),
    (; chain=:Hipparcos_IAD_iad_pmdec,  μ=hip_sol.pm_de,  σ=hip_sol.e_pmde, label="μδ [mas/yr]"),
)

fig = Figure(size=(1080, 720))
ax = nothing
j = i = 1
for prop in comparisons
    global i, j, ax
    ax = Axis(fig[j, i], xlabel=prop.label)
    i += 1
    if i > 3
        j += 1
        i = 1
    end
    n = Normal(prop.μ, prop.σ)
    n0, n1 = quantile.(n, (1e-4, 1 - 1e-4))
    nxs = range(n0, n1, length=200)
    h = fit(Histogram, chain[prop.chain][:], nbins=55)
    h = normalize(h, mode=:pdf)
    barplot!(ax, (h.edges[1][1:end-1] .+ h.edges[1][2:end]) ./ 2, h.weights,
             gap=0, color=:red, label="posterior")
    lines!(ax, nxs, pdf.(n, nxs), label="Hipparcos Catalog", color=:black, linewidth=2)
end
Legend(fig[i-1, j+1], ax, tellwidth=false)
fig
```

Every panel should sit on top of the catalog curve. A 500-sample run of this model recovers
HIP 21547's published five-parameter solution to well within 1σ on all five:

| | posterior | Hipparcos catalog |
|---|---|---|
| `plx` [mas] | 33.99 ± 0.36 | 33.98 ± 0.34 |
| `iad_Δra` [mas] | −0.00 ± 0.28 | 0 ± 0.29 |
| `iad_Δdec` [mas] | 0.01 ± 0.19 | 0 ± 0.19 |
| `iad_pmra` [mas/yr] | 44.21 ± 0.33 | 44.22 ± 0.34 |
| `iad_pmdec` [mas/yr] | −64.39 ± 0.26 | −64.39 ± 0.27 |

!!! note 
    `HipparcosIADObs` defaults to `recalibrate=false`, so it reproduces the
    published catalog data as-is. `G23HObs` applies the Brandt, Michalik & Brandt shift
    (+0.140 mas on the residuals, 2.25 mas of extra dispersion) unconditionally; pass
    `recalibrate=true` here to match its `:iad_hip` channel exactly.

## Constrain Planet Mass

We now allow the planet to have a non-zero mass and a free orbit. We start by specifying
relative astrometry data on the planet, collated by Jason Wang and co. on
[whereistheplanet.com](http://whereistheplanet.com).

```@example 1
astrom_dat = Table(;
    epoch = [57009.1, 57052.1, 57053.1, 57054.3, 57266.4, 57332.2, 57374.2, 57376.2, 57415.0, 57649.4, 57652.4, 57739.1, 58068.3, 58442.2],
    sep   = [454.24, 451.81, 456.8, 461.5, 455.1, 452.88, 455.91, 455.01, 454.46, 454.81, 451.43, 449.39, 447.54, 434.22],
    σ_sep = [1.88, 2.06, 2.57, 23.9, 2.23, 5.41, 6.23, 3.03, 6.03, 2.02, 2.67, 2.15, 3.02, 2.01],
    pa    = [2.98835, 2.96723, 2.97038, 2.97404, 2.91994, 2.89934, 2.89131, 2.89184, 2.8962, 2.82394, 2.82272, 2.79357, 2.70927, 2.61171],
    σ_pa  = [0.00401426, 0.00453786, 0.00523599, 0.0523599, 0.00453786, 0.00994838, 0.00994838, 0.00750492, 0.00890118, 0.00453786, 0.00541052, 0.00471239, 0.00680678, 0.00401426]
)

astrom_obs = RelAstromObs(
    astrom_dat,
    target = :b,      # the body the measurement is *of*
    ref    = :A,      # …and the body it is measured *against*
    name = "VLT/SPHERE",
    variables = @variables begin
        # Fixed values for this example - could be free variables:
        jitter = 0        # mas [could use: jitter ~ Uniform(0, 10)]
        northangle = 0    # radians [could use: northangle ~ Normal(0, deg2rad(1))]
        platescale = 1    # relative [could use: platescale ~ truncated(Normal(1, 0.01), lower=0)]
    end
)
nothing # hide
```

We specify our full model.

```@example 1
A_mass = Body(
    name="A",
    variables=@variables begin
        mass ~ truncated(Normal(1.75, 0.05), lower=0.03)   # Msol
    end
)

planet_b_mass = Body(
    name="b",
    about=A_mass,
    variables=@variables begin
        mass ~ LogUniform(0.1mjup, 100mjup)  # Msol (mjup is a constant, not a unit system)
        a ~ truncated(Normal(10, 1), lower=0.1)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        θ ~ Uniform(0, 2pi)
        epoch = 58442.2
    end
)

hip_obs_mass = Octofitter.HipparcosIADObs(
    hip_id=21547,
    host=A_mass,
    companions=(planet_b_mass,),
    ref=Barycentre,
)

sys_mass = System(
    name="cEri",
    bodies=[A_mass, planet_b_mass],
    observations=[hip_obs_mass, astrom_obs],
    variables=@variables begin
        plx ~ Uniform(20, 40)
    end
)

model = Octofitter.LogDensityModel(sys_mass)
```

Initialize the starting points, and confirm the data are entered correctly:
```@example 1
init_chain = initialize!(model, (;
    plx = 34.0,
    observations = (;
        Hipparcos_IAD = (;
            iad_Δra = 0.0,
            iad_Δdec = 0.0,
            iad_Δpmra = 0.0,
            iad_Δpmdec = 0.0,
        ),
    ),
))
```

Now we sample:
```@example 1
chain, pt = octofit_pigeons(model, n_rounds=8, explorer=SliceSampler())
chain
```

```@example 1
octoplot(model, chain)
```

We see that we constrained both the orbit and the parallax. The mass is not strongly
constrained by Hipparcos alone; `chain["b_mass"]` is in solar masses, so divide by `mjup`
to read it in Jupiter masses:

```@example 1
using StatsBase
mass_mjup = chain["b_mass"][:] ./ mjup
println("b mass [Mjup]: ", round.(quantile(mass_mjup, (0.16, 0.5, 0.84)), digits=2))
```

## See Also

- [Joint Gaia-Hipparcos (G23H)](@ref fit-g23h) — the same abscissae as one channel of a
  joint Gaia + Hipparcos fit.
- [Proper Motion Anomaly](@ref fit-pma) — HGCA-style proper motion anomaly.
