# Resolved Gaia DR3 Companion

This page describes [`GaiaDR3CompanionObs`](@ref): a likelihood for a **resolved,
gravitationally bound stellar companion** that has its own published Gaia DR3
five-parameter astrometric solution (position, parallax, proper motion, and
optionally radial velocity).

The typical use cases are:

1. **Weighing a dark central object.** You have a visible star with a published
   Gaia DR3 entry that is bound to an unseen (dark) massive object. Folding in
   the companion's Gaia solution, together with an assumed companion mass,
   constrains the **mass and location of the central object** — even when it
   emits no light.
2. **Adding a known companion's astrometry/RV** to a wider model, accounting for
   the acceleration/parallax/proper-motion it induces.

## Physical picture

The model treats the system as a host star (the possibly-dark central object,
whose mass is `system.M` and whose barycentric astrometry is the system-level
`ra, dec, plx, pmra, pmdec`) orbited by a `Planet` that represents the resolved,
bound companion. Because the companion is bound, its measured Gaia astrometry is
the **sum of the barycentre's space motion and the companion's reflex orbit**
about the barycentre. The likelihood predicts the companion's absolute
five-parameter solution at the Gaia DR3 reference epoch (J2016.0) and compares it
to the catalogue values using their full 5×5 covariance.

!!! note "AbsoluteVisual orbits are required"
    The companion `Planet` must use an `AbsoluteVisual{...}` orbit basis. This
    engages the rigorous 3D non-linear stellar-motion propagation (perspective
    acceleration, changing parallax, light-travel time) in
    [PlanetOrbits.jl](https://github.com/sefffal/PlanetOrbits.jl). A plain
    `Visual{...}` orbit will raise an error.

!!! note "How system-level variables reach the companion"
    In Octofitter each planet's orbit is built by merging the system and planet
    variables, so `ra, dec, plx, pmra, pmdec, ref_epoch, rv` are baked into the
    companion's `AbsoluteVisual` orbit automatically. The observation is simply
    attached to the `Planet` representing the companion.

## Specifying the data

You can either query the Gaia archive by DR3 `source_id` (the default,
convenient path) or supply the five parameters yourself.

By source id (catalogue values and full covariance fetched and cached):
```julia
GaiaDR3CompanionObs(gaia_id_dr3 = 1234567890123456789)
```

Or specifying the solution manually (useful offline, or to override values):
```julia
GaiaDR3CompanionObs(
    ra = 150.0, dec = -30.0,          # degrees
    parallax = 50.0,                  # mas
    pmra = 80.0, pmdec = -40.0,       # mas/yr
    ra_error = 0.04, dec_error = 0.04,        # mas (σ of α*=α·cosδ and δ)
    parallax_error = 0.03,                    # mas
    pmra_error = 0.04, pmdec_error = 0.04,    # mas/yr
    # correlation coefficients optional, default 0:
    corr = (; ra_dec_corr = 0.1, parallax_pmra_corr = -0.05),
)
```

### Error inflation

Published Gaia uncertainties are multiplied by `inflation_factor` (default
`1.37`, following Brandt 2019) before forming the likelihood. Set
`inflation_factor = 1.0` to use the catalogue uncertainties unchanged.

### Optional radial velocity

If `include_rv = true`, the companion's predicted barycentric radial velocity is
additionally constrained against the Gaia DR3 `radial_velocity` (km/s). The RV is
taken from the catalogue query, or you can supply `radial_velocity` and
`radial_velocity_error` directly.

## Example: weighing a dark central object

```julia
using Octofitter, Distributions

companion = Planet(
    name = "B",
    basis = AbsoluteVisual{KepOrbit},
    observations = [
        GaiaDR3CompanionObs(gaia_id_dr3 = 1234567890123456789, inflation_factor = 1.37),
    ],
    variables = @variables begin
        mass = system.M_companion                     # assumed companion mass [Mjup]
        a ~ truncated(Normal(20, 0.3), lower=1)       # from an existing orbit fit
        e ~ truncated(Normal(0.3, 0.02), lower=0, upper=0.99)
        i ~ truncated(Normal(0.6, 0.02), lower=0, upper=pi)
        ω ~ Normal(0.4, 0.02)
        Ω ~ Normal(1.0, 0.02)
        M = system.M
        tp = 50000.0
    end
)

sys = System(
    name = "darkmass",
    companions = [companion],
    observations = [],
    variables = @variables begin
        M_central ~ LogUniform(0.3, 5.0)              # the dark central mass [Msol]
        M_companion ~ truncated(Normal(314, 30), lower=0)  # assumed companion mass [Mjup]
        M = M_central + M_companion*Octofitter.mjup2msol
        # The barycentric (systemic) motion — constrain it from whatever you know
        # about the central object (here, informatively):
        plx  ~ Normal(50.0, 0.05)
        pmra ~ Normal(80.0, 0.2)
        pmdec ~ Normal(-40.0, 0.2)
        rv = 10_000.0
        ra = 150.0
        dec = -30.0
        ref_epoch = Octofitter.meta_gaia_DR3.ref_epoch_mjd   # J2016.0
    end
)

model = Octofitter.LogDensityModel(sys)
initialize!(model)
chain = octofit(model)
```

!!! tip "Identifiability"
    A single companion's five-parameter solution does not by itself break the
    degeneracy between the barycentric motion and the companion's reflex. To
    constrain the central mass you need an additional anchor on the barycentre's
    motion and/or the companion's orbital geometry — e.g. an informative prior on
    the systemic proper motion (as above), a relative-astrometry/RV orbit, or the
    central object's own astrometry. With those in hand, the Gaia companion
    solution pins the central mass.

## Modelling notes

- The predicted proper motion is the instantaneous value at the Gaia DR3
  reference epoch, evaluated by a central finite difference of the companion's
  absolute sky position (step `fd_step_days`, default 180 days). This is accurate
  for companions whose orbital motion is approximately linear over the ~34-month
  Gaia DR3 baseline (i.e. wide / long-period companions — the regime this
  likelihood targets). Very short-period companions are smeared within the Gaia
  baseline and would typically be flagged as non-single-star solutions instead.
- An optional `astrometric_jitter` (mas) observation variable is added in
  quadrature to the position uncertainties if present in the `@variables` block.

## API

```@docs
GaiaDR3CompanionObs
```
