# API Documentation

<!--
Coverage list only — the prose belongs on the topic pages. Every name below
resolves and carries a docstring against the v2 surface (checked during the
port integration pass). Names retired with v1 (`Planet`, `PlanetRelAstromObs`,
`StarAbsoluteRVObs`, `PlanetRelativeRVObs`, `Visual{KepOrbit}`,
`AbsoluteVisual{KepOrbit}`, `RadialVelocityOrbit`, the `*_nss` helpers) were
removed rather than aliased; see `v2-migration.md`. Orbit and system *element*
types now live in PlanetOrbits and are documented there.
-->

## The model

```@docs
@variables
System
Body
Orbit
UniformCircular
Sine
KDEDist
```

`Barycentre` and `Photocentre` are the reference grammar's two derived points;
they are singletons and are described in [`System`](@ref)'s docstring.

## Observations

```@docs
RelAstromObs
RadialVelocityObs
MarginalizedRVObs
PhotometryObs
GaiaDR4AstromObs
G23HObs
HGCAObs
HipparcosIADObs
```

These four live in the two unregistered subpackages, which `make.jl` loads so
that their docstrings resolve here:

```@docs
ImageObs
LogLikelihoodMapObs
InterferometryObs
GRAVITYWideKPObs
```

## Prior terms

```@docs
ObsPriorONeil2019
OrbitOrderPrior
NonCrossingPrior
LimitClosestApproachAUPrior
HillStabilityPrior
```

## Fitting

```@docs
octofit
octofit_rejection
initialize!
startingpoints!
Octofitter.advancedhmc
octofit_pigeons
```

`octofit_pigeons` needs `using Pigeons`: its methods live in a package
extension, so without that import the function exists but has no methods.

## Analysis

```@docs
construct_system
drawfrompriors
generate_from_params
prior_only_model
Octofitter.pointwise_like
Octofitter.calibrationhmc
Octofitter.sbctrial
completeness_jobs
run_completeness_trial
assemble_completeness
completeness_map
CompletenessJob
CompletenessResult
CompletenessMap
```

## Plotting

```@docs
octoplot
octocorner
rvpostplot
rvpostplot_animated
dotplot
gaiastarplot
gaiatimeplot
skytrackplot
hipparcosplot
completenessplot
completenessplot!
PlotChannel
plotchannels
defaultpanels
ObservableQuery
PosteriorSeries
OctoPlotResult
obscontext
modelcurves
mapcurve
timeseriespanel!
skypanel!
phasefoldpanel!
photometrypanel!
likemappanel!
Octofitter.residuals
```

`orbitlines!`, `plot_epochs`, `orbit_track_epochs`, `orbit_theme`,
`add_mjd_axis!`, `MJDConversion` and `paraminfo` are re-exported from
PlanetOrbits and documented there.

## Loading and saving

```@docs
Octofitter.savechain
Octofitter.loadchain
Octofitter.checkchain
Octofitter.savehdf5
Octofitter.loadhdf5
Octofitter.Whereistheplanet_astrom
```

## Substellar models and catalog helpers

```@docs
sonora_photometry_interpolator
sonora_cooling_interpolator
Octofitter.bhac15_mass_age_interpolator
gaia_plx
gaia_dr3_solution
mjd
years2mjd
mjd2date
```

## Extending Octofitter

The observation interface a custom likelihood implements, and the hooks a
Gaussian-process backend plugs into.

```@docs
ObsContext
refspecs
epochs
solutionat
resolverefs
likelihoodname
Octofitter.likeobj_from_epoch_subset
sky_offset
sky_offset!
sky_calibration
FrameOffset
frame_offset
frame_offset_alongscan
Octofitter.gp_condition
Octofitter.gp_ln_like
Octofitter.gp_predict
```

## Retired names

**Aliases.** These resolve to their v2 replacements and exist only so v1
scripts keep parsing: `PhotometryLikelihood` → [`PhotometryObs`](@ref),
`PlanetOrderPrior` → [`OrbitOrderPrior`](@ref),
`ObsPriorAstromONeil2019` → [`ObsPriorONeil2019`](@ref),
`InterferometryLikelihood` → [`InterferometryObs`](@ref),
`GRAVITYWideKPLikelihood` → [`GRAVITYWideKPObs`](@ref).

**Error stubs.** These are defined, exported, and raise an error naming the
replacement, because in each case the v2 spelling takes different arguments
or lives somewhere else in the model — an alias would not have helped, and a
bare `UndefVarError` says nothing:

| v1 name | v2 |
|---|---|
| `Planet` | [`Body`](@ref) + the observation moves to [`System`](@ref)'s `observations=` |
| `θ_at_epoch_to_tperi` | declare `θ` and `epoch` as orbital elements |
| `PlanetRelAstromObs` | [`RelAstromObs`](@ref)`(tab; target, ref)` |
| `StarAbsoluteRVObs` | [`RadialVelocityObs`](@ref)`(tab; target=A, ref=Barycentre)` |
| `PlanetRelativeRVObs` | [`RadialVelocityObs`](@ref)`(tab; target=b, ref=A)` |
| `MarginalizedStarAbsoluteRVObs` | [`MarginalizedRVObs`](@ref)`(tab; target, ref)` |
| `masspostplot` | [`octocorner`](@ref), or `hist(vec(chain[:b_mass]))` |
| `HGCAInstantaneousObs`, `GaiaCatalogFitObs` | [`HGCAObs`](@ref) / [`G23HObs`](@ref) |

The HGCA pair is the only case where the *modelling code* was not ported at
all — it was subsumed by [`G23HObs`](@ref). The rest are pure renames whose
call signature changed.

Orbit element types, observables (`raoff`, `decoff`, `radvel`, `posx`,
`projectedseparation`, …), `orbitsolve`, `Trajectory`, `WeightedPoint`,
`barycentre`, `photocentre` and `fluxes` come from PlanetOrbits, which
Octofitter re-exports; they are documented in the PlanetOrbits manual.
