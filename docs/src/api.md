# API Documentation

<!--
Coverage list only — the prose belongs on the topic pages. Orbit and system
*element* types live in PlanetOrbits and are documented there.
-->

## The model

```@docs
@variables
System
Body
Orbit
Barycentre
Photocentre
UniformCircular
Sine
KDEDist
```

[`Barycentre`](@ref) and [`Photocentre`](@ref) are the reference grammar's two
derived points. They are singletons, and both are also callable to name a
subsystem or a blended subset.

### The anchored frame

A `variables=` block whose absolute frame is parameterized by an anchor
source's catalogue solution rather than by the barycentre's — the
parameterization several Gaia sources on one frame need. See
[Anchoring the frame to a source](@ref g23h-anchored).

```@docs
AnchoredFrame
anchored_frame
anchor_offsets
barycentre_parallax
```

### Sampling-coordinate Jacobians

```@docs
logjac_cartesian_to_campbell
```

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
FGSEpochAstromObs
```

These four live in the two unregistered subpackages, which `make.jl` loads so
that their docstrings resolve here:

```@docs
ImageObs
LogLikelihoodMapObs
InterferometryObs
GRAVITYWideKPObs
ClosurePhases
KernelPhases
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
Octofitter.CorrectionReport
Octofitter.CorrectionDecision
Octofitter.ObsImpact
Octofitter.correction_impact
Octofitter.has_correction_impact
Octofitter.reduced_lighttime_free
Octofitter.recheck_corrections
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
rvplot
rvplot_animated
dotplot
gaiastarplot
gaiastarplot!
gaiatimeplot
skytrackplot
hipparcosplot
completenessplot
completenessplot!
PlotChannel
plotchannels
plotobs
defaultpanels
sharepanel
datacalibration
noisemodel
default_queries
predictedchannels
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
Octofitter.phasebinmeans
Octofitter.residuals
```

The pre-v9 names still work and forward with a deprecation warning:

```@docs
rvpostplot
rvpostplot_animated
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
g23h_scan_uncertainty
gaia_dr4_transit_template
Octofitter.GOST_forecast
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

**Error stubs.** Two names are still defined, and raise an error naming the
replacement rather than a bare `UndefVarError`. They are the two an old script
hits first, and in both cases the replacement takes different arguments, so an
alias would not have helped:

| Retired name | Use instead |
|---|---|
| `Planet` | [`Body`](@ref) + the observation moves to [`System`](@ref)'s `observations=` |
| `θ_at_epoch_to_tperi` | declare `θ` and `epoch` as orbital elements |
| `HGCAInstantaneousObs`, `GaiaCatalogFitObs` | [`HGCAObs`](@ref) / [`G23HObs`](@ref) |

The HGCA pair is the only case where the *modelling code* was not ported at
all — it was subsumed by [`G23HObs`](@ref).

**Everything else is simply gone**, so the name is free again:
`PlanetRelAstromObs` → [`RelAstromObs`](@ref), `StarAbsoluteRVObs` and
`PlanetRelativeRVObs` → [`RadialVelocityObs`](@ref),
`MarginalizedStarAbsoluteRVObs` → [`MarginalizedRVObs`](@ref),
`masspostplot` → [`octocorner`](@ref), `PhotometryLikelihood` →
[`PhotometryObs`](@ref), `PlanetOrderPrior` → [`OrbitOrderPrior`](@ref),
`ObsPriorAstromONeil2019` → [`ObsPriorONeil2019`](@ref),
`InterferometryLikelihood` → [`InterferometryObs`](@ref),
`GRAVITYWideKPLikelihood` → [`GRAVITYWideKPObs`](@ref). See
[Migrating to Octofitter v9](@ref v9-migration) for the call signatures.

Orbit element types, observables (`raoff`, `decoff`, `radvel`, `posx`,
`projectedseparation`, …), `orbitsolve`, `Trajectory`, `WeightedPoint`,
`barycentre`, `photocentre` and `fluxes` come from PlanetOrbits, which
Octofitter re-exports; they are documented in the PlanetOrbits manual.
