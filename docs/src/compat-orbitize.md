# Compatibility with Orbitize!

The [Orbitize!](https://orbitize.readthedocs.io/en/latest/) python library is a popular package for fitting astrometric and radial velocity orbits.

Octofitter has support for loading and saving posteriors in HDF5 format--the same format used by Orbitize!. This is useful if you want to load an Orbitize! posterior into Octofitter for plotting or to compare results.
Similarily, you can export Octofitter chains to use with Orbitize! analysis tools, including the popular [whereistheplanet.com](http://whereistheplanet.com) website for predicting planet locations from stored posteriors.

!!! warning
    The Orbitize! import/export functionality only works with simple models with visual orbits and only one companion.

In addition, it is possible to load a orbit posterior and/or astrometry data directly from whereistheplanet.com by target name.


## Loading an Orbitize! posterior

```julia
chain = Octofitter.loadhdf5("fname.h5")
```

Both tools use the same orbital-element conventions, so the import is a rename plus one change of phase variable: orbitize!'s `tau` — a fraction of a period past a reference epoch — becomes Octofitter's `tp`, which is added as an extra column alongside the original `tau`.

Column names are derived from **v2 body names**, which you can choose:

```julia
# m0 -> A_mass, sma1 -> b_a, ecc1 -> b_e, …   (the defaults)
chain = Octofitter.loadhdf5("fname.h5")

# …or name the bodies yourself
chain = Octofitter.loadhdf5("fname.h5"; host=:Aa, bodynames=("Ab", "B"))
```

If the file stores the standard eight-column basis with a total mass (`mtot`) rather than individual masses, `mtot` keeps its own name and is *not* renamed onto a body — it is a property of the pair, not of either component.

If the file's `parameter_labels` attribute is missing, pass the column names yourself:

```julia
chain = Octofitter.loadhdf5("fname.h5";
    colnames=["sma1", "ecc1", "inc1", "aop1", "pan1", "tau1", "plx", "mtot"])
```

A file holding several concatenated chains can be split back apart:

```julia
chain = Octofitter.loadhdf5("fname.h5", 4)   # 4 chains
```

!!! note "Changed from v1"
    * **Masses are solar masses.** v1 divided imported companion masses by
      `mjup2msol` on the way in, because a v1 `Planet.mass` was in Jupiter masses.
      v2 has one mass unit throughout. Anything in your analysis that compares a
      mass against a literal threshold needs updating (`mass > 10` meaning
      10 M_jup becomes `mass > 10mjup`).
    * **There is no synthesised `<planet>_M` column.** v1 built a per-planet total
      mass by summing every body with a smaller semi-major axis — the epicyclic
      superposition selection this release removes. A body's dynamical mass now
      comes from the hierarchy.
    * `numchains > 1` works. In v1 it always errored, because the derived `tp`
      columns were `hcat`ed onto an already-3-D array.

## Save a posterior in Orbitize! format

```julia
Octofitter.savehdf5("fname.h5", model, chain)
```

Note the **three** positional arguments: `savehdf5` needs the model in order to know which body is the companion and which is its host. A fourth argument selects the companion explicitly when the model has more than one:

```julia
Octofitter.savehdf5("fname.h5", model, chain, :c)
```

Only the eight columns of orbitize!'s standard basis are written — `sma`, `ecc`, `inc`, `aop`, `pan`, `tau`, `plx`, `mtot` — so this exports **one companion's** visual orbit and nothing else. No data are exported.

A few constraints follow from that basis:

* The chain needs a `<body>_tp` column, since `tau` is derived from it. A model
  parametrized on `θ` + `epoch` does not have one.
* The chain needs `<host>_mass` and `<body>_mass`, whose sum becomes `mtot`.
* A companion that orbits a **barycentre** rather than a single body is rejected
  rather than guessed at — orbitize!'s two-body standard basis cannot express it.

`savehdf5` accepts either a `LogDensityModel` or a bare `System`.

!!! note
    v2 also writes the `parameter_labels` HDF5 attribute, not just `col_names`.
    v1 wrote only the latter while its own `loadhdf5` read the former, so
    Octofitter could not round-trip its own export without falling through to a
    guess plus a warning.

## Loading an Orbitize! posterior saved to Whereistheplanet.com

```julia
chain = Octofitter.loadhdf5("51erib")
```

Passing a name rather than a filename looks the target up on whereistheplanet.com. `Octofitter.Whereistheplanet_search("51erib")` returns the matching file path, and reports similar names if there is no exact match.

## Loading Astrometry Data saved to Whereistheplanet.com

```julia
astrom_obs_seppa, astrom_obs_radec = Octofitter.Whereistheplanet_astrom(
    "51erib"; object=1, target=:b, ref=:A)
```

Two different astrometry likelihood objects are returned since orbitize supports both PA/sep and RA/DEC formats. Octofitter also supports both formats, but they must be placed into separate likelihood objects. Simply add both to your `System`'s `observations=` list to include all the data.

`target` and `ref` say which model references the astrometry measures — the companion and the host in the usual case. They take the full v2 reference grammar, so a `Body` model node, a `Symbol`, `Barycentre(...)` or `Photocentre(...)` all work; whereistheplanet's "object 1" is a companion only by convention.

The two objects are named `"<name>_seppa"` and `"<name>_radec"` (default `name="whereistheplanet"`), which is what lets both go into one `System` — observation names must be unique.

!!! note "Changed from v1"
    v1's `return out` was commented out, so the documented
    `astrom1, astrom2 = Whereistheplanet_astrom(...)` actually destructured `nothing`.
    It returns its vector now.

    Note also that only formats *present in the file* come back: a target with only
    sep/PA data yields a one-element vector, and destructuring it into two names
    will error. Assign the result and check its length if you are not sure.
