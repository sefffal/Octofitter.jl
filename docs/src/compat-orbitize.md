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

## Save a posterior in Orbitize! format

```julia
Octofitter.savehdf5("fname.h5", model, chain)
```


## Loading an Orbitize! posterior saved to Whereistheplanet.com

```julia
chain = Octofitter.loadhdf5("51erib",)
```

## Loading Astrometry Data saved to Whereistheplanet.com

```julia
astrom_like1, astro_like2 = Octofitter.Whereistheplanet_astrom("51erib"; object=1)
```

Two different astrometry likelihood objects are returned since orbitize supports both PA/sep and RA/DEC formats. Octofitter also supports both formats, but they must be placed into separate likelihood objects. Simply add both to the model to include all data.