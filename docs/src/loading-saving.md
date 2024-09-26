# [Loading and Saving Data](@id loading-saving)


## Loading Observations
For models with lots of data points, it becomes cumbersome to write all your data in your model script.
Intead, you can load your observations (astrometry, proper motion anomaly, radial velocity, etc) from any Tables.jl
compatible source. These could include a TypedTable, a DataFrame, a CSV file, an Arrow file, Excel, etc.

Here is an example of loading data from a CSV file:
```julia
using CSV
astrom = PlanetRelAstromLikelihood(CSV.File("astrom.csv"))
# Or equivalently
astrom = CSV.read("astrom.csv", PlanetRelAstromLikelihood)
```

The list of columns necessary for each type of observation are listed in the API documentation for e.g. PlanetRelAstromLikelihood.

This works for other observation types too:
```julia
pma = CSV.read("pma.csv", PropMotionAnom)
```

This pattern also allows you to load data directly from remote databases using any Tables.jl compatible library.

Once loaded, you can access the underlying table using e.g. `astrom.table`.

## Saving Chains

There are two ways you can save chains for later analysis. The first is a built in function that stores the chain and metadata into a FITS table. The second is converting the chain to a Table and saving it using any Tables.jl compatible package (CSV, Arrow, SQL, etc.)

#### Example: Saving chains and metadata to FITS Table
The default, and recommended way to save your chains is to a FITS table:
```julia
Octofitter.savechain("mychain.fits", chain)

chain = Octofitter.loadchain(fname)
```

#### Example: Saving chains to Orbitize format
For compatbility purposes, orbit posteriors can be exported and loaded from the Orbitize! HDF5 format. This only works for basic two-object orbits. FITS format (above) should be preferred.
```julia
Octofitter.savehdf5("mychain.h5", chain)

chain = Octofitter.loadhdf5(fname)
```


#### Example: Saving to CSV

Converting chain to a TypedTables.jl Table (re-exported by this package)
```julia
tbl = Table(chain)
```

Converting chain to a DataFrames.jl DataFrame:
```julia
df = DataFrame(Chain)
```

Saving chains:
```julia
using CSV
CSV.write("chains.csv", tbl) # or df

using Arrow
Arrow.write("chains.arrow", tbl) # or df
```

You can also convert a chain object to general Array which you can save in any format you wish:
```julia
arr = Array(chains)
```


## Saving and Restoring Models
We recommend that you save each model as a script that generates the model, e.g. in a julia file called `model-systemname.jl`. 

For convenience, it is also possible to save and restore the full model. This is not garuanteed to work across Julia versions or between computers, but is very fast for interactive work etc.

Saving model:
```julia
using Serialization
serialize("mymodel-systemname.jls", model)
```

Restoring model:
```julia
using Octofitter # must load all previously used dependencies
using Serialization 
model = deserialize("mymodel-systemname.jls")
```