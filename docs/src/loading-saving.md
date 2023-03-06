# [Loading and Saving Data](@id loading-saving)


## Loading Observations
For models with lots of data points, it becomes cumbersome to write all your data in your model script.
Intead, you can load your observations (astrometry, proper motion anomaly, radial velocity, etc) from any Tables.jl
compatible source. These could include a TypedTable, a DataFrame, a CSV file, an Arrow file, Excel, etc.

Here is an example of loading data from a CSV file:
```julia
using CSV
astrom = AstrometryLikelihood(CSV.File("astrom.csv"))
# Or equivalently
astrom = CSV.read("astrom.csv", AstrometryLikelihood)
```

The list of columns necessary for each type of observation are listed in the API documentation for e.g. AstrometryLikelihood.

This works for other observation types too:
```julia
pma = CSV.read("pma.csv", PropMotionAnom)
```

This pattern also allows you to load data directly from remote databases using any Tables.jl compatible library.

Once loaded, you can access the underlying table using e.g. `astrom.table`.

## Saving Chains

You can save your chains using any Tables.jl compatible data package. I recommend Arrow.jl since it is a binary table format that takes up less space but is still cross compatible with other languages like Python. Another good choice is CSV.jl.


### Examples

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
