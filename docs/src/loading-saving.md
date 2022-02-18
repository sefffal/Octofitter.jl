# [Loading and Saving Data](@id loading-saving)


## Loading Observations
For models with lots of data points, it becomes cumbersome to write all your data in your model script.
Intead, you can load your observations (astrometry, proper motion anomaly, radial velocity, etc) from any Tables.jl
compatible source. These could include a TypedTable, a DataFrame, a CSV file, an Arrow file, Excel, etc.

Here is an example of loading data from a CSV file:
```julia
using CSV
astrom = Astrometry(CSV.File("astrom.csv"))
# Or equivalently
astrom = CSV.read("astrom.csv", Astrometry)
```

The list of columns necessary for each type of observation are listed in the API documentation for e.g. Astrometry.

This works for other observation types too:
```julia
pma = CSV.read("pma.csv", PropMotionAnom)
```

This pattern also allows you to load data directly from remote databases using any Tables.jl compatible library.

Once loaded, you can access the underlying table using e.g. `astrom.table`.

## Saving Chains
The easiest way to save your chains is using the `JLD2` package. This saves the chains as well as the model for future reference, but is only compatible with the same Julia version. That is, chains saved with JLD2 in Julia 1.7 may not be loadable in Julia 1.8.

Saving chains with JLD2:
```julia
using JLD2
JLD2.jldsave("chains.jld2"; chains)
```

Other forward-compatible ways to save your chains are to convert them into a table using DataFrames:
```julia
using DataFrames
df = DataFrame(chains)
```
which can then be saved to any Tables.jl compatible format:
```julia
using CSV
CSV.write("chains.csv", df)

using Arrow
Arrow.write("chains.arrow", df)
```

Or a general Array which you can save in any format you wish:
```julia
arr = Array(chains)
```