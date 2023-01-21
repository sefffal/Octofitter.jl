# [Chains](@id chains)

This page describes the format of the Monte Carlo chains created by Octofitter.jl


The output of the samplers in Octofitter is an MCMCChains.Chains object

A column will be present for each variable in your model, both defined in the Priors blocks or as Derived variables. 

Variables defined for the System as a whole can be accessed directly. For example:
```julia
chain["M"]
```
This will return an array of μ values from the posterior. The format is a matrix of N-samples by N-walkers.

Variables defined for an individual Planet are grouped according to the standard MCMCChains format. For example:
```julia
chain["b_a"]
```
This returns an array of semi-major axis values (`a`) for the planet `b` sampled from the posterior.

## Diagnostics
Printing the chains will display a number of useful summaries for each quantity, like the mean, 0.25, 0.5, and 0.75 quantiles, and convergence metrics. See MCMCChains documentation for more details.

## Exporting Chains
There are two useful ways to export chains. One is with the JLD2 library which preserves all the information and structure of the chain (but is only easy to open again in Julia) and the other is converting them to tables.


### As a table
You can convert your chains to any Tables.jl compatible table. `TypedTables.Table` is included with this package, but `DataFrames.DataFrame` works well too.
```julia
tbl = Table(chain)

using DataFrames
df = DataFrame(chain)
```

You can then use a wide variety of Tables.jl source or sink libraries to persist your data to a file or database. The easiest is probably Arrow.jl:

```julia
Arrow.write("mychain.csv", tbl)
```

Other useful formats could be CSV.jl or SQLite.jl.
In these formats, the data can be archived and imported easily into other programs; however, there is not yet an automatic way to return the data into the `MCMCChains.Chains` format it originated in.
