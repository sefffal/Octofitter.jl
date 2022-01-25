# [Chains](@id chains)

This page describes the format of the Monte Carlo chains created by DirectDetections.jl


The output of the samplers in DirectDetections is an MCMCChains.Chains object

A column will be present for each variable in your model, both defined in the Priors blocks or as Derived variables. 

Variables defined for the System as a whole can be accessed directly. For example:
```julia
chain["μ"]
```
This will return an array of μ values from the posterior. The format is a matrix of N-samples by N-walkers.

Variables defined for an individual Planet are grouped according to the standard MCMCChains format. For example:
```julia
chain["b[a]"]
```
This returns an array of semi-major axis values (`a`) for the planet `b` sampled from the posterior.

## Diagnostics
Printing the chains will display a number of useful summaries for each quantity, like the mean, 0.25, 0.5, and 0.75 quantiles, and convergence metrics. See MCMCChains documentation for more details.

## Exporting Chains
There are two useful ways to export chains. One is with the JLD2 library which preserves all the information and structure of the chain (but is only easy to open again in Julia) and the other is converting them to tables.

### Using JLD2
Using the JLD2 library, you can save your chain like so:
```julia
using JLD2
@save "mychain.jld2" chain
```
And open it again (if using the same version of Julia and DirectDetections):
```julia
@load "mychain.jld2"
```

### As a table
To convert your chains into a table, first use the DataFrames.jl library:
```julia
using DataFrames
df = DataFrame(chain)
```

You can then use a wide variety of Tables.jl source or sink libraries to persist your data to a file or database. The easiest is probably CSV.jl:

```julia
CSV.write("mychain.csv", df)
```

Other useful formats could be Arrow.jl or SQLite.jl which may be more efficient.
In these formats, the data can be archived and imported easily into other programs; however, there is not yet an automatic way to return the data into the `MCMCChains.Chains` format it originated in.