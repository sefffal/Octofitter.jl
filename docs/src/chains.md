# [Chains](@id chains)

This page describes the format of the Monte Carlo chains created by Octofitter.jl


The output of the samplers in Octofitter is an MCMCChains.Chains object

A column will be present for each variable in your model, whether it was defined as a prior (`~`) or as a derived variable (`=`).

Column names follow the namespace the variable was defined in:

| Defined on | Column name | Example |
|---|---|---|
| the `System` | the bare variable name | `plx` |
| a `Body` (or `Orbit`) node | `<node>_<variable>` | `A_mass`, `b_a`, `b_e` |
| an observation | `<observation>_<variable>` | `GPI_jitter` |

Variables defined for the System as a whole can be accessed directly. For example:
```julia
chain["plx"]
```
This will return an array of parallax values from the posterior. The format is a matrix of N-samples by N-chains.

Variables belonging to a body are prefixed with that body's name:
```julia
chain["b_a"]
```
This returns an array of semi-major axis values (`a`) for the planet `b` sampled from the posterior.

!!! note "How columns are named"
    Body variables are prefixed with the body's name (`b_a`, `b_e`, `b_flux_H`);
    observation variables are prefixed with the observation's name (`GPI_jitter`);
    system variables have no prefix (`plx`). 

## Reconstructing orbits from a chain

To evaluate orbits from the posterior, rebuild the whole system for a draw with
[`construct_system`](@ref) and query it. Every observable takes
`(solution, target, reference)`:

```julia
posys = construct_system(model, chain, 1)         # draw #1
traj  = orbitsolve(posys, [mjd("2030-01-01")])
raoff(traj[1], :b, :A)                            # [mas]

# The barycentre and photocentre are resolved against the sample's masses/fluxes:
bc = PlanetOrbits.barycentre(posys)
radvel(traj[1], :A, bc)                           # [m/s] stellar reflex
```

!!! note "`Barycentre` vs `barycentre`"
    `Barycentre` and `Photocentre` (capitalized) are *model declaration* specs — they are
    what you write in an observation's `target=`/`ref=`, or in an
    [`ObservableQuery`](@ref). To call an observable directly you need a reference
    resolved against a particular sample: `PlanetOrbits.barycentre(posys)` or
    `PlanetOrbits.photocentre(posys; band=:H)`. Passing the capitalized spec straight to
    `radvel` is a `MethodError`.

There is no per-planet orbit object: a draw is one `PlanetOrbits.System` containing
every body and every orbit.

## Diagnostics
Printing the chains will display a summary of the chains size and columns. Running `describe(chain)` will output a number of useful summaries for each quantity, like the mean, 0.25, 0.5, and 0.75 quantiles, and convergence metrics. See MCMCChains documentation for more details.

Sampler diagnostics (`tree_depth`, `numerical_error`, `step_size`, …) are stored in the
chain's `:internals` section rather than as model parameters, so they do not show up in
`describe(chain)` or corner plots. That section survives a `savechain`/`loadchain`
round-trip.

## Exporting Chains

See [Loading and Saving Data](@ref loading-saving) for the FITS and HDF5 formats, which
preserve the model metadata. This page covers plain tabular export.

### As a table
You can convert your chains to any Tables.jl compatible table. `TypedTables.Table` is included with this package, but `DataFrames.DataFrame` works well too.
```julia
tbl = Table(chain)

using DataFrames
df = DataFrame(chain)
```

You can then use a wide variety of Tables.jl source or sink libraries to persist your data to a file or database. The easiest is probably Arrow.jl:

```julia
Arrow.write("mychain.arrow", tbl)
```

Other useful formats could be CSV.jl or SQLite.jl.
In these formats, the data can be archived and imported easily into other programs; however, there is not yet an automatic way to return the data into the `MCMCChains.Chains` format it originated in.
