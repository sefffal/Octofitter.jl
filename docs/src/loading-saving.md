# [Loading and Saving Data](@id loading-saving)


## Loading Observations
For models with lots of data points, it becomes cumbersome to write all your data in your model script.
Intead, you can load your observations (astrometry, radial velocity, epoch astrometry, etc) from any Tables.jl
compatible source. These could include a TypedTable, a DataFrame, a CSV file, an Arrow file, Excel, etc.

Here is an example of loading data from a CSV file:
```julia
using CSV
astrom_dat = CSV.read("astrom.csv", Table)
astrom = RelAstromObs(astrom_dat; target=b, ref=A, name="GPI")
```

The list of columns necessary for each type of observation is listed in the API documentation for e.g. [`RelAstromObs`](@ref).

This works for other observation types too:
```julia
rv_dat = CSV.read("rvs.csv", Table)
rvs = RadialVelocityObs(rv_dat; target=A, ref=Barycentre, name="HARPS",
    variables=@variables begin
        offset ~ Normal(0, 100)
        jitter ~ LogUniform(0.1, 100)
    end)
```

This pattern also allows you to load data directly from remote databases using any Tables.jl compatible library.

Once loaded, you can access the underlying table using e.g. `astrom.table`.


## Saving Chains

There are two ways you can save chains for later analysis. The first is a built in function that stores the chain and metadata into a FITS table. The second is converting the chain to a Table and saving it using any Tables.jl compatible package (CSV, Arrow, SQL, etc.)

#### Example: Saving chains and metadata to FITS Table
The default, and recommended way to save your chains is to a FITS table:
```julia
Octofitter.savechain("mychain.fits", chain)

chain = Octofitter.loadchain("mychain.fits"; model)
```

Passing `model=` is the recommended spelling. It runs [`Octofitter.checkchain`](@ref), which asserts that the chain carries every free parameter the model expects. That check is worth doing because the failure it catches is silent: `mcmcchain2result` looks each parameter up by name and yields `missing` for anything absent, so a chain from a *different* model flows on into plotting and post-prediction producing nonsense rather than an error.

You can also run the check by hand on a chain you already have in memory:
```julia
Octofitter.checkchain(model, chain)              # errors on a mismatch
Octofitter.checkchain(model, chain; strict=false) # warns instead
```

!!! note "Loading chains written by Octofitter v8 or earlier"
    They load, with a warning. Their samples are unchanged, but their **column names**
    follow the older model surface: an observation attached to a companion was named
    `<planet>_<observation>_<variable>` (`b_GPI_jitter`), where it is now
    `<observation>_<variable>` (`GPI_jitter`), and the host star's mass moves from the
    system-level `M` to `A_mass`.

    [`Octofitter.checkchain`](@ref) recognizes that pattern and prints the suggested
    renames rather than leaving you to discover the mismatch downstream. See
    [Migrating to Octofitter v9](@ref v9-migration).

#### Example: Saving chains to Orbitize format
For compatbility purposes, orbit posteriors can be exported and loaded from the Orbitize! HDF5 format. This only works for basic two-object orbits. FITS format (above) should be preferred.
```julia
Octofitter.savehdf5("mychain.h5", model, chain)

chain = Octofitter.loadhdf5("mychain.h5")
```

Note the **three** positional arguments to `savehdf5`: it needs the model in order to know which body is the companion and which is its host. See [Compatibility with Orbitize!](@ref compat-orbitize) for the details and for the additional keywords.

!!! note
    orbitize!'s standard basis stores the epoch of periastron, so `savehdf5` requires a
    `<body>_tp` column in the chain. A model parametrized on `θ` + `epoch` (position angle
    at a reference epoch) does not have one — either add `tp` explicitly, or export from a
    model that samples `tp`.


#### Example: Saving to CSV

Converting chain to a TypedTables.jl Table (re-exported by this package)
```julia
tbl = Table(chain)
```

Converting chain to a DataFrames.jl DataFrame:
```julia
df = DataFrame(chain)
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
arr = Array(chain)
```

Note that these plain-table formats keep only the numbers. The run metadata, the
`:parameters`/`:internals` split, and the format stamp that `loadchain` uses are all lost;
use the FITS format if you want the chain to survive as a chain.


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

## Rebuilding orbits from a chain

Neither of the formats above stores orbits — they store the parameters an orbit is built from. To go back the other way, use [`construct_system`](@ref), which turns one chain row into a full `PlanetOrbits.System`:

```julia
posys = construct_system(model, chain, 1)          # row 1 of the chain
traj  = orbitsolve(posys, [59000.0, 59100.0])
raoff(traj[1], :b, :A)
```

