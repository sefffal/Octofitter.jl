# [Full G23H Example Script](@id g23h-example)

This page provides a complete, production-ready script for fitting orbital models using the G23H catalog. This script supports:

!!! warning "v2 changes this script depends on"
    * `G23HObs` now declares **source membership** — `host=`, `companions=`, `ref=` — and
      reads companion flux ratios from the bodies' own `flux_G` / `flux_Hp` variables.
    * `Planet` / `companions=` are replaced by `Body` / `bodies=`, and **all masses are in
      solar masses** (`mjup` is a plain multiplicative constant).
    * `StarAbsoluteRVObs` is gone: it is `RadialVelocityObs(…; target=A, ref=Barycentre)`,
      and `offset`/`jitter` are no longer auto-injected.
    * `PlanetOrderPrior` is renamed [`OrbitOrderPrior`](@ref) and takes `Body` nodes or
      `Symbol`s.
    * The incremental parallel-tempering loop this script is built around
      ([`octofit_pigeons`](@ref) + `increment_n_rounds!` + `Pigeons.stepping_stone_pair`)
      is unchanged: `octofit_pigeons` lives in a package extension, so `using Pigeons`
      must be in scope, and it returns `(;chain, pt)` as before.
    * `Octofitter.dotplot` and `Octofitter.rvpostplot` both work, with the same
      arguments; `rvpostplot` is now sugar over [`octoplot`](@ref)'s panels and lives in
      core rather than in an RV-specific Makie extension.

- Single or multi-planet systems
- Optional Gaussian Process modeling of stellar activity in RV data
- RV data from DACE or custom RDB files
- Flexible selection of which astrometric observation types to include
- Incremental sampling with checkpoint saving

## Prerequisites

Before running this script, ensure you have installed the necessary packages:

2. **Required packages**:
   ```julia
   using Pkg
   Pkg.add([
       "Octofitter",
       "OctofitterRadialVelocity",
       "Pigeons",
       "Arrow",
       "CSV",
       "DataFrames",
       "Distributions",
       "CairoMakie",
       "PairPlots",
       "MCMCChains",
       "ArgParse",
       "DelimitedFiles",
       "FiniteDiff",
       "DifferentiationInterface",
   ])
   ```

3. **Optional for DACE RV data**:
   ```julia
   using CondaPkg
   # Then in Pkg mode: conda pip_add dace-query
   ```

## Script Overview

The script is organized into several sections:

1. **Argument Parsing** - Command-line interface for specifying targets and options
2. **Observation Setup** - Create G23H astrometry likelihood object (catalog downloaded automatically)
3. **RV Data Loading** - Load RV data from DACE or RDB files
4. **Model Definition** - Define planet and system models with optional GP
5. **Sampling** - Run MCMC
6. **Analysis** - Generate plots and save results

## Complete Script

```julia
#=
G23H Joint Astrometry + RV Fitting Script
=========================================

This script fits orbital models using the G23H catalog for joint
Gaia-Hipparcos astrometry modeling, optionally combined with RV data.

Usage:
    julia script.jl --hip 12345 --host-mass 1.0 --n-rounds 10
    julia script.jl --gaia 756291174721509376 --host-mass 1.0 --dace-rv
    julia script.jl --hip 12345 --host-mass 1.0 --rdb-rv mydata.rdb --gp

=#

# First time only: install dace-query
# using CondaPkg;  Then type "]" to go into Pkg-mode and paste: conda pip_add dace-query

using ArgParse
using Arrow
using CSV
using CairoMakie
using DataFrames
using Dates
using Distributions
using Octofitter, OctofitterRadialVelocity
using Pigeons
using LinearAlgebra
using PairPlots
using Printf
using Statistics
using StatsBase
using MCMCChains
using DelimitedFiles
using FiniteDiff, DifferentiationInterface

## ============================================================================
## ARGUMENT PARSING
## ============================================================================

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--host-mass"
            help = "Host star mass in solar masses"
            arg_type = Float64
            required = true
        "--n-planets"
            help = "Number of planets in the system"
            arg_type = Int
            default = 1
        "--hip"
            help = "Hipparcos ID (mutually exclusive with --gaia)"
            arg_type = Int
        "--gaia"
            help = "Gaia source ID (mutually exclusive with --hip)"
            arg_type = Int64
        "--gaia-rv"
            help = "Include Gaia DR3 RV in astrometry likelihood"
            action = :store_true
        "--dace-rv"
            help = "Query DACE for RV data"
            action = :store_true
        "--rdb-rv"
            help = "Load RV data from RDB file (can be specified multiple times)"
            action = :append_arg
        "--gp"
            help = "Enable Gaussian Process kernel for stellar activity"
            action = :store_true
        "--n-rounds"
            help = "Number of MCMC rounds"
            arg_type = Int
            default = 12
        "--obs-types"
            help = "Comma-separated list of observation types to include (default: all)"
            arg_type = String
        "--circular"
            help = "Force a circular (e=0) orbit for all companions"
            action = :store_true
    end

    return parse_args(s)
end

args = parse_commandline()

# Validate mutually exclusive arguments
if isnothing(args["hip"]) && isnothing(args["gaia"])
    error("Either --hip or --gaia must be specified")
end
if !isnothing(args["hip"]) && !isnothing(args["gaia"])
    error("Cannot specify both --hip and --gaia")
end

# Extract configuration
host_mass = args["host-mass"]
num_planets = args["n-planets"]
n_rounds = args["n-rounds"]
use_gaia_rv = args["gaia-rv"]
use_dace_rv = args["dace-rv"]
rdb_files = isnothing(args["rdb-rv"]) ? String[] : args["rdb-rv"]
use_gp = args["gp"]
circular = args["circular"]

# Generate planet names
planet_names = Symbol.(Char.(Int('b') .+ (0:num_planets-1)))

output_dir = "results"
mkpath(output_dir)

## ============================================================================
## CREATE G23H OBSERVATION OBJECT
## ============================================================================

# G23HObs accepts either gaia_id or hip_id directly.
# The G23H catalog is automatically downloaded on first use via DataDeps.
# `host` / `companions` name the bodies this catalog source is a measurement of. They
# accept `Body` model nodes or `Symbol`s; we use Symbols here because the bodies are
# built further down. `ref=Barycentre` says the host's reflex is measured against the
# system barycentre, which is what a catalog proper motion is.
if !isnothing(args["hip"])
    hip_id = args["hip"]
    sysname′ = "HIP$hip_id"
    absastrom = Octofitter.G23HObs(;
        hip_id=hip_id,
        host=:A,
        companions=Tuple(planet_names),
        ref=Barycentre,
        include_rv=use_gaia_rv,
    )
else
    gaia_id = args["gaia"]
    absastrom = Octofitter.G23HObs(;
        gaia_id=gaia_id,
        host=:A,
        companions=Tuple(planet_names),
        ref=Barycentre,
        include_rv=use_gaia_rv,
    )
    # Determine system name from catalog
    if isfinite(absastrom.catalog.hip_id)
        hip_id = round(Int, absastrom.catalog.hip_id)
        sysname′ = "HIP$hip_id"
    else
        sysname′ = "GDR3$gaia_id"
    end
end

ra = absastrom.catalog.ra
dec = absastrom.catalog.dec
println("Found target: RA=$ra, Dec=$dec")

# Determine which observation types to include
if isnothing(args["obs-types"])
    # Default: include all available types
    obs_to_include = [
        :iad_hip,
        :ra_hip,  :dec_hip,
        :ra_hg,   :dec_hg,
        :ra_dr2,  :dec_dr2,
        :ra_dr32, :dec_dr32,
        :ra_dr3,  :dec_dr3,
        :ueva_dr3,
    ]
    if use_gaia_rv
        push!(obs_to_include, :rv_dr3)
    end
else
    # Parse user-specified observation types
    obs_to_include = Symbol.(split(args["obs-types"], ","))
    # Add rv_dr3 if --gaia-rv is specified
    if use_gaia_rv && !(:rv_dr3 in obs_to_include)
        push!(obs_to_include, :rv_dr3)
    end
end

# `channels=` on the constructor does the same thing, filtering the same `table.kind`
# rows; the post-hoc route below is kept because the channel list is built from
# command-line arguments.
indices = findall([
    available_obs_kind ∈ obs_to_include
    for available_obs_kind in absastrom.table.kind
])
absastrom = Octofitter.likeobj_from_epoch_subset(absastrom, indices)

sysname = "$sysname′-" * join(string.(absastrom.table.kind), "-")

## ============================================================================
## GAUSSIAN PROCESS KERNEL (OPTIONAL)
## ============================================================================

if use_gp
    function gauss_proc_approx_quasi_periodic(θ_system)
        B = θ_system.gp_B
        C = θ_system.gp_C
        L = θ_system.gp_L
        Prot = θ_system.gp_Prot

        kernel = OctofitterRadialVelocity.Celerite.RealTerm(
            log(B * (1 + C) / (2 + C)),  # log_a
            log(1 / L)                    # log_c
        ) + OctofitterRadialVelocity.Celerite.ComplexTerm(
            log(B / (2 + C)),            # log_a
            -Inf,                         # log_b
            log(1 / L),                  # log_c
            log(2π / Prot)               # log_d
        )

        return OctofitterRadialVelocity.Celerite.CeleriteGP(kernel)
    end
end

## ============================================================================
## RV HELPER FUNCTIONS
## ============================================================================

"""
Create RV likelihood for a single instrument.
Handles both GP and non-GP cases using composable variable blocks.
"""
function create_rv_likelihood(dat, name, mean_epoch, use_gp)
    # Base variables always present
    base_vars = @variables begin
        offset ~ Uniform(-1000, 1000)
        jitter ~ LogUniform(0.1, 100)
        m = system.m
    end

    # GP variables (conditional)
    gp_vars = @variables begin
        gp_B = system.gp_B
        gp_C = system.gp_C
        gp_L = system.gp_L
        gp_Prot = system.gp_Prot
    end

    # Combine variables based on whether GP is enabled
    vars = use_gp ? vcat(base_vars, gp_vars) : base_vars

    # v1 had `StarAbsoluteRVObs`; v2 has one RV type, and `target`/`ref` say which
    # kind of RV it is. `offset` and `jitter` are declared explicitly above because v2
    # never invents a prior for you.
    return RadialVelocityObs(
        dat;
        target = :A,
        ref = Barycentre,
        name,
        gaussian_process = use_gp ? gauss_proc_approx_quasi_periodic : nothing,
        variables = vars,
        trend_function = (θ_obs, epoch) -> θ_obs.m * (epoch - mean_epoch)
    )
end

"""
Perform outlier rejection on RV data using MAD (Median Absolute Deviation).
"""
function reject_rv_outliers(d, inst_names; clip_threshold=3.0)
    d_filt_parts = DataFrame[]
    total_original = 0
    total_filtered = 0

    for inst in inst_names
        d_inst = d[d.ins_name .== inst, :]
        n_inst_original = nrow(d_inst)
        total_original += n_inst_original

        # Outlier rejection using MAD within this instrument
        sigma = 1.4826 * mad(d_inst.rv)  # MAD to std conversion
        median_rv = median(d_inst.rv)

        # Create mask for non-outliers
        rv_mask = abs.(d_inst.rv .- median_rv) .< (clip_threshold * sigma)
        d_inst_filt = d_inst[rv_mask, :]

        n_inst_filtered = nrow(d_inst_filt)
        n_inst_rejected = n_inst_original - n_inst_filtered

        println("  $inst: kept $n_inst_filtered/$n_inst_original (rejected $n_inst_rejected)")

        if n_inst_filtered < 3
            @warn "Too few measurements for $inst after filtering ($n_inst_filtered < 3); keeping all for this instrument"
            d_inst_filt = d_inst
        end

        d_inst_filt.rv .-= median(d_inst_filt.rv)

        push!(d_filt_parts, d_inst_filt)
        total_filtered += nrow(d_inst_filt)
    end

    d_filt = vcat(d_filt_parts...)
    println("Total: kept $total_filtered/$total_original measurements across all instruments")

    return d_filt
end

"""
Process RV data for multiple instruments: outlier rejection, plotting, and likelihood creation.
Reusable for both DACE and RDB data sources.
"""
function process_rv_data(d_filt, inst_names; sysname, use_gp, source_name="")
    # Plot the filtered data
    f, a, p = errorbars(
        d_filt.rjd,
        d_filt.rv,
        d_filt.rv_err,
        color = map(i -> findfirst(==(i), inst_names), d_filt.ins_name),
        colormap = :turbo,
        axis = (; xlabel="time", ylabel="rv [m/s]", title=source_name)
    )

    plot_filename = isempty(source_name) ?
        "$sysname-rvs.png" :
        "$sysname-rvs-$(lowercase(source_name)).png"
    save(plot_filename, f, px_per_unit=4)

    # Create observations for each instrument
    observations = map(inst_names) do name
        dat = d_filt[d_filt.ins_name .== name, :]
        dat.epoch .= dat.rjd
        dat.rv .= dat.rv
        dat.σ_rv .= dat.rv_err
        dat.rv .-= median(dat.rv)

        mean_epoch = mean(dat.epoch)

        create_rv_likelihood(dat, name, mean_epoch, use_gp)
    end

    return observations
end

## ============================================================================
## LOAD RV DATA
## ============================================================================

rvlikes = []

# Load RV data from DACE
if use_dace_rv
    @info "Loading RV data from DACE"
    using PythonCall
    dace_spectroscopy = pyimport("dace_query.spectroscopy")

    hip_for_dace = !isnothing(args["hip"]) ? args["hip"] : absastrom.catalog.hip_id
    result = dace_spectroscopy.Spectroscopy.get_timeseries(
        target="HIP $(Int(hip_for_dace))",
        sorted_by_instrument=false,
        output_format="pandas"
    )
    d = DataFrame(pyconvert(PyTable, result))

    if isempty(d)
        @warn "NO RVS Found from DACE!"
    else
        # Perform outlier rejection per instrument
        println("Performing per-instrument outlier rejection for DACE data...")
        d_insts = unique(d.ins_name)

        d_filt = reject_rv_outliers(d, d_insts)

        # Process and create observations using helper function
        dace_likes = process_rv_data(d_filt, d_insts; sysname, use_gp, source_name="dace")
        append!(rvlikes, dace_likes)
    end
end

# Load RV data from RDB files
for rdb_file in rdb_files
    @info "Loading RV data from RDB file: $rdb_file"

    if !isfile(rdb_file)
        @warn "RDB file not found: $rdb_file"
        continue
    end

    # Read RDB file (tab-separated, skip first 2 header lines)
    d_raw = readdlm(rdb_file, '\t', skipstart=2)

    # Create DataFrame from RDB data
    # Columns: rjd, vrad, svrad, s_mw, sig_s, ha, sig_ha, weight
    d = DataFrame(
        rjd = Float64.(d_raw[:, 1]),
        rv = Float64.(d_raw[:, 2]),      # vrad
        rv_err = Float64.(d_raw[:, 3]),  # svrad
    )

    # Use filename as instrument name
    inst_name = replace(basename(rdb_file), ".rdb" => "")
    d.ins_name .= inst_name

    # Basic outlier rejection using MAD
    println("Performing outlier rejection for RDB file: $inst_name...")
    d_insts = [inst_name]

    d_filt = reject_rv_outliers(d, d_insts)

    # Process and create observations using helper function
    rdb_likes = process_rv_data(d_filt, d_insts; sysname, use_gp, source_name=inst_name)
    append!(rvlikes, rdb_likes)
end

if isempty(rvlikes)
    @info "No RV data loaded"
end

has_rv_data = !isempty(rvlikes)

## ============================================================================
## DEFINE PLANET MODEL(S)
## ============================================================================

# We should always use RUWE mode
ueva_mode = :RUWE

# There is no `basis=` in v2: which frame the model is in is decided by which frame
# variables the *system* block supplies. Declaring the full set
# `plx, ra, dec, pmra, pmdec, rv, ref_epoch` (below) is what v1 spelled
# `AbsoluteVisual{KepOrbit}`.

ref_planet_pos = mean(absastrom.gaia_table.epoch)

# The host star is a model node like any other. `flux_G` / `flux_Hp` = 1.0 make every
# other body's flux a contrast ratio against it, which is how G23H forms the Gaia
# photocentre and modulates the Hipparcos abscissa (v1 used a positional `fluxratio`
# vector on the observation).
host = Body(
    name="A",
    variables=@variables begin
        mass = $host_mass    # Msol
        flux_G  = 1.0
        flux_Hp = 1.0
    end
)

planets = Body[]
for planet_i in 1:num_planets
    planet = Body(
        name=planet_names[planet_i],
        about=host,
        variables=@variables begin
            planet_i = $planet_i
            a ~ LogUniform(0.001, 1024)
            planet_present = system.n_planets >= planet_i
            # Masses are M⊙ throughout; `mjup` is a plain multiplicative constant.
            # NB: `$`-interpolation works on the right of `=` (derived) but NOT on
            # the right of `~` — a prior expression closes over the enclosing scope
            # directly, so write `host_mass`, not `$host_mass`.
            mass′ ~ LogUniform(0.01mjup, host_mass)
            mass = mass′ * planet_present
            e ~ Uniform(0, 0.9)
            ω ~ Uniform(0, 2pi)
            i = system.i         # coplanar: shared inclination…
            Ω = system.Ω         # …and node
            # `M0` (mean anomaly at `epoch`) is the direct replacement for v1's
            # `τ ~ Uniform(0,1)` plus a hand-written `tp = τ*P*365.25 + …`; the orbit's
            # total mass comes from the hierarchy, so `M = …` is gone too.
            M0 ~ Uniform(0, 2pi)
            epoch = $ref_planet_pos
            flux_G  = 0.0        # dark companion
            flux_Hp = 0.0
        end
    )
    push!(planets, planet)
end

if num_planets == 0
    planets = [Body(name=:b, about=host, variables=@variables begin
        a = 1.0
        mass = 0.0
        e = 0.0
        ω = 0.0
        i = 0.0
        Ω = 0.0
        tp = 0.0
        flux_G  = 0.0
        flux_Hp = 0.0
    end)]
end

## ============================================================================
## DEFINE SYSTEM MODEL
## ============================================================================

rv = isnan(absastrom.catalog.radial_velocity) ? 0 : absastrom.catalog.radial_velocity * 1e3

observations = []
push!(observations, absastrom)
append!(observations, rvlikes)

if num_planets > 1
    # `PlanetOrderPrior` is renamed `OrbitOrderPrior` and takes `Body` nodes (or
    # `Symbol`s). Both terms go in the system's `observations=` list, not on a planet.
    push!(observations, OrbitOrderPrior(planets...))
    # With no `bodies=` list this covers *every* hierarchy row, which is v1's behaviour.
    # Naming the bodies restricts it, which matters for hierarchical systems.
    push!(observations, NonCrossingPrior(bodies=Tuple(planets)))
end

ref_epoch = Octofitter.meta_gaia_DR3.ref_epoch_mjd

# Base system variables. The host mass now lives on the host body, not here.
base_system_vars = @variables begin
    n_planets ~ truncated(NegativeBinomial(), upper=num_planets)

    ref_epoch = $ref_epoch
    plx ~ truncated(Normal(absastrom.catalog.parallax, absastrom.catalog.parallax_error),
                    lower=absastrom.catalog.parallax / 2)

    pmra ~ Uniform((Float64(absastrom.catalog.pmra_dr3) .+ (-10, 10))...)
    pmdec ~ Uniform((Float64(absastrom.catalog.pmdec_dr3) .+ (-10, 10))...)

    dec = $dec
    ra = $ra
    rv = $rv

    i ~ Sine()
    Ω ~ Uniform(0, 2pi)
end

# RV trend variable (conditional on having RV data)
rv_trend_vars = @variables begin
    m ~ Normal(0, 0.1)  # m/s/d trend
end

# GP hyperparameters (conditional on GP being enabled)
gp_system_vars = @variables begin
    gp_B ~ Uniform(0.00001, 2000000)
    gp_C ~ Uniform(0.00001, 200)
    gp_L ~ Uniform(40, 2000)
    gp_Prot ~ Uniform(35, 45)
end

# Compose system variables based on configuration
system_vars = base_system_vars
if has_rv_data
    system_vars = vcat(system_vars, rv_trend_vars)
end
if use_gp
    system_vars = vcat(system_vars, gp_system_vars)
end

sys = System(;
    name=sysname,
    bodies=[host; planets],
    observations,
    variables=system_vars
)

model = Octofitter.LogDensityModel(sys; verbosity=4, autodiff=AutoFiniteDiff())

## ============================================================================
## INITIALIZE AND SAMPLE
## ============================================================================

# Enable threaded Kepler solving for performance
Octofitter._kepsolve_use_threads[] = true

# Find good starting positions
initial_θ = collect(Octofitter.guess_starting_position(model, 10000)[1])
model.starting_points = fill(collect(model.link(initial_θ)), 100)
@time model.ℓπcallback(model.starting_points[1])

# Run initial sampling. G23H posteriors are frequently multi-modal, which is why this
# is tempered rather than HMC; seed the sampler well (above) regardless.
#
# The variational leg is switched off here: a Gaussian reference is a poor fit to these
# posteriors, and dropping it makes the ladder's behaviour easier to read off `Λ`.
explorer = SliceSampler()

chain, pt = octofit_pigeons(
    model,
    n_chains=32,
    n_chains_variational=0,
    variational=nothing,
    n_rounds=1,
    explorer=explorer,
    multithreaded=true,
)

## ============================================================================
## INCREMENTAL SAMPLING WITH CHECKPOINTS
## ============================================================================
#
# Pigeons can resume an existing run, so rather than one long call we add a round at a
# time and write a checkpoint after each. A long G23H fit can then be inspected (and
# survive an interruption) without starting over.

while length(pt.shared.reports.summary.last_round_max_time) < n_rounds
    global chain, pt
    increment_n_rounds!(pt, 1)
    chain, pt = octofit_pigeons(pt)
    display(chain)

    # Save checkpoint, tagged with the log evidence ratio for this round
    chain = MCMCChains.setinfo(chain, (; chain.info...,
        logevidence_ratio=Pigeons.stepping_stone_pair(pt)[2]))
    chain_file = joinpath(output_dir, "$sysname-post-round$(pt.inputs.n_rounds).fits")
    Octofitter.savechain(chain_file, chain)

    # Save sampling summary -- watch `Λ` and the log evidence ratio settle across rounds
    CSV.write(joinpath(output_dir, "$(sysname)-pigeons-summary.csv"),
              pt.shared.reports.summary)
end

Octofitter.savechain(joinpath(output_dir, "$sysname-post.fits"), chain)

## ============================================================================
## ANALYSIS AND VISUALIZATION
## ============================================================================

# Final chain summary
println("\n" * "="^60)
println("FINAL RESULTS")
println("="^60)
display(chain)

# Corner plot
corner_fig = octocorner(model, chain, small=true)
save(joinpath(output_dir, "$(sysname)-corner.png"), corner_fig, px_per_unit=2)

# Orbit plot. `octoplot` returns an `OctoPlotResult` (figure + named axes + the shared
# PosteriorSeries), so save via its `fname=` keyword or its `.figure` field.
orbit_res = octoplot(model, chain; fname=joinpath(output_dir, "$(sysname)-orbits.png"))

# Mass vs. semi-major axis summary, coloured by eccentricity.
dot_fig = Octofitter.dotplot(model, chain)
save(joinpath(output_dir, "$(sysname)-mass-sep.png"), dot_fig, px_per_unit=2)

# The same two columns, fully custom.
mass_fig = pairplot(
    (; a = chain["b_a"][:], mass = chain["b_mass"][:] ./ mjup) => (
        PairPlots.Scatter(markersize=4),
        PairPlots.MarginHist(),
        PairPlots.MarginQuantileText(),
    ),
    labels = Dict(:a => "semi-major axis [au]", :mass => "mass [Mⱼᵤₚ]"),
)
save(joinpath(output_dir, "$(sysname)-mass-sma.png"), mass_fig, px_per_unit=2)

# RV panels come from `octoplot` now: `RadialVelocityObs` declares its plot channels, so
# the orbit figure above already carries the RV time series, its residual strip, and a
# phase-folded panel per hierarchy row. `rvpostplot` is the RV-only slice of that same
# figure -- `octoplot` with `channels=radvel` and no sky panel.
rv_res = rvpostplot(model, chain; fname=joinpath(output_dir, "$(sysname)-rvpostplot.png"))

println("\nResults saved to: $output_dir")
```

## Command Line Usage

Run the script from the command line with various options:

```bash
# Basic usage with Hipparcos ID
julia script.jl --hip 12345 --host-mass 1.0 --n-rounds 10

# Using Gaia source ID
julia script.jl --gaia 756291174721509376 --host-mass 1.0

# With RV data from DACE
julia script.jl --hip 12345 --host-mass 1.0 --dace-rv

# With RV data from an RDB file
julia script.jl --hip 12345 --host-mass 1.0 --rdb-rv harps_data.rdb

# Multiple RDB files
julia script.jl --hip 12345 --host-mass 1.0 --rdb-rv harps.rdb --rdb-rv hires.rdb

# With Gaussian Process for stellar activity
julia script.jl --hip 12345 --host-mass 1.0 --dace-rv --gp

# Multi-planet fit
julia script.jl --hip 12345 --host-mass 1.0 --n-planets 2 --n-rounds 15

# Force circular orbits
julia script.jl --hip 12345 --host-mass 1.0 --circular

# Include Gaia DR3 RV variability constraint
julia script.jl --hip 12345 --host-mass 1.0 --gaia-rv

# Select specific observation types
julia script.jl --hip 12345 --host-mass 1.0 --obs-types "ra_hip,dec_hip,ra_hg,dec_hg,ueva_dr3"
```

## Command Line Arguments

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `--hip` | Int | Yes* | Hipparcos ID |
| `--gaia` | Int64 | Yes* | Gaia DR3 source ID |
| `--host-mass` | Float64 | Yes | Host star mass in solar masses |
| `--n-planets` | Int | No | Number of planets (default: 1) |
| `--n-rounds` | Int | No | Number of MCMC rounds (default: 12) |
| `--dace-rv` | Flag | No | Query DACE for RV data |
| `--rdb-rv` | String | No | Path to RDB file (can repeat) |
| `--gp` | Flag | No | Enable GP for stellar activity |
| `--gaia-rv` | Flag | No | Include Gaia RV variability |
| `--obs-types` | String | No | Comma-separated observation types |
| `--circular` | Flag | No | Force circular orbits |

*One of `--hip` or `--gaia` must be specified.

## Performance Tips

1. **Start with fewer rounds** for initial exploration, then increase for publication-quality results.

2. **Use `freeze_epochs=true`** for faster sampling during testing:
   ```julia
   absastrom = G23HObs(;
       gaia_id=gaia_id,
       host=:A, companions=(:b,), ref=Barycentre,
       freeze_epochs=true,  # Faster but approximate
   )
   ```

3. **Run Julia with multiple threads**:
   ```bash
   julia --threads=auto script.jl --hip 12345 --host-mass 1.0
   ```

4. **Consider using a subset of observations** for initial exploration, then add more data for final fits.

## See Also

- [G23H Overview](@ref fit-g23h) - Conceptual overview of the G23H method
- [Cross Validation](@ref cross-validation) - Techniques for model validation
- [Samplers](@ref samplers) - More on MCMC sampling options
