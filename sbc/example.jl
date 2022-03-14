# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

# Input/Output
using JSON
using JLD2

## Plotting
using Plots
using PairPlots

## Orbit fitting
using Distributions
using DirectOrbits
using DirectDetections

## Simulation-based calibration
include("sbc.jl")

# ----------------------------------------------------------------------------------------------------------------------
# Example System
# ----------------------------------------------------------------------------------------------------------------------

## HD 984 B
@named b = Planet{KeplerianElements}(
    Variables(
        a = LogUniform(1, 1_000),
        e = Beta(0.867, 3.03),
        i = Sine(),
        ω = Uniform(-π, π),
        Ω = Uniform(-π, π),
        τ = Uniform(0, 1),
        mass = LogUniform(1, 200)
    ),
    Astrometry(
        (epoch=56126.0, ra=179.86335942198232, dec=-61.23048209379709, σ_ra=19.202512762873308, σ_dec=11.411998215929994),
        (epoch=56128.0, ra=196.7857546361294, dec=-67.37482298521509, σ_ra=22.063191571806936, σ_dec=12.99482437830172),
        (epoch=56910.0, ra=201.45140394838123, dec=-7.738982312044018, σ_ra=0.40537050357694315, σ_dec=1.758062198846047),
        (epoch=57264.0, ra=214.82281149517453, dec=25.235880434585788, σ_ra=1.001921920191633, σ_dec=1.1308442515245785),
        (epoch=57264.0, ra=216.54202958509396, dec=24.289080327758366, σ_ra=0.7007853226883667, σ_dec=0.7598909993473824),
        (epoch=58671.0, ra=197.49128022878907, dec=125.13846024941299, σ_ra=1.647113247206374, σ_dec=1.388300437337799),
        (epoch=59060.0, ra=190.38547085349813, dec=150.8435695941057, σ_ra=1.4980010280192033, σ_dec=1.3641618232112809)
    )
)

## HD 984
@named HD984 = System(
    Variables(
        plx = gaia_plx(gaia_id=2431157720981843200),
        M = TruncatedNormal(1.2, 0.1, 0, Inf),
        v_pma = Normal(0, 100),
        θ_pma = Uniform(-π, π),
        rv = Normal(0, 100),
        jitter = Gamma(2, 2),
        pmra = sys -> sys.v_pma*sin(sys.θ_pma),
        pmdec = sys -> sys.v_pma*cos(sys.θ_pma)
    ),
    # ProperMotionAnomHGCA(gaia_id=2431157720981843200),
    RadialVelocity(
        (epoch=58737.342, rv=-30, σ_rv=110),
        (epoch=58766.250, rv=-10, σ_rv=90),
        (epoch=58771.241, rv=-30, σ_rv=50),
        (epoch=58819.119, rv=-40, σ_rv=90),
        (epoch=58820.107, rv=10, σ_rv=3),
        (epoch=58821.096, rv=-10, σ_rv=80),
        (epoch=58839.050, rv=10, σ_rv=60),
        (epoch=59063.443, rv=140, σ_rv=90),
        (epoch=59072.417, rv=20, σ_rv=60),
        (epoch=59085.378, rv=-45, σ_rv=11),
        (epoch=59129.254, rv=-40, σ_rv=40),
        (epoch=59145.232, rv=10, σ_rv=80),
        (epoch=59159.164, rv=20, σ_rv=100)
    ),
    b
)

# ----------------------------------------------------------------------------------------------------------------------
# Calibrate System
# ----------------------------------------------------------------------------------------------------------------------

## Calibrate
chainparams = Dict(
    "acceptance" => 0.65,
    "num_chains" => 1,
    "adaptation" => 1000,
    "iterations" => 5000,
    "tree_depth" => 13,
    "verbosity" => 2
)

# Put relevant path here
saveas = "C:/Users/Jensen/Documents/undergrad/coop/haa/sbc/accept=0.650/test_output/test"

calibrate(HD984, chainparams, saveas)

# ----------------------------------------------------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------------------------------------------------

## Get chain parameters
paramsdict = open(saveas * "_chain_params.json", "r") do f
    JSON.parse(f)
end
display(paramsdict)

## Get sampled values 
sampleddict = open(saveas * "_sampled_vals.json", "r") do f
    JSON.parse(f)
end
display(sampleddict)

## Get CDF values 
cdfdict = open(saveas * "_cdf_vals.json", "r") do f
    JSON.parse(f)
end
display(cdfdict)

## Get chains 
chains = load(saveas * "_chains.jld2", "chains")
display(chains)

## Model plot
plotmodel(chains, alpha=0.1)

## Corner plot: just the highlights
table = (
    M_pri = chains["M"],
    M_sec = chains["b[mass]"],
    a = chains["b[a]"],
    e = chains["b[e]"],
    i = rad2deg.(chains["b[i]"]),
)
corner(table)

## Full corner plot
corner(chains)

## Parameters against time
timeplotgrid(chains)

## Trace plot
plot(chains, seriestype=:traceplot)

# ----------------------------------------------------------------------------------------------------------------------