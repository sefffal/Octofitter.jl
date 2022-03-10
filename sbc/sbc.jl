# ----------------------------------------------------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------------------------------------------------

## General
using Plots
using StatsBase
using Distributions
using StaticArrays
using TypedTables
using LinearAlgebra

## Orbit fitting
using DirectOrbits
using DirectDetections

# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------

## Draw from priors
function drawfrompriors(system::System)
    θ = DirectDetections.sample_priors(system)
    arr2nt = DirectDetections.make_arr2nt(system)
    θnt = arr2nt(θ)
    return θnt
end

## New data
# Generate new astrometry observations
function newobs(obs::Astrometry, elem::KeplerianElements)
    epochs = obs.table.epoch
    σ_ras = obs.table.σ_ra 
    σ_decs = obs.table.σ_dec
    ras = DirectOrbits.raoff.(elem, epochs)
    decs = DirectOrbits.decoff.(elem, epochs)
    astrometry_table = Table(epoch=epochs, ra=ras, dec=decs, σ_ra=σ_ras, σ_dec=σ_decs)
    return Astrometry(astrometry_table)
end

# Generate new radial velocity observations for a planet
function newobs(obs::RadialVelocity, elem::KeplerianElements)
    epochs = obs.table.epoch 
    σ_rvs = obs.table.σ_rv 
    rvs = DirectOribts.radvel.(elem, epochs)
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)
    return RadialVelocity(radvel_table)
end

# Generate new radial velocity observations for a star
function newobs(obs::RadialVelocity, elems::Vector{<:KeplerianElements}, θ_system)
    epochs = obs.table.epoch 
    σ_rvs = obs.table.σ_rv 
    planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun
    rvs = DirectOrbits.radvel.(reshape(elems, :, 1), epochs, transpose(planet_masses))
    noise = randn(length(epochs)) .* θ_system.jitter
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv .+ noise
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)
    return RadialVelocity(radvel_table)
end

# Generate new images
# Need images function

## Generate calibration data
function calibrationsystem(system::System)
    θ_newsystem = drawfrompriors(system)

    elements = map(eachindex(system.planets)) do i
        planet = system.planets[i]
        θ_newplanet = θ_newsystem.planets[i]

        if (hasproperty(θ_newplanet, :a) && θ_newplanet.a <= 0) ||
            (hasproperty(θ_newplanet, :e) && !(0 <= θ_newplanet.e < 1))
            out_of_bounds[] = true
        end

        neworbit = DirectDetections.construct_elements(DirectDetections.orbittype(planet), θ_newsystem, θ_newplanet)

        return neworbit
    end

    newplanets = map(1:length(system.planets)) do i
        planet = system.planets[i]
        elem = elements[i]

        newplanet_obs = map(planet.observations) do obs
            return newobs(obs, elem)
        end
        newplanet = Planet{DirectDetections.orbittype(planet)}(planet.priors, planet.derived, newplanet_obs..., name=planet.name)
    end

    newstar_obs = map(system.observations) do obs
        return newobs(obs, collect(elements), θ_newsystem)
    end

    newsystem = System(system.priors, system.derived, newstar_obs..., newplanets..., name=system.name)

    return θ_newsystem, newsystem
end

## Run chains on new model
function calibrationhmc(system::System; accept=0.65, adapt=1000, iter=25000, td=13, v=2)
    θ_newsystem, newsystem = calibrationsystem(system)
    θ_array = DirectDetections.result2mcmcchain(newsystem, [θ_newsystem])

    @time chains = DirectDetections.hmc(
        newsystem, accept,
        adaptation = adapt,
        iterations = iter,
        tree_depth = td,
        verbosity = v
    )

    chainkeys = string.(keys(chains))
    cdfdict = Dict()

    for key in chainkeys
        paramdist = vec(chains[key])
        paramcdf = ecdf(paramdist)
        trueval = θ_array[key][1]
        cdfval = paramcdf(trueval)
        cdfdict[key] = cdfval
    end

    return chains, cdfdict
end

# ----------------------------------------------------------------------------------------------------------------------
# Test
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

## Draw from priors
chains, cdfdict = calibrationhmc(
    HD984,
    accept = 0.80,
    adapt = 2500,
    iter = 10000
)

## Results
display(cdfdict)