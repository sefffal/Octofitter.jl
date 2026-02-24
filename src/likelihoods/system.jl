
# Threading coordination flag for initialization parallelism
const _kepsolve_use_threads = Ref(false)

# Generate calibration data
"""
    generate_from_params(system, θ=drawfrompriors(system); add_noise=false)

Generate a new system and observations from an existing model.
"""
function generate_from_params(system::System, θ_newsystem = drawfrompriors(system); add_noise=false)

    # Generate new orbits for each planet in the system
    orbits = map(eachindex(system.planets)) do i
        planet = system.planets[i]
        neworbit = Octofitter.construct_elements(θ_newsystem, planet.name)
        return neworbit
    end

    # Helper: solve all orbits at a given set of epochs
    function _solve_at_epochs(obs)
        if hasproperty(obs, :table) && hasproperty(obs.table, :epoch)
            epochs = obs.table.epoch
            return ntuple(length(orbits)) do i
                [orbitsolve(orbits[i], ep) for ep in epochs]
            end
        else
            return ntuple(_ -> [], length(orbits))
        end
    end

    # Generate new observations for each planet in the system
    newplanets = map(1:length(system.planets)) do i
        planet = system.planets[i]
        orbit = orbits[i]
        θ_newplanet = θ_newsystem.planets[i]

        newplanet_obs = map(planet.observations) do obs
            # Get the observation-specific variables if they exist
            obs_name = normalizename(likelihoodname(obs))
            θ_obs = hasproperty(θ_newplanet.observations, obs_name) ?
                    getproperty(θ_newplanet.observations, obs_name) :
                    (;)

            # Solve all orbits at this observation's epochs
            orbit_solutions = _solve_at_epochs(obs)

            # Construct context object
            ctx = PlanetObservationContext(
                θ_newsystem,
                θ_newplanet,
                θ_obs,
                orbits,
                orbit_solutions,
                i,  # planet index
            )

            # Call with context
            return generate_from_params(obs, ctx; add_noise)
        end

        newplanet = Planet(
            variables=(planet.priors, planet.derived),
            basis=Octofitter.orbittype(planet),
            observations=newplanet_obs,
            name=planet.name
        )
        return newplanet
    end

    # Generate new observations for the star
    newstar_obs = map(system.observations) do obs
        # Get the observation-specific variables if they exist
        obs_name = normalizename(likelihoodname(obs))
        θ_obs = hasproperty(θ_newsystem.observations, obs_name) ?
                getproperty(θ_newsystem.observations, obs_name) :
                (;)

        # Solve all orbits at this observation's epochs
        orbit_solutions = _solve_at_epochs(obs)

        # Construct context object
        ctx = SystemObservationContext(
            θ_newsystem,
            θ_obs,
            orbits,
            orbit_solutions,
        )

        # Call with context
        return generate_from_params(obs, ctx; add_noise)
    end

    # Generate new system
    newsystem = System(
        variables=(system.priors, system.derived),
        observations=newstar_obs,
        companions=newplanets,
        name=system.name
    )

    return newsystem
end
export generate_from_params

"""
    make_ln_like(system, θ_system_example)

Build a log-likelihood function `(system, θ_system) -> ll::Float64` for post-processing
(chain diagnostics, SBC, cross-validation). Not used in the hot sampling path.
"""
function make_ln_like(system::System, θ_system_example)
    n_planets = length(system.planets)
    val_n = Val(n_planets)

    construct_orbits = _make_construct_orbits(system)

    let construct_orbits=construct_orbits, val_n=val_n
        function (system::System, θ_system)
            T = _system_number_type(θ_system)
            ll = zero(T)

            orbits = construct_orbits(θ_system)

            # Helper: solve orbits at an observation's epochs
            function _solve_at(obs)
                if hasproperty(obs, :table) && hasproperty(obs.table, :epoch)
                    epochs = obs.table.epoch
                    return ntuple(val_n) do i
                        [orbitsolve(orbits[i], ep) for ep in epochs]
                    end
                else
                    return ntuple(_ -> [], val_n)
                end
            end

            # Planet observations
            for i in 1:n_planets
                planet = system.planets[i]
                θ_planet = θ_system.planets[i]
                for obs in planet.observations
                    orbit_solutions = _solve_at(obs)
                    obs_name = hasproperty(obs, :name) ? normalizename(likelihoodname(obs)) : nothing
                    θ_obs = if !isnothing(obs_name) && hasproperty(θ_planet.observations, obs_name)
                        getproperty(θ_planet.observations, obs_name)
                    else
                        (;)
                    end
                    ctx = PlanetObservationContext(θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i)
                    ll += ln_like(obs, ctx)
                end
            end

            # System observations
            for obs in system.observations
                orbit_solutions = _solve_at(obs)
                obs_name = normalizename(likelihoodname(obs))
                θ_obs = hasproperty(θ_system.observations, obs_name) ?
                        getproperty(θ_system.observations, obs_name) : (;)
                ctx = SystemObservationContext(θ_system, θ_obs, orbits, orbit_solutions)
                ll += ln_like(obs, ctx)
            end

            return ll
        end
    end
end
