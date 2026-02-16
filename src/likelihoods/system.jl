using Bumper


# Helper to generate orbit-solving expressions for a single observation term.
# Returns (solve_exprs, sol_keys) where sol_keys are symbols for per-planet solution arrays.
function _make_per_term_solve_exprs(term_id::String, like, planet_keys::Vector{Symbol}, n_planets::Int)
    epochs = if hasproperty(like, :table) && hasproperty(like.table, :epoch)
        collect(Float64, like.table.epoch)
    else
        Float64[]
    end
    n_epochs = length(epochs)

    sol_exprs = Expr[]
    sol_keys = Symbol[]

    if n_epochs == 0 || n_planets == 0
        # No epochs or no planets: empty solution tuples
        for ip in 1:n_planets
            sk = Symbol("sols_$(term_id)_p$(ip)")
            push!(sol_keys, sk)
            push!(sol_exprs, :($sk = ()))
        end
    else
        # Allocate epochs array (shared by all planets for this term)
        epochs_sym = Symbol("epochs_$(term_id)")
        push!(sol_exprs, quote
            $epochs_sym = @alloc(Float64, $n_epochs)
            $((:($(epochs_sym)[$je] = $(epochs[je])) for je in 1:n_epochs)...)
        end)

        # For each planet, solve at this term's epochs
        for ip in 1:n_planets
            sk = Symbol("sols_$(term_id)_p$(ip)")
            push!(sol_keys, sk)
            pk = planet_keys[ip]
            push!(sol_exprs, quote
                let _orbit = $pk, _eps = $epochs_sym
                    _sol0 = orbitsolve(_orbit, _eps[1])
                    $sk = @alloc(typeof(_sol0), $n_epochs)
                    $(sk)[1] = _sol0
                    for _j in 2:$n_epochs
                        $(sk)[_j] = orbitsolve(_orbit, _eps[_j])
                    end
                end
            end)
        end
    end

    return sol_exprs, sol_keys, n_epochs
end

function make_ln_like(system::System, θ_system)

    # --- Planet orbit construction (shared across all terms) ---
    planet_keys = Symbol[]
    planet_construction_exprs = Expr[]
    planet_declarations = Expr[]
    n_planets = length(system.planets)

    for i in 1:n_planets
        planet = system.planets[i]
        OrbitType = _planet_orbit_type(planet)
        key = Symbol("planet_$i")
        push!(planet_declarations, :($key = nothing))
        push!(planet_keys, key)
        push!(planet_construction_exprs, quote
            $key = $(OrbitType)(;merge(θ_system, θ_system.planets[$i])...)
        end)
    end

    # --- Per-term evaluation expressions ---
    # Each observation term gets its own @no_escape block that solves orbits
    # at only its own epochs. This enables per-term differentiation in Stage 5.
    j = 0  # running ll variable counter
    term_exprs = Expr[]

    # Planet observations
    for i in 1:n_planets
        planet = system.planets[i]
        for (i_like, like) in enumerate(planet.observations)
            term_id = "p$(i)_o$(i_like)"
            sol_exprs, sol_keys, n_epochs = _make_per_term_solve_exprs(term_id, like, planet_keys, n_planets)
            solutions_tuple = :(tuple($(sol_keys...)))

            # Build context expression
            obs_name = hasproperty(like, :name) ? normalizename(likelihoodname(like)) : nothing
            if !isnothing(obs_name)
                ctx_expr = :(PlanetObservationContext(
                    θ_system,
                    θ_system.planets[$i],
                    hasproperty(θ_system.planets[$i].observations, $(Meta.quot(obs_name))) ?
                        θ_system.planets[$i].observations.$(obs_name) :
                        (;),
                    elems,
                    $solutions_tuple,
                    $i,
                ))
            else
                ctx_expr = :(PlanetObservationContext(
                    θ_system,
                    θ_system.planets[$i],
                    (;),
                    elems,
                    $solutions_tuple,
                    $i,
                ))
            end

            ll_call = :(ln_like(system.planets[$i].observations[$i_like], $ctx_expr))

            if n_epochs > 0
                push!(term_exprs, :(
                    $(Symbol("ll$(j+1)")) = $(Symbol("ll$j")) + @no_escape begin
                        $(sol_exprs...)
                        $ll_call
                    end
                ))
            else
                push!(term_exprs, :(
                    $(Symbol("ll$(j+1)")) = $(Symbol("ll$j")) + $ll_call
                ))
            end
            j += 1
        end
    end

    # System observations
    for (i_obs, like) in enumerate(system.observations)
        term_id = "sys_o$(i_obs)"
        sol_exprs, sol_keys, n_epochs = _make_per_term_solve_exprs(term_id, like, planet_keys, n_planets)
        solutions_tuple = :(tuple($(sol_keys...)))

        obs_name = normalizename(likelihoodname(like))
        ctx_expr = :(SystemObservationContext(
            θ_system,
            hasproperty(θ_system.observations, $(Meta.quot(obs_name))) ?
                θ_system.observations.$(obs_name) :
                (;),
            elems,
            $solutions_tuple,
        ))

        ll_call = :(ln_like(system.observations[$i_obs], $ctx_expr))

        if n_epochs > 0
            push!(term_exprs, :(
                $(Symbol("ll$(j+1)")) = $(Symbol("ll$j")) + @no_escape begin
                    $(sol_exprs...)
                    $ll_call
                end
            ))
        else
            push!(term_exprs, :(
                $(Symbol("ll$(j+1)")) = $(Symbol("ll$j")) + $ll_call
            ))
        end
        j += 1
    end

    return @RuntimeGeneratedFunction(:(function (system::System, θ_system)
        T = _system_number_type(θ_system)
        ll0 = zero(T)

        # Declare all planet variables before the try block
        $(planet_declarations...)

        # Try-catch only for planet construction
        try
            $(planet_construction_exprs...)
        catch err
            @warn "Failed to constructor orbit:" exception=(err, catch_backtrace())
            return convert(T, -Inf)
        end

        # Construct a tuple of existing planet orbital elements
        elems = tuple($(planet_keys...))

        # Evaluate each term with per-term orbit solving
        $(term_exprs...)

        $(Symbol("ll$j"))
    end))
end

# Threading flag used by sampling/initialization for outer-loop parallelism
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

