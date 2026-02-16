
using LogDensityProblems

# Define the target distribution using the `LogDensityProblem` interface
mutable struct LogDensityModel{D,Tℓπ,T∇ℓπ,TSys,TLink,TInvLink,TArr2nt,TPriSamp,ADType}
    # Dimensionality
    const D::Int
    # Calculate the log-posterior density given transformed parameters
    const ℓπcallback::Tℓπ
    # Calculate the log-posterior density and gradient given transformed parameters
    const ∇ℓπcallback::T∇ℓπ
    # The underlying System object
    const system::TSys
    # Convert flat parameter vector into transformed domain
    const link::TLink
    # Convert flat transformed parameter vector back to natural domain
    const invlink::TInvLink
    # Convert a flat parameter vector into a nested named tuple structure,
    # matching the variable definitions in the System and Planet blocks
    const arr2nt::TArr2nt
    # Sample IID from the model's priors
    const sample_priors::TPriSamp
    # A set of starting points that can be sampled from to initialize a sampler, or nothing
    starting_points::Union{Nothing,Vector} 
    function LogDensityModel(system::System; autodiff=nothing, verbosity=2, chunk_sizes=nothing)
        verbosity >= 1 && @info "Preparing model"

        sample_priors = make_prior_sampler(system)

        # Choose parameter dimensionality and initial parameter value
        initial_θ_0 = sample_priors(Random.default_rng())
        D = length(initial_θ_0)
        verbosity >= 2 && @info "Determined number of free variables" D

        # We support models with discrete or mixed variables, but in these
        # cases we can't support autodiff.
        # Detect this case, warn the user, and skip over defining ∇ℓπcallback
        contains_discrete_variables = autodiff === false || any(isa.(sample_priors(Random.default_rng(), system),Integer))
        if contains_discrete_variables && verbosity >= 1
            @info "Model contains discrete variables; model gradients not supported."
        end

        # autodiff is kept as-is (nothing, false, or a specific backend).
        # Per-term backend resolution happens in resolve_ad_backend().

        ln_prior_transformed = make_ln_prior_transformed(system)
        # ln_prior = make_ln_prior(system)
        arr2nt = Octofitter.make_arr2nt(system) 
        ln_like_generated = make_ln_like(system, arr2nt(initial_θ_0))

        # Check number type
        θ_nt = arr2nt(initial_θ_0)
        T = _system_number_type(θ_nt)
        verbosity >= 2 && @info "Determined number type" T
        if !(T <: Real)
            error("Error: inferred that you wanted to use $(T) as the number type, which is not supported. It must be a floating point number or similar. Check that all the variables you provided in your model are promotable to a float (e.g. not `nothing` or `missing`)")
        end


        priors_vec = _list_priors(system)
        Bijector_invlinkvec = make_Bijector_invlinkvec(priors_vec)
        # TODO : this could be unrolled like invlink if used anywhere performance sensitive. Currently
        # it just transforms things back after all sampling is done.
        Bijector_linkvec = let priors_vec=priors_vec
            # (θ) -> Bijectors.link.(priors_vec, θ)
            function (θ)
                i = 0
                out = zeros(eltype(θ), length(θ))
                for prior in priors_vec
                    if length(prior) == 1
                        i += 1
                        out[i] = Bijectors.link(prior, θ[i])
                    else
                        param_outs = Bijectors.link(prior, θ[i+1:i+length(prior)])
                        for param_out in param_outs
                            i += 1
                            out[i] = param_out
                        end
                    end
                end
                return out
            end
        end
        initial_θ_0_t = Bijector_linkvec(initial_θ_0)

        # Test out model likelihood and prior computations. This way, if they throw
        # an error, we'll see it right away instead of burried in some deep stack
        # trace from the sampler, autodiff, etc.
        ln_like_generated(system, arr2nt(initial_θ_0))

        ln_prior_transformed(initial_θ_0,false)


        # We use let blocks to prevent type instabilities from closures
        # A function barrier would also work.
        ℓπcallback, ∇ℓπcallback = (function(
            arr2nt,
            system,
            Bijector_invlinkvec,
            ln_prior_transformed,
            ln_like_generated,
            D,
        )

            # Capture these variables in a let binding to improve performance
            # We also set up temporary storage to reduce allocations
            # ForwardDiff is used to compute the likelihood gradients using the in-place 
            # API. This ensures type stability.
            function ℓπcallback(
                θ_transformed,
                system=system,
                arr2nt=arr2nt,
                Bijector_invlinkvec=Bijector_invlinkvec,
                ln_prior_transformed=ln_prior_transformed,
                ln_like_generated=ln_like_generated;sampled=true)#::eltype(θ_transformed)

                lpost = zero(eltype(θ_transformed))
                # Stop right away if we are given non-finite arguments
                if any(!isfinite, θ_transformed)
                    # @warn "non finite parameters encountered (maxlog=1)" θ_transformed maxlog=1
                    lpost = convert(eltype(θ_transformed), -Inf)
                    return lpost
                end
                # Transform back from the unconstrained support to constrained support for the likelihood function
                θ_natural = @inline Bijector_invlinkvec(θ_transformed)
                θ_structured = @inline arr2nt(θ_natural)
                lpost += @inline ln_prior_transformed(θ_natural,sampled)
                # Don't compute likelihood if we fell outside the prior bounds
                if !isfinite(lpost)
                    # @warn "non finite log prior (maxlog=1)" lpost maxlog=1
                    return lpost
                end
                lpost += @inline ln_like_generated(system, θ_structured)
                # if !isfinite(lpost)
                #     # TODO: check for performance impact here
                #     # Display parameters that caused an invalid log-likelihood to be calculated
                #     # Strip off any forward diff Dual tags, as these make it impossible to see
                #     # what's going on.
                #     θ_transformed_primals = ForwardDiff.value.(θ_transformed)
                #     θ_structured = arr2nt(Bijector_invlinkvec(θ_transformed_primals))
                #     llike = ln_like_generated(system, θ_structured)
                #     @warn "Invalid log likelihood encountered. (maxlog=1)" θ=θ_structured llike θ_transformed=θ_transformed_primals
                # end
                return lpost
            end

            args = (
                system,
                arr2nt,
                Bijector_invlinkvec,
                ln_prior_transformed,
                ln_like_generated,
            )
            # Test likelihood function immediately to give user a clean error
            # if it fails for some reason.
            # Also Display their run time. If something is egregiously wrong we'll notice something
            # here in the output logs.
            if verbosity >= 1
                (function(ℓπcallback, θ)
                    @showtime ℓπcallback(θ)
                end)(ℓπcallback, initial_θ_0_t)
            end

            if contains_discrete_variables
                return ℓπcallback, nothing
            end

            # --- Per-term gradient composition ---
            # Each observation term gets its own AD closure so different terms
            # can use different backends (e.g. ForwardDiff vs Enzyme).

            n_planets = length(system.planets)

            # Build orbit constructors: one closure per planet capturing its OrbitType
            orbit_constructors = ntuple(n_planets) do i
                OT = _planet_orbit_type(system.planets[i])
                let OT=OT, i=i
                    (θ_system) -> OT(;merge(θ_system, θ_system.planets[i])...)
                end
            end

            # Helper: construct all orbits from structured parameters
            function _build_orbits(orbit_constructors, θ_system, n_planets)
                ntuple(i -> orbit_constructors[i](θ_system), n_planets)
            end

            # Prior closure: θ_transformed → scalar prior log-density
            prior_closure = let invlink=Bijector_invlinkvec, ln_prior=ln_prior_transformed
                function(θ_transformed)
                    if any(!isfinite, θ_transformed)
                        return convert(eltype(θ_transformed), -Inf)
                    end
                    θ_natural = invlink(θ_transformed)
                    ln_prior(θ_natural, true)
                end
            end

            # Resolve prior backend (prior has no ad_backend, always use default)
            prior_backend = isnothing(autodiff) ? AutoForwardDiff(chunksize=D) : autodiff

            # Build per-term closures and resolve backends
            term_closures = Any[]
            term_backends = Any[]

            # Planet observation closures
            for i in 1:n_planets
                planet = system.planets[i]
                for (i_like, obs) in enumerate(planet.observations)
                    obs_name = hasproperty(obs, :name) ? normalizename(likelihoodname(obs)) : nothing
                    epochs = _get_obs_epochs(obs)

                    closure = let i=i, obs=obs, obs_name=obs_name, epochs=epochs,
                                  invlink=Bijector_invlinkvec, arr2nt=arr2nt,
                                  orbit_constructors=orbit_constructors, n_planets=n_planets
                        function(θ_transformed)
                            θ_natural = invlink(θ_transformed)
                            θ_system = arr2nt(θ_natural)
                            local orbits
                            try
                                orbits = _build_orbits(orbit_constructors, θ_system, n_planets)
                            catch
                                return convert(eltype(θ_transformed), -Inf)
                            end
                            solutions = _solve_all_orbits(orbits, epochs)
                            θ_planet = θ_system.planets[i]
                            θ_obs = if obs_name !== nothing && hasproperty(θ_planet.observations, obs_name)
                                getproperty(θ_planet.observations, obs_name)
                            else
                                (;)
                            end
                            ctx = PlanetObservationContext(θ_system, θ_planet, θ_obs, orbits, solutions, i)
                            ln_like(obs, ctx)
                        end
                    end
                    push!(term_closures, closure)
                    push!(term_backends, resolve_ad_backend(obs, autodiff, D))
                end
            end

            # System observation closures
            for (i_obs, obs) in enumerate(system.observations)
                obs_name = normalizename(likelihoodname(obs))
                epochs = _get_obs_epochs(obs)

                closure = let obs=obs, obs_name=obs_name, epochs=epochs,
                              invlink=Bijector_invlinkvec, arr2nt=arr2nt,
                              orbit_constructors=orbit_constructors, n_planets=n_planets
                    function(θ_transformed)
                        θ_natural = invlink(θ_transformed)
                        θ_system = arr2nt(θ_natural)
                        local orbits
                        try
                            orbits = _build_orbits(orbit_constructors, θ_system, n_planets)
                        catch
                            return convert(eltype(θ_transformed), -Inf)
                        end
                        solutions = _solve_all_orbits(orbits, epochs)
                        θ_obs = if hasproperty(θ_system.observations, obs_name)
                            getproperty(θ_system.observations, obs_name)
                        else
                            (;)
                        end
                        ctx = SystemObservationContext(θ_system, θ_obs, orbits, solutions)
                        ln_like(obs, ctx)
                    end
                end
                push!(term_closures, closure)
                push!(term_backends, resolve_ad_backend(obs, autodiff, D))
            end

            # Prepare gradient infrastructure for each term
            prior_prep = prepare_gradient(prior_closure, prior_backend, zero(initial_θ_0_t))
            prior_grad = similar(initial_θ_0_t)

            n_terms = length(term_closures)
            term_preps = [prepare_gradient(term_closures[i], term_backends[i], zero(initial_θ_0_t)) for i in 1:n_terms]
            term_grads = [similar(initial_θ_0_t) for _ in 1:n_terms]
            total_grad = similar(initial_θ_0_t)

            ∇ℓπcallback =
            let prior_closure=prior_closure, prior_backend=prior_backend,
                prior_prep=prior_prep, prior_grad=prior_grad,
                term_closures=term_closures, term_backends=term_backends,
                term_preps=term_preps, term_grads=term_grads,
                total_grad=total_grad, n_terms=n_terms
                function(θ_transformed)
                    # Prior gradient
                    ll_p, _ = value_and_gradient!(prior_closure, prior_grad, prior_prep, prior_backend, θ_transformed)

                    if !isfinite(ll_p)
                        fill!(total_grad, zero(eltype(θ_transformed)))
                        return ll_p, total_grad
                    end

                    ll = ll_p
                    copyto!(total_grad, prior_grad)

                    # Per-term gradients
                    for i in 1:n_terms
                        ll_i, _ = value_and_gradient!(term_closures[i], term_grads[i], term_preps[i], term_backends[i], θ_transformed)
                        ll += ll_i
                        total_grad .+= term_grads[i]
                    end

                    return ll, total_grad
                end
            end

            # Run the callback once right away. If there is a coding error in the users
            # model, we want to surface it ASAP.
            if verbosity >= 1
                (function(∇ℓπcallback, θ)
                    @showtime ∇ℓπcallback(θ)
                end)(∇ℓπcallback, initial_θ_0_t)
            else
                ∇ℓπcallback(initial_θ_0_t)
            end


            ℓπcallback, ∇ℓπcallback
        end)(
            arr2nt,
            system,
            Bijector_invlinkvec,
            ln_prior_transformed,
            ln_like_generated,
            D
        )
        
        # Perform some quick diagnostic checks to warn users for performance-gtochas
        out_type_model = Core.Compiler.return_type(ℓπcallback, typeof((initial_θ_0_t,)))
        if !contains_discrete_variables
            out_type_model_grad = Core.Compiler.return_type(∇ℓπcallback, typeof((initial_θ_0_t,)))
        end
        out_type_arr2nt = Core.Compiler.return_type(arr2nt, typeof((initial_θ_0_t,)))
        out_type_prior = Core.Compiler.return_type(ln_prior_transformed, typeof((initial_θ_0,false,)))
        out_type_like = Core.Compiler.return_type(ln_like_generated, typeof((system,arr2nt(initial_θ_0),)))

        if isconcretetype(out_type_prior) &&
            isconcretetype(out_type_like) &&
            isconcretetype(out_type_arr2nt) &&
            !isconcretetype(out_type_model)
            @warn "\nThis model's log density function is not type stable, but all of its components are. \nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub." out_type_prior out_type_like out_type_arr2nt out_type_model
        end
        if !contains_discrete_variables && !isconcretetype(out_type_model_grad)
            @warn "\nThis model's log density gradient is not type stable, but all of its components are. \nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub."
            end
        if !isconcretetype(out_type_prior)
            @warn "\nThis model's prior sampler does not appear to be type stable, which will likely hurt sampling performance.\nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub." out_type_prior out_type_like out_type_arr2nt out_type_model
        end
        if !isconcretetype(out_type_like)
            @warn "\nThis model's likelihood function does not appear to be type stable, which will likely hurt sampling performance.\nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub." out_type_prior out_type_like out_type_arr2nt out_type_model
        end
        if !isconcretetype(out_type_arr2nt)
            @warn "\nThis model specification (arr2nt) is not type-stable, which will likely hurt sampling performance.\nCheck for global variables used within your model definition, and prepend these with `\$`.\nIf that doesn't work, you could trying running:\n`Cthulhu.@descend model.arr2nt(model.sample_priors(Random.Xoshiro(0)))` for more information.\nFor assistance, please file an issue on GitHub." out_type_prior out_type_like out_type_arr2nt out_type_model
        end

        # Return fully concrete type wrapping all these functions
        new{
            D,
            typeof(ℓπcallback),
            typeof(∇ℓπcallback),
            typeof(system),
            typeof(Bijector_linkvec),
            typeof(Bijector_invlinkvec),
            typeof(arr2nt),
            typeof(sample_priors),
            autodiff # an ADType
        }(
            D,
            ℓπcallback,
            ∇ℓπcallback,
            system,
            Bijector_linkvec,
            Bijector_invlinkvec,
            arr2nt,
            sample_priors,
            nothing # no starting points set
        )
    end
end
LogDensityProblems.logdensity(p::LogDensityModel, θ) = p.ℓπcallback(θ)
LogDensityProblems.logdensity_and_gradient(p::LogDensityModel, θ) = p.∇ℓπcallback(θ)
LogDensityProblems.dimension(p::LogDensityModel{D}) where D = D
LogDensityProblems.capabilities(::Type{<:LogDensityModel{D,Tℓπ,Nothing}}) where {D,Tℓπ} = LogDensityProblems.LogDensityOrder{0}()
LogDensityProblems.capabilities(::Type{<:LogDensityModel{D,Tℓπ,T∇ℓπ}}) where {D,Tℓπ,T∇ℓπ} = LogDensityProblems.LogDensityOrder{1}()

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::LogDensityModel)
    L = _count_epochs(p.system)
    println(io, "LogDensityModel for System $(p.system.name) of dimension $(p.D) and $(L) epochs with fields .ℓπcallback and .∇ℓπcallback")
end
function Base.show(io::IO, @nospecialize p::LogDensityModel)
    L = _count_epochs(p.system)
    println(io, "LogDensityModel for System $(p.system.name) of dimension $(p.D) and $(L) epochs with fields .ℓπcallback and .∇ℓπcallback")
end
