
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
        arr2nt = Octofitter.make_arr2nt(system)

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

        # --- Build model configuration structs ---

        n_planets = length(system.planets)

        # Compute per-term active parameter indices for gradient sparsity
        active_indices_all, _D_check = _compute_active_indices(system)
        @assert _D_check == D "Active index computation disagrees on D: $_D_check vs $D"

        # Build orbit constructor: a RuntimeGeneratedFunction with orbit types
        # and planet indices interpolated as compile-time constants.
        construct_orbits = _make_construct_orbits(system)

        # Extract per-planet Kepler solvers. Specifying e.g. Markley() avoids
        # Enzyme having to shadow unused solver branches (like the e>=1 Roots.jl fallback).
        planet_solvers = ntuple(i -> system.planets[i].solver, Val(n_planets))

        mcfg = ModelEvalConfig(Bijector_invlinkvec, arr2nt, construct_orbits, planet_solvers, Val(n_planets))

        # Pre-compute orbit solution types for buffer pre-allocation.
        _θ_nat_0 = Bijector_invlinkvec(initial_θ_0_t)
        _θ_sys_0 = arr2nt(_θ_nat_0)
        _orbits_0 = construct_orbits(_θ_sys_0)

        # --- Build per-term configs ---
        tcfgs = ()
        term_idx = 0

        # Planet observation configs
        for i in 1:n_planets
            planet = system.planets[i]
            for (i_like, obs) in enumerate(planet.observations)
                term_idx += 1
                active = active_indices_all[term_idx]
                obs_name = hasproperty(obs, :name) ? normalizename(likelihoodname(obs)) : nothing
                epochs = _get_obs_epochs(obs)
                tcfg = TermEvalConfig(obs, obs_name, epochs, active, i)
                tcfgs = (tcfgs..., tcfg)
            end
        end

        # System observation configs
        for (i_obs, obs) in enumerate(system.observations)
            term_idx += 1
            active = active_indices_all[term_idx]
            obs_name = normalizename(likelihoodname(obs))
            epochs = _get_obs_epochs(obs)
            tcfg = TermEvalConfig(obs, obs_name, epochs, active, 0)
            tcfgs = (tcfgs..., tcfg)
        end

        # --- Build primal workspace and callback ---
        primal_tws = _build_term_workspaces(mcfg, tcfgs, _θ_sys_0, Float64, D)

        # Test primal evaluation
        ℓπcallback = let tcfgs=tcfgs, primal_tws=primal_tws,
                        invlink=Bijector_invlinkvec, arr2nt=arr2nt,
                        ln_prior_transformed=ln_prior_transformed,
                        construct_orbits=construct_orbits,
                        planet_solvers=planet_solvers
            @inline function(θ_transformed)
                lpost = zero(eltype(θ_transformed))
                @inbounds for k in eachindex(θ_transformed)
                    isfinite(θ_transformed[k]) || return convert(eltype(θ_transformed), -Inf)
                end
                θ_natural = invlink(θ_transformed)
                θ_structured = arr2nt(θ_natural)
                lpost += ln_prior_transformed(θ_natural, true)
                if !isfinite(lpost)
                    return lpost
                end
                # Sum term likelihoods using pre-allocated workspaces
                orbits = construct_orbits(θ_structured)
                lpost += _sum_term_likelihoods(θ_structured, orbits, tcfgs, primal_tws, planet_solvers)
                return lpost
            end
        end

        # Test likelihood function immediately
        if verbosity >= 1
            (function(ℓπcallback, θ)
                @showtime ℓπcallback(θ)
            end)(ℓπcallback, initial_θ_0_t)
        end

        if contains_discrete_variables
            # Return fully concrete type wrapping all these functions
            return new{
                D,
                typeof(ℓπcallback),
                Nothing,
                typeof(system),
                typeof(Bijector_linkvec),
                typeof(Bijector_invlinkvec),
                typeof(arr2nt),
                typeof(sample_priors),
                autodiff
            }(
                D,
                ℓπcallback,
                nothing,
                system,
                Bijector_linkvec,
                Bijector_invlinkvec,
                arr2nt,
                sample_priors,
                nothing
            )
        end

        # --- Build gradient workspace and callback ---

        # Resolve per-term backends, evaluators, and workspaces.
        # Evaluators hold only mcfg (with zero-size RuntimeGeneratedFunctions) so
        # Enzyme can prove them readonly for Const annotation. tcfg (with heap data)
        # is passed via DI.Constant, workspace via DI.Cache (→ Duplicated).
        # Iterates over tcfgs (built above) to guarantee identical term ordering.
        term_backends = ()
        term_evaluators = ()
        term_workspaces = ()
        for term_idx in 1:length(tcfgs)
            tcfg_i = tcfgs[term_idx]
            active = tcfg_i.active_indices
            D_active = length(active)
            obs_backend = resolve_ad_backend(tcfg_i.obs, autodiff, D_active)
            all_active = D_active == D

            tw = _make_term_workspace(mcfg, tcfg_i, _θ_sys_0, Float64, D)

            if all_active
                evaluator = TermEvaluator(mcfg)
            else
                evaluator = SparseTermEvaluator(mcfg)
            end

            term_evaluators = (term_evaluators..., evaluator)
            term_backends = (term_backends..., obs_backend)
            term_workspaces = (term_workspaces..., tw)
        end

        # Prior evaluator and gradient prep
        prior_evaluator = PriorEvaluator(Bijector_invlinkvec, ln_prior_transformed)
        # Prior evaluator has no mutable state → Const annotation is fine
        prior_backend = isnothing(autodiff) ? AutoForwardDiff() : autodiff
        θ_zero = zero(initial_θ_0_t)
        prior_prep = prepare_gradient(prior_evaluator, prior_backend, θ_zero)
        prior_grad = similar(initial_θ_0_t)
        total_grad = similar(initial_θ_0_t)
        θ_full_base = Vector{Float64}(undef, D)

        # Prepare gradient infrastructure: type-stable heterogeneous tuple
        # DI context wrappers (Constant, Cache) are pre-created inside each TermGradSpec
        # and reused every gradient call — θ_full_base is shared with the closure below.
        term_specs = _prepare_term_specs(term_evaluators, term_backends, active_indices_all,
                                         term_workspaces, tcfgs, initial_θ_0_t, θ_full_base)

        # Wrap term_specs in a Ref so the closure captures an 8-byte pointer instead
        # of storing the large heterogeneous Tuple inline (~2400 bytes). Without this,
        # sizeof(closure) ≈ 2440 bytes and Julia heap-allocates a copy on every call.
        term_specs_ref = Ref(term_specs)

        ∇ℓπcallback =
        let prior_evaluator=prior_evaluator, prior_backend=prior_backend,
            prior_prep=prior_prep, prior_grad=prior_grad,
            total_grad=total_grad, term_specs_ref=term_specs_ref,
            θ_full_base=θ_full_base
            function(θ_transformed)
                copyto!(θ_full_base, θ_transformed)

                # Prior gradient (full-D — prior touches all params)
                ll_p, _ = value_and_gradient!(prior_evaluator, prior_grad, prior_prep, prior_backend, θ_transformed)

                if !isfinite(ll_p)
                    fill!(total_grad, zero(eltype(θ_transformed)))
                    return ll_p, total_grad
                end

                copyto!(total_grad, prior_grad)

                # Per-term sparse gradients (type-stable recursion over heterogeneous tuple)
                ll = _accumulate_term_gradients!(total_grad, ll_p, term_specs_ref[], θ_transformed)

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

        # Perform some quick diagnostic checks to warn users for performance-gotchas
        out_type_model = Core.Compiler.return_type(ℓπcallback, typeof((initial_θ_0_t,)))
        out_type_model_grad = Core.Compiler.return_type(∇ℓπcallback, typeof((initial_θ_0_t,)))
        out_type_arr2nt = Core.Compiler.return_type(arr2nt, typeof((initial_θ_0_t,)))
        out_type_prior = Core.Compiler.return_type(ln_prior_transformed, typeof((initial_θ_0,false,)))

        if isconcretetype(out_type_prior) &&
            isconcretetype(out_type_arr2nt) &&
            !isconcretetype(out_type_model)
            @warn "\nThis model's log density function is not type stable, but all of its components are. \nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub." out_type_prior out_type_arr2nt out_type_model
        end
        if !isconcretetype(out_type_model_grad)
            @warn "\nThis model's log density gradient is not type stable, but all of its components are. \nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub."
            end
        if !isconcretetype(out_type_prior)
            @warn "\nThis model's prior sampler does not appear to be type stable, which will likely hurt sampling performance.\nThis may indicate a performance bug in Octofitter; please consider filing an issue on GitHub." out_type_prior out_type_arr2nt out_type_model
        end
        if !isconcretetype(out_type_arr2nt)
            @warn "\nThis model specification (arr2nt) is not type-stable, which will likely hurt sampling performance.\nCheck for global variables used within your model definition, and prepend these with `\$`.\nIf that doesn't work, you could trying running:\n`Cthulhu.@descend model.arr2nt(model.sample_priors(Random.Xoshiro(0)))` for more information.\nFor assistance, please file an issue on GitHub." out_type_prior out_type_arr2nt out_type_model
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
