
using LogDensityProblems

# ---------------------------------------------------
# ForwardDiff chunk width
#
# Total gradient cost goes as ⌈D/N⌉ × c(N) for chunk width N. While c is
# roughly linear in N — one Dual multiply is a product rule, 2N+1 flops
# against 1 — fewer, wider chunks win, and a single chunk of width D wins
# outright. That holds further than ForwardDiff's own default assumes
# (`pickchunksize` caps at 12), but not forever: past ~24 partials the
# `Partials{N}` tuple stops living in registers and c(N) turns superlinear.
#
# Measured on this codebase (∇ℓπ, µs, N-companion relative astrometry):
#
#   D=23 | chunk  12:178.7   16:217.2   20:250.5   23:149.0   <- single chunk
#   D=30 | chunk  10:410.7   12:414.6   15:346.1   18:377.5   20:388.6
#        |        24:459.8   30:443.7               ^ best, and 1.28x under D
#   D=37 | chunk  10:782     13:672     19:569     37:809     ^ 1.42x under D
#
# The rule below — the *balanced* split into as few chunks as keeps every
# chunk at or under `CHUNK_WIDTH_MAX` — reproduces the measured optimum at all
# three: 23 -> 23, 30 -> 15, 37 -> 19. Balance matters as much as width; at
# D=30, an unbalanced 24+6 (459.8) is worse than 15+15 (346.1) even though 24
# is the more efficient width, because the 6-wide chunk pays nearly a full
# traversal for a quarter of the partials.
#
# `CHUNK_WIDTH_MAX` is a machine property (register file, vector width), not a
# mathematical one, and 24 is where it sits on the AVX2 hardware this was
# measured on. Users on very different hardware, or with a model whose shape
# moves the knee, can measure instead: `LogDensityModel(sys; chunk_sizes=[…])`
# times each candidate and picks the best.
# ---------------------------------------------------

const CHUNK_WIDTH_MAX = Ref(24)

"""
    _chunk_heuristic(D) -> Int

Chunk width for `D` free parameters: the balanced split into the fewest chunks
whose width does not exceed `CHUNK_WIDTH_MAX[]`. Reduces to a single chunk of
width `D` for the great majority of models.
"""
function _chunk_heuristic(D::Int)
    D <= 0 && return 1
    nmax = max(1, CHUNK_WIDTH_MAX[])
    D <= nmax && return D
    return cld(D, cld(D, nmax))
end

"""
    _probe_chunk_size(ℓπ, θ, D, candidates, verbosity) -> Int

Time a gradient at each candidate chunk width and return the fastest.

Opt-in (`LogDensityModel(sys; chunk_sizes=[…])`) rather than the default,
because each candidate compiles a fresh gradient at that `Dual` width — which
for a wide chunk is by far the most expensive thing model construction does.
The heuristic is free and matches the measured optimum on everything tried;
this is for hardware or model shapes where it might not.
"""
function _probe_chunk_size(ℓπ, θ, D::Int, candidates, verbosity::Int)
    cands = sort(unique(Int[c for c in candidates if 1 <= c <= D]))
    if isempty(cands)
        @warn "No candidate in `chunk_sizes` is between 1 and $D; using the default." candidates
        return _chunk_heuristic(D)
    end
    verbosity >= 1 && @info "Timing ForwardDiff chunk widths (this compiles one gradient per candidate)" candidates=cands
    best, best_t = first(cands), Inf
    for c in cands
        t = try
            ad = AutoForwardDiff(chunksize=c)
            prep = prepare_gradient(ℓπ, ad, zero(θ))
            grad = similar(θ)
            value_and_gradient!(ℓπ, grad, prep, ad, θ)   # warm up / compile
            # Best of a few, not a mean: we are after the achievable cost, and
            # a model build shares the machine with whatever else is running.
            minimum(1:5) do _
                t0 = time_ns()
                value_and_gradient!(ℓπ, grad, prep, ad, θ)
                time_ns() - t0
            end
        catch err
            err isa InterruptException && rethrow()
            @warn "Chunk width $c failed; skipping it." exception=err
            continue
        end
        verbosity >= 2 && @info "  chunk $c: $(round(t/1e3, digits=2)) µs ($(cld(D, c)) chunks)"
        if t < best_t
            best, best_t = c, t
        end
    end
    verbosity >= 1 && @info "Selected ForwardDiff chunk width" chunksize=best time_µs=round(best_t/1e3, digits=2)
    return best
end

"""
    LogDensityModel(system; autodiff=nothing, verbosity=2, chunk_sizes=nothing)

Compile `system` into a log-posterior and its gradient, ready for a sampler.

# Keywords
  - `autodiff` — an ADType (e.g. `AutoForwardDiff(chunksize=8)`) to use
    verbatim, or `false` for no gradient. Defaults to ForwardDiff at a chunk
    width chosen by [`_chunk_heuristic`](@ref).
  - `chunk_sizes` — candidate ForwardDiff chunk widths to **measure** rather
    than assume, e.g. `chunk_sizes=[8, 12, 16, 24, D]`. Each candidate costs a
    gradient compile, so this is worth it only for a long run on a large model,
    or to check the default on unfamiliar hardware. Ignored if `autodiff` is
    given.
  - `verbosity` — 0 silences the timing and selection reports.
"""
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

        # Deferred: with `chunk_sizes` given we time each candidate, which needs
        # `ℓπcallback` to exist first. See `_chunk_heuristic`.
        autodiff_requested = autodiff

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
        ℓπcallback, ∇ℓπcallback, autodiff = (function(
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
                # No gradient, but `ADType` is still a type parameter of the
                # model and `autodiff_type` reads it — keep what was asked for.
                return ℓπcallback, nothing,
                    isnothing(autodiff_requested) ?
                        AutoForwardDiff(chunksize=_chunk_heuristic(D)) : autodiff_requested
            end

            autodiff = if !isnothing(autodiff_requested)
                autodiff_requested
            elseif isnothing(chunk_sizes)
                AutoForwardDiff(chunksize=_chunk_heuristic(D))
            else
                AutoForwardDiff(chunksize=_probe_chunk_size(
                    ℓπcallback, initial_θ_0_t, D, chunk_sizes, verbosity))
            end

            ∇ℓπcallback =
            let ℓπcallback=ℓπcallback,
                autodiff=autodiff,
                prep = prepare_gradient(ℓπcallback, autodiff, zero(initial_θ_0_t)),
                grad = similar(initial_θ_0_t)
                function(initial_θ_0_t)
                    return value_and_gradient!(ℓπcallback, grad, prep, autodiff, initial_θ_0_t)
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


            ℓπcallback, ∇ℓπcallback, autodiff
        end)(
            arr2nt,
            system,
            Bijector_invlinkvec,
            ln_prior_transformed,
            ln_like_generated,
            D
        )
        # `autodiff` is chosen inside the closure (the measured probe needs
        # `ℓπcallback`), and it is a type parameter of the model, so it has to
        # come back out. Discrete models return `nothing` for it.
        
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
