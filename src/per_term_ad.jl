#=
Per-term AD backend dispatch and resolution.

Each observation type can declare a preferred AD backend via `ad_backend(obs)`.
The `resolve_ad_backend` function resolves the final backend to use, considering
model-level overrides.
=#

using ADTypes: AutoForwardDiff, AutoFiniteDiff

"""
Backends compatible with Bumper.jl's `@no_escape`/`@alloc` pattern.
Non-Bumper backends (Enzyme, Reactant, Mooncake, etc.) use plain heap allocation.
"""
_uses_bumper(::AutoForwardDiff) = true
_uses_bumper(::AutoFiniteDiff) = true
_uses_bumper(::Nothing) = true
_uses_bumper(_) = false

"""
    ad_backend(obs::AbstractObs)

Return the preferred AD backend for this observation type, or `nothing` to use
the default (ForwardDiff). External packages can specialize this to request
e.g. `AutoEnzyme()` for specific observation types.

This function is public but not exported. Access as `Octofitter.ad_backend(obs)`.
"""
ad_backend(::AbstractObs) = nothing

"""
    resolve_ad_backend(obs, model_override, D::Int)

Resolve which AD backend to use for a given observation term.

- If `model_override` is an AD backend type (from DifferentiationInterface) -> use it
- If `model_override === false` -> return `nothing` (no gradients)
- If `model_override === nothing` -> check `ad_backend(obs)`, fall back to `AutoForwardDiff(chunksize=D)`
"""
function resolve_ad_backend(obs, model_override, D::Int)
    if model_override === false
        return nothing
    elseif !isnothing(model_override)
        # Model-level override takes precedence
        return model_override
    else
        # Check per-type preference
        backend = ad_backend(obs)
        if !isnothing(backend)
            return backend
        end
        # Default to ForwardDiff
        return AutoForwardDiff(chunksize=D)
    end
end

"""
    _get_obs_epochs(obs)

Get the epochs from an observation object, or an empty Float64 vector if none.
"""
function _get_obs_epochs(obs)
    if hasproperty(obs, :table) && hasproperty(obs.table, :epoch)
        return collect(Float64, obs.table.epoch)
    end
    return Float64[]
end

"""
    _solve_all_orbits(orbits::NTuple{N}, epochs) where N

Solve all orbits at given epochs. Returns a tuple of solution vectors.
Uses heap allocation — for use in `generate_from_params` and other
non-hot-path code.
"""
function _solve_all_orbits(orbits::NTuple{N}, epochs) where N
    if isempty(epochs)
        return ntuple(Returns(()), Val(N))
    end
    ntuple(Val(N)) do i
        map(ep -> orbitsolve(orbits[i], ep), epochs)
    end
end

"""
    _solve_all_orbits_bumper(orbits::NTuple{N}, epochs) where N

Solve all orbits at given epochs using Bumper allocation.
**Must be called from within a `@no_escape` block.**
Uses `Bumper.alloc!` (the function) rather than `@alloc` (the macro)
so it can be called from a separate function scope.
"""
function _solve_all_orbits_bumper(orbits::NTuple{N}, epochs) where N
    if isempty(epochs)
        return ntuple(Returns(()), Val(N))
    end
    buf = Bumper.default_buffer()
    ntuple(Val(N)) do i
        orbit = orbits[i]
        sol0 = orbitsolve(orbit, first(epochs))
        sols = Bumper.alloc!(buf, typeof(sol0), length(epochs))
        sols[1] = sol0
        for j in 2:length(epochs)
            sols[j] = orbitsolve(orbit, epochs[j])
        end
        sols
    end
end

#=
Type-stable per-term gradient infrastructure.

Each term's closure, AD backend, and factories for prep/gradient objects are
bundled into a TermGradSpec and stored in a heterogeneous Tuple. Actual mutable
buffers are obtained via task-local storage for parallel MCMC chain safety.
=#

"""
    _get_task_local(key, factory)

Get or create a task-local value. Uses `Base.task_local_storage` so that each
MCMC chain (running on its own Task) gets independent mutable buffers.
"""
function _get_task_local(key, factory)
    get!(factory, task_local_storage(), key)
end

"""
    TermGradSpec{TC,TB,TPF,TGF}

Bundles a single term's gradient infrastructure for type-stable storage
in a heterogeneous Tuple. Stores factories rather than mutable buffers
so that task-local copies can be created on demand.
"""
struct TermGradSpec{TC,TB,TPF,TGF}
    closure::TC
    backend::TB
    prep_factory::TPF
    grad_factory::TGF
    tls_prep_key::Symbol
    tls_grad_key::Symbol
end

"""
    _prepare_term_specs(closures::Tuple, backends::Tuple, θ_example) -> Tuple{TermGradSpec...}

Recursively build a type-stable heterogeneous tuple of TermGradSpec, storing
factories for `prepare_gradient` and gradient buffers instead of mutable state.
"""
_prepare_term_specs(::Tuple{}, ::Tuple{}, θ_example) = ()
function _prepare_term_specs(closures::Tuple, backends::Tuple, θ_example)
    c = first(closures)
    b = first(backends)
    θ_zero = zero(θ_example)
    prep_factory = let c=c, b=b, θ_zero=θ_zero
        () -> prepare_gradient(c, b, θ_zero)
    end
    grad_factory = let θ_example=θ_example
        () -> similar(θ_example)
    end
    key_id = gensym(:term)
    spec = TermGradSpec(c, b, prep_factory, grad_factory,
                        Symbol(key_id, :_prep), Symbol(key_id, :_grad))
    rest = _prepare_term_specs(Base.tail(closures), Base.tail(backends), θ_example)
    return (spec, rest...)
end

"""
    _accumulate_term_gradients!(total_grad, ll, specs::Tuple, θ) -> ll_total

Recursively evaluate each term's gradient and accumulate into `total_grad`.
Type-stable because each recursive call peels off a concrete element from
the heterogeneous tuple. Gradient buffers and prep objects are task-local.
"""
_accumulate_term_gradients!(total_grad, ll, ::Tuple{}, θ) = ll
function _accumulate_term_gradients!(total_grad, ll, specs::Tuple, θ)
    spec = first(specs)
    prep = _get_task_local(spec.tls_prep_key, spec.prep_factory)
    grad = _get_task_local(spec.tls_grad_key, spec.grad_factory)
    ll_i, _ = value_and_gradient!(spec.closure, grad, prep, spec.backend, θ)
    ll += ll_i
    total_grad .+= grad
    return _accumulate_term_gradients!(total_grad, ll, Base.tail(specs), θ)
end
