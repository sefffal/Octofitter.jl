#=
Per-term AD backend dispatch and resolution.

Each observation type can declare a preferred AD backend via `ad_backend(obs)`.
The `resolve_ad_backend` function resolves the final backend to use, considering
model-level overrides.
=#

using ADTypes: AutoForwardDiff, AutoFiniteDiff, AutoEnzyme
using Bumper: @no_escape
using Enzyme: Duplicated

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
    _ensure_duplicated(backend::AutoEnzyme)

Override function_annotation to Duplicated for Enzyme backends.
Required when the differentiated closure captures pre-allocated mutable
buffers (orbit solution buffers, θ_embedded buffer) that Enzyme must
create shadows for.
"""
_ensure_duplicated(b::AutoEnzyme) = AutoEnzyme(mode=b.mode, function_annotation=Duplicated)
_ensure_duplicated(b) = b  # non-Enzyme backends unchanged

# --- Config and workspace types for per-term evaluation ---

"""
    ModelEvalConfig{TInvlink, TArr2nt, TOCs}

Immutable model-level configuration captured once and shared by all per-term
evaluation closures. Replaces per-closure captures of invlink, arr2nt, etc.
"""
struct ModelEvalConfig{TInvlink, TArr2nt, TOCs}
    invlink::TInvlink
    arr2nt::TArr2nt
    orbit_constructors::TOCs   # NTuple of per-planet closures
    n_planets::Int
end

"""
    TermEvalConfig{TObs, IsPlanet}

Per-observation-term configuration. Each term gets its own `TermEvalConfig`
capturing the observation object, its name, epochs, active parameter indices,
and which planet it belongs to (0 for system-level observations).

The `IsPlanet` type parameter (`Val{true}` or `Val{false}`) enables dispatch-based
specialization of `_get_obs_params` and `_make_obs_context`, eliminating union return
types that force Enzyme to generate shadow code for dead branches.
"""
struct TermEvalConfig{TObs, IsPlanet}
    obs::TObs
    obs_name::Union{Nothing, Symbol}
    epochs::Vector{Float64}
    active_indices::Vector{Int}
    planet_index::Int          # 0 for system obs
    _is_planet::IsPlanet       # Val{true} or Val{false}
end
function TermEvalConfig(obs, obs_name, epochs, active_indices, planet_index)
    is_planet = planet_index > 0 ? Val(true) : Val(false)
    TermEvalConfig(obs, obs_name, epochs, active_indices, planet_index, is_planet)
end

"""Zero-size sentinel for Bumper-based workspace (ForwardDiff/FiniteDiff)."""
struct BumperTermWS end

"""Pre-allocated workspace for non-Bumper backends (Enzyme/Mooncake/Reactant)."""
struct PreallocTermWS{TSolBufs}
    sol_bufs::TSolBufs         # NTuple{N, Vector{OrbitSolution{Float64}}}
    θ_embedded::Vector{Float64}
end

@inline function _get_obs_params(tcfg::TermEvalConfig{<:Any, Val{true}}, θ_system)
    θ_planet = θ_system.planets[tcfg.planet_index]
    θ_obs = if tcfg.obs_name !== nothing && hasproperty(θ_planet.observations, tcfg.obs_name)
        getproperty(θ_planet.observations, tcfg.obs_name)
    else; (;) end
    return θ_planet, θ_obs
end
@inline function _get_obs_params(tcfg::TermEvalConfig{<:Any, Val{false}}, θ_system)
    θ_obs = if tcfg.obs_name !== nothing && hasproperty(θ_system.observations, tcfg.obs_name)
        getproperty(θ_system.observations, tcfg.obs_name)
    else; (;) end
    return (;), θ_obs
end

@inline function _make_obs_context(tcfg::TermEvalConfig{<:Any, Val{true}}, θ_system, θ_planet, θ_obs, orbits, solutions)
    PlanetObservationContext(θ_system, θ_planet, θ_obs, orbits, solutions, tcfg.planet_index)
end
@inline function _make_obs_context(tcfg::TermEvalConfig{<:Any, Val{false}}, θ_system, θ_planet, θ_obs, orbits, solutions)
    SystemObservationContext(θ_system, θ_obs, orbits, solutions)
end

function eval_term(θ, mcfg::ModelEvalConfig, tcfg::TermEvalConfig, ::BumperTermWS)
    θ_natural = mcfg.invlink(θ)
    θ_system  = mcfg.arr2nt(θ_natural)
    orbits    = ntuple(i -> mcfg.orbit_constructors[i](θ_system), mcfg.n_planets)
    θ_planet, θ_obs = _get_obs_params(tcfg, θ_system)
    @no_escape begin
        solutions = _solve_all_orbits_bumper(orbits, tcfg.epochs)
        ctx = _make_obs_context(tcfg, θ_system, θ_planet, θ_obs, orbits, solutions)
        ln_like(tcfg.obs, ctx)
    end
end

function eval_term(θ, mcfg::ModelEvalConfig, tcfg::TermEvalConfig, ws::PreallocTermWS)
    θ_natural = mcfg.invlink(θ)
    θ_system  = mcfg.arr2nt(θ_natural)
    orbits    = ntuple(i -> mcfg.orbit_constructors[i](θ_system), mcfg.n_planets)
    θ_planet, θ_obs = _get_obs_params(tcfg, θ_system)
    solutions = _solve_all_orbits!(ws.sol_bufs, orbits, tcfg.epochs)
    ctx = _make_obs_context(tcfg, θ_system, θ_planet, θ_obs, orbits, solutions)
    ln_like(tcfg.obs, ctx)
end

function eval_term_sparse(θ_active, mcfg::ModelEvalConfig, tcfg::TermEvalConfig,
                          ::BumperTermWS, θ_full_base)
    @no_escape begin
        buf = Bumper.default_buffer()
        θ = Bumper.alloc!(buf, eltype(θ_active), length(θ_full_base))
        _embed_active!(θ, θ_full_base, θ_active, tcfg.active_indices)
        eval_term(θ, mcfg, tcfg, BumperTermWS())
    end
end

function eval_term_sparse(θ_active, mcfg::ModelEvalConfig, tcfg::TermEvalConfig,
                          ws::PreallocTermWS, θ_full_base)
    _embed_active!(ws.θ_embedded, θ_full_base, θ_active, tcfg.active_indices)
    eval_term(ws.θ_embedded, mcfg, tcfg, ws)
end

"""
    _make_term_closure(mcfg, tcfg, obs_backend, use_bumper, D_active, D,
                       sol_type_examples, epochs)

Build a per-term evaluation closure and (possibly mutated) backend.
Replaces the 12-path closure matrix with 4 paths: {Bumper,Prealloc} × {all-active,sparse}.
"""
function _make_term_closure(mcfg, tcfg, obs_backend, use_bumper, D_active, D,
                            sol_type_examples, epochs)
    all_active = D_active == D
    n_planets = mcfg.n_planets

    if use_bumper
        if all_active
            closure = let mcfg=mcfg, tcfg=tcfg
                (θ,) -> eval_term(θ, mcfg, tcfg, BumperTermWS())
            end
        else
            closure = let mcfg=mcfg, tcfg=tcfg
                (θ_active, θ_fb) -> eval_term_sparse(θ_active, mcfg, tcfg, BumperTermWS(), θ_fb)
            end
        end
    else
        sol_bufs = ntuple(n_planets) do ip
            Vector{typeof(sol_type_examples[ip])}(undef, length(epochs))
        end
        ws = PreallocTermWS(sol_bufs, Vector{Float64}(undef, D))
        obs_backend = _ensure_duplicated(obs_backend)
        if all_active
            closure = let mcfg=mcfg, tcfg=tcfg, ws=ws
                (θ,) -> eval_term(θ, mcfg, tcfg, ws)
            end
        else
            closure = let mcfg=mcfg, tcfg=tcfg, ws=ws
                (θ_active, θ_fb) -> eval_term_sparse(θ_active, mcfg, tcfg, ws, θ_fb)
            end
        end
    end
    return closure, obs_backend
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

"""
    _solve_all_orbits!(sol_bufs::NTuple{N}, orbits::NTuple{N}, epochs) where N

Solve all orbits at given epochs, writing into pre-allocated buffers.
Returns `sol_bufs` (as a Tuple) for use as `solutions` in observation contexts.

Pre-allocated buffers avoid heap allocation inside the differentiated closure,
which is critical for Enzyme performance (each heap alloc inside the AD boundary
requires shadow memory and bookkeeping).
"""
function _solve_all_orbits!(sol_bufs::NTuple{N}, orbits::NTuple{N}, epochs) where N
    if isempty(epochs)
        return sol_bufs
    end
    for i in 1:N
        orbit = orbits[i]
        buf = sol_bufs[i]
        for j in eachindex(epochs)
            buf[j] = orbitsolve(orbit, epochs[j])
        end
    end
    return sol_bufs
end

"""
    _embed_active!(θ_embedded, θ_full_base, θ_active, active_indices)

Embed active parameters into a full-length vector. Inactive positions
get Float64 base values (with zero AD partials), active positions get
the AD-tracked values from θ_active.
"""
@inline function _embed_active!(θ_embedded, θ_full_base, θ_active, active_indices)
    @inbounds for k in eachindex(θ_embedded)
        θ_embedded[k] = θ_full_base[k]
    end
    @inbounds for (j, idx) in enumerate(active_indices)
        θ_embedded[idx] = θ_active[j]
    end
    return θ_embedded
end

"""
    _compute_active_indices(system) -> (Tuple of Vector{Int}, D)

Compute per-term active parameter indices. Returns a tuple of Vector{Int}
(one per observation term, in closure construction order: planet obs first,
then system obs) and the total parameter count D.

Each term's active indices include:
- All system prior indices (always needed for orbit construction)
- All planet prior indices (always needed — epicycle approximation)
- This observation's own prior indices

Parameters from other observations are excluded (structurally zero gradient).
"""
function _compute_active_indices(system)
    idx = 0

    # System priors
    sys_prior_start = idx + 1
    for key in keys(system.priors.priors)
        idx += length(system.priors.priors[key])
    end
    sys_prior_end = idx

    # System observation priors — track each obs separately
    sys_obs_ranges = UnitRange{Int}[]
    for obs in system.observations
        obs_start = idx + 1
        if hasproperty(obs, :priors)
            for key in keys(obs.priors.priors)
                idx += length(obs.priors.priors[key])
            end
        end
        push!(sys_obs_ranges, obs_start:idx)
    end

    # Planet priors and planet observation priors
    planet_prior_ranges = UnitRange{Int}[]
    planet_obs_ranges = Vector{UnitRange{Int}}[]
    for planet in system.planets
        planet_start = idx + 1
        for key in keys(planet.priors.priors)
            idx += length(planet.priors.priors[key])
        end
        push!(planet_prior_ranges, planet_start:idx)

        obs_ranges = UnitRange{Int}[]
        for obs in planet.observations
            obs_start = idx + 1
            if hasproperty(obs, :priors)
                for key in keys(obs.priors.priors)
                    idx += length(obs.priors.priors[key])
                end
            end
            push!(obs_ranges, obs_start:idx)
        end
        push!(planet_obs_ranges, obs_ranges)
    end

    D = idx

    # "Always active" = system priors + all planet priors
    always_active = Int[]
    for i in sys_prior_start:sys_prior_end
        push!(always_active, i)
    end
    for r in planet_prior_ranges
        for i in r
            push!(always_active, i)
        end
    end
    sort!(always_active)

    # Build active indices per term, in closure order:
    # Planet obs first, then system obs
    active_indices_list = Vector{Int}[]

    for (ip, planet) in enumerate(system.planets)
        for (io, obs) in enumerate(planet.observations)
            obs_range = planet_obs_ranges[ip][io]
            if isempty(obs_range)
                push!(active_indices_list, copy(always_active))
            else
                active = sort!(union(always_active, collect(obs_range)))
                push!(active_indices_list, active)
            end
        end
    end

    for (io, obs) in enumerate(system.observations)
        obs_range = sys_obs_ranges[io]
        if isempty(obs_range)
            push!(active_indices_list, copy(always_active))
        else
            active = sort!(union(always_active, collect(obs_range)))
            push!(active_indices_list, active)
        end
    end

    return Tuple(active_indices_list), D
end

#=
Type-stable per-term gradient infrastructure.

Each term's closure, AD backend, and factories for prep/gradient objects are
bundled into a TermGradSpec and stored in a heterogeneous Tuple. Actual mutable
buffers are obtained via task-local storage for parallel MCMC chain safety.
=#

"""
    _get_task_local(::Type{T}, key, factory) where T

Get or create a task-local value. Uses `Base.task_local_storage` so that each
MCMC chain (running on its own Task) gets independent mutable buffers.
The type parameter `T` is required for type stability, since `task_local_storage()`
returns `Dict{Symbol, Any}`.
"""
@inline function _get_task_local(::Type{T}, key, factory) where T
    get!(factory, task_local_storage(), key)::T
end

"""
    TermGradSpec{TC,TB,TPF,TGF,TAF,TP,TG,TA}

Bundles a single term's gradient infrastructure for type-stable storage
in a heterogeneous Tuple. Stores factories rather than mutable buffers
so that task-local copies can be created on demand.

Non-all-active terms always receive `θ_full_base` as a second argument
passed via `DifferentiationInterface.Constant` context.
"""
struct TermGradSpec{TC,TB,TPF,TGF,TAF,TP,TG,TA}
    closure::TC
    backend::TB
    prep_factory::TPF
    grad_factory::TGF
    θ_active_factory::TAF
    tls_prep_key::Symbol
    tls_grad_key::Symbol
    tls_θ_active_key::Symbol
    active_indices::Vector{Int}
    all_active::Bool  # true when active_indices == 1:D (skip copy)
    prep_type::Type{TP}
    grad_type::Type{TG}
    θ_active_type::Type{TA}
end

"""
    _prepare_term_specs(closures, backends, active_indices_tuple, θ_example, θ_full_example)

Recursively build a type-stable heterogeneous tuple of TermGradSpec, storing
factories for `prepare_gradient` and gradient buffers instead of mutable state.
Each term uses D_active-length buffers based on its active parameter indices.

Non-all-active terms automatically get `prepare_gradient` called with a
`DifferentiationInterface.Constant` context to match `_accumulate_term_gradients!`.
"""
_prepare_term_specs(::Tuple{}, ::Tuple{}, ::Tuple{}, θ_example, θ_full_example) = ()
function _prepare_term_specs(closures::Tuple, backends::Tuple, active_indices_tuple::Tuple,
                             θ_example, θ_full_example)
    c = first(closures)
    b = first(backends)
    active = first(active_indices_tuple)
    θ_active_example = θ_example[active]
    θ_zero = zero(θ_active_example)
    D_total = length(θ_example)
    D_active = length(active)
    all_active = D_active == D_total

    # Create prep factory — includes Constant context for sparse terms
    prep_factory = if !all_active
        let c=c, b=b, θ_zero=θ_zero, ctx=DifferentiationInterface.Constant(zero(θ_full_example))
            () -> prepare_gradient(c, b, θ_zero, ctx)
        end
    else
        let c=c, b=b, θ_zero=θ_zero
            () -> prepare_gradient(c, b, θ_zero)
        end
    end

    grad_factory = let θ_ex=θ_active_example
        () -> similar(θ_ex)
    end
    θ_active_factory = let θ_ex=θ_active_example
        () -> similar(θ_ex)
    end
    # Compute concrete types for type-stable task-local retrieval
    prep_example = prep_factory()
    grad_example = grad_factory()
    θ_active_example_buf = θ_active_factory()
    key_id = gensym(:term)
    spec = TermGradSpec(c, b, prep_factory, grad_factory, θ_active_factory,
                        Symbol(key_id, :_prep), Symbol(key_id, :_grad),
                        Symbol(key_id, :_θ_active),
                        active, all_active,
                        typeof(prep_example), typeof(grad_example),
                        typeof(θ_active_example_buf))
    rest = _prepare_term_specs(Base.tail(closures), Base.tail(backends),
                               Base.tail(active_indices_tuple), θ_example,
                               θ_full_example)
    return (spec, rest...)
end

"""
    _accumulate_term_gradients!(total_grad, ll, specs::Tuple, θ, θ_full_base) -> ll_total

Recursively evaluate each term's gradient and accumulate into `total_grad`.
Type-stable because each recursive call peels off a concrete element from
the heterogeneous tuple. Gradient buffers and prep objects are task-local.

Each term operates on a reduced-dimension active subset of θ, then scatters
its gradient contributions back to the full-D total_grad vector.

`θ_full_base` is a `Vector{Float64}` with the current (untransformed) parameter
values, passed as `Constant` context to closures that use the sparse embedding
pattern (avoids `task_local_storage` inside the AD boundary).
"""
_accumulate_term_gradients!(total_grad, ll, ::Tuple{}, θ, θ_full_base) = ll
function _accumulate_term_gradients!(total_grad, ll, specs::Tuple, θ, θ_full_base)
    spec = first(specs)
    prep = _get_task_local(spec.prep_type, spec.tls_prep_key, spec.prep_factory)
    active_grad = _get_task_local(spec.grad_type, spec.tls_grad_key, spec.grad_factory)

    if spec.all_active
        # Fast path: all parameters active, no subsetting needed
        ll_i, _ = value_and_gradient!(spec.closure, active_grad, prep, spec.backend, θ)
        ll += ll_i
        @inbounds for k in eachindex(total_grad)
            total_grad[k] += active_grad[k]
        end
    else
        # Extract active subset into preallocated buffer
        θ_active = _get_task_local(spec.θ_active_type, spec.tls_θ_active_key, spec.θ_active_factory)
        @inbounds for (j, idx) in enumerate(spec.active_indices)
            θ_active[j] = θ[idx]
        end

        # Gradient on reduced-dimension input — always pass θ_full_base via Constant
        ll_i, _ = value_and_gradient!(spec.closure, active_grad, prep, spec.backend,
                                      θ_active, DifferentiationInterface.Constant(θ_full_base))
        ll += ll_i

        # Scatter back to full gradient
        @inbounds for (j, idx) in enumerate(spec.active_indices)
            total_grad[idx] += active_grad[j]
        end
    end

    return _accumulate_term_gradients!(total_grad, ll, Base.tail(specs), θ, θ_full_base)
end
