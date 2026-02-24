#=
Per-term AD backend dispatch and resolution.

Each observation type can declare a preferred AD backend via `ad_backend(obs)`.
The `resolve_ad_backend` function resolves the final backend to use, considering
model-level overrides.

Workspace model:
- Model = immutable configuration (no closures, no mutable state)
- Workspace = mutable state (pre-allocated buffers, gradient infrastructure)
- Callbacks = plain functions taking (model, θ, workspace) instead of stored closures

All backends are Enzyme-based (reverse-mode, source-to-source AD). Enzyme keeps
the element type as Float64 and tracks derivatives via shadow memory, so all
pre-allocated workspace buffers can be plain Float64 arrays.
=#

using ADTypes: AutoForwardDiff, AutoFiniteDiff, AutoEnzyme
using Enzyme: Duplicated, Reverse, Const, set_runtime_activity

"""Default Enzyme backend for all observation terms. Uses reverse-mode with
runtime activity for maximum compatibility."""
const _default_enzyme_backend = AutoEnzyme(mode=set_runtime_activity(Reverse), function_annotation=Const)

"""
    alloc_obs_workspace(obs::AbstractObs, ::Type{T}) where T

Allocate observation-specific workspace buffers typed with element type `T`.
Returns `nothing` by default. Observation types that need pre-allocated buffers
should specialize this method to return a named tuple of typed arrays.
"""
alloc_obs_workspace(::AbstractObs, ::Type{T}) where T = nothing

"""
    ad_backend(obs::AbstractObs)

Return the preferred AD backend for this observation type, or `nothing` to use
the default (Enzyme reverse-mode). External packages can specialize this to request
a different backend for specific observation types.

This function is public but not exported. Access as `Octofitter.ad_backend(obs)`.
"""
ad_backend(::AbstractObs) = nothing

"""
    ReactantBackend{TD}

Backend sentinel for Reactant/XLA compilation. Returned by `ad_backend(obs)` when
an observation type opts in to Reactant acceleration (e.g., via a `device` field).

The `device` field carries the target device (e.g., `Reactant.CPU()`).
The type is defined in core Octofitter so that both the main Reactant extension
and per-package extensions can dispatch on it. The actual compilation logic lives
in `OctofitterReactantExt`.
"""
struct ReactantBackend{TD}
    device::TD
end

"""
    resolve_ad_backend(obs, model_override, D::Int)

Resolve which AD backend to use for a given observation term.

- If `model_override` is an AD backend type (from DifferentiationInterface) -> use it
- If `model_override === false` -> return `nothing` (no gradients)
- If `model_override === nothing` -> check `ad_backend(obs)`, fall back to Enzyme reverse-mode
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
        # Default to Enzyme reverse-mode
        return _default_enzyme_backend
    end
end

# --- Config and workspace types for per-term evaluation ---

"""
    _make_construct_orbits(system::System)

Build a `RuntimeGeneratedFunction` that constructs all planet orbits from
structured parameters `θ_system`. Orbit types and planet indices are
interpolated as compile-time constants, guaranteeing type stability
regardless of the caller's closure complexity.

Returns `θ_system -> (orbit_1, orbit_2, ...)`.
"""
function _make_construct_orbits(system::System)
    n_planets = length(system.planets)
    planet_exprs = map(1:n_planets) do i
        OT = _planet_orbit_type(system.planets[i])
        :($(OT)(;merge(θ_system, θ_system.planets[$i])...))
    end
    @RuntimeGeneratedFunction(:(function(θ_system)
        tuple($(planet_exprs...))
    end))
end

"""
    ModelEvalConfig{TInvlink, TArr2nt, N, TConstructOrbits}

Immutable model-level configuration captured once and shared by all per-term
evaluation closures. Replaces per-closure captures of invlink, arr2nt, etc.

The number of planets `N` is encoded as a type parameter for `_n_planets`.
`construct_orbits` is a `RuntimeGeneratedFunction` that takes `θ_system`
and returns an `NTuple{N}` of orbits with all types and indices as
compile-time constants, ensuring type stability.
"""
struct ModelEvalConfig{TInvlink, TArr2nt, N, TConstructOrbits}
    invlink::TInvlink
    arr2nt::TArr2nt
    construct_orbits::TConstructOrbits
    function ModelEvalConfig(invlink::TI, arr2nt::TA, construct_orbits::TCO, ::Val{N}) where {TI, TA, TCO, N}
        new{TI, TA, N, TCO}(invlink, arr2nt, construct_orbits)
    end
end
@inline _n_planets(::ModelEvalConfig{<:Any, <:Any, N}) where N = Val(N)

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

"""
    TermWorkspace{T, TSolBufs, TθEmbed, TObsWS}

Pre-allocated workspace for per-term evaluation. Replaces BumperTermWS and
PreallocTermWS with a single unified type.

# Fields
- `sol_bufs` — NTuple{N, Vector{OrbitSolution{...}}} for orbit solutions
- `θ_embedded` — Vector{T} for sparse parameter embedding, or nothing
- `obs_workspace` — observation-specific buffers from `alloc_obs_workspace`
"""
struct TermWorkspace{T, TSolBufs, TθEmbed, TObsWS}
    sol_bufs::TSolBufs
    θ_embedded::TθEmbed
    obs_workspace::TObsWS
end

"""
    TermEvaluator{TMcfg}

Callable struct replacing closures for per-term evaluation.
Contains only `mcfg` (with zero-size RuntimeGeneratedFunction fields) so that
Enzyme can prove the function argument is readonly when annotated as `Const`.

`tcfg` (which contains heap-allocated `Vector` fields and observation data) is
passed separately as `DI.Constant` — inactive for differentiation but not
subject to Enzyme's readonly proof requirement for the function argument.
`workspace` is passed as `DI.Cache` → `Duplicated` for Enzyme.

Call signatures:
- `(θ, tcfg, workspace)`: Enzyme path — pre-allocated workspace (zero alloc)
- `(θ, tcfg)`: ForwardDiff path — allocating (dual numbers change element type)
"""
struct TermEvaluator{TMcfg}
    mcfg::TMcfg
end
(te::TermEvaluator)(θ, tcfg, workspace) = eval_term(θ, te.mcfg, tcfg, workspace)
(te::TermEvaluator)(θ, tcfg) = eval_term_alloc(θ, te.mcfg, tcfg)

"""
    SparseTermEvaluator{TMcfg}

Callable struct for sparse (active parameter subset) evaluation.
Like `TermEvaluator`, contains only `mcfg` for Enzyme `Const` compatibility.
`tcfg`, `θ_full_base`, and `workspace` are passed as DI contexts.

Call signatures:
- `(θ_active, θ_full_base, tcfg, workspace)`: Enzyme path
- `(θ_active, θ_full_base, tcfg)`: ForwardDiff path — allocating
"""
struct SparseTermEvaluator{TMcfg}
    mcfg::TMcfg
end
(te::SparseTermEvaluator)(θ_active, θ_full_base, tcfg, workspace) =
    eval_term_sparse(θ_active, te.mcfg, tcfg, workspace, θ_full_base)
(te::SparseTermEvaluator)(θ_active, θ_full_base, tcfg) =
    eval_term_sparse_alloc(θ_active, te.mcfg, tcfg, θ_full_base)

"""
    ReactantTermEvaluator{TMcfg, TTcfg}

Callable struct for Reactant-compiled term evaluation.
No mutable workspace — Reactant allocates arrays via broadcasting inside the
compiled XLA graph. This struct is what gets passed to `@compile`.
"""
struct ReactantTermEvaluator{TMcfg, TTcfg}
    mcfg::TMcfg
    tcfg::TTcfg
end
# Callable defined in OctofitterReactantExt (uses @allowscalar from Reactant)

"""
    PriorEvaluator{TInvlink, TLnPrior}

Callable struct replacing the prior closure. Combines inverse link transform
and prior evaluation into a single callable.
"""
struct PriorEvaluator{TInvlink, TLnPrior}
    invlink::TInvlink
    ln_prior::TLnPrior
end
function (pe::PriorEvaluator)(θ_transformed)
    if any(!isfinite, θ_transformed)
        return convert(eltype(θ_transformed), -Inf)
    end
    θ_natural = pe.invlink(θ_transformed)
    pe.ln_prior(θ_natural, true)
end

"""
    TermGradSpec{TEval, TB, TP, TG, TA, TWs}

Bundles a single term's gradient infrastructure for type-stable storage
in a heterogeneous Tuple. Stores direct buffer references (no factories,
no task-local storage).

Backend dispatch determines the gradient computation path:
- **Enzyme** (`AutoEnzyme`): DI's `value_and_gradient!` with workspace passed as
  `DI.Cache`. DI wraps this as `Duplicated(workspace, shadow)` for Enzyme, ensuring
  mutable workspace buffers are properly tracked (not hidden inside a `Const` struct).
- **ForwardDiff** (`AutoForwardDiff`): DI's `value_and_gradient!` with allocating eval.
  ForwardDiff dual numbers change element type, so Float64 workspace isn't used.
  `workspace` is `nothing`.
"""
struct TermGradSpec{TEval, TB, TP, TG, TA, TWs, TTcfg, TDiTcfg, TDiCache, TDiθFull}
    evaluator::TEval
    backend::TB
    prep::TP
    grad_buf::TG
    θ_active_buf::TA
    active_indices::Vector{Int}
    all_active::Bool
    workspace::TWs
    tcfg::TTcfg
    di_tcfg::TDiTcfg        # pre-created DI.Constant(tcfg) — reused every gradient call
    di_cache::TDiCache       # pre-created DI.Cache(workspace) for Enzyme, nothing otherwise
    di_θ_full_base::TDiθFull # pre-created DI.Constant(θ_full_base) for sparse terms, nothing for all_active
end

# --- Observation parameter extraction and context construction ---

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

@inline function _make_obs_context(tcfg::TermEvalConfig{<:Any, Val{true}}, θ_system, θ_planet, θ_obs, orbits, solutions, obs_workspace)
    PlanetObservationContext(θ_system, θ_planet, θ_obs, orbits, solutions, tcfg.planet_index, obs_workspace)
end
@inline function _make_obs_context(tcfg::TermEvalConfig{<:Any, Val{false}}, θ_system, θ_planet, θ_obs, orbits, solutions, obs_workspace)
    SystemObservationContext(θ_system, θ_obs, orbits, solutions, obs_workspace)
end

# Backward-compatible method without workspace
@inline function _make_obs_context(tcfg::TermEvalConfig{<:Any, Val{true}}, θ_system, θ_planet, θ_obs, orbits, solutions)
    PlanetObservationContext(θ_system, θ_planet, θ_obs, orbits, solutions, tcfg.planet_index)
end
@inline function _make_obs_context(tcfg::TermEvalConfig{<:Any, Val{false}}, θ_system, θ_planet, θ_obs, orbits, solutions)
    SystemObservationContext(θ_system, θ_obs, orbits, solutions)
end

# --- Unified eval_term (no Bumper, no branching) ---

"""
    eval_term(θ, mcfg, tcfg, tw::TermWorkspace)

Evaluate a single observation term with pre-allocated workspace.
Single unified method — no branching on backend type.

Enzyme keeps the element type as Float64, so all pre-allocated buffers
(sol_bufs, obs_workspace) are directly usable without type conversion.
"""
function eval_term(θ, mcfg::ModelEvalConfig, tcfg::TermEvalConfig, tw::TermWorkspace)
    θ_natural = mcfg.invlink(θ)
    θ_system  = mcfg.arr2nt(θ_natural)
    orbits    = mcfg.construct_orbits(θ_system)
    θ_planet, θ_obs = _get_obs_params(tcfg, θ_system)
    _solve_all_orbits!(tw.sol_bufs, orbits, tcfg.epochs)
    ctx = _make_obs_context(tcfg, θ_system, θ_planet, θ_obs, orbits, tw.sol_bufs, tw.obs_workspace)
    ln_like(tcfg.obs, ctx)
end

"""
    eval_term_sparse(θ_active, mcfg, tcfg, tw::TermWorkspace, θ_full_base)

Evaluate with sparse active parameter subsetting. Embeds active parameters
into the full-length vector, then calls eval_term.
"""
function eval_term_sparse(θ_active, mcfg::ModelEvalConfig, tcfg::TermEvalConfig,
                          tw::TermWorkspace, θ_full_base)
    _embed_active!(tw.θ_embedded, θ_full_base, θ_active, tcfg.active_indices)
    eval_term(tw.θ_embedded, mcfg, tcfg, tw)
end

"""
    eval_term_alloc(θ, mcfg, tcfg)

Allocating version of eval_term for use with ForwardDiff.
ForwardDiff changes the element type to Dual numbers, so pre-allocated
Float64 workspace buffers cannot be used. Allocates fresh buffers instead.
"""
function eval_term_alloc(θ, mcfg::ModelEvalConfig, tcfg::TermEvalConfig)
    θ_natural = mcfg.invlink(θ)
    θ_system  = mcfg.arr2nt(θ_natural)
    orbits    = mcfg.construct_orbits(θ_system)
    θ_planet, θ_obs = _get_obs_params(tcfg, θ_system)
    solutions = _solve_all_orbits(orbits, tcfg.epochs)
    ctx = _make_obs_context(tcfg, θ_system, θ_planet, θ_obs, orbits, solutions)
    ln_like(tcfg.obs, ctx)
end

"""
    eval_term_sparse_alloc(θ_active, mcfg, tcfg, θ_full_base)

Allocating sparse eval for ForwardDiff. Creates a fresh embedded vector
since ForwardDiff dual numbers change the element type.
"""
function eval_term_sparse_alloc(θ_active, mcfg::ModelEvalConfig, tcfg::TermEvalConfig, θ_full_base)
    θ_embedded = copy(θ_full_base)
    @inbounds for (j, idx) in enumerate(tcfg.active_indices)
        θ_embedded[idx] = θ_active[j]
    end
    eval_term_alloc(θ_embedded, mcfg, tcfg)
end

# --- Helper functions ---

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
    _solve_all_orbits_bulk(orbits::NTuple{N}, epochs) where N

Solve all orbits at given epochs using `orbitsolve_bulk`, returning
`NTuple{N, OrbitSolutionBulk}`. Fully vectorized for Reactant/XLA traceability.
"""
function _solve_all_orbits_bulk(orbits::NTuple{N}, epochs) where N
    if isempty(epochs)
        return ntuple(Returns(()), Val(N))
    end
    ntuple(i -> orbitsolve_bulk(orbits[i], epochs), Val(N))
end

"""
    _embed_active_broadcast(θ_base, θ_active, active_indices)

Embed active parameters into a copy of the base vector using broadcasting.
Unlike `_embed_active!`, this creates a new array (no mutation), making it
suitable for use inside Reactant's traced computation graph.
"""
function _embed_active_broadcast(θ_base, θ_active, active_indices)
    θ_full = copy(θ_base)
    for (j, idx) in enumerate(active_indices)
        θ_full[idx] = θ_active[j]
    end
    θ_full
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

# --- Workspace construction helpers ---

"""
    _make_typed_sol_bufs(construct_orbits, ::Val{N}, θ_system, epochs, ::Type{T}) where {N, T}

Create pre-allocated orbit solution buffers with the correct element type.
Runs `orbitsolve` once to determine the concrete `OrbitSolution{T}` type,
then pre-allocates vectors of that type.
"""
function _make_typed_sol_bufs(construct_orbits, ::Val{N}, θ_system, epochs, ::Type{T}) where {N, T}
    if isempty(epochs) || N == 0
        return ntuple(i -> Vector{Nothing}(undef, 0), Val(N))
    end
    orbits_example = construct_orbits(θ_system)
    ntuple(Val(N)) do ip
        sol0 = orbitsolve(orbits_example[ip], first(epochs))
        Vector{typeof(sol0)}(undef, length(epochs))
    end
end

"""
    _make_term_workspace(mcfg, tcfg, θ_system_example, ::Type{T}, D) where T

Create a TermWorkspace{T} with pre-allocated buffers for a single term.
"""
function _make_term_workspace(mcfg, tcfg, θ_system_example, ::Type{T}, D) where T
    sol_bufs = _make_typed_sol_bufs(mcfg.construct_orbits, _n_planets(mcfg), θ_system_example, tcfg.epochs, T)
    θ_embedded = Vector{T}(undef, D)
    obs_ws = alloc_obs_workspace(tcfg.obs, T)
    TermWorkspace{T, typeof(sol_bufs), typeof(θ_embedded), typeof(obs_ws)}(sol_bufs, θ_embedded, obs_ws)
end

# --- Primal and gradient workspace types ---

"""
    PrimalWorkspace{TTermWS}

Pre-allocated workspace for primal (non-gradient) evaluation.
Contains a heterogeneous Tuple of TermWorkspace{Float64}, one per observation term.
"""
struct PrimalWorkspace{TTermWS}
    term_workspaces::TTermWS
end

"""
    GradWorkspace{TTermSpecs, TPriorEval, TPriorPrep}

Pre-allocated workspace for gradient evaluation.
Contains per-term gradient specs (with typed buffers), prior evaluator,
and shared gradient accumulation buffers.
"""
struct GradWorkspace{TTermSpecs, TPriorEval, TPriorPrep}
    term_specs::TTermSpecs
    prior_evaluator::TPriorEval
    prior_prep::TPriorPrep
    prior_grad::Vector{Float64}
    total_grad::Vector{Float64}
    θ_full_base::Vector{Float64}
end

# --- Prepare functions ---

"""
    _build_term_workspaces(mcfg, tcfgs, θ_system_example, ::Type{T}, D)

Recursively build a type-stable heterogeneous tuple of TermWorkspace{T}.
"""
_build_term_workspaces(mcfg, ::Tuple{}, θ_system_example, ::Type{T}, D) where T = ()
function _build_term_workspaces(mcfg, tcfgs::Tuple, θ_system_example, ::Type{T}, D) where T
    tw = _make_term_workspace(mcfg, first(tcfgs), θ_system_example, T, D)
    (tw, _build_term_workspaces(mcfg, Base.tail(tcfgs), θ_system_example, T, D)...)
end

"""
    _sum_term_likelihoods(θ_sys, orbits, tcfgs, tws)

Recursive type-stable tuple peeling for summing term likelihoods.
Replaces `ln_like_generated` RuntimeGeneratedFunction.
"""
_sum_term_likelihoods(θ_sys, orbits, ::Tuple{}, ::Tuple{}) = 0.0
function _sum_term_likelihoods(θ_sys, orbits, tcfgs::Tuple, tws::Tuple)
    tcfg = first(tcfgs)
    tw = first(tws)
    θ_planet, θ_obs = _get_obs_params(tcfg, θ_sys)
    _solve_all_orbits!(tw.sol_bufs, orbits, tcfg.epochs)
    ctx = _make_obs_context(tcfg, θ_sys, θ_planet, θ_obs, orbits, tw.sol_bufs, tw.obs_workspace)
    ll = ln_like(tcfg.obs, ctx)
    ll + _sum_term_likelihoods(θ_sys, orbits, Base.tail(tcfgs), Base.tail(tws))
end

#=
Type-stable per-term gradient infrastructure.

Each term's evaluator, AD backend, prep, and gradient buffers are
bundled into a TermGradSpec and stored in a heterogeneous Tuple.
Buffers are stored directly (no factories, no task-local storage).
=#

"""
    _prepare_term_specs(evaluators, backends, active_indices_tuple, workspaces, tcfgs, θ_example, θ_full_base)

Recursively build a type-stable heterogeneous tuple of term gradient specs, storing
direct buffer references and gradient buffers.
Each term uses D_active-length buffers based on its active parameter indices.

`θ_full_base` is the runtime mutable buffer that the ∇ℓπcallback closure writes
`θ_transformed` into before each gradient call. DI context wrappers (`DI.Constant`,
`DI.Cache`) are pre-created here and stored in each `TermGradSpec`, eliminating
per-call wrapper allocation in the hot gradient loop.

Also triggers compilation eagerly at model construction time.

Per-term work is delegated to `_make_term_grad_spec`, which dispatches on
backend type. Extensions (e.g., OctofitterReactantExt) can add methods for
custom backends to provide alternative compilation strategies.
"""
_prepare_term_specs(::Tuple{}, ::Tuple{}, ::Tuple{}, ::Tuple{}, ::Tuple{}, θ_example, θ_full_base) = ()
function _prepare_term_specs(evaluators::Tuple, backends::Tuple, active_indices_tuple::Tuple,
                             workspaces::Tuple, tcfgs::Tuple, θ_example, θ_full_base)
    eval = first(evaluators)
    b = first(backends)
    active = first(active_indices_tuple)
    workspace = first(workspaces)
    tcfg = first(tcfgs)
    θ_active_example = θ_example[active]
    θ_zero = zero(θ_active_example)
    D_total = length(θ_example)
    D_active = length(active)
    all_active = D_active == D_total

    grad_buf = similar(θ_active_example)
    θ_active_buf = similar(θ_active_example)
    spec = _make_term_grad_spec(eval, b, workspace, tcfg, θ_zero, θ_full_base,
                                 grad_buf, θ_active_buf, active, all_active)

    rest = _prepare_term_specs(Base.tail(evaluators), Base.tail(backends),
                               Base.tail(active_indices_tuple), Base.tail(workspaces),
                               Base.tail(tcfgs),
                               θ_example, θ_full_base)
    return (spec, rest...)
end

# Enzyme backend: DI with tcfg as DI.Constant, workspace as DI.Cache (→ Duplicated)
# DI wrappers are pre-created here and reused every gradient call to avoid per-call allocation.
function _make_term_grad_spec(eval, b::AutoEnzyme, workspace, tcfg, θ_zero, θ_full_base,
                               grad_buf, θ_active_buf, active, all_active)
    di_tcfg = DifferentiationInterface.Constant(tcfg)
    di_cache = DifferentiationInterface.Cache(workspace)
    if !all_active
        di_θ_full_base = DifferentiationInterface.Constant(θ_full_base)
        # Use zero vector for prep to ensure safe values during compilation
        prep = prepare_gradient(eval, b, θ_zero,
            DifferentiationInterface.Constant(zero(θ_full_base)),
            di_tcfg, di_cache)
    else
        di_θ_full_base = nothing
        prep = prepare_gradient(eval, b, θ_zero, di_tcfg, di_cache)
    end
    TermGradSpec(eval, b, prep, grad_buf, θ_active_buf, active, all_active,
                 workspace, tcfg, di_tcfg, di_cache, di_θ_full_base)
end

# Non-Enzyme backend (ForwardDiff, FiniteDiff, etc.): DI with allocating eval (no workspace)
# DI wrappers are pre-created here and reused every gradient call to avoid per-call allocation.
function _make_term_grad_spec(eval, b, workspace, tcfg, θ_zero, θ_full_base,
                               grad_buf, θ_active_buf, active, all_active)
    di_tcfg = DifferentiationInterface.Constant(tcfg)
    if !all_active
        di_θ_full_base = DifferentiationInterface.Constant(θ_full_base)
        prep = prepare_gradient(eval, b, θ_zero,
            DifferentiationInterface.Constant(zero(θ_full_base)),
            di_tcfg)
    else
        di_θ_full_base = nothing
        prep = prepare_gradient(eval, b, θ_zero, di_tcfg)
    end
    TermGradSpec(eval, b, prep, grad_buf, θ_active_buf, active, all_active,
                 nothing, tcfg, di_tcfg, nothing, di_θ_full_base)
end

"""
    _accumulate_term_gradients!(total_grad, ll, specs::Tuple, θ) -> ll_total

Recursively evaluate each term's gradient and accumulate into `total_grad`.
Type-stable because each recursive call peels off a concrete element from
the heterogeneous tuple.

DI context wrappers (`di_tcfg`, `di_cache`, `di_θ_full_base`) are pre-created
in each `TermGradSpec` at model construction time, eliminating per-call allocation.
"""
_accumulate_term_gradients!(total_grad, ll, ::Tuple{}, θ) = ll
function _accumulate_term_gradients!(total_grad, ll, specs::Tuple, θ)
    spec = first(specs)
    ll = _eval_term_gradient!(total_grad, ll, spec, θ)
    return _accumulate_term_gradients!(total_grad, ll, Base.tail(specs), θ)
end

# --- Enzyme path: DI value_and_gradient! with pre-created DI wrappers ---
function _eval_term_gradient!(total_grad, ll, spec::TermGradSpec{<:Any, <:AutoEnzyme}, θ)
    fill!(spec.grad_buf, zero(eltype(spec.grad_buf)))

    if spec.all_active
        ll_i, _ = value_and_gradient!(spec.evaluator, spec.grad_buf, spec.prep, spec.backend,
                                       θ,
                                       spec.di_tcfg,
                                       spec.di_cache)
        ll += ll_i
        @inbounds for k in eachindex(total_grad)
            total_grad[k] += spec.grad_buf[k]
        end
    else
        @inbounds for (j, idx) in enumerate(spec.active_indices)
            spec.θ_active_buf[j] = θ[idx]
        end

        ll_i, _ = value_and_gradient!(spec.evaluator, spec.grad_buf, spec.prep, spec.backend,
                                       spec.θ_active_buf,
                                       spec.di_θ_full_base,
                                       spec.di_tcfg,
                                       spec.di_cache)
        ll += ll_i

        @inbounds for (j, idx) in enumerate(spec.active_indices)
            total_grad[idx] += spec.grad_buf[j]
        end
    end
    return ll
end

# --- Non-Enzyme path (ForwardDiff, FiniteDiff, etc.): DI value_and_gradient! with pre-created wrappers ---
function _eval_term_gradient!(total_grad, ll, spec::TermGradSpec, θ)
    fill!(spec.grad_buf, zero(eltype(spec.grad_buf)))

    if spec.all_active
        ll_i, _ = value_and_gradient!(spec.evaluator, spec.grad_buf, spec.prep, spec.backend,
                                       θ,
                                       spec.di_tcfg)
        ll += ll_i
        @inbounds for k in eachindex(total_grad)
            total_grad[k] += spec.grad_buf[k]
        end
    else
        @inbounds for (j, idx) in enumerate(spec.active_indices)
            spec.θ_active_buf[j] = θ[idx]
        end

        ll_i, _ = value_and_gradient!(spec.evaluator, spec.grad_buf, spec.prep, spec.backend,
                                       spec.θ_active_buf,
                                       spec.di_θ_full_base,
                                       spec.di_tcfg)
        ll += ll_i

        @inbounds for (j, idx) in enumerate(spec.active_indices)
            total_grad[idx] += spec.grad_buf[j]
        end
    end
    return ll
end

