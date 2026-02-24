"""
    OctofitterReactantExt

Generic Reactant/XLA compilation extension for Octofitter. Knows nothing about
specific observation types — works for any obs that provides a vectorized `ln_like`.

Provides:
- `ReactantGradSpec` — compiled gradient spec using `@compile`d Enzyme.gradient
- Override of `_make_term_grad_spec` for `ReactantBackend`
- Override of `_accumulate_term_gradients!` for `ReactantGradSpec`
- Callable for `ReactantTermEvaluator` using `@allowscalar arr2nt`
"""
module OctofitterReactantExt

using Octofitter
using Octofitter: ReactantBackend, ReactantTermEvaluator,
                  _n_planets, _get_obs_params, _make_obs_context,
                  _embed_active_broadcast, ln_like,
                  _make_term_grad_spec
using Reactant
using Reactant: @compile, @allowscalar, ConcreteRArray
using Enzyme
using Distributions: Distribution

# ── Reactant tracing overrides ────────────────────────────────────────
# Tell Reactant to treat these types as opaque constants (not traced).
# Without these, Reactant fails with "Unhandled type Distribution" when
# trying to trace the ReactantTermEvaluator struct hierarchy, which contains
# observation objects with Distribution priors.

function Reactant.traced_type_inner(::Type{T}, seen::Dict{Type,Type},
        mode::Reactant.TraceMode, track_numbers::Type, sharding, runtime) where {T<:Distribution}
    T
end
function Reactant.traced_type_inner(::Type{T}, seen::Dict{Type,Type},
        mode::Reactant.TraceMode, track_numbers::Type, sharding, runtime) where {T<:Octofitter.AbstractObs}
    T
end
function Reactant.traced_type_inner(::Type{T}, seen::Dict{Type,Type},
        mode::Reactant.TraceMode, track_numbers::Type, sharding, runtime) where {T<:Octofitter.Priors}
    T
end
function Reactant.traced_type_inner(::Type{T}, seen::Dict{Type,Type},
        mode::Reactant.TraceMode, track_numbers::Type, sharding, runtime) where {T<:Octofitter.Derived}
    T
end
# TypedTables.Table — observational data (constants, not parameters to differentiate)
using TypedTables: Table
function Reactant.traced_type_inner(::Type{T}, seen::Dict{Type,Type},
        mode::Reactant.TraceMode, track_numbers::Type, sharding, runtime) where {T<:Table}
    T
end

# ── ReactantTermEvaluator callable ────────────────────────────────────

"""
Evaluate a single term using Reactant-traceable operations.
Uses `@allowscalar` around `arr2nt` for parameter extraction,
and `_solve_all_orbits_bulk` for vectorized orbit solving.
"""
function (te::ReactantTermEvaluator)(θ)
    mcfg = te.mcfg
    tcfg = te.tcfg
    # Both invlink and arr2nt use scalar indexing on the parameter vector
    # (RuntimeGeneratedFunctions that access arr[i]). @allowscalar permits this
    # during Reactant tracing — the MLIR compiler optimizes away the overhead.
    θ_natural = @allowscalar mcfg.invlink(θ)
    θ_system  = @allowscalar mcfg.arr2nt(θ_natural)
    orbits    = @allowscalar mcfg.construct_orbits(θ_system)
    θ_planet, θ_obs = _get_obs_params(tcfg, θ_system)
    # Pass orbits to ln_like via context — orbit solving is done inside ln_like
    # itself (e.g., via orbitsolve_bulk for vectorized Reactant-traceable evaluation).
    # We pass empty orbit_solutions; the Reactant ln_like method must not rely on them.
    orbit_solutions = ntuple(Returns(()), _n_planets(mcfg))
    ctx = _make_obs_context(tcfg, θ_system, θ_planet, θ_obs, orbits, orbit_solutions)
    ln_like(tcfg.obs, ctx)
end

# ── ReactantGradSpec ──────────────────────────────────────────────────

"""
    ReactantGradSpec

Pre-compiled gradient spec for a Reactant-compiled term.
Stores a single compiled thunk (from `@compile Enzyme.gradient(ReverseWithPrimal, ...)`)
that returns both the primal value and gradient in one pass.

The `eval_fn` field stores the callable that was compiled (either a
`ReactantTermEvaluator` for all-active, or a sparse closure for sparse terms).
This must be re-passed to the compiled gradient thunk at call time per the
`@compile` calling convention.
"""
struct ReactantGradSpec{TCompiledGrad, TEvalFn}
    compiled_grad::TCompiledGrad
    eval_fn::TEvalFn
    grad_buf::Vector{Float64}
    active_indices::Vector{Int}
    all_active::Bool
end

# ── _make_term_grad_spec override for ReactantBackend ─────────────────

function Octofitter._make_term_grad_spec(eval, b::ReactantBackend, workspace, tcfg, θ_zero, θ_full_base,
                                          grad_buf, θ_active_buf, active, all_active)
    D_active = length(active)
    D_total = length(θ_full_base)

    # Extract mcfg from the evaluator; tcfg is passed separately
    mcfg = eval.mcfg

    # Build ReactantTermEvaluator (no mutable workspace needed)
    rt_eval = ReactantTermEvaluator(mcfg, tcfg)

    if all_active
        # Example input as ConcreteRArray
        θ_ra = ConcreteRArray(zeros(Float64, D_total))

        # Compile gradient+primal: θ → (gradient(ll), ll) in a single pass
        compiled_grad = @compile Enzyme.gradient(Enzyme.ReverseWithPrimal, rt_eval, θ_ra)

        eval_fn = rt_eval
    else
        # Sparse case: build a callable struct that embeds active params into full vector
        sparse_eval = SparseReactantEval(rt_eval, copy(θ_full_base), collect(Int, active))

        θ_ra = ConcreteRArray(zeros(Float64, D_active))

        compiled_grad = @compile Enzyme.gradient(Enzyme.ReverseWithPrimal, sparse_eval, θ_ra)

        eval_fn = sparse_eval
    end

    grad_buf = zeros(Float64, D_active)
    ReactantGradSpec(compiled_grad, eval_fn, grad_buf, collect(Int, active), all_active)
end

"""
    SparseReactantEval

Callable struct for sparse Reactant evaluation. Embeds active parameters into
the full parameter vector before calling the ReactantTermEvaluator.
Uses a struct instead of a closure so it can be re-passed to compiled thunks.
"""
struct SparseReactantEval{TEval, TBase}
    rt_eval::TEval
    θ_base::TBase
    active_indices::Vector{Int}
end
function (se::SparseReactantEval)(θ_active)
    θ_full = @allowscalar _embed_active_broadcast(se.θ_base, θ_active, se.active_indices)
    se.rt_eval(θ_full)
end

# ── _accumulate_term_gradients! override for ReactantGradSpec ─────────

function Octofitter._accumulate_term_gradients!(
    total_grad, ll, specs::Tuple{<:ReactantGradSpec, Vararg}, θ
)
    spec = first(specs)

    # Convert active params to ConcreteRArray
    θ_input = if spec.all_active
        ConcreteRArray(collect(Float64, θ))
    else
        ConcreteRArray(collect(Float64, θ[spec.active_indices]))
    end

    # Call compiled gradient+primal — returns (gradient, primal) in one pass
    result = spec.compiled_grad(Enzyme.ReverseWithPrimal, spec.eval_fn, θ_input)

    # Extract gradient and primal value
    # Enzyme.gradient(ReverseWithPrimal, f, x) returns ((df/dx,), f(x))
    # result[1] is a tuple of gradients (one per positional arg), result[2] is the primal
    grad_jl = Array(result[1][1])
    ll += Float64(result[2])

    # Scatter gradient to total
    if spec.all_active
        @inbounds for k in eachindex(total_grad)
            total_grad[k] += grad_jl[k]
        end
    else
        @inbounds for (j, idx) in enumerate(spec.active_indices)
            total_grad[idx] += grad_jl[j]
        end
    end

    return Octofitter._accumulate_term_gradients!(total_grad, ll, Base.tail(specs), θ)
end

end # module OctofitterReactantExt
