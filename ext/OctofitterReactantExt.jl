"""
    OctofitterReactantExt

Generic Reactant/XLA compilation extension for Octofitter. Knows nothing about
specific observation types — works for any obs that provides a vectorized `ln_like`.

Provides:
- `ReactantGradSpec` — compiled gradient spec using `@compile`d Enzyme.gradient
- Override of `_prepare_one_term` for `ReactantBackend`
- Override of `_accumulate_term_gradients!` for `ReactantGradSpec`
- Callable for `ReactantTermEvaluator` using `@allowscalar arr2nt`
"""
module OctofitterReactantExt

using Octofitter
using Octofitter: ReactantBackend, ReactantTermEvaluator,
                  _n_planets, _get_obs_params, _make_obs_context,
                  _solve_all_orbits_bulk, _embed_active_broadcast, ln_like
using Reactant
using Reactant: @compile, @allowscalar, ConcreteRArray
using Enzyme
using Distributions: Distribution

# ── Reactant math overrides ──────────────────────────────────────────
# LogExpFunctions functions used by Bijectors.jl invlink transforms.
# Reactant TracedRNumber doesn't have methods for these, so we provide
# fallbacks using basic arithmetic that Reactant can trace.
# Access LogExpFunctions through Distributions (which loads it transitively).
const _LogExpFunctions = let
    pkgid = Base.PkgId(Base.UUID("2ab3a3ac-af41-5b50-aa03-7779005ae688"), "LogExpFunctions")
    Base.loaded_modules[pkgid]
end

_LogExpFunctions.logistic(x::Reactant.TracedRNumber) = one(x) / (one(x) + exp(-x))
_LogExpFunctions.logit(x::Reactant.TracedRNumber) = log(x / (one(x) - x))

# Workaround for Reactant bug: stablehlo.atan2 gradient is sign-flipped.
# Use the half-angle identity: atan(y,x) = 2*atan(y/(sqrt(x²+y²)+x))
# which decomposes into 1-arg atan (correct gradient) + standard arithmetic.
# Valid for all (x, y) except the negative real axis (x < 0, y = 0), which
# does not occur in our orbital mechanics use case (true anomaly computation).
Base.atan(y::Reactant.TracedRNumber, x::Reactant.TracedRNumber) =
    2 * atan(y / (sqrt(x^2 + y^2) + x))

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
    orbits    = ntuple(i -> mcfg.orbit_constructors[i](θ_system), _n_planets(mcfg))
    θ_planet, θ_obs = _get_obs_params(tcfg, θ_system)
    orbit_solutions = _solve_all_orbits_bulk(orbits, tcfg.epochs)
    ctx = _make_obs_context(tcfg, θ_system, θ_planet, θ_obs, orbits, orbit_solutions)
    ln_like(tcfg.obs, ctx)
end

# ── ReactantGradSpec ──────────────────────────────────────────────────

"""
    ReactantGradSpec

Pre-compiled gradient spec for a Reactant-compiled term.
Stores a compiled gradient thunk (from `@compile Enzyme.gradient(...)`)
and a compiled forward evaluation thunk.

The `eval_fn` field stores the callable that was compiled (either a
`ReactantTermEvaluator` for all-active, or a sparse closure for sparse terms).
This must be re-passed to the compiled gradient thunk at call time per the
`@compile` calling convention.
"""
struct ReactantGradSpec{TCompiledGrad, TCompiledFwd, TEvalFn}
    compiled_grad::TCompiledGrad
    compiled_fwd::TCompiledFwd
    eval_fn::TEvalFn
    grad_buf::Vector{Float64}
    active_indices::Vector{Int}
    all_active::Bool
end

# ── _prepare_one_term override for ReactantBackend ────────────────────

function Octofitter._prepare_one_term(eval, backend::ReactantBackend, active, θ_example, θ_full_example)
    D_active = length(active)
    D_total = length(θ_example)
    all_active = D_active == D_total

    # Extract mcfg and tcfg from the evaluator
    mcfg = eval.mcfg
    tcfg = eval.tcfg

    # Build ReactantTermEvaluator (no mutable workspace needed)
    rt_eval = ReactantTermEvaluator(mcfg, tcfg)

    if all_active
        # Example input as ConcreteRArray
        θ_ra = ConcreteRArray(zeros(Float64, D_total))

        # Compile gradient: θ → gradient(ll) w.r.t. θ
        # @compile calling convention: compiled(Reverse, f, args...)
        compiled_grad = @compile Enzyme.gradient(Enzyme.Reverse, rt_eval, θ_ra)

        # Compile forward pass: θ → ll value
        compiled_fwd = @compile rt_eval(θ_ra)

        eval_fn = rt_eval
    else
        # Sparse case: build a callable struct that embeds active params into full vector
        sparse_eval = SparseReactantEval(rt_eval, copy(θ_full_example), collect(Int, active))

        θ_ra = ConcreteRArray(zeros(Float64, D_active))

        compiled_grad = @compile Enzyme.gradient(Enzyme.Reverse, sparse_eval, θ_ra)
        compiled_fwd = @compile sparse_eval(θ_ra)

        eval_fn = sparse_eval
    end

    grad_buf = zeros(Float64, D_active)
    ReactantGradSpec(compiled_grad, compiled_fwd, eval_fn, grad_buf, collect(Int, active), all_active)
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
    total_grad, ll, specs::Tuple{<:ReactantGradSpec, Vararg}, θ, θ_full_base
)
    spec = first(specs)

    # Convert active params to ConcreteRArray
    θ_input = if spec.all_active
        ConcreteRArray(collect(Float64, θ))
    else
        ConcreteRArray(collect(Float64, θ[spec.active_indices]))
    end

    # Call compiled gradient — must pass Reverse + eval_fn + data
    grad_result = spec.compiled_grad(Enzyme.Reverse, spec.eval_fn, θ_input)

    # Extract gradient as Julia array
    grad_jl = Array(grad_result[1])

    # Forward pass for ll value
    # @compile f(args) calling convention: compiled(args) — f is NOT re-passed
    ll_ra = spec.compiled_fwd(θ_input)
    # ll_ra is a ConcreteRNumber (scalar) — convert directly to Float64
    ll_val = Float64(ll_ra)
    ll += ll_val

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

    return Octofitter._accumulate_term_gradients!(total_grad, ll, Base.tail(specs), θ, θ_full_base)
end

end # module OctofitterReactantExt
