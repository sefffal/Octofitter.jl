# ---------------------------------------------------
# Chains ↔ parameter NamedTuples
#
# The v2 parameter structure is exactly two levels deep:
#
#     (; <system vars>…, bodies = (; <node> = (; <var>…)…),
#                        observations = (; <obs> = (; <var>…)…))
#
# v1's was three, because observations could hang off a companion
# (`planets.b.observations.GPI.jitter`), and every conversion here carried a
# hand-written case for each depth plus prefix-matching to tell a planet
# variable from a planet-observation variable. Flattening the model flattened
# these functions too: one generic pass over the two groups replaces all of
# it.
#
# Chain column names stay `<owner>_<var>`, so a chain from a v1 fit of the
# same model still reads the same.
# ---------------------------------------------------

const _NT_GROUPS = (:bodies, :observations)

"""
    flatten_named_tuple(nt)

Flatten one parameter sample into the flat `<owner>_<var>` NamedTuple used
as a chain row. Vector-valued variables expand to `<name>_1`, `<name>_2`, ….
"""
function flatten_named_tuple(nt)
    pairs = Pair{Symbol,Float64}[]
    _flatten_into!(pairs, nt, nothing)
    for group in _NT_GROUPS
        for owner in keys(get(nt, group, (;)))
            _flatten_into!(pairs, nt[group][owner], owner)
        end
    end
    return namedtuple(pairs)
end

function _flatten_into!(pairs, ns, prefix)
    for key in keys(ns)
        key in _NT_GROUPS && continue
        name = prefix === nothing ? key : Symbol(prefix, '_', key)
        v = ns[key]
        if v isa Number
            push!(pairs, name => v)
        elseif v isa AbstractArray || v isa Tuple
            for i in eachindex(v)
                v[i] isa Number && push!(pairs, Symbol(name, '_', i) => v[i])
            end
        end
    end
    return pairs
end

"""
    result2mcmcchain(samples, sectionmap=Dict())

Convert a vector of parameter NamedTuples into an `MCMCChains.Chains`.
"""
function result2mcmcchain(chain_in, sectionmap=Dict())
    rows = map(flatten_named_tuple, chain_in)
    labels = keys(first(rows))
    data = zeros(length(rows), length(labels))
    for (i, row) in enumerate(rows), (j, l) in enumerate(labels)
        data[i, j] = row[l]
    end
    return Chains(data, [string(l) for l in labels], sectionmap)
end

# MCMCChains v7 no longer defines `haskey` for `Chains`, which Octofitter
# relies on to test whether a parameter is present. Restore it with a single
# method rather than clobbering `Base.haskey` for Dicts/NamedTuples.
Base.haskey(chain::MCMCChains.Chains, key) = key ∈ names(chain)

"""
    mcmcchain2result(model, chain, ii=(:))

The inverse of [`result2mcmcchain`](@ref): rebuild the nested parameter
NamedTuple(s) for the requested sample index (or all of them).
"""
function mcmcchain2result(model, chain, ii=(:))
    # A template sample tells us the shape: which names exist in each
    # namespace, and which of them are vector-valued.
    template = model.arr2nt(collect(model.sample_priors(Random.default_rng())))

    # (owner-or-nothing, variable name, chain column(s))
    plan = Tuple{Union{Nothing,Symbol},Symbol,Vector{Symbol}}[]
    _plan_into!(plan, template, nothing)
    for group in _NT_GROUPS, owner in keys(get(template, group, (;)))
        _plan_into!(plan, template[group][owner], owner)
    end

    read(i, j, col) = haskey(chain, col) ? chain[i, col, j] : missing
    getvalue(i, j, cols) = length(cols) == 1 ? read(i, j, cols[1]) :
                           [read(i, j, c) for c in cols]

    function reform((i, j))
        sysvals = Pair{Symbol,Any}[]
        groups = Dict{Symbol,Vector{Pair{Symbol,Any}}}(
            g => Pair{Symbol,Any}[] for g in _NT_GROUPS)
        owners = Dict{Symbol,Dict{Symbol,Vector{Pair{Symbol,Any}}}}(
            g => Dict{Symbol,Vector{Pair{Symbol,Any}}}() for g in _NT_GROUPS)
        for (owner, key, cols) in plan
            v = getvalue(i, j, cols)
            if owner === nothing
                push!(sysvals, key => v)
            else
                g = haskey(get(template, :bodies, (;)), owner) ? :bodies : :observations
                push!(get!(owners[g], owner, Pair{Symbol,Any}[]), key => v)
            end
        end
        result = namedtuple(sysvals)
        for g in _NT_GROUPS
            sub = (; (o => namedtuple(vs) for (o, vs) in owners[g])...)
            result = merge(result, NamedTuple{(g,)}((sub,)))
        end
        return result
    end

    IIs = broadcast((i, j) -> (i, j), 1:size(chain, 1), (1:size(chain, 3))')
    return ii isa Number ? reform(IIs[ii]) : broadcast(reform, IIs[ii])
end

function _plan_into!(plan, ns, owner)
    for key in keys(ns)
        key in _NT_GROUPS && continue
        name = owner === nothing ? key : Symbol(owner, '_', key)
        v = ns[key]
        if v isa Number
            push!(plan, (owner, key, [name]))
        elseif v isa AbstractArray || v isa Tuple
            push!(plan, (owner, key, [Symbol(name, '_', i) for i in eachindex(v)]))
        end
    end
    return plan
end

# ---------------------------------------------------
# Rebuilding the orbital system from a posterior sample
# ---------------------------------------------------

"""
    construct_system(model, chain, i)
    construct_system(model, θ)

The `PlanetOrbits.System` for posterior sample `i` — every body, every
hierarchy row, and the frame, exactly as the likelihood saw them.

This replaces v1's `construct_elements(chain, :b, i)`. There is no
per-planet orbit object to hand back any more: a companion's motion is a
property of the *system*, which is the whole point of the change. Query it
the same way a likelihood does:

    sys  = construct_system(model, chain, 1)
    traj = orbitsolve(sys, epochs)
    raoff.(traj, :b, :A)
"""
construct_system(model, θ::NamedTuple) = _build_po_system(model.system, θ)
construct_system(model, chain::Chains, i::Number) =
    construct_system(model, mcmcchain2result(model, chain, i))
construct_system(model, chain::Chains, ii=(:)) =
    map(θ -> construct_system(model, θ), mcmcchain2result(model, chain, ii))
export construct_system

# Non-generated counterpart of the builder `make_ln_like` compiles, for use
# outside the hot loop. Kept as a separate, obvious implementation rather
# than reusing the RGF: this one runs once per plot, and being readable
# matters more than being fast.
function _build_po_system(sys::System, θ)
    T = _system_number_type(θ)
    pobodies = map(bodynodes(sys)) do n
        nsp = θ.bodies[n.name]
        mass = hasproperty(nsp, :mass) ? nsp.mass : zero(T)
        fl = _flux_vars(n)
        flux = (; (_flux_band(f) => nsp[f] for f in fl)...)
        return PlanetOrbits.Body(; mass, flux, name=n.name)
    end
    byname = Dict(zip(sys.bodynames, pobodies))
    spec(members) = length(members) == 1 ? byname[members[1]] :
                    Tuple(byname[m] for m in members)
    poorbits = map(sys.rows) do (owner, ext, int)
        node = _node(sys, owner)
        nsp = θ.bodies[owner]
        kws = (; (e => nsp[e] for e in _element_vars(node))...)
        return PlanetOrbits.Orbit(spec(ext); about=spec(int), kws...)
    end
    frame = (; (f => θ[f] for f in sys.framevars)...)
    return PlanetOrbits.System(Tuple(pobodies), Tuple(poorbits); frame...)
end

# ---------------------------------------------------

sample_priors(sys::System) = sample_priors(Random.default_rng(), sys)
sample_priors(rng::Random.AbstractRNG, sys::System) =
    collect(make_prior_sampler(sys)(rng))
sample_priors(rng::Random.AbstractRNG, sys::System, N::Number) =
    [sample_priors(rng, sys) for _ in 1:N]
sample_priors(model, N::Number) = sample_priors(Random.default_rng(), model, N)
sample_priors(rng::Random.AbstractRNG, model, N::Number) =
    [collect(model.sample_priors(rng)) for _ in 1:N]
export sample_priors

# Helper for displaying nested named tuples compactly.
stringify_nested_named_tuple(val::Any) = string(val)
stringify_nested_named_tuple(num::Number) = string(round(num, digits=1))
stringify_nested_named_tuple(nt::NamedTuple) =
    "(;" * join(("$k=" * stringify_nested_named_tuple(nt[k]) for k in keys(nt)), ", ") * ")"
