using MCMCChains: MCMCChains
using MCMCChains: Chains


export sample_priors
sample_priors(planet::Planet) = rand.(ComponentArray(planet.priors.priors))
# sample_priors(planet::Planet,N) = rand.(planet.priors.priors,N)
sample_priors(planet::Planet,N) = [sample_priors(planet) for _ in 1:N]

# function sample_priors(planet::ReparameterizedPlanet)
#     sample = NamedTuple(sample_priors(planet.planet))

#     plus = sample.ω + sample.Ω
#     minus = sample.ω - sample.Ω
#     reparameterized = merge(
#         delete(sample,  (:ω, :Ω)),
#         (ωΩ⁺ = plus, ωΩ⁻=minus)
#     )
    
#     # Φ = sample.ω + sample.Ω + 2π*sample.τ
#     # Φω⁻ = Φ - sample.ω
#     # ΦΩ⁻ = Φ - sample.Ω
#     # reparameterized = merge(
#     #     delete(sample,  (:τ, :ω, :Ω)),
#     #     (;Φ, Φω⁻, ΦΩ⁻)
#     # )

#     return ComponentVector(reparameterized)
# end
# function sample_priors(planet::ReparameterizedPlanet, N)
#     sample = NamedTuple(sample_priors(planet.planet, N))
#     plus = sample.ω .+ sample.Ω
#     minus = sample.ω .- sample.Ω
#     reparameterized = merge(
#         delete(sample,  (:ω, :Ω)),
#         (ωΩ⁺ = plus, ωΩ⁻=minus)
#     )

#     # Φ = sample.ω .+ sample.Ω .+ 2π*sample.τ
#     # Φω⁻ = Φ .- sample.ω
#     # ΦΩ⁻ = Φ .- sample.Ω
#     # reparameterized = merge(
#     #     delete(sample,  (:τ, :ω, :Ω)),
#     #     (;Φ, Φω⁻, ΦΩ⁻)
#     # )
#     return ComponentVector(reparameterized)
# end


function sample_priors(system::System)
    sampled = ComponentVector(
        merge(NamedTuple(rand.(ComponentArray(system.priors.priors))),
        # (;planets=[sample_priors(planet) for planet in system.planets])
        (;planets=namedtuple(collect(keys(system.planets)), [
            ComponentArray(NamedTuple([k=>v for (k,v) in pairs(NamedTuple(sample_priors(planet)))]))
            for planet in system.planets
        ]))
    ))
    return getdata(sampled)
end
# function sample_priors(system::System,N)
#     sampled = ComponentVector(
#         merge(NamedTuple(rand.(system.priors.priors,N)),
#         # (;planets=[sample_priors(planet,N) for planet in system.planets])
#         (;planets=[
#             ComponentArray(NamedTuple([k=>v for (k,v) in pairs(NamedTuple(sample_priors(planet,N)))]))
#             for planet in system.planets
#         ])
#     ))
#     return sampled
# end

sample_priors(system::System,N) = [sample_priors(system) for _ in 1:N]



# TODO: there should be a more efficient way to do this.
# Worst case, can unroll using (runtime) generated function?
# It appears that the bottleneck is actually converting back and forth to NamedTuples
# Can we do this another way? Do we really need component arrays?
# They do make the API quite nice but we're using them in fewer and fewer spots.
# They are basically an implementation detail at this point.
# We could write our own function to go from input array to named tuple?




function resolve_deterministic(system::System, θ)
    if isnothing(system.deterministic)
        θ_resolved = NamedTuple(θ)
    else
        resolved_values = map(values(system.deterministic.variables)) do func
            val = func(θ)
        end
        # for (key, func) in pairs(system.deterministic.variables)
        #     val = func(θ)
        #     θ_resolved = merge(θ_resolved, NamedTuple{(key,), Tuple{eltype(θ)}}(val))
        # end
        θ_resolved = merge(
            NamedTuple(θ),
            namedtuple(keys(system.deterministic.variables), resolved_values)
        )
    end
    # return θ_resolved

    resolved_planets = map(keys(system.planets)) do key
        θ_planet = θ.planets[key]
        planet_model = system.planets[key]

        θ_planet_resolved = NamedTuple()
        if !isnothing(planet_model.deterministic)
            for (key, func) in pairs(planet_model.deterministic.variables)
                val = func(θ_resolved, θ_planet)
                θ_planet_resolved = merge(θ_planet_resolved, NamedTuple{(key,), Tuple{eltype(θ)}}(val))
            end
        end
        return merge(θ_planet_resolved, NamedTuple(θ_planet))
    end
    resolved_planets_nt = namedtuple([pl.name for pl in system.planets], resolved_planets)
    θ_resolved_complete = merge(θ_resolved, (;planets=resolved_planets_nt))
    return ComponentVector(θ_resolved_complete)
    # return θ_resolved_complete
end

# # _ntkeys(nt::Type{NamedTuple{TKeys,<:Any}}) where TKeys = TKeys

# _determinekeysvals(::Type{Deterministic{NamedTuple{Keys,Vals}}}) where {Keys,Vals} = Keys, Vals

# function _planetkeysvals(::Type{NamedTuple{Keys,Vals}}) where {Keys,Vals} 
#     Keys, Vals, findfirst(==(:planets), Keys)
# end

# @generated function resolve_deterministic(system::System{TDetermine}, θ) where TDetermine


#     if TDetermine != Nothing
#         keys, funcs = _determinekeysvals(TDetermine)
        
#         body = Expr[]
#         for key in keys
#             # ex = :($key = system.deterministic.variables.$key(θ))
#             ex = :(system.deterministic.variables.$key(θ))
#             push!(body, ex)
#         end

#         # Now need to expand out each planet as well.
#         # Need to pass in (planet, θ_sys, θ_pl)
#         planets = 

#         NTType = NamedTuple{keys}
#         ex = :(
#             merge(
#                 $NTType($(body...)),
#                 NamedTuple(θ),
#                 planets = 
#             )
#         )

#         return ex

#     end
#     return θ
#     # resolved_planets = map(zip(θ.planets, system.planets)) do (θ_planet, planet_model)
#     # resolved_planets = map(keys(system.planets)) do key
#     #     θ_planet = θ.planets[key]
#     #     planet_model = system.planets[key]

#     #     θ_planet_resolved = NamedTuple()
#     #     if !isnothing(planet_model.deterministic)
#     #         for (key, func) in pairs(planet_model.deterministic.variables)
#     #             val = func(θ_resolved, θ_planet)
#     #             θ_planet_resolved = merge(θ_planet_resolved, NamedTuple{(key,), Tuple{eltype(θ)}}(val))
#     #         end
#     #     end
#     #     return merge(θ_planet_resolved, NamedTuple(θ_planet))
#     # end
#     # resolved_planets_nt = namedtuple([pl.name for pl in system.planets], resolved_planets)
#     # θ_resolved = merge(θ_resolved, (;planets=resolved_planets_nt))
#     # return ComponentVector(θ)
# end


@generated function resolve_deterministic(planet::Planet{TDetermine}, θ_sys, θ_pl) where TDetermine
    if TDetermine != Nothing
        keys, funcs = _determinekeysvals(TDetermine)
        body = Expr[]
        for key in keys
            ex = :(planet.deterministic.variables.$key(θ_sys,θ_pl))
            push!(body, ex)
        end
        NTType = NamedTuple{keys}
        ex = :(
            merge($NTType($(body...)), NamedTuple(θ))
        )
        return ex
    end
    return θ
end


function make_resolve_deterministic(planet::Planet)

    body = Expr[]

    dkeys = keys(planet.deterministic.variables)
    funcs = values(planet.deterministic.variables)
    for (key,func) in zip(dkeys, funcs)
        ex = :($key = $func(θ_sys, θ_pl))
        push!(body, ex)
    end

    ex = :(function (θ_sys, θ_pl)
        return merge(
            (;$(body...)),
            NamedTuple(θ_pl)
        )
    end)

    ln_prior = @RuntimeGeneratedFunction(ex)
    return ln_prior
end
function make_resolve_deterministic(system::System)

    body = Expr[]

    dkeys = keys(system.deterministic.variables)
    funcs = values(system.deterministic.variables)
    for (key,func) in zip(dkeys, funcs)
        ex = :($key = $func(θ))
        push!(body, ex)
    end

    ex = :(function (θ)
        return merge(
            (;$(body...)),
            NamedTuple(θ)
        )
    end)

    ln_prior = @RuntimeGeneratedFunction(ex)
    return ln_prior
end




# Instead of just calling mean for the distributions, we sample and then take the mean of the data.
# This does add a little jitter, but some distributions do not directly define the mean function!
# Specifically, truncated(InverseGamma()) does not work, and this is very useful.
# mean_priors(planet::Planet) = Statistics.mean.(Statistics.rand.(planet.priors.priors,1000))
function mean_priors(system::System)
    priors_all = ComponentVector(;
        NamedTuple(system.priors.priors)...,
        planets=[planet.priors.priors for planet in system.planets]
    )
    # return Statistics.mean.(Statistics.rand.(priors_all,1000))
    return Statistics.mean.(priors_all)
end
# mean_priors(planet::Planet) = Statistics.mean.(Statistics.rand.(planet.priors.priors,1000))
# function mean_priors(system::System)
#     N = 5000
#     sampled = ComponentVector(
#         merge(NamedTuple(mean.(rand.(system.priors.priors,N))),
#         # (;planets=[mean.(sample_priors(planet,N)) for planet in system.planets])
#         (;planets=[
#             ComponentArray(NamedTuple([k=>mean(v) for (k,v) in pairs(NamedTuple(sample_priors(planet,N)))]))
#             for planet in system.planets
#         ])
#     ))
#     return sampled
# end


function guess_starting_position(system, N=500_000)

    @info "Guessing a good starting location by sampling from priors" N
    # TODO: this shouldn't have to allocate anything, we can just loop keeping the best.
    θ = sample_priors(system, N)
    arr2nt = DirectDetections.make_arr2nt(system) 
    θr = arr2nt.(θ)

    posts = zeros(N)
    ln_prior = make_ln_prior(system)
    arr2nt = DirectDetections.make_arr2nt(system) 
    Threads.@threads for i in eachindex(posts)
        θ_res = arr2nt(θ[i])
        posts[i] = ln_prior(θ[i]) + ln_like(θ_res, system)
    end
    # posts = map(eachrow(A)) do c
    #     DirectDetections.ln_post(ComponentVector(c, ax), system)
    # end
    mapv,mapi = findmax(posts)
    best = θ[mapi]
    
    # @info "Found good location" mapv best=NamedTuple(best)

    return best
end



function get_ωΩτ(θ_system, θ_planet)
    if haskey(θ_planet, :ωΩ⁺)
        
        # Derivation:
        # ωΩ⁺ = ω+Ω
        # ωΩ⁻ = ω-Ω
        # ωΩ⁺ + ωΩ⁻ = ω+Ω + ω-Ω = ω
        # ωΩ⁺ - ωΩ⁻ = ω+Ω - ω+Ω = Ω
        ω = θ_planet.ωΩ⁺ + θ_planet.ωΩ⁻ 
        Ω = θ_planet.ωΩ⁺ - θ_planet.ωΩ⁻
        τ = θ_planet.τ
    elseif haskey(θ_planet, :Φ)
    
        # Φ = sample.ω .+ sample.Ω .+ 2π*sample.τ
        # Φω⁻ = Φ .- sample.ω
        # ΦΩ⁻ = Φ .- sample.Ω
    
        τ2pi = θ_planet.Φ + θ_planet.Φω⁻ + θ_planet.ΦΩ⁻
        ω = θ_planet.Φ - τ2pi + θ_planet.ΦΩ⁻
        Ω = θ_planet.Φ - τ2pi + θ_planet.Φω⁻
        τ = τ2pi/2π

    else
        ω = θ_planet.ω
        τ = θ_planet.τ

        
        if hasproperty(θ_planet,:Ω)
            Ω = θ_planet.Ω
        elseif hasproperty(θ_system,:Ω)
            Ω = θ_system.Ω
        else
            error("property `Ω` not specified for planet or system")
        end
    end
    return ω, Ω, τ
end


function construct_elements(θ_system, θ_planet)
    return KeplerianElements((;
        θ_system.μ,
        θ_system.plx,
        θ_planet.i,
        θ_planet.Ω,
        θ_planet.ω,
        θ_planet.e,
        θ_planet.τ,
        θ_planet.a,
    ))
end
function construct_elements(θ_system, θ_planet, i)
    return KeplerianElements((;
        μ=θ_system.μ[i],
        plx=θ_system.plx[i],
        i=θ_planet.i[i],
        Ω=θ_planet.Ω[i],
        ω=θ_planet.ω[i],
        e=θ_planet.e[i],
        τ=θ_planet.τ[i],
        a=θ_planet.a[i],
    ))
end

function construct_elements(chain::Chains, planet_key::Union{String,Symbol}, i::Integer)
    pk = string(planet_key)
    return KeplerianElements((;
        μ=chain["μ"][i],
        plx=chain["plx"][i],
        i=chain[pk*"[i]"][i],
        Ω=chain[pk*"[Ω]"][i],
        ω=chain[pk*"[ω]"][i],
        e=chain[pk*"[e]"][i],
        τ=chain[pk*"[τ]"][i],
        a=chain[pk*"[a]"][i],
    ))
end

function construct_elements(chain::Chains, planet_key::Union{String,Symbol}, ii::AbstractArray{<:Integer})
    pk = string(planet_key)
    μs=chain["μ"]
    plxs=chain["plx"]
    is=chain[pk*"[i]"]
    Ωs=chain[pk*"[Ω]"]
    ωs=chain[pk*"[ω]"]
    es=chain[pk*"[e]"]
    τs=chain[pk*"[τ]"]
    as=chain[pk*"[a]"]
    return map(ii) do i
        KeplerianElements((;
            μ=μs[i],
            plx=plxs[i],
            i=is[i],
            Ω=Ωs[i],
            ω=ωs[i],
            e=es[i],
            τ=τs[i],
            a=as[i],
        ))
    end
end




function mcmc(
    system::System;
    burnin,
    numwalkers,
    numsamples_perwalker,
    thinning = 1,
    squash = true,
)
    ln_post(θ) = ln_prior(θ, system) + ln_like(θ, system)
 
    # ln_prior_system_specialized = make_ln_prior(θ, system)
    # ln_post(θ) = ln_prior_system_specialized(θ, system) + ln_like(θ, system)

    @info "Finding starting point"
    initial_walkers = [sample_priors(system) for _ in 1:numwalkers]

    # Convert the initial walkers into static arrays for stack allocation.
    # This messy line should have no impact on the semantics of the code.
    initial_walkers_static = [
        ComponentVector{SVector{length(cv)}}(; NamedTuple(cv)...) for cv in initial_walkers
    ]

    # Run the MCMC
    thetase, _accept_ratioe = KissMCMC.emcee(
        ln_post,
        initial_walkers_static;
        nburnin = burnin * numwalkers,
        use_progress_meter = true,
        nthin = thinning,
        niter = numsamples_perwalker * numwalkers,
    )

    # Convert the output into an MCMCChains.Chain.
    # Use reinterpret to avoid re-allocating all that memory
    if squash
        thetase′, _ = KissMCMC.squash_walkers(thetase, _accept_ratioe)
        reinterptted = reinterpret(reshape, eltype(first(thetase′)), thetase′)
        chains = ComponentArray(collect(eachrow(reinterptted)), getaxes(thetase′[1]))
    else
        matrix_paramxstep_per_walker = [reinterpret(reshape, eltype(first(θ)), θ) for θ in thetase]
        A = reshape(
            mapreduce(θ->reinterpret(reshape, eltype(first(θ)), θ), hcat, thetase),
            (length(thetase[1][1]), :, numwalkers,)
        )
        ax = getaxes(thetase[1][1])
        chains = ComponentArray(collect(eachslice(A,dims=1)), ax)
    end

    # return Chains(reinterptted, column_names)
    return chains
end

"""
Given a component vector and a matching component vector that is a boolean mask,
return only the values that are true
"""
function select_cv(cv, mask)
    dat = getdata(cv)
    ax = getaxes(cv)
    dat, ax
    # ComponentArray(dat[mask],ax[1][mask])
    ax[1][mask]
end

# using Zygote
# https://github.com/FluxML/Zygote.jl/issues/570
# @Zygote.adjoint (T::Type{<:SArray})(x::Number...) = T(x...), y->(nothing, y...)
function hmc(
    system::System, target_accept=0.9;
    adaptation,
    iterations,
    tree_depth=10,
    include_adapatation=false,
    initial_samples=100_000,
    initial_parameters=nothing,
)

    # Choose parameter dimensionality and initial parameter value
    initial_θ_0 = sample_priors(system)
    D = length(initial_θ_0)

    ln_prior = make_ln_prior(system)
    arr2nt = DirectDetections.make_arr2nt(system) 

    # Capture these variables in a let binding to improve performance
    ℓπ = let system=system, ln_prior=ln_prior, arr2nt=arr2nt
        function (θ)
            θ_res = arr2nt(θ)
            ll = ln_prior(θ) + ln_like(θ_res, system)
            return ll
        end
    end

    if isnothing(initial_parameters)
        initial_θ_mut = guess_starting_position(system,initial_samples)
    else
        initial_θ_mut = initial_parameters
    end

    # Use a static comopnent array for efficiency
    # initial_θ = SVector{length(initial_θ_mut)}(initial_θ_mut)
    initial_θ = initial_θ_mut

    # Define a Hamiltonian system
    metric = DenseEuclideanMetric(D)
    # metric = DiagEuclideanMetric(D)
    hamiltonian = Hamiltonian(metric, ℓπ, ForwardDiff)
    # hamiltonian = Hamiltonian(metric, ℓπ, ℓπ_grad)

    initial_ϵ = find_good_stepsize(hamiltonian, getdata(initial_θ))
    # initial_ϵ = 0.005

    integrator = Leapfrog(initial_ϵ)
    # 1.05 improves the sampling over standard leapfrog, but 4.0 is too much. It just stays stuck.
    # 1.5 seems better but seems to favour certain trajectories.
    # integrator = TemperedLeapfrog(initial_ϵ, 1.05)
    proposal = NUTS(integrator, max_depth=tree_depth) 


    # # We have change some parameters when running with image data
    if !isnothing(system.images) && target_accept > 0.6
        target_accept = 0.6
        @info "Sampling from images, lowering target_accept to 0.6"
    end

    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(target_accept, integrator)) 
    # adaptor = MassMatrixAdaptor(metric)

    logger = SimpleLogger(stdout, Logging.Error)
    start_time = time()
    samples, stat = with_logger(logger) do
        drop_warmup=!(adaptor isa AdvancedHMC.NoAdaptation) || include_adapatation
        sample(
            hamiltonian,
            proposal,
            getdata(initial_θ),
            iterations,
            adaptor,
            adaptation;
            progress=true,
            drop_warmup=drop_warmup
        )
    end
    stop_time = time()
    @info "Resolving deterministic variables"
    chain_res = arr2nt.(samples)
    @info "Constructing chains"
    mcmcchain = DirectDetections.result2mcmcchain(system, chain_res)
    mcmcchain_with_info = MCMCChains.setinfo(
        mcmcchain,
        (;start_time, stop_time, model=system)
    )
    return mcmcchain_with_info, stat
end


"""
Convert a vector of component arrays returned from sampling into an MCMCChains.Chains
object.
"""
function result2mcmcchain(system, chains_in_0)
    chains_in = ComponentArray.(chains_in_0)
    # `system` not currently used, but a more efficient/robust mapping in future might require it.

    # There is a specific column name convention used by MCMCChains to indicate
    # that multiple parameters form a group. Instead of planets.X.a, we adapt our to X[a] 
    # accordingly
    flattened_labels = replace.(labels(first(chains_in)), r"planets\.([^\.]+).([^\.]+)" => s"\1[\2]")
    c = Chains(
        reduce(vcat, getdata.(chains_in)'),
        flattened_labels
    )
end

"""
    planet_keys(system, key)

For planet `key`, return the column names used in the output chains.

Examples:
```julia-repl
julia> planet_keys(system, :X)
7-element Vector{String}:
 "f[a]"
 "f[e]"
 "f[τ]"
 "f[ω]"
 "f[i]"
 "f[Ω]"
 "f[Keck_L′]"
```
"""
function planet_keys(system, key)
    # There is a specific column name convention used by MCMCChains to indicate
    # that multiple parameters form a group. Instead of planets.X.a, we adapt our to X[a] 
    # accordingly
    
    θ = sample_priors(system)
    θ_r = resolve_deterministic(system, θ)

    flattened_labels = replace.(labels(θ_r), r"planets\.([^\.]+).([^\.]+)" => s"\1[\2]")

    planet_key_start = string(key) * "["
    return filter(startswith(planet_key_start), flattened_labels)
end
export planet_keys