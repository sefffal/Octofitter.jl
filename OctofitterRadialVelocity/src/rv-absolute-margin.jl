# This observation object is the same as StarAbsoluteRVObs,
# but instead of specifying RV offsets using variables in the model,
# the RV zero point is marginalized out.
# This object only supports one instrument, but you can add as many
# observation objects as you want (one per instrument).

using Octofitter: SystemObservationContext
using ADTypes: AutoEnzyme
using Enzyme: Const, Reverse, set_runtime_activity

"""
    MarginalizedStarAbsoluteRVObs(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11);

        name="inst name",
        variables=@variables begin
            jitter ~ LogUniform(0.1, 100.0)  # RV jitter (m/s)
        end
    )

Represents a likelihood function of absolute radial velocity of a host star with analytical
marginalization over the RV zero point.
`:epoch` (mjd), `:rv` (m/s), and `:σ_rv` (m/s) are all required.
In addition to the example above, any Tables.jl compatible source can be provided.

The `jitter` variable should be defined in the variables block and represents additional
uncertainty to be added in quadrature to the formal measurement errors.

Unlike `StarAbsoluteRVLikelihood`, this version analytically marginalizes over the instrument
RV zero point. This makes it faster in most cases.
That said, if you have a specific prior you want to apply for the RV zero point, correlations between
zero points, hierarchical models, etc, you should use `StarAbsoluteRVLikelihood`.

Additionally, a gaussian process is not supported with the analytical marginalization.

!!! note
    If you don't supply a `variables` argument, the detault priors are `jitter ~ LogUniform(0.001, 100)`

"""
struct MarginalizedStarAbsoluteRVObs{TTable<:Table,TF,TDevice} <: Octofitter.AbstractObs
    table::TTable
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    name::String
    trend_function::TF
    device::TDevice
end
function MarginalizedStarAbsoluteRVObs(
    observations;
    variables::Union{Nothing,Tuple{Octofitter.Priors,Octofitter.Derived}}=nothing,
    name::String,
    trend_function=(θ_obs, epoch)->zero(typeof(epoch)),
    device=nothing,
)
    if isnothing(variables)
        variables = @variables begin
            jitter ~ LogUniform(0.001, 100)
        end
        @info "Added the following observation variables:"
        display(variables[1])
        display(variables[2])
    end
    (priors,derived)=variables
    table = Table(observations)[:,:,1]
    if !Octofitter.equal_length_cols(table)
        error("The columns in the input data do not all have the same length")
    end
    if !issubset(rv_cols, Tables.columnnames(table))
        error("Expected columns $rv_cols")
    end
    rows = map(eachrow(table)) do row′
        row = (;row′[1]..., rv=float(row′[1].rv[1]))
        return row
    end
    ii = sortperm([r.epoch for r in rows])
    table = Table(rows[ii])

    if any(>=(mjd("2050")),  table.epoch) || any(<=(mjd("1950")),  table.epoch)
        @warn "The data you entered fell outside the range year 1950 to year 2050. The expected input format is MJD (modified julian date). We suggest you double check your input data!"
    end

    return MarginalizedStarAbsoluteRVObs{typeof(table),typeof(trend_function),typeof(device)}(
        table, priors, derived, name, trend_function, device
    )
end
function Octofitter.likeobj_from_epoch_subset(obs::MarginalizedStarAbsoluteRVObs, obs_inds)
    error("""
    Data subsetting is not supported for MarginalizedStarAbsoluteRVObs.

    MarginalizedStarAbsoluteRVObs analytically marginalizes over the RV zero point,
    which couples all observations together. This makes it incompatible with
    cross-validation and PSIS-LOO methods that require independent pointwise likelihoods.

    To use cross-validation or PSIS-LOO with absolute RV data, use `StarAbsoluteRVObs`
    instead, which explicitly models the RV zero point as a parameter.
    """)
end

# Backwards compatibility aliases
const MarginalizedStarAbsoluteRVLikelihood = MarginalizedStarAbsoluteRVObs
const StarAbsoluteRVMarginLikelihood = MarginalizedStarAbsoluteRVObs

export MarginalizedStarAbsoluteRVObs, MarginalizedStarAbsoluteRVLikelihood, StarAbsoluteRVMarginLikelihood

# Use Enzyme reverse-mode AD for this observation type (when no device is set).
# MarginalizedStarAbsoluteRVObs is a good candidate because its ln_like
# is a simple loop without complex control flow.
# The explicit Nothing constraint avoids ambiguity with extension methods
# that dispatch on specific device types (e.g., Reactant.XLA.AbstractDevice).
Octofitter.ad_backend(::MarginalizedStarAbsoluteRVObs{<:Any, <:Any, Nothing}) = AutoEnzyme(mode=set_runtime_activity(Reverse), function_annotation=Const)


# In-place simulation logic for MarginalizedStarAbsoluteRVObs (performance-critical)
function Octofitter.simulate!(rv_model_buf, rvlike::MarginalizedStarAbsoluteRVObs, θ_system, θ_obs, planet_orbits::Tuple, orbit_solutions)
    L = length(rvlike.table.epoch)
    T = Octofitter._system_number_type(θ_system)
    
    # Start with trend function (no offset for marginalized version)
    rv_model_buf .= rvlike.trend_function.(Ref(θ_obs), rvlike.table.epoch)

    # Add RV contribution from all planets:
    for planet_i in eachindex(planet_orbits)
        orbit = planet_orbits[planet_i]
        planet_mass = θ_system.planets[planet_i].mass
        for epoch_i in eachindex(rvlike.table.epoch)
            rv_model_buf[epoch_i] += radvel(
                orbit_solutions[planet_i][epoch_i],
                planet_mass*Octofitter.mjup2msol
            )
        end
    end

    return (rv_model = rv_model_buf, epochs = rvlike.table.epoch)
end

# Allocating simulation logic for MarginalizedStarAbsoluteRVObs (convenience method)
function Octofitter.simulate(rvlike::MarginalizedStarAbsoluteRVObs, θ_system, θ_obs, planet_orbits::Tuple, orbit_solutions)
    T = Octofitter._system_number_type(θ_system)
    L = length(rvlike.table.epoch)
    rv_model_buf = Vector{T}(undef, L)
    return Octofitter.simulate!(rv_model_buf, rvlike, θ_system, θ_obs, planet_orbits, orbit_solutions)
end


"""
Absolute radial velocity likelihood (for a star) with analytical marginalization over RV zero point.

Uses `orbitsolve_bulk` for vectorized orbit solving. All array operations use
broadcasting for Reactant/XLA traceability.
"""
function Octofitter.ln_like(
    rvlike::MarginalizedStarAbsoluteRVObs,
    ctx::SystemObservationContext
)
    (; θ_system, θ_obs, orbits) = ctx
    epochs = rvlike.table.epoch
    jitter = θ_obs.jitter

    # Bulk solve + RV for each planet (vectorized)
    # Initialize rv_model from the first planet to establish the correct element type
    # (TracedRNumber during Reactant tracing, Float64 for Enzyme)
    rv_model = let
        sol = orbitsolve_bulk(orbits[1], epochs)
        planet_mass = θ_system.planets[1].mass
        radvel(sol, planet_mass * Octofitter.mjup2msol)
    end
    for planet_i in 2:length(orbits)
        sol = orbitsolve_bulk(orbits[planet_i], epochs)
        planet_mass = θ_system.planets[planet_i].mass
        rv_model = rv_model .+ radvel(sol, planet_mass * Octofitter.mjup2msol)
    end

    # Residuals (vectorized)
    resid = rvlike.table.rv .- rv_model

    # Variance per observation (vectorized)
    var = rvlike.table.σ_rv .^ 2 .+ jitter^2

    # Marginalize out the instrument zero point (Orvara paper)
    A = sum(1.0 ./ var)
    B = -2.0 * sum(resid ./ var)
    C = sum(resid .^ 2 ./ var)
    ll = -sum(log.(2π .* var))
    ll -= -B^2 / (4A) + C + log(A)

    return ll
end



# Generate new radial velocity observations for a star
function Octofitter.generate_from_params(like::MarginalizedStarAbsoluteRVObs, ctx::SystemObservationContext; add_noise)
   (; θ_system, θ_obs, orbits, orbit_solutions) = ctx
    # Get epochs and uncertainties from observations
    epochs = like.table.epoch
    σ_rvs = like.table.σ_rv

    # Use the same simulation method as ln_like to generate model RV values
    sim = Octofitter.simulate(like, θ_system, θ_obs, orbits, orbit_solutions)

    # For generate_from_params, we add a zero offset (since marginalized version doesn't have offset parameter)
    rvs = sim.rv_model .+ 0.0

    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    if add_noise
        jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : 0.0
        radvel_table.rv .+= randn.() .* hypot.(σ_rvs, jitter)
    end

    return MarginalizedStarAbsoluteRVObs(
        radvel_table;
        name=likelihoodname(like),
        trend_function=like.trend_function,
        variables=(like.priors, like.derived),
        device=like.device,
    )
end

