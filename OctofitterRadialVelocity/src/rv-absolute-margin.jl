# This likelihood object is the same as MarginalizedStarAbsoluteRVLikelihood,
# but instead of specifying RV offsets using variables in the model,
# the RV zero point is marginalized out.
# This object only supports one instrument, but you can add as many 
# likelihood objects as you want (one per instrument).

using Bumper

"""
    MarginalizedStarAbsoluteRVLikelihood(
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
struct MarginalizedStarAbsoluteRVLikelihood{TTable<:Table,TF} <: Octofitter.AbstractLikelihood
    table::TTable
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    name::String
    trend_function::TF
end
function MarginalizedStarAbsoluteRVLikelihood(
    observations;
    variables::Union{Nothing,Tuple{Octofitter.Priors,Octofitter.Derived}}=nothing,
    name::String,
    trend_function=(θ_obs, epoch)->zero(Octofitter._system_number_type(θ_obs)),
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
    
    return MarginalizedStarAbsoluteRVLikelihood{typeof(table),typeof(trend_function)}(
        table, priors, derived, name, trend_function
    )
end
function Octofitter.likeobj_from_epoch_subset(obs::MarginalizedStarAbsoluteRVLikelihood, obs_inds)
    table = Table(first(eachcol(obs.table[obs_inds])))
    return MarginalizedStarAbsoluteRVLikelihood{
        typeof(table),typeof(obs.trend_function)
    }(
        table, obs.priors, obs.derived, likelihoodname(obs), obs.trend_function
    )
end
export MarginalizedStarAbsoluteRVLikelihood


# In-place simulation logic for MarginalizedStarAbsoluteRVLikelihood (performance-critical)
function Octofitter.simulate!(rv_model_buf, rvlike::MarginalizedStarAbsoluteRVLikelihood, θ_system, θ_obs, planet_orbits::Tuple, orbit_solutions, orbit_solutions_i_epoch_start)
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
                orbit_solutions[planet_i][epoch_i+orbit_solutions_i_epoch_start],
                planet_mass*Octofitter.mjup2msol
            )
        end
    end        
    
    return (rv_model = rv_model_buf, epochs = rvlike.table.epoch)
end

# Allocating simulation logic for MarginalizedStarAbsoluteRVLikelihood (convenience method)
function Octofitter.simulate(rvlike::MarginalizedStarAbsoluteRVLikelihood, θ_system, θ_obs, planet_orbits::Tuple, orbit_solutions, orbit_solutions_i_epoch_start)
    T = Octofitter._system_number_type(θ_system)
    L = length(rvlike.table.epoch)
    rv_model_buf = Vector{T}(undef, L)
    return Octofitter.simulate!(rv_model_buf, rvlike, θ_system, θ_obs, planet_orbits, orbit_solutions, orbit_solutions_i_epoch_start)
end


"""
Absolute radial velocity likelihood (for a star) with analytical marginalization over RV zero point.
"""
function Octofitter.ln_like(
    rvlike::MarginalizedStarAbsoluteRVLikelihood,
    θ_system,
    θ_obs,
    planet_orbits::Tuple,
    orbit_solutions,
    orbit_solutions_i_epoch_start
)
    L = length(rvlike.table.epoch)
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)
    
    jitter = θ_obs.jitter

    @no_escape begin
        # Allocate buffers using bump allocator
        rv_model_buf = @alloc(T, L)
        resid = @alloc(T, L)

        # Use in-place simulation method to get model values
        sim = Octofitter.simulate!(rv_model_buf, rvlike, θ_system, θ_obs, planet_orbits, orbit_solutions, orbit_solutions_i_epoch_start)

        # Data for this instrument:
        epochs = rvlike.table.epoch
        σ_rvs = rvlike.table.σ_rv
        rvs = rvlike.table.rv

        # RV residual calculation: measured RV - model
        resid .= rvs .- rv_model_buf
        
        # Marginalize out the instrument zero point using math from the Orvara paper
        A = zero(T)
        B = zero(T)
        C = zero(T)
        for i_epoch in eachindex(epochs)
            # The noise variance per observation is the measurement noise and the jitter added
            # in quadrature
            var = σ_rvs[i_epoch]^2 + jitter^2
            A += 1/var
            B -= 2resid[i_epoch]/var
            C += resid[i_epoch]^2/var
            # Penalize RV likelihood by total variance 
            ll -= log(2π*var) #  ll += log(1/var)
        end
        ll -= -B^2 / (4A) + C + log(A)
    end

    return ll
end



# Generate new radial velocity observations for a star
function Octofitter.generate_from_params(like::MarginalizedStarAbsoluteRVLikelihood, θ_system,  θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start; add_noise)

    # Get epochs and uncertainties from observations
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 

    # Use the same simulation method as ln_like to generate model RV values
    sim = Octofitter.simulate(like, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    
    # For generate_from_params, we add a zero offset (since marginalized version doesn't have offset parameter)
    rvs = sim.rv_model .+ 0.0
    
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    if add_noise
        jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : 0.0
        radvel_table.rv .+= randn.() .* hypot.(σ_rvs, jitter)
    end

    return MarginalizedStarAbsoluteRVLikelihood(
        radvel_table;
        name=likelihoodname(like),
        trend_function=like.trend_function,
        variables=(like.priors, like.derived)
    )
end

