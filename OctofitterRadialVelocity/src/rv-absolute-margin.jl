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
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11),

        jitter=:jitter_1,
        instrument_name="inst name",
    )

Represents a likelihood function of relative astometry between a host star and a secondary body.
`:epoch` (mjd), `:rv` (m/s), and `:σ_rv` (m/s)  are all required.
In addition to the example above, any Tables.jl compatible source can be provided.

The`jitter` parameters specify which variables should be read from the model for the 
    jitter of this instrument.
    
Unlike `StarAbsoluteRVLikelihood`, this version analytically marginalizes over the instrument
RV zero point. This makes it faster in most cases.
That said, if you have a specific prior you want to apply for the RV zero point, correlations between
zero points, hierarchical models, etc, you should use `StarAbsoluteRVLikelihood`.

Additionally, a gaussian and trend fit are not supported with the analytically marginalization.
"""
struct MarginalizedStarAbsoluteRVLikelihood{TTable<:Table,TF,jitter_symbol} <: Octofitter.AbstractLikelihood
    table::TTable
    instrument_name::String
    trend_function::TF
    jitter_symbol::Symbol
    function MarginalizedStarAbsoluteRVLikelihood(
        observations...;
        jitter,
        instrument_name=nothing,
        trend_function=(θ_system, epoch)->zero(Octofitter._system_number_type(θ_system)),
    )
        table = Table(observations...)
        if !issubset(rv_cols, Tables.columnnames(table))
            error("Expected columns $rv_cols")
        end
        rows = map(eachrow(table)) do row′
            row = (;inst_idx=1, row′[1]..., rv=float(row′[1].rv[1]))
            return row
        end
        ii = sortperm([r.epoch for r in rows])
        table = Table(rows[ii])
        # Check instrument indexes are contiguous
        if isnothing(instrument_name)
            instrument_name = ""
        end
        return new{typeof(table),typeof(trend_function),jitter}(table, instrument_name, trend_function, jitter)
    end
end
MarginalizedStarAbsoluteRVLikelihood(observations::NamedTuple...;kwargs...) = MarginalizedStarAbsoluteRVLikelihood(observations; kwargs...)
function Octofitter.likeobj_from_epoch_subset(obs::MarginalizedStarAbsoluteRVLikelihood, obs_inds)
    return MarginalizedStarAbsoluteRVLikelihood(
        obs.table[obs_inds,:,1]...;
        jitter=obs.jitter_symbol,
        trend_function=obs.trend_function,
        instrument_name=obs.instrument_name,
    )
end
export MarginalizedStarAbsoluteRVLikelihood


function _getjitter(::MarginalizedStarAbsoluteRVLikelihood{TTable,TF,jitter_symbol}, θ_system) where {TTable,TF,jitter_symbol}
    jitter = getproperty(θ_system, jitter_symbol)
    return jitter
end

"""
Absolute radial velocity likelihood (for a star).
"""
function Octofitter.ln_like(
    rvlike::MarginalizedStarAbsoluteRVLikelihood,
    θ_system,
    planet_orbits::Tuple,
    orbit_solutions,
    orbit_solutions_i_epoch_start
)
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    @no_escape begin

        
        jitter = _getjitter(rvlike, θ_system)

        # Data for this instrument:
        epochs = rvlike.table.epoch
        σ_rvs = rvlike.table.σ_rv
        rvs = rvlike.table.rv

        # RV residual calculation: measured RV - model
        resid = @alloc(T, length(rvs))
        fill!(resid, 0)
        resid .+= rvs
        # Start with model ^

        # Subtract trend function
        resid .-= rvlike.trend_function(θ_system, rvlike.table.epoch)

        # Go through all planets and subtract their modelled influence on the RV signal:
        for planet_i in eachindex(planet_orbits)
            planet_mass = θ_system.planets[planet_i].mass
            for i_epoch in eachindex(epochs)
                # We've already solved each orbit at each epoch; grab the solution
                sol_i = i_epoch+orbit_solutions_i_epoch_start
                if sol_i > length(orbit_solutions)
                    # @show length(orbit_solutions)
                    # @show size(orbit_solutions[planet_i])
                    # @show sol_i i_epoch orbit_solutions_i_epoch_start
                end
                sol = orbit_solutions[planet_i][sol_i]
                resid[i_epoch] -= radvel(sol, planet_mass*Octofitter.mjup2msol)
            end
        end
        
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
