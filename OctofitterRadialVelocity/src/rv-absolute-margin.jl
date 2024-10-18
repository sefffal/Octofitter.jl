# This likelihood object is the same as MarginalizedStarAbsoluteRVLikelihood,
# but instead of specifying RV offsets using variables in the model,
# the RV zero point is marginalized out.
# This object only supports one instrument, but you can add as many 
# likelihood objects as you want (one per instrument).

# Gaussian processes not supported.
struct MarginalizedStarAbsoluteRVLikelihood{TTable<:Table,jitter_symbol} <: Octofitter.AbstractLikelihood
    table::TTable
    instrument_name::String
    jitter_symbol::Symbol
    function MarginalizedStarAbsoluteRVLikelihood(observations...; jitter, instrument_name=nothing)
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
        return new{typeof(table),jitter}(table, instrument_name, jitter)
    end
end
MarginalizedStarAbsoluteRVLikelihood(observations::NamedTuple...;kwargs...) = MarginalizedStarAbsoluteRVLikelihood(observations; kwargs...)
export MarginalizedStarAbsoluteRVLikelihood


function _getjitter(::MarginalizedStarAbsoluteRVLikelihood{TTable,jitter_symbol}, θ_system) where {TTable,jitter_symbol}
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

    
    jitter = _getjitter(rvlike, θ_system)

    # Data for this instrument:
    epochs = rvlike.table.epoch
    σ_rvs = rvlike.table.σ_rv
    rvs = rvlike.table.rv

    # RV residual calculation: measured RV - model
    resid = zeros(T, length(rvs))
    resid .+= rvs
    # Start with model ^

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

    return ll
end
