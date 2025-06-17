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
    variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin;end),
    name="",
    trend_function=(θ_obs, epoch)->zero(Octofitter._system_number_type(θ_obs)),
)
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
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    @no_escape begin

        
        jitter = θ_obs.jitter

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
        resid .-= rvlike.trend_function(θ_obs, rvlike.table.epoch)

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
