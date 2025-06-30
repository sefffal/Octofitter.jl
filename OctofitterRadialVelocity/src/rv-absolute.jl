# Radial Velocity data type
const rv_cols = (:epoch, :rv, :σ_rv)

"""
    StarAbsoluteRVLikelihood(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11);
        
        name="inst name",
        variables=@variables begin
            offset ~ Normal(0, 100)           # RV zero-point (m/s)
            jitter ~ LogUniform(0.1, 100.0)  # RV jitter (m/s)
        end
    )

    # Example with trend function and Gaussian Process:
    StarAbsoluteRVLikelihood(
        (;epoch=5000.0,  rv=−6.54, σ_rv=1.30),
        (;epoch=5050.1,  rv=−3.33, σ_rv=1.09),
        (;epoch=5100.2,  rv=7.90,  σ_rv=.11);
        
        name="inst name",
        trend_function = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000),  # Linear trend
        gaussian_process = θ_obs -> GP(θ_obs.gp_η₁^2 * SqExponentialKernel() ∘ ScaleTransform(1/θ_obs.gp_η₂)),
        variables=@variables begin
            offset ~ Normal(0, 100)             # RV zero-point (m/s)
            jitter ~ LogUniform(0.1, 100.0)    # RV jitter (m/s)
            trend_slope ~ Normal(0, 1)          # Linear trend slope (m/s/day)
            gp_η₁ ~ LogUniform(1.0, 100.0)      # GP amplitude
            gp_η₂ ~ LogUniform(1.0, 100.0)      # GP length scale
        end
    )

Represents a likelihood function of absolute radial velocity of a host star.
`:epoch` (mjd), `:rv` (m/s), and `:σ_rv` (m/s) are all required.

In addition to the example above, any Tables.jl compatible source can be provided.

The `offset` and `jitter` variables should be defined in the variables block and represent the 
RV zero-point and additional uncertainty to be added in quadrature to the formal measurement errors.

When using a trend function, it should be a function that takes `θ_obs` (observation parameters) 
and `epoch` and returns an RV offset. Trend parameters should be defined in the variables block.

When using a Gaussian process, the `gaussian_process` parameter should be a function that takes
`θ_obs` (observation parameters) and returns a GP kernel. GP hyperparameters should be defined
in the variables block and accessed via `θ_obs.parameter_name`.
"""
struct StarAbsoluteRVLikelihood{TTable<:Table,GP,TF} <: Octofitter.AbstractLikelihood
    table::TTable
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    held_out_table::TTable
    name::String
    gaussian_process::GP
    trend_function::TF
end
function StarAbsoluteRVLikelihood(
    observations;
    variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin;end),
    trend_function=(θ_obs, epoch)->zero(Octofitter._system_number_type(θ_obs)),
    name::String,
    gaussian_process=nothing
)
    (priors,derived)=variables
    table = Table(observations)[:,:,1]
    if !Octofitter.equal_length_cols(table)
        error("The columns in the input data do not all have the same length")
    end
    if !issubset(rv_cols, Tables.columnnames(table))
        error("Expected columns $rv_cols")
    end
    if hasproperty(table, :inst_idx) && length(unique(table.inst_idx)) > 1
        error("Deprecated: data from separate RV instruments should now be placed into different StarAbsoluteRVLikelihood likelihood objects, rather than specified by an inst_idx parameter.")
    end
    rows = map(eachrow(table)) do row′
        row = (;row′[1]..., rv=float(row′[1].rv[1]))
        return row
    end
    # We sort the data first by instrument index then by epoch to make some later code faster
    ii = sortperm([r.epoch for r in rows])
    table = Table(rows[ii])

    if any(>=(mjd("2050")),  table.epoch) || any(<=(mjd("1950")),  table.epoch)
        @warn "The data you entered fell outside the range year 1950 to year 2050. The expected input format is MJD (modified julian date). We suggest you double check your input data!"
    end

    # We need special book keeping for computing cross-validataion scores
    # We keep a table of "held out" data if needed for that purpose.
    # Here we leave it empty.
    held_out_table = empty(table)

    return StarAbsoluteRVLikelihood{typeof(table),typeof(gaussian_process),typeof(trend_function)}(
        table, priors, derived, held_out_table, name, gaussian_process, trend_function
    )
end
# StarAbsoluteRVLikelihood(observations::NamedTuple...;kwargs...) = StarAbsoluteRVLikelihood(observations; kwargs...)
function Octofitter.likeobj_from_epoch_subset(obs::StarAbsoluteRVLikelihood, obs_inds)
    # Due to TypedTables bug, the line below creates a "matrix" table that isn't the same type as the input.
    # table = typeof(obs.table)(obs.table[setdiff(1:size(obs.table,1), obs_inds),:,1])
    # table = Table(collect(eachrow(obs.table))[setdiff(1:size(obs.table,1), obs_inds)]...)
    table = Table(first(eachcol(obs.table[setdiff(1:size(obs.table,1), obs_inds)])))
    if obs_inds isa Number
        held_out_table = obs.table[obs_inds,:,1]
    else
        held_out_table = Table(first(eachcol(obs.table[obs_inds])))
    end
    return StarAbsoluteRVLikelihood{
        typeof(table),typeof(obs.gaussian_process),typeof(obs.trend_function)
    }(
        table, obs.priors, obs.derived, held_out_table, likelihoodname(obs), obs.gaussian_process, obs.trend_function
    )
end
export StarAbsoluteRVLikelihood



"""
Absolute radial velocity likelihood (for a star).
"""
function Octofitter.ln_like(
    rvlike::StarAbsoluteRVLikelihood,
    θ_system,
    θ_obs,
    planet_orbits::Tuple,
    orbit_solutions,
    orbit_solutions_i_epoch_start
)
    L = length(rvlike.table.epoch)
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    offset = hasproperty(θ_obs, :offset) ? θ_obs.offset : zero(T)
    jitter = hasproperty(θ_obs, :jitter) ? θ_obs.jitter : zero(T)
    

    @no_escape begin


        # Vector of radial velocity of the star at each epoch. Go through and sum up the influence of
        # each planet and put it into here. 
        rv_buf =  @alloc(T, L)
        rv_var_buf =   @alloc(T, L)

        # RV "data" calculation: measured RV + our barycentric rv calculation
        rv_buf .= rvlike.table.rv .- offset .- rvlike.trend_function(θ_obs, rvlike.table.epoch)

        # Go through all planets and subtract their modelled influence on the RV signal:
        # You could consider `rv_star` as the residuals after subtracting these.
        
        for planet_i in eachindex(planet_orbits)
            orbit = planet_orbits[planet_i]
            planet_mass = θ_system.planets[planet_i].mass
            for epoch_i in eachindex(rvlike.table.epoch)
                rv_buf[epoch_i] -= radvel(
                    orbit_solutions[planet_i][epoch_i+orbit_solutions_i_epoch_start],
                    planet_mass*Octofitter.mjup2msol
                )
            end
        end        

        # The noise variance per observation is the measurement noise and the jitter added
        # in quadrature
        rv_var_buf .= rvlike.table.σ_rv.^2 .+ jitter^2

        # Two code paths, depending on if we are modelling the residuals by 
        # a Gaussian process or not.
        if isnothing(rvlike.gaussian_process)
            # Don't fit a GP
            fx = MvNormal(Diagonal((rv_var_buf)))
            ll += logpdf(fx, rv_buf)
        else
            # Fit a GP
            local gp
            try
                gp = @inline rvlike.gaussian_process(θ_obs)
            catch err
                if err isa DomainError
                    ll = convert(T,-Inf)
                else
                    rethrow(err)
                end
            end

            local gp, fx
            try
                gp = @inline rvlike.gaussian_process(θ_obs)
                if gp isa Celerite.CeleriteGP
                    Celerite.compute!(gp, rvlike.table.epoch, sqrt.(rv_var_buf))# TODO: is this std or var?
                else
                    fx = gp(rvlike.table.epoch, rv_var_buf)
                end

            catch err
                if err isa DomainError
                    ll = convert(T,-Inf)
                elseif err isa PosDefException
                    ll = convert(T,-Inf)
                elseif err isa ArgumentError
                    ll = convert(T,-Inf)
                else
                    rethrow(err)
                end
            end 

            # early return not allowed with Bumper
            if isfinite(ll)
                try
                    # Normal path: evaluate likelihood against all data
                    if isempty(rvlike.held_out_table)
                        if gp isa Celerite.CeleriteGP
                            ll += Celerite.log_likelihood(gp, rv_buf)
                        else
                            ll += logpdf(fx, rv_buf)
                        end
                    # Cross validation path: condition against rvlike.table, but evaluate against
                    # rvlike.held_out_table
                    else
                        # If we have held out data, that means we are doing cross-validataion ----
                        # we are conditioning on a subset of data, and computing the likelihood for
                        # the held out data. 


                        # Vector of radial velocity of the star at each epoch. Go through and sum up the influence of
                        # each planet and put it into here. 
                        rv_buf_held_out =  @alloc(T, length(rvlike.held_out_table.epoch))
                        rv_var_buf_held_out =  @alloc(T, length(rvlike.held_out_table.epoch))

                        # RV "data" calculation: measured RV + our barycentric rv calculation
                        rv_buf_held_out .= rvlike.held_out_table.rv .- offset .- rvlike.trend_function(θ_obs, rvlike.held_out_table.epoch)

                        # Go through all planets and subtract their modelled influence on the RV signal:
                        # You could consider `rv_star` as the residuals after subtracting these.
                        
                        for planet_i in eachindex(planet_orbits)
                            orbit = planet_orbits[planet_i]
                            planet_mass = θ_system.planets[planet_i].mass
                            for epoch_i in eachindex(rvlike.held_out_table.epoch)
                                rv_buf_held_out[epoch_i] -= radvel(
                                    # We can't look into the pre-populated orbit solutions here, since these 
                                    # are only generated for entries in an likelihood objects `table`.
                                    # We have to solve it ourselves as we go. This should have negligible
                                    # performance impact unless we are holding out many data points and
                                    # we would miss the multi-threaded solve.
                                    # orbit_solutions[planet_i][epoch_i+orbit_solutions_i_epoch_start],
                                    orbitsolve(orbit, rvlike.held_out_table.epoch[epoch_i]),
                                    planet_mass*Octofitter.mjup2msol
                                )
                            end
                        end        

                        # The noise variance per observation is the measurement noise and the jitter added
                        # in quadrature
                        rv_var_buf_held_out .= rvlike.held_out_table.σ_rv.^2 .+ jitter^2
                        
                        # Compute GP model
                        if gp isa Celerite.CeleriteGP
                            pred, var = Main.Celerite.predict(gp, rv_buf, rvlike.held_out_table.epoch; return_var=true)
                            for i_held_out in 1:size(rvlike.held_out_table.epoch,1)
                                ll += logpdf(Normal(pred[i_held_out], sqrt(var[i_held_out] + rv_var_buf_held_out[i_held_out])), rv_buf_held_out[i_held_out])
                            end
                        else
                            # TODO: need to implement the prediction for the AbstractGPs case
                            throw(NotImplementedException())
                        end
                    end

                catch err
                    if err isa PosDefException || err isa DomainError
                        @warn "err" exception=(err, catch_backtrace()) θ_system
                        ll = convert(T,-Inf)
                    else
                        rethrow(err)
                    end
                end
            end
        end
    end
    return ll
end




# Generate new radial velocity observations for a star
function Octofitter.generate_from_params(like::StarAbsoluteRVLikelihood, θ_system,  θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)

    # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = like.table.epoch 
    σ_rvs = like.table.σ_rv 
    planet_masses = [θ_planet.mass for θ_planet in θ_system.planets] .* 0.000954588 # Mjup -> Msun

    # Generate new star radial velocity data
    rvs = radvel.(reshape(collect(orbits), :, 1), epochs, transpose(planet_masses))
    rvs = sum(rvs, dims=2)[:,1] .+ θ_system.rv
    radvel_table = Table(epoch=epochs, rv=rvs, σ_rv=σ_rvs)

    return StarAbsoluteRVLikelihood(
        radvel_table;
        name=likelihoodname(like),
        gaussian_process=like.gaussian_process,
        trend_function=like.trend_function,
        variables=(like.priors, like.derived)
    )
end

