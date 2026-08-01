# ---------------------------------------------------
# Radial velocity
#
# One likelihood covers what v1 split between the stellar-reflex path (a
# hand-summed superposition over companions) and the relative-RV path: both
# are `radvel(sol, target, ref)`, differing only in which references they
# name.
#   star vs the system barycentre → the classic reflex signal
#   companion vs the host star    → relative RV
# ---------------------------------------------------

"""
    RadialVelocityObs(data; target, ref=Barycentre, name, variables=…)

Radial velocities [m/s] of `target` measured against `ref`, positive
receding.

    rvs = RadialVelocityObs(tab; target=A, ref=Barycentre, name="HARPS",
        variables=@variables begin
            offset ~ Normal(0, 100)      # m/s, instrument zero point
            jitter ~ LogUniform(0.1, 100) # m/s, added in quadrature
        end)

`data` needs `:epoch` [MJD], `:rv` [m/s] and `:σ_rv` [m/s].

One `RadialVelocityObs` per instrument: the offset and jitter are per
instrument, and that is exactly the granularity this object carries.
"""
struct RadialVelocityObs{TTable<:Table,TT,TR} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    target::TT
    ref::TR
    name::String
end

const rv_cols = (:epoch, :rv, :σ_rv)

function RadialVelocityObs(observations;
                           target, ref=Barycentre,
                           name,
                           variables::Tuple{Priors,Derived}=(Priors(), Derived()))
    (priors, derived) = variables
    table = Table(observations)
    equal_length_cols(table) ||
        error("The columns in the input data do not all have the same length")
    issubset(rv_cols, Tables.columnnames(table)) ||
        error("Expected columns $rv_cols")
    table = table[sortperm(vec(table.epoch))]
    t, r = refspec(target), refspec(ref)
    return RadialVelocityObs{typeof(table),typeof(t),typeof(r)}(
        table, priors, derived, t, r, String(name))
end

export RadialVelocityObs

refspecs(obs::RadialVelocityObs) = (obs.target, obs.ref)

likeobj_from_epoch_subset(obs::RadialVelocityObs, inds) = RadialVelocityObs(
    obs.table[inds, :, 1]; target=obs.target, ref=obs.ref, obs.name,
    variables=(obs.priors, obs.derived))

function simulate!(rv_model, obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    offset = hasproperty(ctx.θ_obs, :offset) ? ctx.θ_obs.offset : zero(T)
    target = ref(ctx, obs.target)
    reference = ref(ctx, obs.ref)
    @inbounds for i in eachindex(obs.table.epoch)
        rv_model[i] = radvel(solutionat(ctx, i), target, reference) + offset
    end
    return (; rv_model, epochs=obs.table.epoch)
end
function simulate(obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    return simulate!(Vector{T}(undef, length(obs.table.epoch)), obs, ctx)
end

function ln_like(obs::RadialVelocityObs, ctx::ObsContext)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    ll = zero(T)
    L = length(obs.table.epoch)
    @no_escape ctx.buf begin
        rv_model = @alloc(T, L)
        simulate!(rv_model, obs, ctx)
        for i in eachindex(obs.table.epoch)
            σ² = obs.table.σ_rv[i]^2 + jitter^2
            resid = obs.table.rv[i] - rv_model[i]
            ll -= (resid^2 / σ² + log(2π * σ²)) / 2
        end
    end
    return ll
end

function generate_from_params(obs::RadialVelocityObs, ctx::ObsContext; add_noise)
    sim = simulate(obs, ctx)
    T = _system_number_type(ctx.θ_system)
    jitter = hasproperty(ctx.θ_obs, :jitter) ? ctx.θ_obs.jitter : zero(T)
    rv = add_noise ? sim.rv_model .+ randn.() .* hypot.(obs.table.σ_rv, jitter) : sim.rv_model
    return RadialVelocityObs(Table(; epoch=obs.table.epoch, rv, σ_rv=obs.table.σ_rv);
        target=obs.target, ref=obs.ref, obs.name, variables=(obs.priors, obs.derived))
end
