# Radial velocity: one type for the reflex signal and for relative RV, plus
# the nuisance machinery v1 attached only to the absolute variant (offset,
# trend, correlated noise, held-out data).
#
# Reuses `reference_system` / `model_system` from `v2/likelihoods.jl`, so this
# file must be included after it.
#
# The Gaussian-process backends themselves are exercised in
# `OctofitterRadialVelocity/test/runtests.jl` — AbstractGPs and Celerite both
# live there. What is tested here is core's side of the contract: that the
# three hooks are called in the right places with the right arguments, that
# the `-Inf` guards fire, and that the held-out path conditions on the right
# rows. A toy backend defined below stands in for a real one, which is
# exactly how `OctofitterRadialVelocity` plugs Celerite in.

const RV_EPOCHS = collect(range(56100.0, 58900.0, length=11))

# Evaluate one observation inside the shared reference model.
function rv_ll(obs)
    sys = model_system(obs=[obs])
    return Octofitter.make_ln_like(sys)(sys, Octofitter.make_arr2nt(sys)(Float64[]))
end

# Modelled RVs of one observation, in table order.
function rv_simulate(obs)
    sys = model_system(obs=[obs])
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep)
    # An observation that declares no variables gets no namespace at all.
    sym = Symbol(Octofitter.likelihoodname(obs))
    θ_obs = hasproperty(nt.observations, sym) ? getproperty(nt.observations, sym) : (;)
    ctx = Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs])
    return Octofitter.simulate(obs, ctx).rv_model
end

# Exact RVs from the same system the model is built from, for either query.
function rv_truth(target, reference)
    posys = reference_system()
    tr = orbitsolve(posys, RV_EPOCHS)
    r = reference === Barycentre ? PlanetOrbits.barycentre(posys) : reference
    return [radvel(tr[k], target, r) for k in eachindex(RV_EPOCHS)]
end

@testset "one type covers the reflex signal and relative RV" begin
    # `target=A, ref=Barycentre` — the classic stellar reflex, which v1 built
    # by summing `radvel(sol_j, mass_j)` over every planet by hand.
    reflex = rv_truth(:A, Barycentre)
    # `target=b, ref=A` — relative RV, v1's separate `PlanetRelativeRVObs`,
    # which reached for the same epicyclic superposition to correct for inner
    # companions.
    rel = rv_truth(:b, :A)
    @test maximum(abs, reflex) > 1.0
    @test maximum(abs, rel) > 100 * maximum(abs, reflex)   # c's reflex is tiny beside b's orbit

    σ = fill(1.0, length(RV_EPOCHS))
    o_reflex = RadialVelocityObs((epoch=RV_EPOCHS, rv=reflex, σ_rv=σ);
        target=:A, ref=Barycentre, name="rv")
    o_rel = RadialVelocityObs((epoch=RV_EPOCHS, rv=rel, σ_rv=σ);
        target=:b, ref=:A, name="rv")
    # Zero residuals both ways: the likelihood is exactly the normalization.
    exact = -length(RV_EPOCHS) * log(2π) / 2
    @test rv_ll(o_reflex) ≈ exact rtol = 1e-12
    @test rv_ll(o_rel) ≈ exact rtol = 1e-12

    # …and a companion measured against *another companion*, which v1's
    # per-planet solutions could not express at all.
    cb = rv_truth(:c, :b)
    o_cb = RadialVelocityObs((epoch=RV_EPOCHS, rv=cb, σ_rv=σ); target=:c, ref=:b, name="rv")
    @test rv_ll(o_cb) ≈ exact rtol = 1e-12
end

@testset "offset and trend enter the model additively" begin
    reflex = rv_truth(:A, Barycentre)
    σ = fill(2.0, length(RV_EPOCHS))
    trend = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000)

    obs = RadialVelocityObs((epoch=RV_EPOCHS, rv=reflex, σ_rv=σ);
        target=:A, ref=Barycentre, name="rv", trend_function=trend,
        variables=@variables begin
            offset = 31.0
            trend_slope = 0.004
            jitter = 3.0
        end)
    model = rv_simulate(obs)
    @test model ≈ reflex .+ 31.0 .+ 0.004 .* (RV_EPOCHS .- 57000)

    # `jitter` in quadrature, and the offset/trend carried through into the
    # residual — spelled out rather than taken from the implementation.
    ll = rv_ll(obs)
    expected = sum(eachindex(RV_EPOCHS)) do i
        σ² = σ[i]^2 + 3.0^2
        logpdf(Normal(0, sqrt(σ²)), reflex[i] - model[i])
    end
    @test ll ≈ expected rtol = 1e-12

    # No trend function at all is the same as a zero one, and stays so under
    # ForwardDiff (the default returns `zero(_system_number_type(θ_obs))`,
    # which must not pin the model to Float64).
    plain = RadialVelocityObs((epoch=RV_EPOCHS, rv=reflex, σ_rv=σ);
        target=:A, ref=Barycentre, name="rv",
        variables=@variables begin
            offset = 31.0
            jitter = 3.0
        end)
    @test rv_simulate(plain) ≈ reflex .+ 31.0
end

@testset "offset is optional, and is no longer injected for you" begin
    # v1's `StarAbsoluteRVObs` added `offset ~ Uniform(-1000, 1000)` and
    # `jitter ~ LogUniform(0.001, 100)` when handed no variables block. v2
    # never invents a prior: a model with no `offset` variable simply has no
    # offset, for absolute and relative RV alike.
    reflex = rv_truth(:A, Barycentre)
    σ = fill(1.0, length(RV_EPOCHS))
    obs = RadialVelocityObs((epoch=RV_EPOCHS, rv=reflex, σ_rv=σ);
        target=:A, ref=Barycentre, name="rv")
    @test isempty(obs.priors.priors)
    @test rv_simulate(obs) ≈ reflex
end

@testset "likeobj_from_epoch_subset keeps the rows it is asked for" begin
    reflex = rv_truth(:A, Barycentre)
    rv = reflex .+ 5.0 .* sinpi.(eachindex(reflex) ./ 3)     # non-zero residuals
    σ = fill(2.0, length(RV_EPOCHS))
    obs = RadialVelocityObs((epoch=RV_EPOCHS, rv=rv, σ_rv=σ);
        target=:A, ref=Barycentre, name="rv",
        variables=@variables begin jitter = 1.0 end)

    lo = Octofitter.likeobj_from_epoch_subset(obs, 1:4)
    hi = Octofitter.likeobj_from_epoch_subset(obs, 5:length(RV_EPOCHS))
    @test length(lo.table.epoch) == 4
    @test isempty(lo.held_out_table)                 # no GP ⇒ rows are independent
    # Uncorrelated pointwise likelihoods add up, which is the property
    # cross-validation and PSIS-LOO are built on.
    @test rv_ll(lo) + rv_ll(hi) ≈ rv_ll(obs) rtol = 1e-12

    # `:` means everything.
    @test length(Octofitter.likeobj_from_epoch_subset(obs, :).table.epoch) == length(RV_EPOCHS)
    # …and the trend/GP/variables survive the round-trip.
    @test lo.trend_function === obs.trend_function
    @test lo.gaussian_process === obs.gaussian_process
    @test lo.target === obs.target && lo.ref === obs.ref
    @test lo.priors === obs.priors
end

# ---------------------------------------------------------------------------
# A toy GP backend.
#
# Core Octofitter depends on neither AbstractGPs nor Celerite; a backend is
# three methods. This one is an exact dense squared-exponential GP, which
# makes the held-out identity below checkable in closed form.
# ---------------------------------------------------------------------------

struct ToyGP
    amp::Float64
    len::Float64
end
_toyk(g::ToyGP, s, t) = g.amp^2 * exp(-(s - t)^2 / (2 * g.len^2))

struct ToyFit
    gp::ToyGP
    x::Vector{Float64}
    Σ::Symmetric{Float64,Matrix{Float64}}    # kernel + per-point white noise
end

function Octofitter.gp_condition(g::ToyGP, x, σ²)
    xs = collect(Float64, x)
    K = [_toyk(g, xs[i], xs[j]) + (i == j ? Float64(σ²[i]) : 0.0)
         for i in eachindex(xs), j in eachindex(xs)]
    return ToyFit(g, xs, Symmetric(K))
end
Octofitter.gp_ln_like(f::ToyFit, r) = logpdf(MvNormal(f.Σ), collect(Float64, r))
function Octofitter.gp_predict(f::ToyFit, r, xs)
    Ks = [_toyk(f.gp, Float64(x), f.x[j]) for x in xs, j in eachindex(f.x)]
    α = f.Σ \ collect(Float64, r)
    μ = Ks * α
    # Latent predictive variance: the held-out point's own white noise is
    # added by the likelihood, not here — the same split Celerite uses.
    v = [_toyk(f.gp, Float64(xs[i]), Float64(xs[i])) - dot(Ks[i, :], f.Σ \ Ks[i, :])
         for i in eachindex(xs)]
    return μ, v
end

# Conditions fine, then fails in the factorization — the second guard.
struct BadGP end
struct BadFit end
Octofitter.gp_condition(::BadGP, x, σ²) = BadFit()
Octofitter.gp_ln_like(::BadFit, r) = throw(PosDefException(1))

function gp_obs(; gp=ToyGP(6.0, 400.0), rows=eachindex(RV_EPOCHS), name="rv")
    reflex = rv_truth(:A, Barycentre)
    rv = reflex .+ 5.0 .* sinpi.(eachindex(reflex) ./ 3)
    σ = fill(2.0, length(RV_EPOCHS))
    return RadialVelocityObs((epoch=RV_EPOCHS[rows], rv=rv[rows], σ_rv=σ[rows]);
        target=:A, ref=Barycentre, name=name, gaussian_process=θ_obs -> gp,
        variables=@variables begin jitter = 1.0 end)
end

@testset "a Gaussian process replaces the diagonal noise model" begin
    obs = gp_obs()
    reflex = rv_truth(:A, Barycentre)
    resid = 5.0 .* sinpi.(eachindex(reflex) ./ 3)
    Σ = Octofitter.gp_condition(ToyGP(6.0, 400.0), RV_EPOCHS, fill(2.0^2 + 1.0^2, length(RV_EPOCHS))).Σ
    @test rv_ll(obs) ≈ logpdf(MvNormal(Σ), resid) rtol = 1e-10

    # A GP is not a decoration: it has to change the answer.
    nogp = RadialVelocityObs((epoch=RV_EPOCHS, rv=reflex .+ resid,
            σ_rv=fill(2.0, length(RV_EPOCHS))); target=:A, ref=Barycentre, name="rv",
        variables=@variables begin jitter = 1.0 end)
    @test !isapprox(rv_ll(obs), rv_ll(nogp); rtol=1e-3)
end

@testset "GP construction failures become -Inf, not crashed chains" begin
    reflex = rv_truth(:A, Barycentre)
    tab = (epoch=RV_EPOCHS, rv=reflex, σ_rv=fill(2.0, length(RV_EPOCHS)))
    mk(f) = RadialVelocityObs(tab; target=:A, ref=Barycentre, name="rv",
        gaussian_process=f, variables=@variables begin jitter = 1.0 end)

    # A sampler that walks a length scale to zero, or into a region where the
    # covariance is not positive definite, must see a rejected sample rather
    # than a stack trace. Both stages are guarded: building the GP…
    @test rv_ll(mk(θ_obs -> throw(DomainError(-1.0, "bad hyperparameter")))) == -Inf
    @test rv_ll(mk(θ_obs -> throw(ArgumentError("bad kernel")))) == -Inf
    @test rv_ll(mk(θ_obs -> throw(PosDefException(1)))) == -Inf
    # …and evaluating it, where the factorization actually happens.
    @test rv_ll(mk(θ_obs -> BadGP())) == -Inf

    # Anything else is a bug in the user's kernel function and must surface.
    @test_throws ErrorException rv_ll(mk(θ_obs -> error("typo in kernel")))
end

@testset "the held-out path conditions on the retained rows" begin
    # With a GP the rows are not independent, so `likeobj_from_epoch_subset`
    # cannot just take a slice: it keeps the rest as the conditioning set.
    obs = gp_obs()
    heldout = 4
    sub = Octofitter.likeobj_from_epoch_subset(obs, heldout)
    @test length(sub.table.epoch) == length(RV_EPOCHS) - 1
    @test length(sub.held_out_table.epoch) == 1
    @test sub.held_out_table.epoch[1] == RV_EPOCHS[heldout]
    # Held-out rows need solved states too, and they are not in `table`;
    # `epochs(obs)` appends them so `solutionat` reaches them. v1 re-ran
    # `orbitsolve` per held-out point instead.
    @test Octofitter.epochs(sub) ==
          vcat(RV_EPOCHS[setdiff(eachindex(RV_EPOCHS), heldout)], RV_EPOCHS[heldout])

    # The chain rule for a joint Gaussian: p(all) = p(train) · p(test | train).
    # This is exact for a *single* held-out point, which is the leave-one-out
    # case cross-validation actually asks for. (For several held-out points at
    # once the likelihood treats them as conditionally independent given the
    # training set, which is v1's behaviour and an approximation.)
    train = gp_obs(rows=setdiff(eachindex(RV_EPOCHS), heldout))
    @test rv_ll(train) + rv_ll(sub) ≈ rv_ll(obs) rtol = 1e-8

    # The AbstractGPs hole is inherited from v1, and is a loud error rather
    # than a wrong number.
    @test_throws r"not implemented" Octofitter.gp_predict(nothing, [1.0], [1.0])
end

@testset "a GP over duplicated epochs is rejected at construction" begin
    # The constructor sorts, but sorted is not strictly increasing, and
    # Celerite's O(N) Cholesky needs the latter. Two rows at one epoch also
    # make any stationary kernel's covariance singular — better an error here
    # than a `PosDefException` the sampler guards quietly turn into `-Inf`.
    dup = [56000.0, 56500.0, 56500.0, 57000.0]
    tab = (epoch=dup, rv=[1.0, 2.0, 3.0, 4.0], σ_rv=fill(1.0, 4))
    @test_throws r"strictly increasing" RadialVelocityObs(tab; target=:A, ref=Barycentre,
        name="rv", gaussian_process=θ_obs -> ToyGP(1.0, 100.0))
    # Without a GP duplicated epochs are perfectly legal — two instruments'
    # worth of simultaneous data, or a night with two exposures.
    @test RadialVelocityObs(tab; target=:A, ref=Barycentre, name="rv") isa RadialVelocityObs
end

@testset "simulate/generate_from_params round-trip" begin
    reflex = rv_truth(:A, Barycentre)
    σ = fill(2.0, length(RV_EPOCHS))
    trend = (θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000)
    obs = RadialVelocityObs((epoch=RV_EPOCHS, rv=fill(0.0, length(RV_EPOCHS)), σ_rv=σ);
        target=:A, ref=Barycentre, name="rv", trend_function=trend,
        variables=@variables begin
            offset = -12.0
            trend_slope = 0.002
            jitter = 0.0
        end)
    sys = model_system(obs=[obs])
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    ctx = Octofitter.ObsContext(nt, nt.observations.rv, posys, orbitsolve(posys, ep), maps[obs])

    gen = Octofitter.generate_from_params(obs, ctx; add_noise=false)
    @test gen isa RadialVelocityObs
    @test gen.trend_function === trend
    @test gen.table.rv ≈ reflex .- 12.0 .+ 0.002 .* (RV_EPOCHS .- 57000)
    # Simulated data must be an exact fit to the parameters that made it.
    @test rv_ll(gen) ≈ sum(logpdf.(Normal.(0, σ), 0.0)) rtol = 1e-12
end
