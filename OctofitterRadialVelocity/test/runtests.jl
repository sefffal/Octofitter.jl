using Test
using Octofitter
using OctofitterRadialVelocity
using OctofitterRadialVelocity: Celerite
using PlanetOrbits
using Distributions
using LinearAlgebra
using AbstractGPs
using TypedTables
using Random

# What this package still owns after the v2 migration: the two Gaussian-process
# backends behind core's `RadialVelocityObs`, the marginalized estimator, and
# the archive loaders. `StarAbsoluteRVObs` / `PlanetRelativeRVObs` are gone —
# both are `RadialVelocityObs` with different `target`/`ref` — so the tests
# that exercised them live in `Octofitter/test/v2/radial-velocity.jl`.

const EPOCHS = collect(range(56100.0, 58900.0, length=12))

# A two-companion system with every parameter fixed, so a likelihood can be
# checked against a closed form rather than against itself.
function fixed_system(; obs)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0 - 6mjup
        flux = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 2mjup
        a = 3.0; e = 0.15; i = 0.9; ω = 1.2; Ω = 2.0; tp = 55000.0
    end)
    c = Octofitter.Body(name="c", about=A, variables=@variables begin
        mass = 4mjup
        a = 9.0; e = 0.05; i = 1.1; ω = 0.3; Ω = 1.0; tp = 56000.0
    end)
    return Octofitter.System(name="ref", bodies=[A, b, c], observations=obs,
        variables=@variables begin plx = 25.0 end)
end

ll_of(obs) = (sys = fixed_system(obs=[obs]);
              Octofitter.make_ln_like(sys)(sys, Octofitter.make_arr2nt(sys)(Float64[])))

# The exact RV the model will produce, from the same system definition.
function truth(target, reference=Barycentre)
    posys = PlanetOrbits.System(
        (PlanetOrbits.Body(mass=1.0 - 6mjup, flux=(; default=1.0), name=:A),
         PlanetOrbits.Body(mass=2mjup, name=:b),
         PlanetOrbits.Body(mass=4mjup, name=:c)),
        (PlanetOrbits.Orbit(PlanetOrbits.Body(mass=2mjup, name=:b),
             about=PlanetOrbits.Body(mass=1.0 - 6mjup, name=:A);
             a=3.0, e=0.15, i=0.9, ω=1.2, Ω=2.0, tp=55000.0),
         PlanetOrbits.Orbit(PlanetOrbits.Body(mass=4mjup, name=:c),
             about=PlanetOrbits.Body(mass=1.0 - 6mjup, name=:A);
             a=9.0, e=0.05, i=1.1, ω=0.3, Ω=1.0, tp=56000.0)); plx=25.0)
    tr = orbitsolve(posys, EPOCHS)
    r = reference === Barycentre ? PlanetOrbits.barycentre(posys) : reference
    return [radvel(tr[k], target, r) for k in eachindex(EPOCHS)]
end

# ---------------------------------------------------------------------------
# Gaussian-process backends
#
# Core Octofitter declares three hooks (`gp_condition`, `gp_ln_like`,
# `gp_predict`) and depends on no GP package. AbstractGPs works through the
# duck-typed defaults; Celerite gets methods in `src/gp-backends.jl`.
# ---------------------------------------------------------------------------

# Residual-generating data: a reflex signal plus a smooth "activity" wiggle.
function gp_data()
    reflex = truth(:A)
    resid = 7.0 .* sinpi.(EPOCHS ./ 900) .+ 2.0 .* cospi.(EPOCHS ./ 310)
    return reflex, resid, fill(2.0, length(EPOCHS))
end

@testset "AbstractGPs backend" begin
    reflex, resid, σ = gp_data()
    η₁, η₂ = 6.0, 400.0
    kern = η₁^2 * SqExponentialKernel() ∘ ScaleTransform(1 / η₂)

    obs = RadialVelocityObs((epoch=EPOCHS, rv=reflex .+ resid, σ_rv=σ);
        target=:A, ref=Barycentre, name="rv",
        gaussian_process=θ_obs -> GP(θ_obs.gp_η₁^2 * SqExponentialKernel() ∘
                                     ScaleTransform(1 / θ_obs.gp_η₂)),
        variables=@variables begin
            jitter = 1.0
            gp_η₁ = $η₁
            gp_η₂ = $η₂
        end)

    # No dispatch needed: `gp(x, σ²)` then `logpdf` is core's default path.
    Σ = kernelmatrix(kern, EPOCHS) + Diagonal(σ .^ 2 .+ 1.0^2)
    @test ll_of(obs) ≈ logpdf(MvNormal(Symmetric(Σ)), resid) rtol = 1e-8

    # Prediction *is* implemented for AbstractGPs — `gp_predict` on a
    # `FiniteGP`, which closes the hole v1 left. Pin it against the textbook
    # call rather than against itself, since both callers (cross-validation,
    # and the correlated-noise curves the plots draw) trust these numbers.
    σ² = σ .^ 2 .+ 1.0^2
    fx = Octofitter.gp_condition(GP(kern), EPOCHS, σ²)
    xs = collect(range(first(EPOCHS), last(EPOCHS), length=7))
    m, v = Octofitter.gp_predict(fx, resid, xs)
    post = AbstractGPs.posterior(GP(kern)(EPOCHS, σ²), resid)
    @test m ≈ AbstractGPs.mean(post, xs) rtol = 1e-12
    @test v ≈ AbstractGPs.var(post, xs) rtol = 1e-12

    # And so cross-validation works, exactly as it does for Celerite:
    # p(all) = p(train) · p(test | train) for one held-out point. The identity
    # is internal to the likelihood — every term uses the same RV model — so
    # it pins the GP conditioning without depending on the model's absolute
    # value.
    mk(rows) = RadialVelocityObs((epoch=EPOCHS[rows], rv=(reflex.+resid)[rows], σ_rv=σ[rows]);
        target=:A, ref=Barycentre, name="rv",
        gaussian_process=θ_obs -> GP(θ_obs.gp_η₁^2 * SqExponentialKernel() ∘
                                     ScaleTransform(1 / θ_obs.gp_η₂)),
        variables=@variables begin
            jitter = 1.0
            gp_η₁ = $η₁
            gp_η₂ = $η₂
        end)
    k = 3
    sub = Octofitter.likeobj_from_epoch_subset(mk(eachindex(EPOCHS)), k)
    @test sub.held_out_table.epoch == [EPOCHS[k]]
    @test ll_of(mk(setdiff(eachindex(EPOCHS), k))) + ll_of(sub) ≈
          ll_of(mk(eachindex(EPOCHS))) rtol = 1e-8
end

@testset "Celerite backend" begin
    reflex, resid, σ = gp_data()
    log_a, log_c = log(36.0), log(1 / 400)
    celerite_gp() = Celerite.CeleriteGP(Celerite.RealTerm(log_a, log_c))

    obs = RadialVelocityObs((epoch=EPOCHS, rv=reflex .+ resid, σ_rv=σ);
        target=:A, ref=Barycentre, name="rv",
        gaussian_process=θ_obs -> celerite_gp(),
        variables=@variables begin jitter = 1.0 end)

    # A `RealTerm` is k(τ) = a·exp(−c|τ|); build that covariance densely and
    # compare, which pins the O(N) Cholesky against a textbook evaluation.
    a, c = exp(log_a), exp(log_c)
    Σ = [a * exp(-c * abs(s - t)) for s in EPOCHS, t in EPOCHS] +
        Diagonal(σ .^ 2 .+ 1.0^2)
    @test ll_of(obs) ≈ logpdf(MvNormal(Symmetric(Σ)), resid) rtol = 1e-8

    # `compute!` takes standard deviations, not variances. Getting that wrong
    # is silent and wrecks the fit, so check the jitter actually lands in the
    # diagonal at the right power.
    obs2 = RadialVelocityObs((epoch=EPOCHS, rv=reflex .+ resid, σ_rv=σ);
        target=:A, ref=Barycentre, name="rv",
        gaussian_process=θ_obs -> celerite_gp(),
        variables=@variables begin jitter = 4.0 end)
    Σ2 = [a * exp(-c * abs(s - t)) for s in EPOCHS, t in EPOCHS] +
         Diagonal(σ .^ 2 .+ 4.0^2)
    @test ll_of(obs2) ≈ logpdf(MvNormal(Symmetric(Σ2)), resid) rtol = 1e-8
end

@testset "Celerite cross-validation: the held-out point is conditioned on the rest" begin
    reflex, resid, σ = gp_data()
    log_a, log_c = log(36.0), log(1 / 400)
    mk(rows) = RadialVelocityObs((epoch=EPOCHS[rows], rv=(reflex.+resid)[rows], σ_rv=σ[rows]);
        target=:A, ref=Barycentre, name="rv",
        gaussian_process=θ_obs -> Celerite.CeleriteGP(Celerite.RealTerm(log_a, log_c)),
        variables=@variables begin jitter = 1.0 end)

    full = mk(eachindex(EPOCHS))
    k = 5
    sub = Octofitter.likeobj_from_epoch_subset(full, k)
    @test length(sub.table.epoch) == length(EPOCHS) - 1
    @test sub.held_out_table.epoch == [EPOCHS[k]]

    # p(all) = p(train) · p(test | train), exactly, for one held-out point.
    # v1 could not run this at all: it called `Main.Celerite.predict`, which
    # resolves only if the user happens to have a `Celerite` binding in `Main`.
    train = mk(setdiff(eachindex(EPOCHS), k))
    @test ll_of(train) + ll_of(sub) ≈ ll_of(full) rtol = 1e-8
end

# ---------------------------------------------------------------------------
# MarginalizedRVObs
# ---------------------------------------------------------------------------

@testset "MarginalizedRVObs construction and references" begin
    tab = Table(epoch=[50000.0, 50100.0], rv=[100.0, 110.0], σ_rv=[10.0, 10.0])
    obs = MarginalizedRVObs(tab; target=:A, ref=Barycentre, name="MargRV",
        variables=@variables begin jitter ~ LogUniform(0.001, 100) end)
    @test length(obs.table) == 2
    @test Octofitter.refspecs(obs) === (Octofitter.refspec(:A), Barycentre)
    # "Star" and "Absolute" were only ever ref choices; relative RV
    # marginalizes an instrument zero point just as happily.
    rel = MarginalizedRVObs(tab; target=:b, ref=:A, name="MargRel",
        variables=@variables begin jitter = 1.0 end)
    @test Octofitter.refspecs(rel) === (Octofitter.refspec(:b), Octofitter.refspec(:A))

    # An explicit offset would be integrated over *and* added.
    @test_throws r"marginalizes the zero point" MarginalizedRVObs(tab;
        target=:A, ref=Barycentre, name="bad",
        variables=@variables begin
            offset ~ Normal(0, 100)
            jitter = 1.0
        end)

    # No pointwise likelihoods ⇒ no cross-validation. Structural, not a gap.
    @test_throws r"cannot be" Octofitter.likeobj_from_epoch_subset(obs, 1:1)
end

@testset "MarginalizedRVObs integrates the zero point out" begin
    reflex = truth(:A)
    σ = fill(3.0, length(EPOCHS))
    jitter = 2.0
    Random.seed!(1234)
    wiggle = 4.0 .* sinpi.(EPOCHS ./ 700)
    data = reflex .+ 61.0 .+ wiggle           # a big unknown instrument offset

    obs = MarginalizedRVObs((epoch=EPOCHS, rv=data, σ_rv=σ);
        target=:A, ref=Barycentre, name="marg",
        variables=@variables begin jitter = $jitter end)
    ll = ll_of(obs)

    # Residuals against the model *without* any offset, which is what the
    # marginalization integrates over.
    r = data .- reflex
    v = σ .^ 2 .+ jitter^2
    A = sum(1 ./ v)
    B = -2sum(r ./ v)
    C = sum(r .^ 2 ./ v)

    # The normalized marginal log-likelihood, ∫ p(data | offset) d(offset)
    # under a flat prior on the offset:
    marginal = -sum(log.(2π .* v)) / 2 + log(2π / A) / 2 - C / 2 + B^2 / (8A)

    # …checked independently by quadrature over the offset.
    peak = -B / (2A)
    width = 1 / sqrt(A)
    grid = range(peak - 12width, peak + 12width, length=20001)
    integrand = [exp(sum(logpdf.(Normal.(0, sqrt.(v)), r .- μ)) - marginal) for μ in grid]
    @test sum(integrand) * step(grid) ≈ 1.0 rtol = 1e-6

    # The likelihood *is* the normalized marginal — this is the v2 fix.
    #
    # v1 computed `−Σ log(2π vᵢ) + B²/(4A) − C − log A`, which equals
    # `2·marginal − log(2π)`. The additive constant was harmless; the factor of
    # two was not, because it is not a constant: it sharpened this instrument's
    # contribution relative to every other term in a joint model and relative to
    # the priors. Corrected on Will's call (2026-08-05).
    @test ll ≈ marginal rtol = 1e-10

    # And explicitly *not* the v1 expression, so a revert cannot pass silently.
    @test !isapprox(ll, 2marginal - log(2π), rtol=1e-6)

    # Whatever the normalization, the estimator must be invariant to the very
    # thing it marginalizes: shifting every datum by a constant cannot change
    # the likelihood.
    shifted = MarginalizedRVObs((epoch=EPOCHS, rv=data .+ 250.0, σ_rv=σ);
        target=:A, ref=Barycentre, name="marg",
        variables=@variables begin jitter = $jitter end)
    @test ll_of(shifted) ≈ ll rtol = 1e-10

    # And the offset the data prefer is the inverse-variance weighted mean
    # residual — what plotting puts the data on.
    sys = fixed_system(obs=[obs])
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    ctx = Octofitter.ObsContext(nt, nt.observations.marg, posys,
        orbitsolve(posys, ep), maps[obs])
    @test OctofitterRadialVelocity.marginalized_offset(obs, ctx) ≈ sum(r ./ v) / A rtol = 1e-10
end

@testset "MarginalizedRVObs models the same trajectory as RadialVelocityObs" begin
    # v1 summed `radvel(sol_j, mass_j)` over every planet by hand here too.
    # The two types must now agree on the forward model, offset aside.
    reflex = truth(:A)
    σ = fill(3.0, length(EPOCHS))
    marg = MarginalizedRVObs((epoch=EPOCHS, rv=reflex, σ_rv=σ);
        target=:A, ref=Barycentre, name="marg",
        variables=@variables begin jitter = 0.0 end)
    plain = RadialVelocityObs((epoch=EPOCHS, rv=reflex, σ_rv=σ);
        target=:A, ref=Barycentre, name="marg",
        variables=@variables begin jitter = 0.0 end)

    sys = fixed_system(obs=[marg])
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    ctx = Octofitter.ObsContext(nt, nt.observations.marg, posys,
        orbitsolve(posys, ep), maps[marg])
    @test Octofitter.simulate(marg, ctx).rv_model == Octofitter.simulate(plain, ctx).rv_model

    # `generate_from_params` round-trips the type and its refs.
    gen = Octofitter.generate_from_params(marg, ctx; add_noise=false)
    @test gen isa MarginalizedRVObs
    @test gen.table.rv ≈ reflex
    @test Octofitter.refspecs(gen) === Octofitter.refspecs(marg)
end

@testset "plotting protocol" begin
    reflex = truth(:A)
    σ = fill(3.0, length(EPOCHS))
    obs = MarginalizedRVObs((epoch=EPOCHS, rv=reflex .+ 61.0, σ_rv=σ);
        target=:A, ref=Barycentre, name="marg",
        variables=@variables begin jitter = 1.0 end)
    sys = fixed_system(obs=[obs])
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    ctx = Octofitter.ObsContext(nt, nt.observations.marg, posys,
        orbitsolve(posys, ep), maps[obs])

    ch = Octofitter.plotchannels(obs)
    @test length(ch) == 1 && ch[1].name === :rv
    res = Octofitter.residuals(obs, ctx)
    # There is no `offset` variable to remove, so the data are calibrated by
    # the marginalized zero point; residuals then have zero weighted mean.
    @test maximum(abs, res.rv.resid) < 1e-8
    @test res.rv.σ_eff ≈ hypot.(σ, 1.0)
end
