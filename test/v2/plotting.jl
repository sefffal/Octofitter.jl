# The backend-agnostic plotting API: queries, the observation plotting
# protocol, and PosteriorSeries. No Makie here — the extension is exercised
# by the smoke script in examples/ and by users; what matters in CI is that
# the protocol agrees with the likelihoods.

using Octofitter
using Distributions
using LinearAlgebra
using Random: Xoshiro
using Statistics: cor, mean, std, var
using DelimitedFiles: readdlm

@testset "ObservableQuery construction" begin
    q = ObservableQuery(radvel, :b, :A)
    @test q.func === radvel
    q2 = ObservableQuery(:raoff, Barycentre(:Aa, :Ab), Barycentre)
    @test q2.func === raoff
    @test Octofitter._query((posangle, :b, :A)) isa ObservableQuery
    @test_throws ErrorException ObservableQuery(:not_an_observable, :b, :A)
    # no self-recursion on already-normalized inputs
    @test ObservableQuery(q.func, q.target, q.ref) isa ObservableQuery
end

# One small model with astrometry (sep/pa), astrometry (ra/dec), and RV.
function _plotting_test_model()
    A = Octofitter.Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.5)
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 0.005
        a ~ Uniform(2.0, 4.0)
        e = 0.25
        i = 0.6
        ω = 1.0
        Ω = 2.0
        tp = 59200.0
    end)
    seppa = RelAstromObs(
        Table(epoch=[59000.0, 59300.0, 59600.0],
            sep=[300.0, 310.0, 295.0], pa=[0.5, 0.7, 0.9],
            σ_sep=[5.0, 5.0, 5.0], σ_pa=[0.01, 0.01, 0.01]);
        target=b, ref=A, name="SEPPA",
        variables=@variables begin
            jitter = 1.0
            platescale = 1.01
            northangle = 0.02
        end)
    radec = RelAstromObs(
        Table(epoch=[59100.0, 59400.0],
            ra=[200.0, 180.0], dec=[150.0, 170.0],
            σ_ra=[4.0, 4.0], σ_dec=[4.0, 4.0]);
        target=b, ref=A, name="RADEC")
    rvs = RadialVelocityObs(
        Table(epoch=[59050.0, 59250.0, 59450.0, 59650.0],
            rv=[30.0, -10.0, 15.0, -25.0], σ_rv=[3.0, 3.0, 3.0, 3.0]);
        target=A, ref=Barycentre, name="RV",
        variables=@variables begin
            offset = 12.0
            jitter = 2.0
        end)
    sys = Octofitter.System(name="plottest", bodies=[A, b],
        observations=[seppa, radec, rvs],
        variables=@variables begin
            plx ~ truncated(Normal(50.0, 0.5), lower=1)
        end)
    return Octofitter.LogDensityModel(sys, verbosity=0), (; seppa, radec, rvs)
end

@testset "plotting protocol" begin
    model, obs = _plotting_test_model()

    # Relative astrometry declares all four channels whatever the table holds:
    # (sep, pa) and (Δα*, Δδ) are one measurement in two bases, so a mixed
    # dataset can show every point on every panel. The two the table does not
    # carry are marked `derived`, and the natives come first.
    chs = plotchannels(obs.seppa)
    @test [ch.name for ch in chs] == [:sep, :pa, :raoff, :decoff]
    @test [ch.derived for ch in chs] == [false, false, true, true]
    @test chs[1].query isa ObservableQuery
    @test chs[2].name === :pa && chs[2].wrap == 360.0 && chs[2].scale ≈ rad2deg(1.0)
    @test [ch.derived for ch in plotchannels(obs.radec)] == [true, true, false, false]
    @test plotchannels(obs.rvs)[1].name === :rv
    @test !plotchannels(obs.rvs)[1].derived

    # Sharing a panel is a per-type policy, not a plot-side guess.
    @test Octofitter.sharepanel(obs.seppa)
    @test !Octofitter.sharepanel(obs.rvs)

    # prior-shaped terms and unknown types are not plottable
    @test isempty(plotchannels(Octofitter.UserLikelihood(:x, :y, "u")))
    @test length(Octofitter.plottable_observations(model.system)) == 3

    # A "chain" of prior draws is enough to exercise the whole path.
    rng = Xoshiro(2)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:40]
    chain = Octofitter.result2mcmcchain(nts)
    series = PosteriorSeries(model, chain; N=10)

    @test length(series) == 10
    @test allunique(series.ii)
    @test issorted(series.ts) && allunique(series.ts)
    @test all(e -> e in series.ts, series.data_epochs)
    @test length(series.trajs) == 10

    # residuals: shapes, calibration, and exact agreement with ln_like
    ctx = obscontext(series, obs.rvs)
    r = Octofitter.residuals(obs.rvs, ctx)
    @test length(r.rv.resid) == 4
    # data is offset-calibrated; model is the pure orbit quantity
    @test r.rv.data ≈ obs.rvs.table.rv .- 12.0
    @test r.rv.resid ≈ r.rv.data .- r.rv.model
    @test all(r.rv.σ_eff .≈ hypot.(obs.rvs.table.σ_rv, 2.0))
    ll_expect = -0.5 * sum(@. r.rv.resid^2 / r.rv.σ_eff^2 + log(2π * r.rv.σ_eff^2))
    @test ll_expect ≈ Octofitter.ln_like(obs.rvs, ctx) rtol = 1e-12

    # sep/pa: residuals in the likelihood's own space
    ctx2 = obscontext(series, obs.seppa)
    r2 = Octofitter.residuals(obs.seppa, ctx2)
    @test keys(r2) == (:sep, :pa, :raoff, :decoff)
    @test r2.sep.data ≈ obs.seppa.table.sep .* 1.01          # platescale applied
    @test all(abs.(r2.pa.resid) .<= 180.0)                    # wrapped, degrees
    @test r2.pa.σ[1] ≈ rad2deg(0.01)

    # ra/dec: no nuisance variables → data pass through
    ctx3 = obscontext(series, obs.radec)
    r3 = Octofitter.residuals(obs.radec, ctx3)
    @test keys(r3) == (:sep, :pa, :raoff, :decoff)
    @test r3.raoff.data ≈ collect(obs.radec.table.ra) rtol = 1e-12
    @test r3.raoff.σ_eff ≈ r3.raoff.σ                        # no jitter

    # ln_like agreement for the diagonal ra/dec case too
    ll3 = -0.5 * sum(@. r3.raoff.resid^2 / r3.raoff.σ^2 + log(2π * r3.raoff.σ^2)) +
          -0.5 * sum(@. r3.decoff.resid^2 / r3.decoff.σ^2 + log(2π * r3.decoff.σ^2))
    @test ll3 ≈ Octofitter.ln_like(obs.radec, ctx3) rtol = 1e-12

    # queries evaluate over the dense grid, per draw
    curves = modelcurves(series, (radvel, :A, Barycentre))
    @test length(curves) == 10 && length(curves[1]) == length(series.ts)
    mc = mapcurve(series, ObservableQuery(:projectedseparation, :b, :A))
    @test length(mc) == length(series.ts)
    @test all(>(0), mc)

    # default sky queries: one per row, exterior about interior
    qs = Octofitter.default_sky_queries(model.system)
    @test length(qs) == 1
    @test qs[1][2] === :b
    @test Octofitter._refstr(qs[1][1].target) == "b"
    @test Octofitter._refstr(qs[1][1].ref) == "A"
end

@testset "the plotting epoch grid can be set" begin
    # Every panel draws over `series.ts` and clips its axis to it, so this is
    # the only way to see the orbit outside the data — `xlims!` past the grid
    # gives an empty axis, not more curve.
    model, obs = _plotting_test_model()
    rng = Xoshiro(21)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:10]
    chain = Octofitter.result2mcmcchain(nts)
    auto = PosteriorSeries(model, chain; N=3)
    tlast = maximum(auto.data_epochs)

    # One end at a time: the other stays exactly where the derivation put it.
    ext = PosteriorSeries(model, chain; N=3, tmax=tlast + 3650)
    @test last(ext.ts) ≈ tlast + 3650
    @test first(ext.ts) ≈ first(auto.ts)
    @test issorted(ext.ts) && allunique(ext.ts)
    # The grid is still allocated per orbital period rather than stretched over
    # a fixed number of points, so a longer span gets more of them.
    @test length(ext.ts) > length(auto.ts)
    @test length(modelcurves(ext, (radvel, :A, Barycentre))[1]) == length(ext.ts)

    # An epoch is an MJD number, or anything `mjd` understands.
    @test last(PosteriorSeries(model, chain; N=2, tmax="2035-01-01").ts) ≈ mjd("2035-01-01")

    # A window is a window: a data epoch outside it is not merged back in,
    # which would stretch every panel's axis past what was asked for.
    win = PosteriorSeries(model, chain; N=3, tmin=59300.0)
    @test first(win.ts) ≈ 59300.0
    @test all(>=(59300.0), win.ts)
    @test 59000.0 in auto.ts          # …and it is there by default

    # `ts=` replaces the grid outright.
    grid = collect(59000.0:25.0:59700.0)
    exact = PosteriorSeries(model, chain; N=3, ts=grid)
    @test exact.ts == grid
    @test length(modelcurves(exact, (radvel, :A, Barycentre))[1]) == length(grid)
    @test exact.ts == PosteriorSeries(model, chain; N=3, ts=reverse(grid)).ts

    @test_throws ErrorException PosteriorSeries(model, chain; N=2, tmin=60000.0, tmax=59000.0)
    @test_throws ErrorException PosteriorSeries(model, chain; N=2, tmin=60000.0, tmax=60000.0)
    @test_throws ErrorException PosteriorSeries(model, chain; N=2, ts=grid, tmax=60000.0)
    @test_throws ErrorException PosteriorSeries(model, chain; N=2, ts=Float64[])
end

@testset "relative astrometry projects between its two bases" begin
    # (sep, pa) and (Δα*, Δδ) are the same measurement rotated, with no free
    # parameter in between, so a fit mixing the two conventions can put every
    # point on every panel. What must hold is that the projection is exact and
    # that it does not touch the residual the likelihood scores.
    model, obs = _plotting_test_model()
    rng = Xoshiro(11)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:10]
    series = PosteriorSeries(model, Octofitter.result2mcmcchain(nts); N=4)

    r = Octofitter.residuals(obs.seppa, obscontext(series, obs.seppa))
    @test r.sep.data ≈ hypot.(r.raoff.data, r.decoff.data) rtol = 1e-12
    @test r.raoff.data ≈ r.sep.data .* sin.(deg2rad.(r.pa.data)) rtol = 1e-12
    @test r.decoff.data ≈ r.sep.data .* cos.(deg2rad.(r.pa.data)) rtol = 1e-12
    # The model side projects the same way, so a derived residual is the
    # derived data minus the derived model and nothing else.
    @test r.raoff.resid ≈ r.raoff.data .- r.raoff.model rtol = 1e-12
    # Derived σ is positive and finite; the jitter still widens σ_eff.
    @test all(>(0), r.raoff.σ) && all(isfinite, r.raoff.σ)
    @test all(r.raoff.σ_eff .>= r.raoff.σ)

    # …and the other way round, from a ra/dec table.
    r2 = Octofitter.residuals(obs.radec, obscontext(series, obs.radec))
    @test r2.sep.data ≈ hypot.(r2.raoff.data, r2.decoff.data) rtol = 1e-12
    @test r2.pa.data ≈ rad2deg.(rem2pi.(atan.(r2.raoff.data, r2.decoff.data), RoundDown)) rtol = 1e-10

    # The native channels are untouched by any of this: `ln_like` still agrees.
    ctx = obscontext(series, obs.radec)
    ll = -0.5 * sum(@. r2.raoff.resid^2 / r2.raoff.σ^2 + log(2π * r2.raoff.σ^2)) +
         -0.5 * sum(@. r2.decoff.resid^2 / r2.decoff.σ^2 + log(2π * r2.decoff.σ^2))
    @test ll ≈ Octofitter.ln_like(obs.radec, ctx) rtol = 1e-12
end

@testset "datacalibration is the exact inverse of what residuals removes" begin
    # A panel that draws raw measurements adds this back to the model curve;
    # `residuals` subtracts it from the data. If the two ever disagreed, the
    # points would sit off their own curve — so they come from one function.
    model, obs = _plotting_test_model()
    rng = Xoshiro(12)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:10]
    series = PosteriorSeries(model, Octofitter.result2mcmcchain(nts); N=3)

    ctx = obscontext(series, obs.rvs)
    ch = first(plotchannels(obs.rvs))
    r = Octofitter.residuals(obs.rvs, ctx).rv
    cal = Octofitter.datacalibration(obs.rvs, ch, ctx, r.epoch)
    @test cal ≈ fill(12.0, length(r.epoch))          # the declared offset, no trend
    @test r.data .+ cal ≈ collect(obs.rvs.table.rv) rtol = 1e-12
    @test r.model .+ cal ≈ Octofitter.simulate(obs.rvs, ctx).rv_model rtol = 1e-12

    # Observations with nothing to calibrate say so, rather than returning zeros
    # a caller would have to recognise.
    @test Octofitter.datacalibration(obs.seppa, first(plotchannels(obs.seppa)),
        obscontext(series, obs.seppa), [59000.0]) === nothing
    # No correlated-noise model here, and none invented.
    @test Octofitter.noisemodel(obs.rvs, ctx, [59000.0]) === nothing
    @test !haskey(r, :gp_mean)
end

# ---------------------------------------------------------------------------
# Correlated noise on the plot side.
#
# An exact dense squared-exponential GP, which core can carry with no GP
# package at all: three methods, and the prediction in closed form so the test
# can check the plotting hooks against the same arithmetic the likelihood runs.
# (`test/v2/radial-velocity.jl` has a sibling of this for the likelihood; the
# names differ because both files are included into one namespace.)
# ---------------------------------------------------------------------------

struct PlotGP
    amp::Float64
    len::Float64
end
_plotk(g::PlotGP, s, t) = g.amp^2 * exp(-(s - t)^2 / (2 * g.len^2))

# Parametrized on the covariance element type, so the model's ForwardDiff
# gradient goes through it like any other likelihood term.
struct PlotGPFit{T}
    gp::PlotGP
    x::Vector{Float64}
    Σ::Matrix{T}
end

function Octofitter.gp_condition(g::PlotGP, x, σ²)
    xs = collect(Float64, x)
    Σ = [_plotk(g, xs[i], xs[j]) + (i == j ? σ²[i] : zero(eltype(σ²)))
         for i in eachindex(xs), j in eachindex(xs)]
    return PlotGPFit(g, xs, Σ)
end
function Octofitter.gp_ln_like(f::PlotGPFit, r)
    rr = collect(r)
    C = cholesky(Symmetric(f.Σ))
    return -(dot(rr, C \ rr) + logdet(C) + length(rr) * log(2π)) / 2
end
function Octofitter.gp_predict(f::PlotGPFit, r, xs)
    Ks = [_plotk(f.gp, Float64(x), f.x[j]) for x in xs, j in eachindex(f.x)]
    C = cholesky(Symmetric(f.Σ))
    μ = Ks * (C \ collect(r))
    # Latent predictive variance: a point's own white noise is added by
    # whoever consumes this, not here.
    v = [_plotk(f.gp, Float64(xs[i]), Float64(xs[i])) - dot(Ks[i, :], C \ Ks[i, :])
         for i in eachindex(xs)]
    return μ, v
end

# A backend that scores fits but cannot predict — the case a figure has to
# survive rather than throw on.
struct MutePlotGP end
struct MutePlotGPFit end
Octofitter.gp_condition(::MutePlotGP, x, σ²) = MutePlotGPFit()
Octofitter.gp_ln_like(::MutePlotGPFit, r) = 0.0

# …and one whose covariance will not factorize for this draw, which `ln_like`
# already turns into `-Inf` rather than a crashed chain.
struct SingularPlotGP end
Octofitter.gp_condition(::SingularPlotGP, x, σ²) = throw(PosDefException(1))

const _GP_EPOCHS = collect(range(59000.0, 59180.0, length=45))
# Smooth stellar activity: two rotation-scale harmonics, the structure a GP is
# for, plus white noise from a fixed seed. Deterministic, so the variance ratio
# below is a number, not a draw.
const _GP_ACTIVITY = @. 18.0 * sinpi(2 * (_GP_EPOCHS - 59000) / 23) +
                        7.0 * sinpi(2 * (_GP_EPOCHS - 59000) / 9 + 0.3)
const _GP_WHITE = 3.0 .* randn(Xoshiro(20250814), length(_GP_EPOCHS))

function _gp_rv_model(rv, gp)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass ~ Uniform(0.004, 0.006)
        a ~ Uniform(0.99, 1.01)
        e = 0.05
        i = 1.2
        ω = 0.4
        Ω = 1.0
        tp = 59000.0
    end)
    rvs = RadialVelocityObs(
        Table(epoch=_GP_EPOCHS, rv=rv, σ_rv=fill(3.0, length(_GP_EPOCHS)));
        target=A, ref=Barycentre, name="HARPS", gaussian_process=θ_obs -> gp,
        variables=@variables begin
            offset = 120.0
            jitter = 1.5
        end)
    sys = Octofitter.System(name="gpplot", bodies=[A, b], observations=[rvs],
        variables=@variables begin
            plx = 25.0
        end)
    return Octofitter.LogDensityModel(sys, verbosity=0), rvs
end

# One draw of a seeded prior "chain". The data never enter the forward model,
# so the same seed gives the same draw for every dataset built below — which
# is what lets the second pass plant activity on top of the first pass' own
# model prediction.
function _gp_series(model; seed=31, ii=[1])
    rng = Xoshiro(seed)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:6]
    return PosteriorSeries(model, Octofitter.result2mcmcchain(nts); ii)
end

# The activity dataset both GP testsets below run on: the model's own
# prediction under the first draw, plus the planted activity and white noise,
# so that draw's residual is exactly `activity + noise`.
function _gp_activity_model(gp)
    m0, o0 = _gp_rv_model(zeros(length(_GP_EPOCHS)), gp)
    truth = collect(Float64,
        Octofitter.simulate(o0, obscontext(_gp_series(m0), o0)).rv_model)
    return _gp_rv_model(truth .+ _GP_ACTIVITY .+ _GP_WHITE, gp)
end

@testset "a Gaussian process is conditioned on the draw's own residuals" begin
    gp = PlotGP(20.0, 6.0)
    model, obs = _gp_activity_model(gp)
    series = _gp_series(model)
    ctx = obscontext(series, obs)
    r = Octofitter.residuals(obs, ctx).rv

    @test r.resid ≈ _GP_ACTIVITY .+ _GP_WHITE rtol = 1e-10
    # `residuals` publishes the noise model alongside the plain residual, and
    # `resid` itself stays `data − model`: subtracting the GP is the consumer's
    # call, and a whitened strip and a phase fold both make it.
    @test haskey(r, :gp_mean) && haskey(r, :gp_var)
    @test length(r.gp_mean) == length(_GP_EPOCHS)
    @test all(>=(0), r.gp_var)
    # σ_eff is still measurement error and jitter; the GP variance is a
    # separate key, added by whoever draws the bar.
    @test r.σ_eff ≈ hypot.(collect(Float64, obs.table.σ_rv), 1.5) rtol = 1e-12

    # The conditioning convention, pinned: on `data − model` with the offset
    # and the trend already in the model, against σ_rv² + jitter². (v8
    # conditioned on the same quantity; where the two could differ — what
    # counts as "the model" — this follows v2's `residuals`.)
    fx = Octofitter.gp_condition(gp, _GP_EPOCHS,
        collect(Float64, obs.table.σ_rv .^ 2 .+ 1.5^2))
    μ, v = Octofitter.gp_predict(fx, r.resid, _GP_EPOCHS)
    @test r.gp_mean ≈ μ rtol = 1e-10
    @test r.gp_var ≈ max.(v, 0.0) rtol = 1e-10

    # And it does the job it is drawn for. Raw, these residuals are the
    # activity: many σ from zero and visibly structured. With the GP's mean
    # removed and its variance in the bar — what a whitened strip and a phase
    # fold both use — they are a fraction of a σ.
    @test var(r.resid .- r.gp_mean) < var(r.resid) / 5
    σ_net = sqrt.(r.σ_eff .^ 2 .+ r.gp_var)
    @test std(r.resid ./ r.σ_eff) > 3
    @test std((r.resid .- r.gp_mean) ./ σ_net) < 1.5

    # Off the data epochs too — the band a panel draws is on the model grid.
    nm = Octofitter.noisemodel(obs, ctx, [59000.5, 59090.0, 59179.0])
    @test nm !== nothing && length(nm.mean) == 3 && all(>=(0), nm.var)
end

@testset "many draws: every curve carries its own conditioned activity" begin
    # What `octoplot` draws over a fitted noise model. Each draw's curve is its
    # own orbit plus its own conditioned GP mean — not the MAP draw's, and not
    # nothing, which is what a many-draw RV panel showed before: Keplerians
    # over residuals the fit does not claim are white.
    gp = PlotGP(20.0, 6.0)
    model, obs = _gp_activity_model(gp)
    series = _gp_series(model; ii=[1, 2, 3])
    nc = Octofitter.noisecurves(series, obs, _GP_EPOCHS)

    @test nc !== nothing
    @test length(nc) == 3
    @test all(length(c.mean) == length(_GP_EPOCHS) for c in nc)
    @test all(all(>=(0), c.var) for c in nc)

    # Per draw, and each one on its *own* residuals: the same quantity
    # `residuals` publishes for that draw, which is what the residual strip and
    # the phase fold subtract. Not one activity model shared by the ensemble.
    for d in 1:3
        rd = Octofitter.residuals(obs, obscontext(series, obs; draw=d)).rv
        @test nc[d].mean ≈ rd.gp_mean rtol = 1e-12
        @test nc[d].var ≈ rd.gp_var rtol = 1e-12
    end
    @test !(nc[1].mean ≈ nc[2].mean)

    # Draw 1 is the draw the dataset was built on, so its residual is exactly
    # the injected activity plus white noise — and the curve's activity term
    # is that activity, not a fit to the noise.
    @test cor(nc[1].mean, _GP_ACTIVITY) > 0.9
    # The curve the panel draws (orbit + offset + trend + activity) sits an
    # order of magnitude closer to the measurements than the Keplerian alone.
    for d in 1:3
        ctx = obscontext(series, obs; draw=d)
        # `rv_model` carries the offset and the trend, so this is the Keplerian
        # curve as a panel draws it.
        orbit = collect(Float64, Octofitter.simulate(obs, ctx).rv_model)
        data = collect(Float64, obs.table.rv)
        @test std(data .- orbit .- nc[d].mean) < std(data .- orbit) / 2
    end

    # Off the data epochs — a curve is drawn on the whole plotting grid, and a
    # stationary kernel's conditional mean has to decay back to the orbit
    # rather than leave a step at the last measurement.
    far = [first(_GP_EPOCHS) - 400.0, last(_GP_EPOCHS) + 400.0]
    ncf = Octofitter.noisecurves(series, obs, far)
    @test all(all(isfinite, c.mean) for c in ncf)
    @test all(maximum(abs, c.mean) < 1e-6 * maximum(abs, nc[1].mean) for c in ncf)

    # `ndraws=` truncates the family the way a panel's does, and an
    # observation with no `gaussian_process` has no family at all.
    @test length(Octofitter.noisecurves(series, obs, _GP_EPOCHS; ndraws=2)) == 2
    model2, obs2 = _plotting_test_model()
    rng2 = Xoshiro(7)
    series2 = PosteriorSeries(model2, Octofitter.result2mcmcchain(
        [model2.arr2nt(collect(model2.sample_priors(rng2))) for _ in 1:4]); N=3)
    @test Octofitter.noisecurves(series2, obs2.rvs, [59000.0]) === nothing
end

@testset "a noise model that cannot be predicted degrades to the orbit alone" begin
    # A figure is not a fit: a backend with no `gp_predict`, or a draw whose
    # covariance will not factorize, costs the band and the GP-corrected
    # residual — not the whole plot. (Both emit a warning; `gp_predict` itself
    # still throws, so cross-validation fails loudly. See
    # `test/v2/radial-velocity.jl`.)
    for gp in (MutePlotGP(), SingularPlotGP())
        model, obs = _gp_rv_model(_GP_ACTIVITY .+ _GP_WHITE .+ 120.0, gp)
        ctx = obscontext(_gp_series(model), obs)
        @test Octofitter.noisemodel(obs, ctx, [59000.0, 59100.0]) === nothing
        r = Octofitter.residuals(obs, ctx).rv
        @test !haskey(r, :gp_mean) && !haskey(r, :gp_var)
        @test length(r.resid) == length(_GP_EPOCHS)
        # And over many draws the whole curve family degrades together, so the
        # panel is Keplerian-only rather than activity on some draws and not
        # others — one warning, not one per draw.
        multi = _gp_series(model; ii=[1, 2, 3])
        @test Octofitter.noisecurves(multi, obs, [59000.0, 59100.0]) === nothing
    end
end

@testset "predicted channels for observables the model has no data for" begin
    model, _ = _plotting_test_model()
    sys = model.system

    # Separations are per hierarchy row, exterior about interior.
    qs = Octofitter.default_queries(sys, PlanetOrbits.projectedseparation)
    @test length(qs) == 1
    @test Octofitter._refstr(qs[1][1].target) == "b"
    @test Octofitter._refstr(qs[1][1].ref) == "A"

    # Reflex signals are the host against the barycentre — not the planet's
    # velocity relative to its star, which is what a per-row rule would give.
    qv = Octofitter.default_queries(sys, radvel)
    @test length(qv) == 1
    @test Octofitter._refstr(qv[1][1].target) == "A"
    @test qv[1][1].ref === Octofitter.refspec(Barycentre)

    # The channel carries display metadata from `plotinfo`, with radians
    # converted to degrees the way every other angular channel is.
    chs = Octofitter.predictedchannels(sys, posangle)
    @test length(chs) == 1
    obs_, ch = chs[1]
    @test obs_ === nothing                       # a model-only channel
    @test ch.name === :posangle
    @test ch.unit == "°" && ch.scale ≈ rad2deg(1.0) && ch.wrap ≈ 360.0
    @test Octofitter.predictedchannels(sys, radvel)[1][2].unit == "m/s"

    # Only an observable's function or name asks for a prediction; a channel
    # name asks for data, and missing data is not a request to invent it.
    @test Octofitter._requestedobservables(radvel) == [radvel]
    @test Octofitter._requestedobservables(:radvel) == [radvel]
    @test isempty(Octofitter._requestedobservables(:sep))
    @test Octofitter._requestedobservables((radvel, :posangle)) == [radvel, posangle]
    @test isempty(Octofitter._requestedobservables(nothing))
end

@testset "row signals and phase folding" begin
    model, obs = _plotting_test_model()
    rng = Xoshiro(4)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:20]
    chain = Octofitter.result2mcmcchain(nts)
    series = PosteriorSeries(model, chain; N=5)
    posys = series.sys_map

    # single row: every signal is the full query, exactly
    q = ObservableQuery(radvel, :A, Barycentre)
    @test Octofitter.foldablerows(posys, q) == [1]
    sig = Octofitter.rowsignal(posys, q, 1)
    @test !sig.scaled
    @test Octofitter.signalcoeff(sig, posys) == 1.0
    traj = series.data_traj_map
    @test Octofitter.evalsignal(sig, posys, traj) ≈ Octofitter.evalquery(q, posys, traj)

    # nonlinear observables fold too when a single row is the whole signal
    qsep = ObservableQuery(projectedseparation, :b, :A)
    @test Octofitter.foldablerows(posys, qsep) == [1]

    # out-of-range and unaffected rows return nothing
    @test Octofitter.rowsignal(posys, q, 2) === nothing

    # fold ephemeris: phase zero is in the cycle containing tmid; radvel
    # convention pins an upward zero crossing of the signal at phase 0
    tmid = sum(extrema(series.data_epochs)) / 2
    P, t0 = Octofitter.foldephemeris(sig, posys, tmid)
    @test P ≈ Octofitter.PlanetOrbits.period(posys, 1)
    @test abs(t0 - tmid) <= P
    v = Octofitter.evalsignal(sig, posys,
        orbitsolve(posys, [t0 - 0.005P, t0 + 0.005P]))
    @test v[1] <= v[2]                       # increasing through the crossing
    @test abs(Octofitter.foldphase(t0, P, t0)) < 1e-12
    @test Octofitter.foldphase(t0 + 0.25P, P, t0) ≈ 0.25
    @test Octofitter.foldphase(t0 + 0.75P, P, t0) ≈ -0.25

    # defaultpanels: the generic default opts in to generic panels
    @test defaultpanels(obs.rvs) === ()
end

@testset "phase-fold binned means and their error bars" begin
    pbm = Octofitter.phasebinmeans
    σ = 4.0
    wof(n) = fill(1 / σ^2, n)
    PW = Octofitter.ProbabilityWeights

    # The two halves of the bar, written out here independently of the
    # implementation: the analytic error on a weighted mean, and the bin's own
    # weighted scatter over √n. `sigma` is their max.
    analytic(w) = 1 / sqrt(sum(w))
    empirical(y, w) = Octofitter.std(y, PW(w); mean=Octofitter.mean(y, PW(w)),
                                    corrected=true) / sqrt(length(y))
    bar(y, w) = max(analytic(w), empirical(y, w))

    # A bin is [lo, hi): 20 bins put the edges on multiples of 0.05.
    b = pbm([-0.05, -0.04], [10.0, 20.0], wof(2), 20)
    @test b.centre ≈ [-0.025]
    @test b.mean ≈ [15.0]
    @test length(b.sigma) == 1 && isfinite(b.sigma[1]) && b.sigma[1] > 0

    # Two 4 m/s measurements 10 m/s apart: the data disagrees with the noise
    # model, so the empirical term wins — 7.071/√2 = 5.0 against a floor of
    # 1/√(2/16) = 2.83. `s_w` is the frequency-weights-corrected weighted
    # scatter, which for weights that are 1/σ² (a variance scale, not a count)
    # means StatsBase's `ProbabilityWeights` correction: Bessel's n/(n−1) on
    # the weighted rms, equivalently `FrequencyWeights` on weights
    # renormalised to sum to n.
    @test analytic(wof(2)) ≈ 4 / sqrt(2)
    @test empirical([10.0, 20.0], wof(2)) ≈ 5.0
    @test b.sigma ≈ [5.0]
    @test b.sigma ≈ [bar([10.0, 20.0], wof(2))]
    # It coincides with the pre-v9 bar here only through the n = 2 accident
    # that √(n/(n−1))/√n = 1: v8 drew the *uncorrected* weighted scatter
    # √(Σw(y−μ)²/Σw), which is a spread of points and not an error on a mean.
    @test Octofitter.std([10.0, 20.0], PW(wof(2))) ≈ 5.0

    # The floor binds when the points agree: two 4 m/s measurements half a m/s
    # apart still cannot make a mean better than 4/√2, however tightly they
    # happen to land. v8 drew 0.25 m/s here — below the mathematical floor.
    let y = [10.0, 10.5], r = pbm([-0.05, -0.04], y, wof(2), 20)
        @test empirical(y, wof(2)) ≈ 0.25
        @test Octofitter.std(y, PW(wof(2))) ≈ 0.25          # the old bar
        @test r.sigma ≈ [4 / sqrt(2)]
        @test r.sigma ≈ [analytic(wof(2))]
        @test r.sigma[1] > Octofitter.std(y, PW(wof(2)))
    end

    # The exact crossover, in the white-noise limit the analytic term is built
    # for: three 4 m/s points whose corrected scatter comes out at exactly
    # 4 m/s (deviations −4, 0, +4 → √(32/2) = 4). Both terms are then σ/√n,
    # so the bar is the plain RadVel/juliet error on the weighted mean.
    let y = [16.0, 20.0, 24.0], w = wof(3), r = pbm([0.01, 0.02, 0.03], y, w, 20)
        @test Octofitter.std(y, PW(w); corrected=true) ≈ σ
        @test empirical(y, w) ≈ σ / sqrt(3)
        @test analytic(w) ≈ σ / sqrt(3)
        @test r.sigma ≈ [σ / sqrt(3)]
    end

    # Unequal weights pull the mean, exactly as Σwy/Σw, and the bar follows
    # the same two-term rule with Σw over the bin's own weights.
    let w = [1 / 1.0^2, 1 / 3.0^2]
        bb = pbm([0.01, 0.02], [0.0, 12.0], w, 20)
        @test bb.mean[1] ≈ sum(w .* [0.0, 12.0]) / sum(w)
        @test bb.sigma ≈ [bar([0.0, 12.0], w)]
    end

    # The sparse limit. One point is not a mean of itself: it is already on
    # the axis at its own phase with its own error bar, and binning it would
    # move it to the bin centre and give it a scatter of zero. So nothing is
    # returned — no marks, no NaN, no zero-length error bars.
    @test isempty(pbm([0.0966], [29.7], wof(1), 20).centre)
    for n in (1, 2, 3, 5, 10)
        x = n == 1 ? [0.0] : collect(range(-0.45, 0.45, length=n))  # ≤1 per bin
        y = 30 .* sin.(2π .* x)
        r = pbm(x, y, wof(n), 20)
        @test isempty(r.centre) && isempty(r.mean) && isempty(r.sigma)
    end

    # Bins that hold a real average still draw, however few there are: three
    # points, two of them sharing a bin.
    let r = pbm([-0.44, -0.43, 0.31], [1.0, 3.0, 9.0], wof(3), 20)
        @test length(r.centre) == 1
        @test r.centre ≈ [-0.425]
        @test r.mean ≈ [2.0]
    end

    # Dense data: 60 points over 20 bins is three per bin. Means are as they
    # were; the bars are the hybrid, and — the point of the change — not one
    # of the twenty now sits below 1/√(3/16) = 2.31 m/s, where under the old
    # convention fourteen of them did (on the rendered 60-point figure it was
    # seven, with the smallest bar at 0.44 m/s).
    let n = 60
        x = collect(range(-0.5, 0.5, length=n + 1))[1:n]
        y = 30 .* sin.(2π .* x) .+ 0.5 .* (-1) .^ (1:n)
        w = wof(n)
        r = pbm(x, y, w, 20)
        @test length(r.centre) == 20
        @test all(isfinite, r.mean) && all(isfinite, r.sigma)
        floor3 = 1 / sqrt(3 / σ^2)
        @test floor3 ≈ σ / sqrt(3)
        @test all(>=(floor3 - 1e-12), r.sigma)
        # …and the old convention really did undershoot it, so this is a
        # change and not a no-op assertion.
        edges = range(-0.5, 0.5, length=21)
        old = [Octofitter.std(y[(edges[k].<=x).&(x.<edges[k+1])],
                              PW(w[(edges[k].<=x).&(x.<edges[k+1])])) for k in 1:20]
        @test count(<(floor3), old) == 14
        @test minimum(old) < 0.5
        for k in eachindex(r.centre)
            @test r.centre[k] ≈ (edges[k] + edges[k+1]) / 2
            m = (edges[k] .<= x) .& (x .< edges[k+1])
            @test count(m) == 3
            @test r.mean[k] ≈ Octofitter.mean(y[m], PW(w[m]))
            @test r.sigma[k] ≈ bar(y[m], w[m])
        end
    end

    # Many points and a scatter far wider than the noise model: the empirical
    # term takes over, by the ratio of the true spread to the quoted σ. Twenty
    # 4 m/s points spread evenly over ±20 m/s have a floor of 4/√20 = 0.894
    # and a corrected scatter of 12.45, so the bar is 2.785 — 3.1× the floor.
    let n = 20
        x = fill(0.31, n)
        y = collect(range(-20.0, 20.0, length=n))
        w = wof(n)
        r = pbm(x, y, w, 20)
        @test length(r.sigma) == 1
        @test r.sigma ≈ [empirical(y, w)]
        @test analytic(w) ≈ σ / sqrt(n)
        @test r.sigma[1] > 3 * analytic(w)
        @test r.sigma[1] ≈ 2.78500138006799
    end

    # Empty bins are skipped rather than drawn at zero — and a bin whose
    # points happen to agree exactly is no longer drawn at zero either: its
    # empirical term is 0 and the analytic floor is the whole bar. Under the
    # old convention both of these marks had error bars of length zero.
    let x = vcat(fill(-0.44, 3), fill(0.31, 4)), y = vcat(fill(2.0, 3), fill(-5.0, 4))
        r = pbm(x, y, wof(7), 20)
        @test r.centre ≈ [-0.425, 0.325]
        @test r.mean ≈ [2.0, -5.0]
        @test all(isfinite, r.sigma)
        @test all(>(0), r.sigma)
        @test r.sigma ≈ [σ / sqrt(3), σ / sqrt(4)]
        @test Octofitter.std(fill(2.0, 3), PW(wof(3))) == 0.0   # the old bar
    end
end

@testset "multi-row signal decomposition" begin
    # Two planets about one star (astrocentric rows): the row signals of a
    # linear observable must sum to the full query.
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 0.01
        a ~ Uniform(2.9, 3.1)
        e = 0.2
        i = 0.6
        ω = 1.0
        Ω = 2.0
        tp = 59200.0
    end)
    c = Octofitter.Body(name="c", about=A, variables=@variables begin
        mass = 0.004
        a = 1.0
        e = 0.05
        i = 0.7
        ω = 0.3
        Ω = 2.0
        tp = 59250.0
    end)
    rvs = RadialVelocityObs(
        Table(epoch=[59050.0, 59250.0, 59450.0], rv=[1.0, -1.0, 0.5], σ_rv=[3.0, 3.0, 3.0]);
        target=A, ref=Barycentre, name="RV")
    sys = Octofitter.System(name="rowsigtest", bodies=[A, b, c], observations=[rvs],
        variables=@variables begin
            plx = 50.0
        end)
    model = Octofitter.LogDensityModel(sys, verbosity=0)
    rng = Xoshiro(5)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:10]
    series = PosteriorSeries(model, Octofitter.result2mcmcchain(nts); N=3)
    posys = series.sys_map
    traj = series.data_traj_map

    q = ObservableQuery(radvel, :A, Barycentre)
    @test Octofitter.foldablerows(posys, q) == [1, 2]
    s1 = Octofitter.rowsignal(posys, q, 1)
    s2 = Octofitter.rowsignal(posys, q, 2)
    @test s1.scaled && s2.scaled
    # The decomposition is of the *kinematic* signal: `radvel`'s Einstein term
    # is quadratic in velocity and 1/r in separation, so it does not telescope
    # into per-row contributions, and a multi-row fold leaves it common-mode
    # (see `_rowfunc`). So the parts reproduce the kinematic total exactly...
    kin = ObservableQuery(Octofitter._kinrv, :A, Barycentre)
    full = Octofitter.evalquery(q, posys, traj)
    fullkin = Octofitter.evalquery(kin, posys, traj)
    parts = Octofitter.evalsignal(s1, posys, traj) .+ Octofitter.evalsignal(s2, posys, traj)
    @test fullkin ≈ parts rtol = 1e-10
    # ...and what they leave behind is precisely the Einstein term, not a
    # decomposition error.
    ein = [radvel(sol, PlanetOrbits.BodyRef(1), PlanetOrbits.barycentre(posys)) -
           Octofitter._kinrv(sol, PlanetOrbits.BodyRef(1), PlanetOrbits.barycentre(posys))
           for sol in traj]
    @test full .- parts ≈ ein rtol = 1e-10
    @test maximum(abs, ein) > 0

    # sep(b, A) is moved only by row 1 under astrocentric rows — exact fold —
    # and row 2's contribution to it cannot be isolated (nonlinear)
    qsep = ObservableQuery(projectedseparation, :b, :A)
    @test Octofitter.foldablerows(posys, qsep) == [1]
    @test !Octofitter.rowsignal(posys, qsep, 1).scaled
end

@testset "RV residuals subtract the trend as well as the offset" begin
    # The `radvel` query the panel draws its curve from knows about neither the
    # instrument offset nor the long-term trend, so `residuals` must remove
    # both or the points will not lie on their own curve. `trend_function`
    # only arrived on the core RV type in v2, so this is the case that was
    # missing.
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 5mjup
        a = 3.0; e = 0.1; i = 0.6; ω = 1.0; Ω = 0.4; tp = 56000.0
    end)
    ep = collect(range(56000.0, 58000.0, length=6))
    rv = RadialVelocityObs(
        Table(epoch=ep, rv=zeros(length(ep)), σ_rv=fill(5.0, length(ep)));
        target=A, ref=Barycentre, name="rvs",
        trend_function=(θ_obs, epoch) -> θ_obs.slope * (epoch - 57000),
        variables=@variables begin
            offset = 30.0
            slope = 0.05          # m/s/day — 50 m/s across this baseline
            jitter = 1.0
        end)
    sys = Octofitter.System(name="trend", bodies=[A, b], observations=[rv],
        variables=@variables begin plx = 25.0 end)

    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    eps, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, eps)
    ctx = Octofitter.ObsContext(nt, nt.observations.rvs, posys, traj, maps[rv])

    r = Octofitter.residuals(rv, ctx)
    tgt = Octofitter.resolveref(posys, rv.target)
    rf = Octofitter.resolveref(posys, rv.ref)
    query = [PlanetOrbits.radvel(traj[k], tgt, rf) for k in maps[rv]]
    @test r.rv.model ≈ query rtol = 1e-12
    # The residual itself is against the *full* model and is unaffected.
    sim = Octofitter.simulate(rv, ctx)
    @test r.rv.resid ≈ rv.table.rv .- sim.rv_model rtol = 1e-12
    @test r.rv.data .- r.rv.model ≈ r.rv.resid rtol = 1e-12
end

@testset "octoplot without Makie errors helpfully" begin
    # (Makie is not loaded in the test env, so the stub must fire.)
    if !haskey(Base.loaded_modules_array() |> ms -> Dict(nameof(m) => m for m in ms), :Makie)
        model, _ = _plotting_test_model()
        @test_throws ErrorException octoplot(model)
    end
end

@testset "a wrapped observation is still plotted as its data" begin
    # `ObsPriorONeil2019` reweights the prior; it does not change what the
    # instrument measured. Before the delegation below it answered the
    # `AbstractObs` defaults — no channels — so `plottable_observations`
    # dropped it and `octoplot` drew the orbit tracks of a wrapped relative
    # astrometry fit with none of its data points on them.
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 0.005
        a ~ Uniform(2.0, 4.0)
        e = 0.25; i = 0.6; ω = 1.0; Ω = 2.0; tp = 59200.0
    end)
    tab = Table(epoch=[59000.0, 59300.0, 59600.0],
        ra=[200.0, 180.0, 160.0], dec=[150.0, 170.0, 165.0],
        σ_ra=[4.0, 4.0, 4.0], σ_dec=[4.0, 4.0, 4.0])
    mkastrom() = RelAstromObs(tab; target=b, ref=A, name="INST",
        variables=@variables begin
            jitter = 1.5
            platescale = 1.02
            northangle = 0.03
        end)
    astrom = mkastrom()
    wrapped = ObsPriorONeil2019(mkastrom())

    # The protocol delegates, so a wrapper answers exactly what it wraps.
    @test Octofitter.plotobs(wrapped) === wrapped.wrapped_like
    @test Octofitter.plotobs(astrom) === astrom
    @test [ch.name for ch in Octofitter.plotchannels(wrapped)] ==
          [ch.name for ch in Octofitter.plotchannels(astrom)]
    @test Octofitter.sharepanel(wrapped) == Octofitter.sharepanel(astrom)

    function series_of(o)
        sys = Octofitter.System(name="w", bodies=[A, b], observations=[o],
            variables=@variables begin
                plx = 50.0
            end)
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        rng = Xoshiro(7)
        nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:20]
        return model, PosteriorSeries(model, Octofitter.result2mcmcchain(nts); N=5)
    end

    mw, sw = series_of(wrapped)
    mp, sp = series_of(astrom)
    @test length(Octofitter.plottable_observations(mw.system)) == 1
    @test isempty(Octofitter.unplottable_observations(mw.system))

    # The *wrapper* is what carries a context: the fitted platescale/northangle
    # /jitter are registered under `obspri_INST`, so unwrapping before building
    # one would silently fall back to the defaults (the v8.3.0 bug, from the
    # other side). The calibrated data must therefore match the unwrapped fit's
    # exactly, not the uncalibrated table.
    rw = Octofitter.residuals(wrapped, obscontext(sw, wrapped))
    rp = Octofitter.residuals(astrom, obscontext(sp, astrom))
    @test keys(rw) == (:sep, :pa, :raoff, :decoff)
    @test rw.raoff.data ≈ rp.raoff.data
    @test rw.decoff.data ≈ rp.decoff.data
    @test rw.sep.data ≈ rp.sep.data
    @test all(rw.raoff.σ_eff .≈ hypot.(tab.σ_ra, 1.5))     # the jitter is live
    @test !(rw.sep.data ≈ collect(hypot.(tab.ra, tab.dec)))  # …and the platescale
end

# ---------------------------------------------------------------------------
# The protocol for the observation types that are not simple epoch series.
#
# These reuse the offline fixtures and model builders from `v2/g23h.jl` and
# `v2/hipparcos.jl`, which `runtests.jl` includes first.
# ---------------------------------------------------------------------------

@testset "every shipped observation type declares plot channels" begin
    # The audit helper is the contract: an observation that carries data and
    # declares no channels is drawn as an empty panel, which reads as a fit
    # with nothing constraining it rather than as a gap in the plotting layer.
    m = hip_model()
    @test isempty(Octofitter.unplottable_observations(m.sys))
    g = g23h_model()
    @test isempty(Octofitter.unplottable_observations(g.sys))
end

@testset "HipparcosIADObs channel agrees with ln_like" begin
    m = hip_model()
    c = hip_context(m.sys, m.obs)
    r = Octofitter.residuals(m.obs, c.ctx).along_scan
    n = length(m.obs.table.epoch)
    @test length(r.resid) == n
    @test r.use == .!collect(m.obs.table.reject)
    # The datum is the catalog abscissa residual, and the residual is still
    # exactly `measured − model` from the likelihood's own simulate.
    sim = Octofitter.simulate(m.obs, c.ctx)
    @test r.data ≈ collect(Float64, m.obs.table.res)
    @test r.resid ≈ collect(Float64, sim.resid)
    @test r.data .- r.model ≈ r.resid
    # σ carries the BINARYS inflation, σ_eff adds the fitted jitter.
    @test r.σ ≈ collect(m.obs.table.sres_renorm .* sim.σ_inflation)
    @test all(r.σ_eff .>= r.σ)
    u = r.use
    ll = -0.5 * sum(@. r.resid[u]^2 / r.σ_eff[u]^2 + log(2π * r.σ_eff[u]^2))
    @test ll ≈ Octofitter.ln_like(m.obs, c.ctx) rtol = 1e-12
end

@testset "G23HObs exposes proper motions and the Hipparcos abscissae" begin
    m = g23h_model()
    c = g23h_context(m.sys, m.obs)
    chs = Octofitter.plotchannels(m.obs)
    @test [ch.name for ch in chs] == [:pmra, :pmdec, :along_scan_hip]
    # The proper-motion channels carry a smooth curve (the host's reflex);
    # the abscissa channel cannot, since the scan angle is per transit.
    @test chs[1].query isa ObservableQuery && chs[1].query.func === PlanetOrbits.pmra
    @test chs[2].query.func === PlanetOrbits.pmdec
    @test chs[3].query === nothing

    r = Octofitter.residuals(m.obs, c.ctx)
    @test keys(r) == (:pmra, :pmdec, :along_scan_hip)
    # Five catalog proper motions per axis: Hip, H−G, DR2, DR3−DR2, DR3.
    @test length(r.pmra.data) == 5
    @test length(r.pmdec.data) == 5
    @test r.pmra.resid ≈ r.pmra.data .- r.pmra.model
    @test all(isfinite, r.pmra.σ) && all(>(0), r.pmra.σ)
    # Each point is an average over a mission window, and says which one.
    @test all(r.pmra.epoch_lo .<= r.pmra.epoch .<= r.pmra.epoch_hi)
    @test r.pmra.epoch_hi[1] - r.pmra.epoch_lo[1] > 365    # the Hipparcos mission

    # The residual is the likelihood's own catalog-minus-model difference: the
    # display calibration is one constant per axis and cancels out of it.
    mom = Octofitter._g23h_catalog_moments(m.obs, c.ctx)
    lut = Dict(l => k for (k, l) in enumerate(mom.labels))
    for (j, k) in enumerate((:ra_hip, :ra_hg, :ra_dr2, :ra_dr32, :ra_dr3))
        @test r.pmra.resid[j] ≈ mom.catalog[lut[k]] - mom.model[lut[k]]
        @test r.pmra.σ[j] ≈ mom.sigma[lut[k]]
    end

    # The Hipparcos channel is the same one `HipparcosIADObs` exposes, from the
    # same transits — and the signed residual squares to the unsigned distance
    # `simulate` publishes.
    sim = Octofitter.simulate(m.obs, c.ctx)
    @test abs.(sim.iad_resid_signed) == sim.iad_resid
    @test r.along_scan_hip.resid ≈ collect(Float64, sim.iad_resid_signed)
    @test length(r.along_scan_hip.data) == m.obs.n_hip
end

@testset "PhotometryObs opts out of the generic time-series panel" begin
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 0.005
        a = 3.0; e = 0.2; i = 0.6; ω = 1.0; Ω = 2.0; tp = 59200.0
        flux_H = 3.8
    end)
    phot = PhotometryObs(Table(phot=[3.9, 3.6, 4.1], σ_phot=[0.3, 0.4, 0.35]);
        target=b, band=:H, name="NIRC2")
    astrom = RelAstromObs(Table(epoch=[59100.0, 59400.0],
            ra=[100.0, 90.0], dec=[80.0, 95.0], σ_ra=[4.0, 4.0], σ_dec=[4.0, 4.0]);
        target=b, ref=A, name="GPI")
    sys = Octofitter.System(name="photplot", bodies=[A, b], observations=[phot, astrom],
        variables=@variables begin plx = 50.0 end)
    @test isempty(Octofitter.unplottable_observations(sys))
    @test [ch.name for ch in Octofitter.plotchannels(phot)] == [:phot]
    # Photometry has no epoch axis, so it declares its own panel rather than
    # being drawn against a meaningless time axis.
    @test !isempty(defaultpanels(phot))
    @test first(first(defaultpanels(phot))) === :phot

    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    eps, maps = Octofitter.epoch_plan(sys)
    ctx = Octofitter.ObsContext(nt, (;), posys, orbitsolve(posys, eps), maps[phot])
    r = Octofitter.residuals(phot, ctx).phot
    @test r.model == fill(3.8, 3)
    @test r.epoch == [1.0, 2.0, 3.0]          # row index: there is no epoch
    @test r.resid ≈ collect(phot.table.phot) .- 3.8
    ll = sum(@. -0.5 * r.resid^2 / r.σ^2 - log(sqrt(2π * r.σ^2)))
    @test ll ≈ Octofitter.ln_like(phot, ctx) rtol = 1e-12
end

# ---------------------------------------------------------------------------
# The parallax ellipse `gaiastarplot` annotates the corner with. The drawing
# lives in the Makie extension; the geometry lives in the package, so this is
# the half of the feature CI can actually hold to account.
# ---------------------------------------------------------------------------

# Ecliptic latitude of an equatorial direction, in degrees.
_eclat(α, δ) = asind(sind(δ) * cosd(Octofitter.obliquity_J2000_deg) -
                     cosd(δ) * sind(Octofitter.obliquity_J2000_deg) * sind(α))

@testset "parallax_ellipse has the semi-axes an ellipse of parallax must have" begin
    # `Δα✶ = A cos λ + B sin λ`, `Δδ = C cos λ + D sin λ` maps the unit circle
    # through `[A B; C D]`, so the semi-axes are that matrix's singular values.
    # They must come out 1 and |sin β| — the whole geometric content of the
    # annotation — and the check is on the matrix rather than on the sampled
    # curve, so it is exact rather than limited by `n`.
    function semiaxes(α, δ)
        sinα, cosα = sincosd(α); sinδ, cosδ = sincosd(δ)
        sinε, cosε = sincosd(Octofitter.obliquity_J2000_deg)
        return svdvals([sinα (-cosε*cosα); (cosα*sinδ) (cosε*sinα*sinδ-sinε*cosδ)])
    end
    for (α, δ) in ((0.0, 0.0), (45.0, 20.0), (158.307, -40.4256),
                   (300.0, -5.0), (12.0, 88.0), (270.0, 66.5607))
        smaj, smin = semiaxes(α, δ)
        @test smaj ≈ 1.0 atol = 1e-14
        @test smin ≈ abs(sind(_eclat(α, δ))) atol = 1e-14
        # …and the sampled curve agrees with the matrix it came from.
        e = parallax_ellipse(α, δ; n=20001)
        r = hypot.(e.raoff, e.decoff)
        @test maximum(r) ≈ smaj rtol = 1e-7
        @test minimum(r) ≈ smin atol = 1e-4
    end

    # The two degenerate cases the formula has to get right, since they are the
    # ones a reader would notice instantly on a figure: a source at the
    # ecliptic pole parallaxes in a circle, one on the ecliptic in a line.
    pole = parallax_ellipse(270.0, 90 - Octofitter.obliquity_J2000_deg; n=2001)
    @test all(≈(1.0; atol=1e-12), hypot.(pole.raoff, pole.decoff))
    onecl = parallax_ellipse(90.0, Octofitter.obliquity_J2000_deg; n=2001)
    @test minimum(hypot.(onecl.raoff, onecl.decoff)) < 1e-12

    # Closed curve, so `lines!`/`poly!` get a loop and not an open arc.
    e = parallax_ellipse(158.307, 40.4256; n=181)
    @test length(e.raoff) == length(e.decoff) == 181
    @test e.raoff[begin] ≈ e.raoff[end] && e.decoff[begin] ≈ e.decoff[end]
    @test_throws ArgumentError parallax_ellipse(0.0, 0.0; n=2)
end

@testset "parallax_ellipse is the locus Gaia's own parallax factors lie on" begin
    # The strong test, and the one that would catch a mirrored or rotated
    # ellipse: a real GOST forecast carries `parallaxFactorAlongScan`, which is
    # the archive's projection of the Earth's *ephemeris* position onto each
    # transit's scan direction. If our ellipse is the right locus in the right
    # frame, every one of those numbers is the projection of some point on it,
    # so none may exceed the ellipse's support function in that direction.
    #
    # Note this is not a restatement of the formula: the numbers come from
    # Gaia, the curve from a circular Earth orbit, and the two only agree if
    # the α/δ convention, the ε rotation and the (sin ψ, cos ψ) scan
    # projection are each right.
    f = joinpath(@__DIR__, "..", "GOST-158.30707896392835-40.42555422701387-dr3.csv")
    raw, hdr = readdlm(f, ',', String, header=true)
    col(name) = findfirst(==(name), strip.(vec(hdr)))
    # From the file's own columns, not its name: the name encodes the query,
    # and a negative declination there is spelled with a second hyphen.
    α = rad2deg(parse(Float64, raw[1, col("ra[rad]")]))
    δ = rad2deg(parse(Float64, raw[1, col("dec[rad]")]))
    ψ = parse.(Float64, raw[:, col("scanAngle[rad]")])
    p = parse.(Float64, raw[:, col("parallaxFactorAlongScan")])
    @test length(p) > 20

    e = parallax_ellipse(α, δ; n=20001)
    # Support function: the largest along-scan parallax factor geometrically
    # available at scan angle ψ. The ellipse is centred, so the lower bound is
    # its negative and |p| ≤ h(ψ) is the whole constraint.
    h(a) = maximum(e.raoff .* sin(a) .+ e.decoff .* cos(a))
    ratio = abs.(p) ./ h.(ψ)
    # 2% of slack for the one thing the circular orbit gets wrong: Earth's
    # eccentricity, 0.0167, which moves it 1.7% in and out over a year.
    @test maximum(ratio) < 1.02
    # …and the bound has to be *tight*, or an ellipse twice too big would pass.
    @test maximum(ratio) > 0.95

    # Teeth. Transposing the ellipse — the classic α/δ swap — keeps both
    # semi-axes and breaks only the orientation, and it must fail.
    ht(a) = maximum(e.decoff .* sin(a) .+ e.raoff .* cos(a))
    @test maximum(abs.(p) ./ ht.(ψ)) > 1.02
    # A source 90° away in RA has a differently-oriented ellipse; it need not
    # violate the bound, but it must stop being tight.
    e2 = parallax_ellipse(α + 90, δ; n=20001)
    h2(a) = maximum(e2.raoff .* sin(a) .+ e2.decoff .* cos(a))
    @test maximum(abs.(p) ./ h2.(ψ)) < 0.95
end
