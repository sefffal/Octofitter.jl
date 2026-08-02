# The backend-agnostic plotting API: queries, the observation plotting
# protocol, and PosteriorSeries. No Makie here — the extension is exercised
# by the smoke script in examples/ and by users; what matters in CI is that
# the protocol agrees with the likelihoods.

using Octofitter
using Distributions
using Random: Xoshiro

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

    chs = plotchannels(obs.seppa)
    @test length(chs) == 2
    @test chs[1].name === :sep && chs[1].query isa ObservableQuery
    @test chs[2].name === :pa && chs[2].wrap == 360.0 && chs[2].scale ≈ rad2deg(1.0)
    @test plotchannels(obs.radec)[1].name === :raoff
    @test plotchannels(obs.rvs)[1].name === :rv

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
    @test haskey(r2, :sep) && haskey(r2, :pa)
    @test r2.sep.data ≈ obs.seppa.table.sep .* 1.01          # platescale applied
    @test all(abs.(r2.pa.resid) .<= 180.0)                    # wrapped, degrees
    @test r2.pa.σ[1] ≈ rad2deg(0.01)

    # ra/dec: no nuisance variables → data pass through
    ctx3 = obscontext(series, obs.radec)
    r3 = Octofitter.residuals(obs.radec, ctx3)
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

@testset "octoplot without Makie errors helpfully" begin
    # (Makie is not loaded in the test env, so the stub must fire.)
    if !haskey(Base.loaded_modules_array() |> ms -> Dict(nameof(m) => m for m in ms), :Makie)
        model, _ = _plotting_test_model()
        @test_throws ErrorException octoplot(model)
    end
end
