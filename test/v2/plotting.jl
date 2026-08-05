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
    full = Octofitter.evalquery(q, posys, traj)
    parts = Octofitter.evalsignal(s1, posys, traj) .+ Octofitter.evalsignal(s2, posys, traj)
    @test full ≈ parts rtol = 1e-10

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
