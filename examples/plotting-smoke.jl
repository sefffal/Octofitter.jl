# ---------------------------------------------------------------------------
# Makie plotting smoke test.
#
# `test/v2/plotting.jl` covers the backend-agnostic protocol — channels,
# residuals, queries, PosteriorSeries — because that is what has to agree with
# the likelihoods. This script covers the *drawing*, which CI does not: it
# builds every panel the extension can produce and writes the figures out, so
# a layout regression is one `julia examples/plotting-smoke.jl` away rather
# than one user report away.
#
#     julia --project=. examples/plotting-smoke.jl [outdir]
#
# The "chains" here are prior draws, not fits, so the curves miss the data by
# design; what is being checked is that every panel is built, labelled, scaled
# and oriented correctly.
# ---------------------------------------------------------------------------

using Octofitter, OctofitterRadialVelocity, Distributions, CairoMakie
using Random: Xoshiro
using PlanetOrbits
# Through OctofitterRadialVelocity rather than directly: AbstractGPs is its
# dependency, not Octofitter's, so this runs under `--project=.` too.
using OctofitterRadialVelocity: AbstractGPs
using .AbstractGPs: GP, Matern52Kernel, ScaleTransform

const OUT = get(ARGS, 1, joinpath(@__DIR__, "plotting-smoke-figs"))
mkpath(OUT)

const MExt = Base.get_extension(Octofitter, :OctofitterMakieExt)

step(msg) = printstyled("\n=== ", msg, " ===\n"; color=:cyan)
ok(msg) = printstyled("  ok  ", msg, "\n"; color=:green)
bad(msg) = printstyled("  !!  ", msg, "\n"; color=:red)
check(cond, msg) = cond ? ok(msg) : bad(msg)

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

# sep/pa astrometry + ra/dec astrometry + two RV instruments with very
# different zero points — the case that shows why RV instruments cannot share
# a panel while relative astrometry can.
function mixed_model()
    A = Octofitter.Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.5)
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 0.005
        a ~ Uniform(2.9, 3.1)
        e = 0.25
        i = 0.6
        ω = 1.0
        Ω = 2.0
        tp = 59200.0
    end)
    seppa = RelAstromObs(
        Table(epoch=[59000.0, 59300.0, 59600.0],
            sep=[150.0, 158.0, 149.0], pa=[0.5, 0.7, 0.9],
            σ_sep=[5.0, 5.0, 5.0], σ_pa=[0.01, 0.01, 0.01]);
        target=b, ref=A, name="SPHERE",
        variables=@variables begin
            jitter = 1.0
            platescale = 1.01
            northangle = 0.02
        end)
    radec = RelAstromObs(
        Table(epoch=[59100.0, 59400.0],
            ra=[80.0, 95.0], dec=[130.0, 125.0],
            σ_ra=[4.0, 4.0], σ_dec=[4.0, 4.0]);
        target=b, ref=A, name="GPI")
    rv1 = RadialVelocityObs(
        Table(epoch=[59050.0, 59250.0, 59450.0, 59650.0],
            rv=[1030.0, 990.0, 1015.0, 975.0], σ_rv=[3.0, 3.0, 3.0, 3.0]);
        target=A, ref=Barycentre, name="HARPS",
        variables=@variables begin
            offset = 1000.0
            jitter = 2.0
        end)
    rv2 = RadialVelocityObs(
        Table(epoch=[59120.0, 59320.0, 59520.0],
            rv=[-500.0, -530.0, -505.0], σ_rv=[4.0, 4.0, 4.0]);
        target=A, ref=Barycentre, name="HIRES",
        variables=@variables begin
            offset = -510.0
            jitter = 3.0
        end)
    sys = Octofitter.System(name="mixed", bodies=[A, b],
        observations=[seppa, radec, rv1, rv2],
        variables=@variables begin
            plx ~ truncated(Normal(50.0, 0.5), lower=1)
        end)
    return Octofitter.LogDensityModel(sys, verbosity=0), (; seppa, radec, rv1, rv2)
end

# Relative astrometry only: the "predict the RV curve before asking for
# spectroscopic time" case.
function astrom_only_model()
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass ~ Uniform(0.001, 0.02)
        a ~ Uniform(2.9, 3.1)
        e = 0.2; i = 0.6; ω = 1.0; Ω = 2.0; tp = 59200.0
    end)
    seppa = RelAstromObs(
        Table(epoch=[59000.0, 59300.0, 59600.0],
            sep=[150.0, 158.0, 149.0], pa=[0.5, 0.7, 0.9],
            σ_sep=[5.0, 5.0, 5.0], σ_pa=[0.01, 0.01, 0.01]);
        target=b, ref=A, name="SPHERE")
    radec = RelAstromObs(
        Table(epoch=[59100.0, 59400.0], ra=[80.0, 95.0], dec=[130.0, 125.0],
            σ_ra=[4.0, 4.0], σ_dec=[4.0, 4.0]);
        target=b, ref=A, name="GPI")
    sys = Octofitter.System(name="astromonly", bodies=[A, b],
        observations=[seppa, radec],
        variables=@variables begin plx = 50.0 end)
    return Octofitter.LogDensityModel(sys, verbosity=0)
end

# An RV fit with a correlated-noise model, for the `noisemodel` band and the
# GP-corrected residual strip.
function gp_model()
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass ~ Uniform(0.001, 0.01)
        a ~ Uniform(0.9, 1.1)
        e = 0.05; i = 1.2; ω = 0.4; Ω = 1.0; tp = 59000.0
    end)
    ep = collect(range(59000.0, 59400.0, length=24))
    rv = RadialVelocityObs(
        Table(epoch=ep,
            rv=30 .* sin.(ep ./ 40) .+ 12 .* sin.(ep ./ 9) .+ 100,
            σ_rv=fill(4.0, length(ep)));
        target=A, ref=Barycentre, name="HARPS-N",
        gaussian_process=θ_obs -> GP(θ_obs.gp_η₁^2 *
                                     Matern52Kernel() ∘ ScaleTransform(1 / θ_obs.gp_η₂)),
        variables=@variables begin
            offset = 100.0
            jitter = 1.0
            gp_η₁ = 25.0
            gp_η₂ = 30.0
        end)
    sys = Octofitter.System(name="gpfit", bodies=[A, b], observations=[rv],
        variables=@variables begin plx = 25.0 end)
    return Octofitter.LogDensityModel(sys, verbosity=0), rv
end

function fakechain(model, n; seed=3)
    rng = Xoshiro(seed)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:n]
    return Octofitter.result2mcmcchain(nts)
end

# ---------------------------------------------------------------------------

model, obs = mixed_model()
chain = fakechain(model, 60)
series = Octofitter.PosteriorSeries(model, chain; N=12)

step("relative astrometry declares all four channels, natives first")
for o in (obs.seppa, obs.radec)
    chs = Octofitter.plotchannels(o)
    println("  ", rpad(likelihoodname(o), 8), [(ch.name, ch.derived) for ch in chs])
end
check([ch.derived for ch in Octofitter.plotchannels(obs.seppa)] == [false, false, true, true],
    "sep/pa table: sep, pa native; ra/dec derived")
check([ch.derived for ch in Octofitter.plotchannels(obs.radec)] == [true, true, false, false],
    "ra/dec table: ra/dec native; sep/pa derived")

step("the projection is exact in both directions")
r = Octofitter.residuals(obs.seppa, Octofitter.obscontext(series, obs.seppa))
check(r.sep.data ≈ hypot.(r.raoff.data, r.decoff.data), "sep/pa -> ra/dec round-trips")
r2 = Octofitter.residuals(obs.radec, Octofitter.obscontext(series, obs.radec))
check(r2.raoff.data ≈ r2.sep.data .* sin.(deg2rad.(r2.pa.data)), "ra/dec -> sep/pa round-trips")
check(keys(r) == (:sep, :pa, :raoff, :decoff), "all four channels in `residuals`")

step("panel policy")
check(Octofitter.sharepanel(obs.seppa), "relative astrometry shares a panel")
check(!Octofitter.sharepanel(obs.rv1), "radial velocity does not")
for (k, entries) in MExt._channelgroups(model.system)
    println("  ", rpad(k, 34), " <- ", [likelihoodname(o) for (o, _) in entries])
end

step("octoplot: the default figure")
res = octoplot(model, chain; N=40, fname=joinpath(OUT, "01-octoplot.png"))
println("  axes: ", keys(res.axes))
check(:rv_HARPS in keys(res.axes) && :rv_HIRES in keys(res.axes),
    "one panel per RV instrument, named for it")
check(count(k -> k in (:sep, :pa, :raoff, :decoff), collect(keys(res.axes))) == 4,
    "one shared panel per relative-astrometry channel")
check(res.axes.sky.sky.xreversed[], "sky panel: RA increases to the left")

step("raw residuals need a single draw")
try
    octoplot(model, chain; N=40, whiten=false)
    bad("many-draw whiten=false was allowed")
catch e
    ok("refused: " * first(split(sprint(showerror, e), "\n")))
end
octoplot(model, chain; N=40, ndraws=1, whiten=false,
    fname=joinpath(OUT, "02-octoplot-one-draw.png"))
ok("single draw with raw residuals draws")

step("rvplot: every instrument on one axis, one draw")
res_rv = rvplot(model, chain; fname=joinpath(OUT, "03-rvplot.png"))
println("  axes: ", keys(res_rv.axes))
check(:rv in keys(res_rv.axes) && :rv_phase_b in keys(res_rv.axes),
    "time panel + phase panel")
try
    rvplot(Octofitter.PosteriorSeries(model, chain; N=10))
    bad("a multi-draw series was accepted")
catch e
    ok("multi-draw series refused")
end
rvpostplot(model, chain; fname=joinpath(OUT, "04-rvpostplot-alias.png"))
ok("the deprecated `rvpostplot` alias still works")

step("channels=: every instrument on the requested panel")
res_sep = octoplot(model, chain; N=40, channels=projectedseparation, show_sky=false,
    fname=joinpath(OUT, "05-separation.png"))
check(keys(res_sep.axes) == (:sep,), "one separation panel, carrying both instruments")

step("show_phase=true folds channels other than RV")
res_pa = octoplot(model, chain; N=40, channels=posangle, show_sky=false, show_phase=true,
    fname=joinpath(OUT, "06-posangle-phase.png"))
check(:pa_phase_b in keys(res_pa.axes), "phase-folded position angle")

step("predicting a quantity the fit has no data for")
m2 = astrom_only_model()
c2 = fakechain(m2, 40; seed=7)
res_pred = octoplot(m2, c2; N=30, channels=radvel, show_sky=false,
    fname=joinpath(OUT, "07-predicted-rv.png"))
check(:radvel in keys(res_pred.axes), "an RV panel from an astrometry-only fit")
octoplot(m2, c2; N=30, fname=joinpath(OUT, "08-astrometry-only.png"))
ok("astrometry-only default figure draws")

step("correlated noise: band, and residuals corrected for it")
m3, rvgp = gp_model()
c3 = fakechain(m3, 30; seed=11)
s3 = Octofitter.PosteriorSeries(m3, c3; N=6)
rgp = Octofitter.residuals(rvgp, Octofitter.obscontext(s3, rvgp)).rv
check(haskey(rgp, :gp_mean) && haskey(rgp, :gp_var), "`residuals` publishes gp_mean/gp_var")
check(Octofitter.noisemodel(rvgp, Octofitter.obscontext(s3, rvgp), [59100.0, 59200.0]) !== nothing,
    "`noisemodel` predicts off the data epochs (AbstractGPs backend)")
rvplot(m3, c3; fname=joinpath(OUT, "09-rvplot-gp.png"))
octoplot(m3, c3; N=20, ndraws=6, fname=joinpath(OUT, "10-octoplot-gp.png"))
ok("GP figures draw")

step("animation")
mktempdir() do d
    f = Octofitter.rvplot_animated(model, chain; N=3, fname=joinpath(d, "rv.gif"))
    check(isfile(f), "rvplot_animated writes a file")
end

printstyled("\nFigures written to $OUT\n"; color=:cyan, bold=true)
