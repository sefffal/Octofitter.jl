# Correction flags (`src/model/corrections.jl`) and the radial-velocity
# provenance declarations (`src/likelihoods/radial-velocity.jl`).
#
# What is under test is not the physics — PlanetOrbits' own suite gates that —
# but the decision procedure: that `:auto` measures the right quantity, ORs
# across observations, is reproducible from its recorded seed, and is recorded
# where a reader of the posterior will find it.

using TypedTables
using Octofitter: MCMCChains

const CORR_EPOCHS = collect(range(56000.0, 59000.0, length=12))

# An observation that does not implement `correction_impact` — i.e. every
# observation type today that has no `simulate`-shaped prediction. The `:auto`
# resolver must treat it as "cannot say" and keep the corrections on.
struct _MuteObs <: Octofitter.AbstractObs
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    name::String
end
_MuteObs() = _MuteObs(Octofitter.Priors(), Octofitter.Derived(), "mute")
Octofitter.epochs(::_MuteObs) = CORR_EPOCHS
Octofitter.refspecs(::_MuteObs) = ()
Octofitter.ln_like(::_MuteObs, ::Octofitter.ObsContext) = 0.0

# A two-body model. `plx` alone unless `absolute`, in which case a Barnard-like
# frame — the hardest case for every correction on this page.
function corr_model(; obs, og=:auto, blt=:auto, ein=:on, absolute=false,
                    verbosity=0, draws=60, epochs=CORR_EPOCHS)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 0.16
        flux = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 2mjup
        a = 1.0; e = 0.2; i = 0.9; ω = 1.1; Ω = 2.2; tp = 57000.0
    end)
    vars = absolute ? @variables(begin
        plx = 546.5
        ra = 45.0; dec = 20.0
        pmra = 7323.0; pmdec = 7323.0; rv = -110e3
        ref_epoch = 57388.0
    end) : @variables(begin
        plx = 25.0
    end)
    return Octofitter.System(name="corr", bodies=[A, b], observations=obs,
        variables=vars, observing_geometry=og, barycentric_lighttime=blt,
        einstein_rv=ein, correction_draws=draws, verbosity=verbosity)
end

# Modelled values of one observation inside a built system, in table order.
function corr_simulate(sys, obs)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep; observing_geometry=sys.observing_geometry,
                      barycentric_lighttime=sys.barycentric_lighttime)
    sym = Symbol(Octofitter.likelihoodname(obs))
    θ_obs = hasproperty(nt.observations, sym) ? getproperty(nt.observations, sym) : (;)
    # Through `_obsctx`, so the system's own solve configuration reaches the
    # likelihood exactly as it does under the sampler.
    return Octofitter.simulate(obs,
        Octofitter._obsctx(sys, nt, θ_obs, posys, traj, maps[obs]))
end

_decision(sys, flag) =
    only(filter(d -> d.flag === flag, sys.corrections.decisions))

# Synthetic relative astrometry of the companion at a chosen precision, over a
# chosen number of epochs on the same baseline.
function corr_astrom(σ; name="astrom", epochs=CORR_EPOCHS)
    n = length(epochs)
    return Octofitter.RelAstromObs(
        Table(; epoch=collect(epochs), ra=fill(100.0, n), dec=fill(50.0, n),
              σ_ra=fill(σ, n), σ_dec=fill(σ, n));
        target=:b, ref=:A, name)
end

function corr_rv(σ; kwargs...)
    n = length(CORR_EPOCHS)
    return Octofitter.RadialVelocityObs(
        Table(; epoch=CORR_EPOCHS, rv=zeros(n), σ_rv=fill(σ, n));
        target=:A, ref=Barycentre, name="rvs", kwargs...)
end

@testset "explicit settings bypass the test and say so" begin
    for (given, want) in ((:on, true), (:off, false), (true, true), (false, false))
        sys = corr_model(obs=[corr_astrom(1.0)], og=given, blt=given)
        @test sys.observing_geometry === want
        @test sys.barycentric_lighttime === want
        d = _decision(sys, :observing_geometry)
        @test d.source === :user
        @test d.resolved === want
        @test isempty(d.impacts)          # nothing was measured
        @test d.ndraws == 0
    end
    @test_throws r"takes `:auto`" corr_model(obs=[corr_astrom(1.0)], og=:sometimes)
end

@testset ":auto decides from the data's own uncertainties" begin
    # Coarse astrometry: the corrections are worth ~1e-2 µas against a 10 mas
    # error bar, so the bias they could introduce is nowhere near 0.1σ.
    coarse = corr_model(obs=[corr_astrom(10.0)])
    @test coarse.observing_geometry === false
    d = _decision(coarse, :observing_geometry)
    @test d.source === :auto
    @test d.ndraws > 0
    w = first(d.impacts)
    @test w.name == "astrom"
    @test 0 <= w.bias < Octofitter.CORRECTION_BIAS_THRESHOLD
    # The bias is the per-point impact accumulated over the data points, and
    # both are recorded so a reader can check the arithmetic.
    @test w.n == 2 * length(CORR_EPOCHS)          # ra and dec per epoch
    @test w.bias ≈ w.impact * sqrt(w.n)

    # The same system at GRAVITY precision keeps them.
    fine = corr_model(obs=[corr_astrom(1e-4)])
    @test fine.observing_geometry === true
    @test first(_decision(fine, :observing_geometry).impacts).bias >=
          Octofitter.CORRECTION_BIAS_THRESHOLD

    # ...and the flip is driven by σ alone: degrading the fine dataset by 1e4
    # is exactly the coarse one.
    @test corr_model(obs=[corr_astrom(1e-4 * 1e4)]).observing_geometry === false
end

@testset "the criterion is on accumulated bias, so it scales with N" begin
    # The point of the √N criterion, isolated: identical per-point precision,
    # identical baseline, identical orbit — only the number of measurements
    # differs. These corrections are common-mode, so their effect on a fitted
    # parameter grows as √N instead of averaging down, and a per-point margin
    # would have given both datasets the same (wrong for one of them) answer.
    σ = 0.15
    few = corr_model(obs=[corr_astrom(σ)], draws=20)
    manyep = range(first(CORR_EPOCHS), last(CORR_EPOCHS), length=1500)
    many = corr_model(obs=[corr_astrom(σ; epochs=manyep)], draws=20)
    @test few.observing_geometry === false
    @test many.observing_geometry === true
    fd = first(_decision(few, :observing_geometry).impacts)
    md = first(_decision(many, :observing_geometry).impacts)
    # Same per-point impact to within the coarser sampling of the extremum...
    @test md.impact ≈ fd.impact rtol = 0.2
    # ...and it is N alone that moved the verdict.
    @test md.n > 100 * fd.n
    @test fd.bias < Octofitter.CORRECTION_BIAS_THRESHOLD <= md.bias
end

@testset "needs are OR'd across observations" begin
    # One precise observation forces the correction on for the whole model,
    # however many coarse ones sit beside it.
    mixed = corr_model(obs=[corr_astrom(10.0; name="coarse"),
                            corr_astrom(1e-4; name="fine")])
    @test mixed.observing_geometry === true
    d = _decision(mixed, :observing_geometry)
    @test length(d.impacts) == 2
    # Worst-first, and it is the precise series that decided it.
    @test first(d.impacts).name == "fine"
    @test first(d.impacts).bias > last(d.impacts).bias
end

@testset "an observation that cannot report predictions keeps corrections on" begin
    # The conservative default: `:auto` can make a model cheaper than the old
    # always-on behaviour, never less accurate. Beside a capable observation
    # the test still runs — the capable one's numbers are worth reporting —
    # and the mute one vetoes the result.
    @test !Octofitter.has_correction_impact(_MuteObs())
    sys = corr_model(obs=[corr_astrom(10.0), _MuteObs()])
    @test sys.observing_geometry === true
    d = _decision(sys, :observing_geometry)
    @test isnan(first(d.impacts).bias)
    @test first(d.impacts).name == "mute"
    @test occursin("does not report comparable predictions", d.note)
end

@testset "a model where nothing can answer takes no draws at all" begin
    # The static capability check. Without it the most expensive models there
    # are (G23H, images, interferometry) would solve three hundred prior draws
    # to reach a conclusion that was fixed before the test started.
    sys = corr_model(obs=[_MuteObs()])
    @test sys.observing_geometry === true
    @test sys.barycentric_lighttime === true
    for flag in (:observing_geometry, :barycentric_lighttime)
        d = _decision(sys, flag)
        @test d.source === :auto
        @test d.ndraws == 0                # nothing was solved
        @test isempty(d.impacts)
        @test occursin("no observation can report comparable predictions", d.note)
        @test occursin("_MuteObs", d.note)  # ...and it names the type
    end
end

@testset "a model with no observations resolves accurate" begin
    sys = corr_model(obs=[])
    @test sys.observing_geometry === true
    @test sys.barycentric_lighttime === true
    @test occursin("no observations", _decision(sys, :observing_geometry).note)
end

@testset "the decision is reproducible from its recorded seed" begin
    a = corr_model(obs=[corr_astrom(1.0)])
    b = corr_model(obs=[corr_astrom(1.0)])
    for flag in (:observing_geometry, :barycentric_lighttime)
        da, db = _decision(a, flag), _decision(b, flag)
        @test da.resolved === db.resolved
        @test da.seed === db.seed
        @test da.impacts == db.impacts
    end
    # A different seed is allowed to differ, but must still be recorded.
    c = Octofitter.System(name="corr", bodies=corr_model(obs=[]).nodes,
        observations=[corr_astrom(1.0)],
        variables=@variables(begin plx = 25.0 end),
        correction_seed=UInt64(12345), correction_draws=60, verbosity=0)
    @test _decision(c, :observing_geometry).seed === UInt64(12345)
end

@testset "the report travels with the model" begin
    sys = corr_model(obs=[corr_astrom(10.0)])
    @test sys.corrections isa Octofitter.CorrectionReport
    @test length(sys.corrections.decisions) == 2
    # It prints, because a reader of a saved chain has only this. The
    # one-line form is what `savechain` puts in a FITS header card, so it has
    # to stay short and ASCII.
    s = sprint(show, MIME"text/plain"(), sys.corrections)
    @test occursin("observing_geometry", s)
    @test occursin("barycentric_lighttime", s)
    one_line = sprint(show, sys.corrections)
    @test !occursin('\n', one_line)
    @test all(isascii, one_line)
    @test occursin("observing_geometry=false(auto)", one_line)
    # A FITS header card silently drops a string value longer than this.
    @test length(one_line) <= 68
end

@testset "the verdict reaches the chain, and a saved one" begin
    sys = corr_model(obs=[corr_astrom(10.0)])
    chain = Octofitter.result2mcmcchain(
        [(; plx=25.0, bodies=(; A=(; mass=0.16, flux=1.0),
                              b=(; mass=2mjup, a=1.0, e=0.2, i=0.9, ω=1.1,
                                 Ω=2.2, tp=57000.0)), observations=(;))])
    chain = MCMCChains.setinfo(chain, (;
        model_name=sys.name,
        corrections=sys.corrections,
        observing_geometry=sys.observing_geometry,
        barycentric_lighttime=sys.barycentric_lighttime))
    # In session, the whole report is there.
    @test chain.info.corrections === sys.corrections
    mktempdir() do dir
        f = joinpath(dir, "c.fits")
        Octofitter.savechain(f, chain)
        back = Octofitter.loadchain(f)
        # The resolved settings round-trip as values, which is the fact a
        # reader of an old chain needs: this posterior was produced with these
        # corrections. The report object itself does not survive a FITS header
        # card (key plus quoted value plus comment exceeds 80 characters), so
        # it is deliberately not asserted here — see the docs.
        @test back.info.observing_geometry == false
        @test back.info.barycentric_lighttime == false
    end
end

# ---------------------------------------------------
# Radial velocity: perspective acceleration and the Einstein hook
# ---------------------------------------------------

@testset "secular_acceleration: what it adds" begin
    # Barnard-like: v_t²/d ≈ 4.5 m/s/yr, so ~45 m/s over the ~8.2 yr baseline.
    modelled = corr_model(obs=[corr_rv(1.0)], absolute=true, og=:on, blt=:on)
    corrected = corr_model(obs=[corr_rv(1.0; secular_acceleration=:data_corrected)],
                           absolute=true, og=:on, blt=:on)
    m = corr_simulate(modelled, only(modelled.observations)).rv_model
    c = corr_simulate(corrected, only(corrected.observations)).rv_model
    drift = m .- c
    yrs = (CORR_EPOCHS[end] - CORR_EPOCHS[1]) / PlanetOrbits.year2day_julian
    @test (maximum(drift) - minimum(drift)) / yrs ≈ 4.5 rtol = 0.15
    # The sign is not hand-derived: it is whatever `_compensate_kinematics`
    # says the propagated frame RV does, differenced against `ref_epoch`.
    fr = Octofitter.make_ln_like(modelled).build(
        Octofitter.make_arr2nt(modelled)(Float64[])).frame
    direct = [PlanetOrbits._compensate_kinematics(fr, t).rv2 - fr.rv
              for t in CORR_EPOCHS]
    @test drift ≈ direct rtol = 1e-3

    # `:data_corrected` is the old behaviour exactly: `radvel` and nothing else.
    posys = Octofitter.make_ln_like(corrected).build(
        Octofitter.make_arr2nt(corrected)(Float64[]))
    traj = orbitsolve(posys, CORR_EPOCHS)
    @test c ≈ [radvel(traj[k], PlanetOrbits.BodyRef(1),
                      PlanetOrbits.barycentre(posys)) for k in eachindex(CORR_EPOCHS)]
end

@testset "secular_acceleration is zero without an absolute frame" begin
    # The frame level chooses the physics; this is not an error.
    a = corr_model(obs=[corr_rv(1.0)], og=:on, blt=:on)
    b = corr_model(obs=[corr_rv(1.0; secular_acceleration=:data_corrected)],
                   og=:on, blt=:on)
    @test corr_simulate(a, only(a.observations)).rv_model ≈
          corr_simulate(b, only(b.observations)).rv_model
end

@testset "secular_acceleration is rejected for a relative series" begin
    n = length(CORR_EPOCHS)
    tab = Table(; epoch=CORR_EPOCHS, rv=zeros(n), σ_rv=fill(1.0, n))
    mk(; kw...) = Octofitter.RadialVelocityObs(tab; target=:b, ref=:A,
                                              name="rel", kw...)
    # Silent by default — the term genuinely does not apply.
    @test mk().secular_acceleration === :not_applicable
    err = try
        mk(secular_acceleration=:model)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    msg = sprint(showerror, err)
    # Both halves of the argument have to be there: the blanket "it cancels"
    # is false, and the half that does not cancel is the large one.
    @test occursin("common-mode", msg)
    @test occursin("differential", msg)
    @test occursin("observing_geometry", msg)
    @test occursin("0.023", msg)
    # An absolute series still validates the value itself.
    @test_throws r"takes `:model`" corr_rv(1.0; secular_acceleration=:maybe)
end

@testset "sampled proper motion moves the RV prediction" begin
    # The point of `:model`: a joint fit's sampled (μ, ϖ) predict both the
    # astrometry and the RV drift, self-consistently.
    function pm_model(pmra)
        A = Octofitter.Body(name="A", variables=@variables begin mass = 0.16 end)
        b = Octofitter.Body(name="b", about=A, variables=@variables begin
            mass = 2mjup
            a = 1.0; e = 0.2; i = 0.9; ω = 1.1; Ω = 2.2; tp = 57000.0
        end)
        Octofitter.System(name="pm", bodies=[A, b], observations=[corr_rv(1.0)],
            variables=@variables(begin
                plx = 546.5
                ra = 45.0; dec = 20.0
                pmra = $pmra; pmdec = 7323.0; rv = -110e3
                ref_epoch = 57388.0
            end), observing_geometry=:on, barycentric_lighttime=:on, verbosity=0)
    end
    lo = pm_model(7323.0)
    hi = pm_model(9323.0)
    a = corr_simulate(lo, only(lo.observations)).rv_model
    b = corr_simulate(hi, only(hi.observations)).rv_model
    @test maximum(abs, a .- b) > 5.0     # m/s
end

@testset "einstein_rv=:off predicts from the kinematic velocity" begin
    # A whole-model switch, not a per-series one: there is no provenance
    # dimension for it, so it lives on `System` beside the other solve
    # settings and reaches the likelihood through the context.
    on = corr_model(obs=[corr_rv(1.0)], absolute=true, og=:on, blt=:on)
    off = corr_model(obs=[corr_rv(1.0)], absolute=true, og=:on, blt=:on, ein=:off)
    @test on.einstein_rv === true
    @test off.einstein_rv === false
    @test_throws r"not `:auto`" corr_model(obs=[corr_rv(1.0)], ein=:auto)
    a = corr_simulate(on, only(on.observations)).rv_model
    b = corr_simulate(off, only(off.observations)).rv_model
    posys = Octofitter.make_ln_like(off).build(
        Octofitter.make_arr2nt(off)(Float64[]))
    traj = orbitsolve(posys, CORR_EPOCHS)
    bary = PlanetOrbits.barycentre(posys)
    k = PlanetOrbits.au2m / PlanetOrbits.year2sec_julian
    drift = [PlanetOrbits.frame_rv(traj[j]) - posys.frame.rv for j in eachindex(CORR_EPOCHS)]
    @test b ≈ [velz(traj[j], PlanetOrbits.BodyRef(1), bary) * k
               for j in eachindex(CORR_EPOCHS)] .+ drift
    @test !(a ≈ b)
end

@testset "einstein_rv reaches every path that re-evaluates a likelihood" begin
    # The trap this guards: `ObsContext`'s convenience constructors default the
    # flag to `true`, so a path that builds one by hand would score a fit
    # against different physics than it sampled. Everything that evaluates a
    # likelihood outside the generated function goes through `_obsctx`.
    off = corr_model(obs=[corr_rv(1.0)], absolute=true, og=:on, blt=:on, ein=:off)
    nt = Octofitter.make_arr2nt(off)(Float64[])
    posys = Octofitter.make_ln_like(off).build(nt)
    ep, maps = Octofitter.epoch_plan(off)
    traj = orbitsolve(posys, ep)
    obs = only(off.observations)
    ctx = Octofitter._obsctx(off, nt, (;), posys, traj, maps[obs])
    @test ctx.einstein_rv === false
    # ...and the generated function agrees with it, which is the property that
    # actually matters: the sampler and every re-read see one physics.
    lnlike = Octofitter.make_ln_like(off)
    @test isfinite(lnlike(off, nt))
    @test Octofitter.ln_like(obs, ctx) ≈ lnlike(off, nt)
end

@testset "the advisories report the Einstein term and the drift" begin
    sys = corr_model(obs=[corr_rv(1.0)], absolute=true, draws=20)
    adv = sys.corrections.advisories
    @test any(a -> occursin("Einstein term", a), adv)
    @test any(a -> occursin("secular-acceleration drift", a), adv)
    # ...and not the drift line for a series that declared it away.
    sys2 = corr_model(obs=[corr_rv(1.0; secular_acceleration=:data_corrected)],
                      absolute=true, draws=20)
    @test !any(a -> occursin("secular-acceleration drift", a), sys2.corrections.advisories)
end

@testset "post-sampling re-check reports both directions" begin
    # A model whose priors are broad enough to force the corrections on, but
    # whose "posterior" sits entirely where they do not matter.
    A = Octofitter.Body(name="A", variables=@variables begin mass = 0.16 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 2mjup
        a ~ Uniform(0.5, 400.0)     # 400 AU at 25 mas is 10" on sky
        e = 0.2; i = 0.9; ω = 1.1; Ω = 2.2; tp = 57000.0
    end)
    # 5 mas. Coarse enough that the tight "posterior" below clears the 10³
    # margin, but the prior's wide tail (ρ² depth scaling at 400 AU is ~470 µas)
    # does not.
    obs = corr_astrom(5.0)
    sys = Octofitter.System(name="recheck", bodies=[A, b], observations=[obs],
        variables=@variables(begin plx = 25.0 end),
        correction_draws=60, verbosity=0)
    @test sys.observing_geometry === true      # the wide tail forces it on
    model = Octofitter.LogDensityModel(sys)
    # A "posterior" concentrated at small separations, where the corrections
    # are far below 50 µas.
    chain = Octofitter.result2mcmcchain(
        [(; plx=25.0, bodies=(; A=(; mass=0.16),
                              b=(; mass=2mjup, a, e=0.2, i=0.9, ω=1.1, Ω=2.2,
                                 tp=57000.0)), observations=(;))
         for a in range(0.5, 0.6, length=40)])
    r = @test_logs (:info,) match_mode = :any recheck_corrections(model, chain)
    d = only(filter(x -> x.flag === :observing_geometry, r.decisions))
    @test d.resolved === false
    @test occursin("de-escalation", d.note)
end
