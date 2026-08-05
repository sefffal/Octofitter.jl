# Observations: reference grammar, epoch planning, and the likelihoods
# themselves.

# A two-companion system used throughout, built from a known PlanetOrbits
# system so the "data" are exact.
function reference_system(; plx=25.0)
    A = PlanetOrbits.Body(mass=1.0 - 6mjup, flux=(; default=1.0), name=:A)
    b = PlanetOrbits.Body(mass=2mjup, name=:b)
    c = PlanetOrbits.Body(mass=4mjup, name=:c)
    return PlanetOrbits.System((A, b, c), (
            PlanetOrbits.Orbit(b, about=A; a=3.0, e=0.15, i=0.9, ω=1.2, Ω=2.0, tp=55000.0),
            PlanetOrbits.Orbit(c, about=A; a=9.0, e=0.05, i=1.1, ω=0.3, Ω=1.0, tp=56000.0),
        ); plx)
end

const EPOCHS_A = collect(range(56000.0, 59000.0, length=9))
const EPOCHS_B = collect(range(56500.0, 58500.0, length=7))

function model_system(; obs)
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

@testset "reference grammar" begin
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    @test Octofitter.refspec(A) === Octofitter.BodyRefSpec{:A}()
    @test Octofitter.refspec(:A) === Octofitter.BodyRefSpec{:A}()
    @test Barycentre === Octofitter.BarycentreSpec{()}()
    @test Barycentre(A, :b) === Octofitter.BarycentreSpec{(:A, :b)}()
    @test Photocentre === Octofitter.PhotocentreSpec{nothing,()}()
    @test Photocentre(:G) === Octofitter.PhotocentreSpec{:G,()}()
    @test_throws r"must be a `Body`" Octofitter.refspec(3.0)

    # Resolution against a solved system produces PlanetOrbits references.
    posys = reference_system()
    @test Octofitter.resolveref(posys, Octofitter.refspec(:b)) === PlanetOrbits.BodyRef(2)
    wp = Octofitter.resolveref(posys, Barycentre)
    @test wp isa PlanetOrbits.WeightedPoint
    @test sum(wp.w) ≈ 1
end

@testset "photocentre subset grammar" begin
    Aa = Octofitter.Body(name="Aa", variables=@variables begin mass = 1.0 end)
    Ab = Octofitter.Body(name="Ab", variables=@variables begin mass = 0.5 end)

    # All three spellings, plus the trailing-argument form.
    @test Photocentre(:G, (Aa, Ab)) === Octofitter.PhotocentreSpec{:G,(:Aa, :Ab)}()
    @test Photocentre(:G, Aa, Ab) === Octofitter.PhotocentreSpec{:G,(:Aa, :Ab)}()
    @test Photocentre(:G, (:Aa, :Ab)) === Octofitter.PhotocentreSpec{:G,(:Aa, :Ab)}()
    @test Photocentre((Aa, Ab)) === Octofitter.PhotocentreSpec{nothing,(:Aa, :Ab)}()
    @test Photocentre([Aa, Ab]) === Octofitter.PhotocentreSpec{nothing,(:Aa, :Ab)}()
    @test Photocentre(:G, (Aa,)) === Octofitter.PhotocentreSpec{:G,(:Aa,)}()
    @test_throws r"empty member list" Photocentre(:G, ())
    @test_throws r"must be a `Body` model node" Photocentre(:G, (3.0,))

    # Every spec is a singleton — the content lives in type parameters, which
    # is what makes `resolveref` constant-fold.
    @test isbits(Photocentre(:G, (Aa, Ab)))
    @test sizeof(Photocentre(:G, (Aa, Ab))) == 0

    # `_refstr` round-trips the spelling; `_refnames` is what build-time
    # validation checks against the body list.
    @test Octofitter._refstr(Photocentre) == "Photocentre"
    @test Octofitter._refstr(Photocentre(:G)) == "Photocentre(:G)"
    @test Octofitter._refstr(Photocentre(:G, (Aa, Ab))) == "Photocentre(:G, (Aa, Ab))"
    @test Octofitter._refstr(Photocentre((Aa, Ab))) == "Photocentre((Aa, Ab))"
    @test Octofitter._refstr(Photocentre(:G, (Aa,))) == "Photocentre(:G, (Aa,))"
    @test Octofitter._refnames(Photocentre(:G, (Aa, Ab))) === (:Aa, :Ab)
    @test Octofitter._refnames(Photocentre(:G)) === ()
    @test Octofitter._refbands(Photocentre(:G, (Aa, Ab))) === (:G,)
    @test Octofitter._refbands(Photocentre) === ()
    @test Octofitter._refbands(Barycentre) === ()

    # Resolution forwards to PlanetOrbits' subset method, and agrees with it.
    posys = PlanetOrbits.System(
        (PlanetOrbits.Body(mass=1.0, flux=(G=1.0,), name=:A),
         PlanetOrbits.Body(mass=0.5, flux=(G=0.25,), name=:b),
         PlanetOrbits.Body(mass=0.2, flux=(G=4.0,), name=:c)),
        (PlanetOrbits.Orbit(PlanetOrbits.Body(mass=0.5, name=:b),
             about=PlanetOrbits.Body(mass=1.0, name=:A); a=1.0, tp=0.0),
         PlanetOrbits.Orbit(PlanetOrbits.Body(mass=0.2, name=:c),
             about=(PlanetOrbits.Body(mass=1.0, name=:A),
                    PlanetOrbits.Body(mass=0.5, name=:b)); a=9.0, tp=0.0)); plx=25.0)
    sub = Octofitter.resolveref(posys, Photocentre(:G, (:A, :b)))
    @test sub isa PlanetOrbits.WeightedPoint
    @test sub.w ≈ PlanetOrbits.photocentre(posys, :A, :b; band=:G).w
    @test sub.w ≈ [1.0 / 1.25, 0.25 / 1.25, 0.0]
    @test Octofitter.resolveref(posys, Photocentre(:G)).w ≈
          PlanetOrbits.photocentre(posys; band=:G).w
    # …and it constant-folds: the spec carries no runtime data at all.
    #
    # Version-gated deliberately, not to make CI green. On 1.11+ this folds to
    # nothing. On 1.10 the same call allocates 64 bytes — the subset weights
    # escape rather than being built in place — so the assertion is kept as a
    # `@test_broken`, which stays visible in the summary and fails loudly if a
    # future 1.10 patch starts folding it. It is a real (if small) cost: this
    # runs once per likelihood evaluation, which is exactly the allocation
    # `resolverefs` exists to avoid. See also the open question about whether
    # Julia 1.10 remains supported at all, given Pigeons cannot be installed
    # there.
    @test (@allocated Octofitter.resolveref(posys, Photocentre(:G, (:A, :b)))) == 0 broken=(VERSION < v"1.11")
end

@testset "photocentre member and band typos are caught at model-build time" begin
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
        flux_G = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 2mjup
        flux_G = 0.1
        a = 3.0; e = 0.1; i = 0.5; ω = 0.0; Ω = 0.0; tp = 55000.0
    end)
    mk(target) = RelAstromObs((epoch=[57000.0], ra=[1.0], dec=[1.0],
        σ_ra=[1.0], σ_dec=[1.0]); target, ref=A, name="astrom")
    build(target) = Octofitter.System(name="s", bodies=[A, b], observations=[mk(target)],
        variables=@variables begin plx = 20.0 end)

    @test build(Photocentre(:G, (A, b))) isa Octofitter.System   # the good case
    @test_throws r"references :zzz" build(Photocentre(:G, (A, :zzz)))
    @test_throws r"asks for the :Ks photocentre" build(Photocentre(:Ks, (A, b)))
    @test_throws r"the bands declared here are G" build(Photocentre(:Ks))
end

@testset "an observation naming a body the system lacks fails at build time" begin
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin a = 3.0 end)
    obs = RelAstromObs((epoch=[57000.0], ra=[1.0], dec=[1.0], σ_ra=[1.0], σ_dec=[1.0]);
        target=:zzz, ref=A, name="astrom")
    @test_throws r"references :zzz" Octofitter.System(
        name="s", bodies=[A, b], observations=[obs],
        variables=@variables begin plx = 20.0 end)
end

@testset "relative astrometry recovers exactly the offsets it was built from" begin
    posys = reference_system()
    tr = orbitsolve(posys, EPOCHS_A)
    ra = [raoff(tr[k], :c, :b) for k in eachindex(EPOCHS_A)]
    dec = [decoff(tr[k], :c, :b) for k in eachindex(EPOCHS_A)]

    # `target=c, ref=b` — a companion measured against *another companion*,
    # which v1's per-planet solutions could not express at all.
    obs = RelAstromObs((epoch=EPOCHS_A, ra=ra, dec=dec,
            σ_ra=fill(1.0, length(ra)), σ_dec=fill(1.0, length(ra)));
        target=:c, ref=:b, name="astrom")
    sys = model_system(obs=[obs])
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    ll = Octofitter.make_ln_like(sys)(sys, nt)
    # Zero residuals: the likelihood is exactly the normalization constant.
    @test ll ≈ length(ra) * logpdf(MvNormal(Diagonal(@SVector[1.0, 1.0])), @SVector[0.0, 0.0]) rtol = 1e-12
end

@testset "radial velocity: reflex against the barycentre" begin
    posys = reference_system()
    tr = orbitsolve(posys, EPOCHS_A)
    rv = [radvel(tr[k], :A, PlanetOrbits.barycentre(posys)) for k in eachindex(EPOCHS_A)]
    obs = RadialVelocityObs((epoch=EPOCHS_A, rv=rv, σ_rv=fill(1.0, length(rv)));
        target=:A, ref=Barycentre, name="rv")
    sys = model_system(obs=[obs])
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    ll = Octofitter.make_ln_like(sys)(sys, nt)
    @test ll ≈ -length(rv) * log(2π) / 2 rtol = 1e-12
end

@testset "epochs are pooled, sorted and deduplicated across observations" begin
    shared = 57000.0
    o1 = RelAstromObs((epoch=[shared, 58000.0], ra=[1.0, 2.0], dec=[1.0, 2.0],
            σ_ra=[1.0, 1.0], σ_dec=[1.0, 1.0]); target=:b, ref=:A, name="a1")
    o2 = RadialVelocityObs((epoch=[56000.0, shared], rv=[1.0, 2.0], σ_rv=[1.0, 1.0]);
        target=:A, ref=Barycentre, name="rv1")
    sys = model_system(obs=[o1, o2])
    ep, maps = Octofitter.epoch_plan(sys)
    @test ep == [56000.0, 57000.0, 58000.0]     # sorted union, shared epoch once
    @test maps[o1] == [2, 3]
    @test maps[o2] == [1, 2]
end

@testset "Gaia DR4 along-scan" begin
    scans, _ = let
        rows, _ = readdlm(joinpath(@__DIR__, "fixtures", "dr4-like-scans.csv"),
            ',', Float64, '\n'; header=true)
        Table(epoch=rows[:, 1], scan_pos_angle=rows[:, 2], parallax_factor_al=rows[:, 3],
            centroid_pos_al=rows[:, 4], centroid_pos_error_al=rows[:, 5],
            outlier_flag=Int.(rows[:, 6])), nothing
    end

    mkobs(; target) = GaiaDR4AstromObs(scans; target, ref=Barycentre, name="GaiaDR4",
        variables=@variables begin
            astrometric_jitter = 0.05
            ra_offset_mas = 0.5
            dec_offset_mas = -0.3
            pmra = 5.3
            pmdec = -24.2
            ref_epoch = 57388.0
        end)

    # With only the host luminous, the flux-weighted photocentre must land
    # exactly on it — the check that the weighted-point path and the direct
    # body path are the same physics.
    sysA = model_system(obs=[mkobs(target=:A)])
    sysP = model_system(obs=[mkobs(target=Photocentre)])
    llA = Octofitter.make_ln_like(sysA)(sysA, Octofitter.make_arr2nt(sysA)(Float64[]))
    llP = Octofitter.make_ln_like(sysP)(sysP, Octofitter.make_arr2nt(sysP)(Float64[]))
    @test llA == llP

    # Detrending removes the linear part of the excursion, so the modelled
    # abscissae differ, but the constant and slope of the *difference* must be
    # what was removed: refitting a line to the detrended excursion gives zero.
    obsd = GaiaDR4AstromObs(scans; target=:A, ref=Barycentre, name="GaiaDR4", detrend=true,
        variables=@variables begin
            astrometric_jitter = 0.05
            ra_offset_mas = 0.0; dec_offset_mas = 0.0; pmra = 0.0; pmdec = 0.0
            ref_epoch = 57388.0
        end)
    sysd = model_system(obs=[obsd])
    ntd = Octofitter.make_arr2nt(sysd)(Float64[])
    lnl = Octofitter.make_ln_like(sysd)
    posys = lnl.build(ntd)
    ep, maps = Octofitter.epoch_plan(sysd)
    traj = orbitsolve(posys, ep)
    ctx = Octofitter.ObsContext(ntd, ntd.observations.GaiaDR4, posys, traj, maps[obsd])
    sim = Octofitter.simulate(obsd, ctx)
    Δt = obsd.detrend_Δt
    @test abs(sum(sim.ra_offset)) < 1e-9
    @test abs(sum(Δt .* sim.ra_offset)) < 1e-9
end

# ---------------------------------------------------------------------------
# Two blended catalog sources in a 2+2 quadruple.
#
# Two tight pairs several arcseconds apart: only intra-pair blending is ever
# possible, so each source's membership is structurally fixed and belongs in
# the spec. What is under test is that the resulting signal carries *both*
# levels of motion — the pair's wide orbit and the intra-pair photocentric
# wobble — from one dot product over absolute body states.
# ---------------------------------------------------------------------------

const QUAD_EPOCHS = collect(range(56000.0, 59650.0, length=48))

function quad_system(; target_A, target_B, flux_Ab=0.6, flux_Bb=0.45)
    Aa = Octofitter.Body(name="Aa", variables=@variables begin
        mass = 0.6
        flux_G = 1.0
    end)
    Ab = Octofitter.Body(name="Ab", about=Aa, variables=@variables begin
        mass = 0.4
        flux_G = $flux_Ab
        a = 0.5; e = 0.1; i = 0.3; ω = 0.2; Ω = 0.1; tp = 55000.0
    end)
    Ba = Octofitter.Body(name="Ba", variables=@variables begin
        mass = 0.8
        flux_G = 1.0
    end)
    Bb = Octofitter.Body(name="Bb", about=Ba, variables=@variables begin
        mass = 0.7
        flux_G = $flux_Bb
        a = 0.6; e = 0.2; i = 0.4; ω = 0.3; Ω = 0.2; tp = 55100.0
    end)
    wide = Octofitter.Orbit(name="AB", exterior=(Ba, Bb), about=(Aa, Ab),
        variables=@variables begin
            a = 50.0; e = 0.3; i = 0.5; ω = 0.4; Ω = 0.3; tp = 50000.0
        end)
    scans = (epoch=QUAD_EPOCHS,
        scan_pos_angle=[0.7 * k for k in eachindex(QUAD_EPOCHS)],
        parallax_factor_al=[0.4 * sinpi(k / 7) for k in eachindex(QUAD_EPOCHS)],
        centroid_pos_al=zeros(length(QUAD_EPOCHS)),
        centroid_pos_error_al=fill(0.05, length(QUAD_EPOCHS)))
    # No reference-point motion: `simulate` then returns the target-vs-ref
    # excursion alone, which is the quantity under test.
    quiet() = @variables begin
        astrometric_jitter = 0.05
        ra_offset_mas = 0.0; dec_offset_mas = 0.0
        pmra = 0.0; pmdec = 0.0; ref_epoch = 57388.0
    end
    obsA = GaiaDR4AstromObs(scans; target=target_A, ref=Barycentre, name="srcA",
        variables=quiet())
    obsB = GaiaDR4AstromObs(scans; target=target_B, ref=Barycentre, name="srcB",
        variables=quiet())
    sys = Octofitter.System(name="quad", bodies=[Aa, Ab, Ba, Bb, wide],
        observations=[obsA, obsB], variables=@variables begin plx = 20.0 end)
    return sys, obsA, obsB
end

# Modelled RA/Dec excursions of one observation, in table order.
function quad_simulate(sys, obs)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep)
    θ_obs = getproperty(nt.observations, Symbol(Octofitter.likelihoodname(obs)))
    ctx = Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs])
    sim = Octofitter.simulate(obs, ctx)
    return sim.ra_offset, sim.dec_offset
end

@testset "2+2 quadruple: two blended sources carry every level of motion" begin
    # (0) The spec survives model construction and shows itself.
    sys, obsA, obsB = quad_system(target_A=Photocentre(:G, (:Aa, :Ab)),
        target_B=Photocentre(:G, (:Ba, :Bb)))
    @test Octofitter.refspecs(obsA) === (Photocentre(:G, (:Aa, :Ab)), Barycentre)
    @test occursin("Photocentre(:G, (Aa, Ab)) vs Barycentre",
        sprint(show, MIME"text/plain"(), obsA))

    raA, decA = quad_simulate(sys, obsA)
    raB, _ = quad_simulate(sys, obsB)

    # (1) Dark secondaries ⇒ each source degrades to its pair's primary,
    # exactly. This is the check that the weighted-point path and the direct
    # body path are the same physics for a subset, not just for the system.
    darksys, darkA, darkB = quad_system(target_A=Photocentre(:G, (:Aa, :Ab)),
        target_B=Photocentre(:G, (:Ba, :Bb)), flux_Ab=0.0, flux_Bb=0.0)
    bodysys, bodyA, bodyB = quad_system(target_A=:Aa, target_B=:Ba,
        flux_Ab=0.0, flux_Bb=0.0)
    @test quad_simulate(darksys, darkA)[1] == quad_simulate(bodysys, bodyA)[1]
    @test quad_simulate(darksys, darkB)[1] == quad_simulate(bodysys, bodyB)[1]

    # (2) The wide-orbit motion is in each source. Against the *pair
    # barycentre* target, the residual is the intra-pair wobble alone.
    barysys, baryA, baryB = quad_system(target_A=Barycentre(:Aa, :Ab),
        target_B=Barycentre(:Ba, :Bb))
    raA_bary, _ = quad_simulate(barysys, baryA)
    raB_bary, _ = quad_simulate(barysys, baryB)
    @test maximum(abs, raA_bary) > 100          # mas — the 50 AU orbit at 20 mas plx
    wobA = raA .- raA_bary
    wobB = raB .- raB_bary
    # the wide orbit dominates each source's track…
    @test maximum(abs, wobA) < 0.2 * maximum(abs, raA_bary)
    # …but the wobble is really there — the sources are not the barycentres.
    # Predicted amplitudes: |f/(1+f) − m₂/M| × pair separation ≈ 0.025 × 10 mas
    # for pair A and 0.156 × 12 mas for pair B. (3) pins them exactly; these
    # only assert the signal is far above numerical noise.
    @test maximum(abs, wobA) > 0.1
    @test maximum(abs, wobB) > 1.0

    # (3) The wobble is exactly the pair's flux-weighted photocentre offset:
    # |f/(1+f) − m₂/M| times the pair's own separation, and its *sign* is
    # load-bearing (here f/(1+f) < m₂/M, so it lands opposite the secondary).
    relsys, relA, _ = quad_system(target_A=:Ab, target_B=:Bb)
    # target=:Ab vs ref=Barycentre isn't the pair separation; build it directly
    posys = Octofitter.make_ln_like(sys).build(Octofitter.make_arr2nt(sys)(Float64[]))
    ep, _ = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep)
    idx = Octofitter.epoch_plan(sys)[2][obsA]
    rel = [raoff(traj[idx[k]], :Ab, :Aa) for k in eachindex(QUAD_EPOCHS)]
    coeffA = 0.6 / 1.6 - 0.4 / 1.0
    @test coeffA < 0
    @test maximum(abs, wobA .- coeffA .* rel) < 1e-9 * maximum(abs, rel)

    # (4) Changing one source's secondary flux moves that source and no other.
    sys2, obsA2, obsB2 = quad_system(target_A=Photocentre(:G, (:Aa, :Ab)),
        target_B=Photocentre(:G, (:Ba, :Bb)), flux_Bb=0.9)
    @test quad_simulate(sys2, obsA2)[1] == raA
    @test quad_simulate(sys2, obsB2)[1] != raB

    # (5) A whole-system photocentre is a *different* model from two blended
    # sources, and must not be silently substituted for one.
    allsys, allA, _ = quad_system(target_A=Photocentre(:G), target_B=Photocentre(:G))
    @test maximum(abs, quad_simulate(allsys, allA)[1] .- raA) > 1.0
end

@testset "northangle sign convention (issues #141/#142)" begin
    # `northangle` must rotate the data the same way on the sky whether the
    # astrometry is given as (sep, pa) or as (ra, dec). Ported from main's
    # regression test for sefffal/Octofitter.jl#141: the v2 branch forked before
    # that fix landed, so the rewritten `ln_like` reintroduced the sign error in
    # the radec branch. The two formats disagreeing is the whole bug.

    posys = reference_system()
    tr = orbitsolve(posys, EPOCHS_A)
    ra_m = [raoff(tr[k], :b, :A) for k in eachindex(EPOCHS_A)]
    dec_m = [decoff(tr[k], :b, :A) for k in eachindex(EPOCHS_A)]

    # Position angle measured North through East, as the `:pa` column uses.
    pa_m = atan.(ra_m, dec_m)
    sep_m = hypot.(ra_m, dec_m)

    # Synthetic data: the model rotated by a known offset, written out both ways.
    ε = 0.05 # radians
    pa_d = pa_m .+ ε
    ra_d = sep_m .* sin.(pa_d)
    dec_d = sep_m .* cos.(pa_d)

    n = length(EPOCHS_A)
    tab_seppa = (epoch=EPOCHS_A, sep=sep_m, pa=pa_d,
                 σ_sep=fill(1.0, n), σ_pa=fill(0.001, n))
    tab_radec = (epoch=EPOCHS_A, ra=ra_d, dec=dec_d,
                 σ_ra=fill(1.0, n), σ_dec=fill(1.0, n))

    # Log-likelihood as a function of the northangle value alone.
    function northangle_ll(tab)
        obs = RelAstromObs(tab; target=:b, ref=:A, name="inst",
            variables=@variables begin
                northangle ~ Normal(0, deg2rad(30))
            end)
        sys = model_system(obs=[obs])
        arr2nt = Octofitter.make_arr2nt(sys)
        lnlike = Octofitter.make_ln_like(sys)
        return δ -> lnlike(sys, arr2nt([δ]))
    end

    ll_seppa = northangle_ll(tab_seppa)
    ll_radec = northangle_ll(tab_radec)

    grid = range(-0.1, 0.1, length=2001)
    best_seppa = grid[argmax(ll_seppa.(grid))]
    best_radec = grid[argmax(ll_radec.(grid))]

    # Data rotated by +ε must be undone by a northangle of -ε, in both formats.
    @test best_seppa ≈ -ε atol = 1e-3
    @test best_radec ≈ -ε atol = 1e-3
    @test sign(best_radec) == sign(best_seppa)

    # A zero northangle is a no-op regardless of format.
    @test ll_seppa(0.0) == ll_seppa(-0.0)
    @test ll_radec(0.0) == ll_radec(-0.0)
end
