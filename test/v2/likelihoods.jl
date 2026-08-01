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
    @test Photocentre(:G) === Octofitter.PhotocentreSpec{:G}()
    @test_throws r"must be a `Body`" Octofitter.refspec(3.0)

    # Resolution against a solved system produces PlanetOrbits references.
    posys = reference_system()
    @test Octofitter.resolveref(posys, Octofitter.refspec(:b)) === PlanetOrbits.BodyRef(2)
    wp = Octofitter.resolveref(posys, Barycentre)
    @test wp isa PlanetOrbits.WeightedPoint
    @test sum(wp.w) ≈ 1
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
