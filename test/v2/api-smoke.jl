# A working example of every ported observation type and prior term, plus the
# public analysis entry points.
#
# This file exists to be *executable documentation*: each testset builds one
# small but complete model, evaluates its log density, and (where the type is
# differentiable in practice) takes a couple of gradient steps. Nothing here
# samples. If a snippet in `docs/src/` disagrees with the corresponding
# testset below, the testset is the one that has been run.
#
# It is deliberately self-contained — no helper from another test file — so it
# can be run on its own:
#
#     julia --project=. -e 'using Test, Octofitter, Distributions, PlanetOrbits,
#                                 LinearAlgebra, Random, Statistics;
#                           include("test/v2/api-smoke.jl")'
#
# The subpackage sections load `OctofitterImages` / `OctofitterRadialVelocity` /
# `OctofitterInterferometry` if they are resolvable in the active environment
# and skip with a message if they are not, so this runs both under core's own
# `Pkg.test()` and in the shared `testenv/` sandbox where all four are
# dev-linked.

using Test
using Octofitter
using Octofitter: Octofitter, Body, System
using Octofitter.TypedTables: Table
using Distributions
using PlanetOrbits
using LinearAlgebra
using Random
using Statistics

# ---------------------------------------------------------------------------
# Shared scaffolding
# ---------------------------------------------------------------------------

const SMOKE_EPOCHS = collect(range(56000.0, 59000.0, length=8))

"""
Host + one or two companions, with everything an observation might read
declared: sampled orbital elements (so the model has free parameters to
differentiate), masses in M⊙, and a fixed `flux_H` per body.

`flux_H` is fixed here only to keep each testset's free parameters to the
thing that testset is about — the `PhotometryObs` testset below samples a
flux for real. (The one configuration that does not work is a flux sampled
with *nothing else* free; see `v2/photometry-and-priors.jl`.)
"""
function smoke_bodies(; two_companions=false, flux_A=1.0, flux_b=0.05, flux_c=0.01)
    A = Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.2)
        flux_H = $flux_A
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass ~ LogUniform(0.5mjup, 50mjup)
        a ~ LogUniform(1.0, 8.0)
        e ~ Uniform(0.0, 0.6)
        i ~ Sine()
        ω ~ Uniform(0, 2pi)
        Ω ~ Uniform(0, 2pi)
        # An offset from a base epoch rather than a raw MJD: a ~5.7e4 absolute
        # value makes any finite-difference check look like an AD bug.
        dtp ~ Uniform(-800, 800)
        tp = 56000.0 + dtp
        flux_H = $flux_b
    end)
    !two_companions && return (A, b)
    c = Body(name="c", about=A, variables=@variables begin
        mass ~ LogUniform(0.5mjup, 50mjup)
        a ~ LogUniform(15.0, 40.0)
        e ~ Uniform(0.0, 0.3)
        i ~ Sine()
        ω ~ Uniform(0, 2pi)
        Ω ~ Uniform(0, 2pi)
        dtp2 ~ Uniform(-3000, 3000)
        tp = 56000.0 + dtp2
        flux_H = $flux_c
    end)
    return (A, b, c)
end

function smoke_system(observations, bodies=smoke_bodies(); name="smoke")
    return System(name=name, bodies=collect(bodies), observations=collect(observations),
        variables=@variables begin
            plx ~ truncated(Normal(25.0, 1.0), lower=1.0)
        end)
end

"""
Build the log-density model, evaluate it at a prior draw, and take `nsteps`
gradient steps. Returns the model so a caller can go on poking at it.
"""
function smoke_eval(sys; nsteps=2, seed=1, gradient=true)
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θ = Octofitter.sample_priors(Random.Xoshiro(seed), sys)
    θt = model.link(θ)
    @test isfinite(model.ℓπcallback(θt))
    if gradient
        for _ in 1:nsteps
            ℓ, g = model.∇ℓπcallback(θt)
            @test isfinite(ℓ)
            @test all(isfinite, g)
            θt = θt .+ 1e-7 .* g
        end
    end
    return model
end

"Exact relative astrometry of `target` about `ref`, from a known trajectory."
function smoke_astrom_table(epochs=SMOKE_EPOCHS; σ=2.0)
    A = PlanetOrbits.Body(mass=1.0 - 5mjup, name=:A)
    b = PlanetOrbits.Body(mass=5mjup, name=:b)
    sys = PlanetOrbits.System((A, b),
        (PlanetOrbits.Orbit(b, about=A; a=3.0, e=0.2, i=0.9, ω=1.2, Ω=2.0, tp=56000.0),);
        plx=25.0)
    tr = orbitsolve(sys, epochs)
    return Table(epoch=collect(epochs),
        ra=[raoff(tr[k], :b, :A) for k in eachindex(epochs)],
        dec=[decoff(tr[k], :b, :A) for k in eachindex(epochs)],
        σ_ra=fill(σ, length(epochs)), σ_dec=fill(σ, length(epochs)))
end

# ---------------------------------------------------------------------------
# Core: relative astrometry
# ---------------------------------------------------------------------------

@testset "RelAstromObs" begin
    A, b = smoke_bodies()

    # RA/Dec spelling.
    astrom = RelAstromObs(smoke_astrom_table();
        target=b, ref=A, name="GPI",
        variables=@variables begin
            jitter ~ LogUniform(0.1, 10.0)     # optional; mas, added in quadrature
        end)
    @test Octofitter.refspecs(astrom) == (Octofitter.refspec(b), Octofitter.refspec(A))
    @test Octofitter.epochs(astrom) == SMOKE_EPOCHS
    smoke_eval(smoke_system([astrom], (A, b)))

    # Sep/PA spelling — the same type, different columns.
    tab = smoke_astrom_table()
    seppa = RelAstromObs((epoch=tab.epoch,
            sep=sqrt.(tab.ra .^ 2 .+ tab.dec .^ 2),
            pa=rem2pi.(atan.(tab.ra, tab.dec), RoundDown),
            σ_sep=fill(2.0, length(tab.epoch)), σ_pa=fill(0.01, length(tab.epoch)));
        target=b, ref=A, name="NIRC2")
    smoke_eval(smoke_system([seppa], (A, b)))

    # A companion measured against another companion, and against the
    # barycentre: `ref` is a real reference, not an implicit host.
    A2, b2, c2 = smoke_bodies(two_companions=true)
    rel = RelAstromObs(smoke_astrom_table(); target=c2, ref=b2, name="c-vs-b")
    bary = RelAstromObs(smoke_astrom_table(); target=b2, ref=Barycentre, name="b-vs-bary")
    smoke_eval(smoke_system([rel, bary], (A2, b2, c2)))

    # Row subsetting, the contract cross-validation is built on.
    sub = Octofitter.likeobj_from_epoch_subset(astrom, [1, 3])
    @test length(Octofitter.epochs(sub)) == 2
end

# ---------------------------------------------------------------------------
# Core: radial velocity
# ---------------------------------------------------------------------------

@testset "RadialVelocityObs" begin
    A, b = smoke_bodies()
    rvtab = Table(epoch=SMOKE_EPOCHS,
        rv=30 .* sin.(SMOKE_EPOCHS ./ 400),
        σ_rv=fill(5.0, length(SMOKE_EPOCHS)))

    # The stellar reflex signal. `offset` and `jitter` are NOT auto-injected in
    # v2 — a v1 model copied across without them fits with no zero point and no
    # jitter.
    reflex = RadialVelocityObs(rvtab;
        target=A, ref=Barycentre, name="HARPS",
        variables=@variables begin
            offset ~ Normal(0, 100)          # m/s
            jitter ~ LogUniform(0.1, 100)    # m/s
        end)
    @test Octofitter.refspecs(reflex) == (Octofitter.refspec(A), Barycentre)
    smoke_eval(smoke_system([reflex], (A, b)))

    # Relative RV: same type, different refs.
    relrv = RadialVelocityObs(rvtab;
        target=b, ref=A, name="CRIRES",
        variables=@variables begin
            jitter ~ LogUniform(0.1, 100)
        end)
    smoke_eval(smoke_system([relrv], (A, b)))

    # A long-term trend. `trend_function(θ_obs, epoch)` reads variables from
    # this observation's own block.
    trended = RadialVelocityObs(rvtab;
        target=A, ref=Barycentre, name="Lick",
        trend_function=(θ_obs, epoch) -> θ_obs.trend_slope * (epoch - 57000),
        variables=@variables begin
            offset ~ Normal(0, 100)
            jitter ~ LogUniform(0.1, 100)
            trend_slope ~ Normal(0, 0.1)     # m/s/day
        end)
    smoke_eval(smoke_system([trended], (A, b)))

    # A Gaussian process, through the three backend hooks. The default hook
    # methods implement AbstractGPs by duck typing; the toy backend below shows
    # the shape a backend has to satisfy, without pulling in a dependency.
    # (AbstractGPs and Celerite are both exercised for real in
    # `OctofitterRadialVelocity/test/runtests.jl`.)
    struct SmokeGP{T}
        amp::T
        len::T
    end
    struct SmokeGPCond{T,V}
        gp::SmokeGP{T}
        Σ::Matrix{V}
    end
    function Octofitter.gp_condition(gp::SmokeGP, epochs, σ²)
        K = [gp.amp^2 * exp(-((s - t) / gp.len)^2 / 2) for s in epochs, t in epochs]
        return SmokeGPCond(gp, K + Diagonal(collect(σ²)))
    end
    Octofitter.gp_ln_like(fx::SmokeGPCond, resid) = logpdf(MvNormal(Symmetric(fx.Σ)), resid)

    gped = RadialVelocityObs(rvtab;
        target=A, ref=Barycentre, name="ESPRESSO",
        gaussian_process=θ_obs -> SmokeGP(θ_obs.gp_η₁, θ_obs.gp_η₂),
        variables=@variables begin
            offset ~ Normal(0, 100)
            jitter ~ LogUniform(0.1, 100)
            gp_η₁ ~ LogUniform(1.0, 100.0)     # amplitude, m/s
            gp_η₂ ~ LogUniform(1.0, 500.0)     # length scale, days
        end)
    smoke_eval(smoke_system([gped], (A, b)))

    # Held-out rows: `ln_like` scores them conditioned on the rest.
    heldout = RadialVelocityObs(rvtab[1:end-1];
        target=A, ref=Barycentre, name="holdout",
        held_out=rvtab[end:end],
        variables=@variables begin
            offset ~ Normal(0, 100)
            jitter ~ LogUniform(0.1, 100)
        end)
    @test length(Octofitter.epochs(heldout)) == length(SMOKE_EPOCHS)
    smoke_eval(smoke_system([heldout], (A, b)))

    # Cross-validation subsetting.
    @test length(Octofitter.epochs(
        Octofitter.likeobj_from_epoch_subset(reflex, [2, 4]))) >= 2
end

# ---------------------------------------------------------------------------
# Core: photometry
# ---------------------------------------------------------------------------

@testset "PhotometryObs" begin
    # The flux parameter lives on the *body*, not on the observation: two
    # instruments in the same band share one `flux_H`.
    A = Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.2)
        flux_H = 1.0                          # sets the contrast scale
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass ~ LogUniform(0.5mjup, 50mjup)
        a ~ LogUniform(1.0, 8.0)
        e ~ Uniform(0.0, 0.6)
        i ~ Sine()
        ω ~ Uniform(0, 2pi)
        Ω ~ Uniform(0, 2pi)
        dtp ~ Uniform(-800, 800)
        tp = 56000.0 + dtp
        flux_H ~ LogUniform(1e-4, 1e-1)       # a sampled contrast ratio
    end)

    phot = PhotometryObs(Table(phot=[0.012, 0.0131], σ_phot=[0.002, 0.0018]);
        target=b, band=:H, name="NIRC2")
    @test Octofitter.refspecs(phot) == (Octofitter.refspec(b),)
    @test isempty(Octofitter.epochs(phot))    # no epoch: forces no trajectory solve
    smoke_eval(smoke_system([phot], (A, b)))

    # `band=:default` (the default) reads a bare `flux` variable instead.
    A2 = Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.2)
        flux = 1.0
    end)
    b2 = Body(name="b", about=A2, variables=@variables begin
        mass = 5mjup
        a ~ LogUniform(1.0, 8.0)
        e = 0.1; i = 0.6; ω = 1.0; Ω = 0.4; tp = 56000.0
        flux ~ LogUniform(1e-4, 1e-1)
    end)
    smoke_eval(smoke_system(
        [PhotometryObs(Table(phot=[0.01], σ_phot=[0.001]); target=b2, name="bare")],
        (A2, b2)))
end

# ---------------------------------------------------------------------------
# Core: the configuration priors (no data, no epochs)
# ---------------------------------------------------------------------------

@testset "dynamical priors" begin
    bodies = smoke_bodies(two_companions=true)
    A, b, c = bodies
    astrom = RelAstromObs(smoke_astrom_table(); target=b, ref=A, name="astrom")

    for prior in (
        OrbitOrderPrior(b, c),
        NonCrossingPrior(),
        NonCrossingPrior(bodies=(b, c)),
        LimitClosestApproachAUPrior(1.0),                 # soft threshold only
        LimitClosestApproachAUPrior(0.2, 1.0),            # hard + soft
        LimitClosestApproachAUPrior(0.2, 1.0; bodies=(b, c)),
        HillStabilityPrior(),
        HillStabilityPrior(bodies=(b, c)),
    )
        @test Octofitter._isprior(prior)
        @test isempty(Octofitter.epochs(prior))
        smoke_eval(smoke_system([astrom, prior], bodies))
    end

end

@testset "ObsPriorONeil2019" begin
    A, b = smoke_bodies()

    # Wrapping relative astrometry: `orbit` defaults to the wrapped
    # observation's own `target`. Attach ONLY the wrapper — listing the
    # wrapped likelihood as well double-counts the data.
    astrom = RelAstromObs(smoke_astrom_table(); target=b, ref=A, name="astrom")
    smoke_eval(smoke_system([ObsPriorONeil2019(astrom)], (A, b)))

    # Wrapping a stellar-reflex RV: the default would be the host, which has no
    # orbit of its own, so `orbit=` is required.
    rvs = RadialVelocityObs(Table(epoch=SMOKE_EPOCHS,
            rv=30 .* sin.(SMOKE_EPOCHS ./ 400), σ_rv=fill(5.0, length(SMOKE_EPOCHS)));
        target=A, ref=Barycentre, name="rvs",
        variables=@variables begin
            offset ~ Normal(0, 100)
            jitter ~ LogUniform(0.1, 100)
        end)
    smoke_eval(smoke_system([ObsPriorONeil2019(rvs; orbit=b)], (A, b)))

    # Several orbits sum.
    A2, b2, c2 = smoke_bodies(two_companions=true)
    rvs2 = RadialVelocityObs(Table(epoch=SMOKE_EPOCHS,
            rv=30 .* sin.(SMOKE_EPOCHS ./ 400), σ_rv=fill(5.0, length(SMOKE_EPOCHS)));
        target=A2, ref=Barycentre, name="rvs",
        variables=@variables begin
            offset ~ Normal(0, 100)
            jitter ~ LogUniform(0.1, 100)
        end)
    smoke_eval(smoke_system([ObsPriorONeil2019(rvs2; orbit=(b2, c2))], (A2, b2, c2)))

end

# ---------------------------------------------------------------------------
# Core: Gaia DR4 along-scan astrometry
# ---------------------------------------------------------------------------

@testset "GaiaDR4AstromObs" begin
    A, b = smoke_bodies()
    n = 60
    ep = collect(range(57389.0, 59000.0, length=n))
    ψ = range(0, 6pi, length=n)
    scans = Table(epoch=ep,
        # `GaiaDR4AstromObs` ingests the scan angle in degrees.
        scan_pos_angle=rad2deg.(ψ),
        parallax_factor_al=cos.(ψ),
        centroid_pos_al=zeros(n),
        centroid_pos_error_al=fill(0.2, n))

    dr4 = GaiaDR4AstromObs(scans;
        target=Photocentre(:G, (A, b)), ref=Barycentre, name="GaiaDR4",
        variables=@variables begin
            ra_offset_mas ~ Normal(0, 100)
            dec_offset_mas ~ Normal(0, 100)
            pmra ~ Normal(0, 100)
            pmdec ~ Normal(0, 100)
            ref_epoch = 57388.5
        end)
    @test Octofitter.refspecs(dr4) == (Photocentre(:G, (A, b)), Barycentre)

    # A photocentre needs the members to declare the band it is weighted by.
    Ag = Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.2)
        flux_G = 1.0
    end)
    bg = Body(name="b", about=Ag, variables=@variables begin
        mass ~ LogUniform(0.5mjup, 50mjup)
        a ~ LogUniform(1.0, 8.0)
        e ~ Uniform(0.0, 0.6)
        i ~ Sine()
        ω ~ Uniform(0, 2pi)
        Ω ~ Uniform(0, 2pi)
        dtp ~ Uniform(-800, 800)
        tp = 56000.0 + dtp
        flux_G = 0.0                         # dark to Gaia
    end)
    dr4g = GaiaDR4AstromObs(scans;
        target=Photocentre(:G, (Ag, bg)), ref=Barycentre, name="GaiaDR4",
        variables=@variables begin
            ra_offset_mas ~ Normal(0, 100)
            dec_offset_mas ~ Normal(0, 100)
            pmra ~ Normal(0, 100)
            pmdec ~ Normal(0, 100)
            ref_epoch = 57388.5
        end)
    smoke_eval(smoke_system([dr4g], (Ag, bg)))
end

# ---------------------------------------------------------------------------
# Core: absolute astrometry (G23H, HGCA, Hipparcos IAD)
#
# These run against the checked-in fixtures in `v2/fixtures/`, so they are
# offline. Construction and one log-density evaluation only: gradients through
# the epoch-marginalization prior are gated in `v2/g23h.jl`, which is where
# they belong.
# ---------------------------------------------------------------------------

@testset "absolute astrometry" begin
    fixdir = joinpath(@__DIR__, "fixtures")
    _fix(name) = Table(Octofitter.CSV.File(joinpath(fixdir, name)))
    function _onerow(t)
        cols = Tuple(Octofitter.Tables.columnnames(t))
        return NamedTuple{cols}(map(c -> getproperty(t, c)[1], cols))
    end
    cat = _onerow(_fix("g23h-catalog-row.csv"))
    hip = (; table=_fix("g23h-hip-iad.csv"), hip_sol=_onerow(_fix("g23h-hip-sol.csv")))

    # Companion flux ratios now default to the bodies' own `flux_G`/`flux_Hp`.
    A = Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.6, 0.1), lower=0.1)
        flux_G = 1.0
        flux_Hp = 1.0
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass ~ LogUniform(1e-4, 0.2)
        a ~ LogUniform(1.0, 50.0)
        e ~ Uniform(0, 0.7)
        i ~ Sine()
        ω ~ Uniform(0, 2pi)
        Ω ~ Uniform(0, 2pi)
        dtp ~ Uniform(-3000, 3000)
        tp = 56000.0 + dtp
        flux_G = 0.0                          # dark
        flux_Hp = 0.0
    end)

    @testset "G23HObs" begin
        # `include_rv=false` — and it is not incidental. The Gaia RV
        # variability channel scores a `NoncentralChisq`, whose log-density
        # does not accept ForwardDiff `Dual`s, so with the default
        # `include_rv=true` the primal log density is finite while the
        # gradient pass collapses the whole model to -Inf. A gradient-based
        # sampler cannot move such a model. Drop the channel (here) or use a
        # derivative-free sampler; `ln_like` warns once if you hit it.
        obs = G23HObs(; target=A, blends=(b,), ref=Barycentre,
            gaia_id=cat.gaia_source_id, catalog=cat,
            forecast_table=_fix("g23h-forecast.csv"), hipparcos=hip,
            dr2_transits_catalog=_fix("g23h-dr2-transits.csv"),
            include_rv=false)
        @test Octofitter._refnames.(Octofitter.refspecs(obs)) == ((:A,), (:b,), ())
        @test :rv_dr3 ∉ obs.table.kind
        # G23H needs the reference point's proper motion. Declaring the *full*
        # absolute frame in the system block is the idiomatic v2 spelling —
        # which of `plx, ra, dec, pmra, pmdec, rv, ref_epoch` are present
        # chooses the frame, and a partial set is rejected at build time. (The
        # alternative, which `v2/g23h.jl` uses to reproduce v1 exactly, is a
        # parallax-only system plus `pmra`/`pmdec` in the observation's block.)
        sys = System(name="g23h", bodies=[A, b], observations=[obs],
            variables=@variables begin
                plx ~ truncated(Normal(cat.parallax, 1.0), lower=1.0)
                ra = $(cat.ra)
                dec = $(cat.dec)
                pmra = $(cat.pmra_dr3)
                pmdec = $(cat.pmdec_dr3)
                rv = 0.0
                ref_epoch = 57388.5
            end)
        # With `:rv_dr3` dropped the model is differentiable, so this takes the
        # gradient steps too — which is the check that the AD-hostile channel
        # was the only obstacle.
        smoke_eval(sys; seed=7)

        # `channels=` restricts the channel set at construction, through the
        # same row filter `likeobj_from_epoch_subset` uses.
        restricted = G23HObs(; target=A, blends=(b,), ref=Barycentre,
            gaia_id=cat.gaia_source_id, catalog=cat,
            forecast_table=_fix("g23h-forecast.csv"), hipparcos=hip,
            dr2_transits_catalog=_fix("g23h-dr2-transits.csv"),
            channels=(:ra_dr3, :dec_dr3), ueva_mode=:none, include_rv=false)
        @test restricted.table.kind == [:ra_dr3, :dec_dr3]
    end

    @testset "HGCAObs" begin
        # A helper over G23HObs, not a type: the six HGCA channels, no UEVA,
        # no IAD, no RV.
        obs = HGCAObs(; target=A, blends=(b,), ref=Barycentre,
            gaia_id=cat.gaia_source_id, catalog=cat,
            forecast_table=_fix("g23h-forecast.csv"), hipparcos=hip,
            dr2_transits_catalog=_fix("g23h-dr2-transits.csv"))
        @test obs isa G23HObs
        @test obs.table.kind ==
              [:ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr3, :dec_dr3]
        sys = System(name="hgca", bodies=[A, b], observations=[obs],
            variables=@variables begin
                plx ~ truncated(Normal(cat.parallax, 1.0), lower=1.0)
                ra = $(cat.ra)
                dec = $(cat.dec)
                pmra = $(cat.pmra_dr3)
                pmdec = $(cat.pmdec_dr3)
                rv = 0.0
                ref_epoch = 57388.5
            end)
        # No `:rv_dr3` channel in the HGCA preset, so this is differentiable.
        smoke_eval(sys; seed=8)

        # The retired types name their replacement rather than failing obscurely.
        @test_throws ErrorException HGCAInstantaneousObs(gaia_id=cat.gaia_source_id)
        @test_throws ErrorException GaiaCatalogFitObs()
    end

    @testset "HipparcosIADObs" begin
        obs = HipparcosIADObs(; target=A, blends=(b,), ref=Barycentre, iad=hip)
        @test Octofitter.refspecs(obs) ==
              (Octofitter.refspec(A), Octofitter.refspec(b), Barycentre)
        sys = System(name="hip", bodies=[A, b], observations=[obs],
            variables=@variables begin
                plx ~ truncated(Normal(cat.parallax, 1.0), lower=1.0)
            end)
        smoke_eval(sys; seed=9)

        # Passing a preloaded IAD twice must not shift it twice.
        obs2 = HipparcosIADObs(; target=A, blends=(b,), iad=hip, recalibrate=true)
        obs3 = HipparcosIADObs(; target=A, blends=(b,), iad=hip, recalibrate=true)
        @test obs2.table.res == obs3.table.res
    end
end

# ---------------------------------------------------------------------------
# The public analysis entry points
# ---------------------------------------------------------------------------

@testset "analysis entry points" begin
    A, b = smoke_bodies()
    astrom = RelAstromObs(smoke_astrom_table(); target=b, ref=A, name="astrom")
    rvtab = Table(epoch=SMOKE_EPOCHS[1:5],
        rv=30 .* sin.(SMOKE_EPOCHS[1:5] ./ 400), σ_rv=fill(5.0, 5))
    rvs = RadialVelocityObs(rvtab; target=A, ref=Barycentre, name="rvs",
        variables=@variables begin
            offset ~ Normal(0, 100)
            jitter ~ LogUniform(0.1, 100)
        end)
    sys = smoke_system([astrom, rvs], (A, b))
    model = Octofitter.LogDensityModel(sys; verbosity=0)

    @testset "priors and construction" begin
        θarr = Octofitter.sample_priors(Random.Xoshiro(2), sys)
        @test all(isfinite, θarr)
        nt = drawfrompriors(sys)
        @test haskey(nt, :bodies)
        # `construct_system` takes the LogDensityModel, and replaces v1's
        # `construct_elements(chain, :b, i)`: there is no per-planet orbit
        # object any more, the trajectory is a property of the system.
        posys = construct_system(model, model.arr2nt(θarr))
        @test posys isa PlanetOrbits.System
        traj = orbitsolve(posys, SMOKE_EPOCHS)
        @test all(isfinite, [raoff(traj[k], :b, :A) for k in eachindex(SMOKE_EPOCHS)])
    end

    @testset "forward simulation" begin
        sim = generate_from_params(sys)                       # noiseless
        @test sim isa Octofitter.System
        simn = generate_from_params(sys, drawfrompriors(sys); add_noise=true)
        @test simn isa Octofitter.System
        @test Octofitter.LogDensityModel(simn; verbosity=0) isa Octofitter.LogDensityModel
    end

    @testset "prior-only models" begin
        # Each data likelihood is replaced by a `BlankLikelihood` carrying the
        # same name and the same `@variables` block, so a prior-only chain
        # lines up column for column with one from the full model.
        p = prior_only_model(sys)
        @test all(o -> o isa Octofitter.BlankLikelihood, p.observations)
        @test Octofitter.likelihoodname.(p.observations) ==
              Octofitter.likelihoodname.(sys.observations)
        pa = prior_only_model(sys; exclude_all=true)
        @test all(o -> o isa Octofitter.BlankLikelihood, pa.observations)
        smoke_eval(p)
        smoke_eval(pa)
    end

    @testset "cross-validation folds" begin
        @test length(Octofitter.generate_systems_per_like(sys)) == 2
        @test length(Octofitter.generate_kfold_systems(sys)) == 2
        filtered = Octofitter.generate_system_filtered_like(o -> Octofitter.likelihoodname(o) == "astrom", sys)
        @test length(filtered.observations) == 1
        per_epoch, eps = Octofitter.generate_system_per_epoch(sys)
        @test length(per_epoch) == length(eps) == length(SMOKE_EPOCHS) + 5
        cum, groups = Octofitter.generate_cumulative_system_per_epoch(sys)
        @test length(cum) == length(groups)
    end

    @testset "pointwise_like" begin
        θs = [model.arr2nt(Octofitter.sample_priors(Random.Xoshiro(s), sys)) for s in (1, 2)]
        chain = Octofitter.result2mcmcchain(θs)
        LL, eps = Octofitter.pointwise_like(model, chain; verbosity=0)
        @test size(LL) == (2, length(SMOKE_EPOCHS) + 5)
        @test all(isfinite, LL)
        # The columns decompose the model's own log-likelihood.
        lnlike = Octofitter.make_ln_like(sys)
        for (i, θ) in enumerate(θs)
            @test sum(LL[i, :]) ≈ lnlike(sys, θ) rtol = 1e-10
        end
    end

    @testset "chain persistence" begin
        θs = [model.arr2nt(Octofitter.sample_priors(Random.Xoshiro(s), sys)) for s in (1, 2)]
        chain = Octofitter.result2mcmcchain(θs)
        mktempdir() do dir
            f = joinpath(dir, "chain.fits")
            Octofitter.savechain(f, chain)
            back = Octofitter.loadchain(f)
            @test Set(String.(keys(back))) == Set(String.(keys(chain)))
            # `model=` checks the chain against the model rather than silently
            # substituting `missing` for a renamed parameter.
            @test Octofitter.loadchain(f; model) isa Octofitter.MCMCChains.Chains
            @test Octofitter.checkchain(model, chain) === chain
        end
    end

    @testset "completeness bookkeeping" begin
        jobs = completeness_jobs(masses=[5mjup, 10mjup], separations=[5.0, 20.0], n_trials=1)
        @test length(jobs) == 4
        @test first(jobs) isa CompletenessJob
        # `assemble_completeness` over an empty result set still produces a map.
        cmap = assemble_completeness(CompletenessResult[], _ -> true;
            masses=[5mjup, 10mjup], separations=[5.0, 20.0])
        @test cmap isa CompletenessMap
        @test size(cmap.completeness) == (2, 2)
    end
end

# ---------------------------------------------------------------------------
# Subpackages, if resolvable in the active environment
# ---------------------------------------------------------------------------

_smoke_have(name::Symbol) = !isnothing(Base.identify_package(String(name)))

if _smoke_have(:OctofitterRadialVelocity)
    @eval using OctofitterRadialVelocity
    @testset "MarginalizedRVObs" begin
        A, b = smoke_bodies()
        rvtab = Table(epoch=SMOKE_EPOCHS,
            rv=30 .* sin.(SMOKE_EPOCHS ./ 400),
            σ_rv=fill(5.0, length(SMOKE_EPOCHS)))
        # The zero point is integrated out analytically, so declaring `offset`
        # here is an error rather than a redundancy.
        marg = MarginalizedRVObs(rvtab;
            target=A, ref=Barycentre, name="HIRES",
            variables=@variables begin
                jitter ~ LogUniform(0.1, 100)
            end)
        @test Octofitter.refspecs(marg) == (Octofitter.refspec(A), Barycentre)
        smoke_eval(smoke_system([marg], (A, b)))
        # Structurally not subsettable: the marginalization couples the rows.
        @test_throws ErrorException Octofitter.likeobj_from_epoch_subset(marg, [1, 2])
    end
else
    @info "api-smoke: OctofitterRadialVelocity not in this environment; skipping."
end

if _smoke_have(:OctofitterImages)
    @eval using OctofitterImages
    # Reached through the subpackage rather than `using AstroImages`, so this
    # file needs no test-only dependency of its own.
    const AstroImages = OctofitterImages.AstroImages
    const AstroImage = AstroImages.AstroImage
    @testset "ImageObs and LogLikelihoodMapObs" begin
        ax = -30.0:1.0:30.0
        img = AstroImages.recenter(AstroImage(
            [exp(-((x - 8.0)^2 + (y + 5.0)^2) / 18) for x in ax, y in ax]))
        lmap = AstroImages.recenter(AstroImage(
            [-(abs(x - 8.0) + abs(y + 5.0)) / 8 for x in ax, y in ax]))

        A, b, c = smoke_bodies(two_companions=true, flux_b=0.9, flux_c=0.3)

        # ONE observation for every companion in the image, not one per planet.
        images = ImageObs(Table(epoch=SMOKE_EPOCHS[1:2],
                image=[img, img], platescale=fill(10.0, 2));
            targets=(b, c), ref=A, band=:H, name="SPHERE",
            variables=@variables begin
                platescale = 1.0
                northangle = 0.0
            end)
        @test Octofitter.refspecs(images) ==
              (Octofitter.refspec(b), Octofitter.refspec(c), Octofitter.refspec(A))
        smoke_eval(smoke_system([images], (A, b, c)))

        likemap = LogLikelihoodMapObs(Table(epoch=SMOKE_EPOCHS[1:2],
                map=[lmap, lmap], platescale=fill(10.0, 2));
            target=b, ref=A, name="GRAVITY",
            variables=@variables begin
                platescale = 1.0
                northangle = 0.0
            end)
        smoke_eval(smoke_system([likemap], (A, b, c)))
    end
else
    @info "api-smoke: OctofitterImages not in this environment; skipping."
end

if _smoke_have(:OctofitterInterferometry)
    @eval using OctofitterInterferometry
    @testset "InterferometryObs" begin
        # A four-telescope array: six baselines and four closure triangles.
        eff_wave = collect(range(2.0e-6, 2.4e-6, length=3))
        baselines = [(46.6, 12.3), (-9.1, 54.8), (33.7, -41.2),
            (-55.7, 42.5), (-12.9, -53.5), (42.8, -96.0)]
        u = [b[1] / eff_wave[l] for b in baselines, l in eachindex(eff_wave)]
        v = [b[2] / eff_wave[l] for b in baselines, l in eachindex(eff_wave)]
        exposure(epoch) = (; epoch, eff_wave, u, v,
            cps_data=fill(3.0, 4, 3), dcps=fill(2.5, 4, 3),
            index_cps1=[1, 1, 2, 4], index_cps2=[4, 5, 6, 6], index_cps3=[2, 3, 3, 5])

        # `flux_K` on every source, host included: there is no privileged
        # primary. Setting the host's to 1.0 makes the rest contrast ratios,
        # which is also the recipe that reproduces v1 bit-for-bit.
        A = Body(name="A", variables=@variables begin
            mass ~ truncated(Normal(1.2, 0.05), lower=0.2)
            flux_K = 1.0
        end)
        b = Body(name="b", about=A, variables=@variables begin
            mass ~ LogUniform(0.5mjup, 50mjup)
            a ~ LogUniform(1.0, 5.0)
            e ~ Uniform(0.0, 0.4)
            i ~ Sine()
            ω ~ Uniform(0, 2pi)
            Ω ~ Uniform(0, 2pi)
            dtp ~ Uniform(-800, 800)
            tp = 58800.0 + dtp
            flux_K = 0.3
        end)

        vis = InterferometryObs(Table([exposure(e) for e in (59000.0, 59200.0)]);
            targets=(A, b), ref=A, band=:K, name="NIRISS-AMI",
            variables=@variables begin
                σ_cp_jitter = 0.0
                platescale = 1.0
                northangle = 0.0
            end)
        @test Octofitter.refspecs(vis) ==
              (Octofitter.refspec(A), Octofitter.refspec(b), Octofitter.refspec(A))
        smoke_eval(System(name="interf", bodies=[A, b], observations=[vis],
            variables=@variables begin
                plx ~ truncated(Normal(100.0, 1.0), lower=1.0)
            end))
    end
else
    @info "api-smoke: OctofitterInterferometry not in this environment; skipping."
end
