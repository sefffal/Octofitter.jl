# G23HObs pointed at a *resolved stellar companion's own* Gaia source.
#
# Will's ask on PR #140 was to be able to "attach a full set of Gaia
# observations to a particular companion without the overhead of G23HObs
# etc." — the case `feature/gaia-dr3-companion` (PR #104) built a separate
# `GaiaDR3CompanionObs` type for. On the v9 surface that needs no new type:
# `target=` names whichever body the catalog source is centred on, so a
# resolved secondary's own DR3 row is `G23HObs(target=B, blends=(), …)`
# restricted to the DR3 channels. See
# `design/gaia-dr3-companion-v2-decision.md`.
#
# This file is the evidence for that claim. It began life on the branch
# `v2-gaia-dr3-companion` as a *characterization* of two gaps in the argument
# surface; both are now fixed, so the same gates assert the fixed behaviour:
#
#   Gap 1 — `blends=()` returned a likelihood that was bit-identical for
#           every orbit, because the fast path was keyed on the declared
#           blend list rather than on whether the target could move relative
#           to `ref`. It must now respond to the orbit.
#   Gap 2 — a resolved pair had no clean spelling: `blends=()` was broken and
#           `blends=(A,)` fell through to the bodies' own `flux_G`, making the
#           bright primary dominate the faint secondary's photocentre. The
#           resolved spelling is now `blends=()` — literally "nothing blends
#           into this source" — and `fluxratio=0.0` says the same thing while
#           still declaring the blend.
#
# Everything runs offline against the `fixtures/` G23H catalog row and scan
# forecast, with `hip_id` set to NaN so the source looks like what it is meant
# to represent: a faint secondary with a Gaia entry and no Hipparcos one.

using Octofitter
using Octofitter: G23HObs, Body, System
using Octofitter.TypedTables: Table
using Octofitter.CSV
using PlanetOrbits
using Distributions, StaticArrays, LinearAlgebra, Random, Logging, Test

const CFIX = joinpath(@__DIR__, "fixtures")
_cfix(name) = Table(CSV.File(joinpath(CFIX, name)))
function _crow(t)
    cols = Tuple(Octofitter.Tables.columnnames(t))
    return NamedTuple{cols}(map(c -> getproperty(t, c)[1], cols))
end

# The catalog row, re-labelled as a secondary: no Hipparcos entry, and the
# DR2 matched-transit count supplied inline so the 300 MB sidecar DataDep is
# not touched.
const COMP_CAT = merge(_crow(_cfix("g23h-catalog-row.csv")),
    (; hip_id=NaN, astrometric_matched_observations_dr2=40.0))
const COMP_FORECAST = _cfix("g23h-forecast.csv")
const COMP_QUIET = ConsoleLogger(stderr, Logging.Warn)

"""
A host `A` and a wide, resolved companion `B` whose own Gaia source is what
the observation models. `flux_G` is declared on both bodies: for a *resolved*
pair those fluxes must not reach the observation, which is what `blends=()`
(and `fluxratio=0.0`) say.
"""
function comp_bodies(; a, mass_A, mass_B=0.30)
    A = Body(name="A", variables=@variables begin
        mass = $mass_A
        flux_G = 1.0
        flux_Hp = 1.0
    end)
    B = Body(name="B", about=:A, variables=@variables begin
        mass = $mass_B
        a = $a
        e = 0.2
        i = 0.9
        ω = 0.4
        Ω = 1.1
        tp = 50000.0
        flux_G = 0.02
        flux_Hp = 0.02
    end)
    return (A, B)
end

"""The DR3-position/proper-motion channels of the secondary's own source."""
function comp_obs(; blends=(), seed=1234, kw...)
    Random.seed!(seed)
    with_logger(COMP_QUIET) do
        G23HObs(; target=:B, blends, ref=Barycentre,
            gaia_id=COMP_CAT.gaia_source_id, catalog=COMP_CAT,
            forecast_table=COMP_FORECAST,
            channels=(:ra_dr3, :dec_dr3), ueva_mode=:none,
            include_rv=false, include_iad=false,
            frame_shift=false, freeze_epochs=true, kw...)
    end
end

function comp_system(obs; a, mass_A)
    A, B = comp_bodies(; a, mass_A)
    vars = @variables begin
        plx = $(COMP_CAT.parallax); ra = $(COMP_CAT.ra); dec = $(COMP_CAT.dec)
        pmra = $(COMP_CAT.pmra_dr3); pmdec = $(COMP_CAT.pmdec_dr3)
        rv = 0.0; ref_epoch = 57388.5
    end
    with_logger(COMP_QUIET) do
        System(name="src", bodies=[A, B], observations=[obs], verbosity=0, variables=vars)
    end
end

comp_lnlike(sys) = (θ = Octofitter.make_arr2nt(sys)(Float64[]);
                    Octofitter.make_ln_like(sys, θ)(sys, θ))

@testset "blends=() is a resolved source, and it moves" begin
    # THE TRIPWIRE. `target=B, blends=()` reads as "this catalog source is
    # just body B", and that is now what it means: B's full reflex about the
    # system barycentre, with nothing blended into it. It used to take the
    # "every per-transit perturbation is exactly zero" branch — true when the
    # target is the body the barycentre sits on, false here — and return the
    # same number for every orbit.
    obs = comp_obs(blends=())
    lls = [comp_lnlike(comp_system(obs; a, mass_A=1.0)) for a in (10.0, 750.0, 3000.0)]
    @test allunique(lls)
    # Wider orbit ⇒ smaller reflex ⇒ closer to the (unperturbed) catalog row.
    @test issorted(lls)

    # Declaring the primary as a blend and gating its light off with the
    # constant keyword is the same physical statement, and gives the same
    # numbers to the last bit: a zero-weight blend is no blend.
    obs0 = comp_obs(blends=(:A,), fluxratio=0.0)
    lls0 = [comp_lnlike(comp_system(obs0; a, mass_A=1.0)) for a in (10.0, 750.0, 3000.0)]
    @test lls0 == lls

    # …and it is not a vacuous agreement: letting the bodies' own fluxes
    # through makes the bright primary dominate B's photocentre, which is a
    # different function.
    blended = [comp_lnlike(comp_system(comp_obs(blends=(:A,)); a, mass_A=1.0))
               for a in (10.0, 750.0, 3000.0)]
    @test all(abs.(blended .- lls) .> 1.0)
end

@testset "recovers an injected orbit and host mass" begin
    # No flux ratios anywhere: `blends=()` needs none, which is the whole
    # point of it being the resolved spelling.
    a_true, m_true = 750.0, 1.0
    obs = comp_obs(blends=())
    sys_truth = comp_system(obs; a=a_true, mass_A=m_true)
    θ_truth = Octofitter.make_arr2nt(sys_truth)(Float64[])
    sys_sim = with_logger(COMP_QUIET) do
        Octofitter.generate_from_params(sys_truth, θ_truth; add_noise=false)
    end
    obs_sim = only(sys_sim.observations)

    # The simulated row is genuinely different from the one we started with:
    # otherwise the profiles below would be testing nothing.
    @test obs_sim.catalog.pmra_dr3 != COMP_CAT.pmra_dr3

    as = [500.0, 650.0, 700.0, 750.0, 800.0, 900.0, 1100.0]
    lls_a = [comp_lnlike(comp_system(obs_sim; a, mass_A=m_true)) for a in as]
    @test as[argmax(lls_a)] == a_true

    ms = [0.6, 0.8, 0.9, 1.0, 1.1, 1.3, 1.7]
    lls_m = [comp_lnlike(comp_system(obs_sim; a=a_true, mass_A=m)) for m in ms]
    @test ms[argmax(lls_m)] == m_true

    # Noiseless data at truth: the DR3 residual is zero, so the two-channel
    # Gaussian sits at its own normalization.
    @test comp_lnlike(comp_system(obs_sim; a=a_true, mass_A=m_true)) >
          maximum(vcat(lls_a[1:3], lls_a[5:end]))
end

@testset "a static source still takes the fast path" begin
    # The fix must not cost the cheap rejection it replaced. With the target
    # *at* the reference point — the ordinary configuration, primary versus
    # barycentre, with every other body massless — every offset really is
    # zero and the fast path is still valid.
    A = Body(name="A", variables=@variables begin
        mass = 1.0
        flux_G = 1.0
    end)
    B = Body(name="B", about=:A, variables=@variables begin
        mass = 0.0
        a = 750.0; e = 0.2; i = 0.9; ω = 0.4; Ω = 1.1; tp = 50000.0
        flux_G = 0.02
    end)
    Random.seed!(1234)
    obs = with_logger(COMP_QUIET) do
        G23HObs(; target=:A, blends=(:B,), ref=Barycentre,
            gaia_id=COMP_CAT.gaia_source_id, catalog=COMP_CAT,
            forecast_table=COMP_FORECAST,
            channels=(:ra_dr3, :dec_dr3), ueva_mode=:none,
            include_rv=false, include_iad=false,
            frame_shift=false, freeze_epochs=true)
    end
    sys = with_logger(COMP_QUIET) do
        System(name="static", bodies=[A, B], observations=[obs], verbosity=0,
            variables=@variables begin
                plx = $(COMP_CAT.parallax); ra = $(COMP_CAT.ra); dec = $(COMP_CAT.dec)
                pmra = $(COMP_CAT.pmra_dr3); pmdec = $(COMP_CAT.pmdec_dr3)
                rv = 0.0; ref_epoch = 57388.5
            end)
    end
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    lnl = Octofitter.make_ln_like(sys, nt)
    posys = lnl.build(nt)
    eps, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, eps; method=sys.method,
        observing_geometry=sys.observing_geometry,
        barycentric_lighttime=sys.barycentric_lighttime)
    θ_obs = getproperty(nt.observations,
        Symbol(Octofitter.normalizename(Octofitter.likelihoodname(obs))))
    ctx = Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs])
    # A massless B leaves the barycentre exactly on A, so the source cannot
    # move: the predicate the fast path is keyed on says so, and the modelled
    # perturbation is exactly zero.
    @test Octofitter._g23h_static_source(Octofitter.ref(ctx, obs.ref),
                                         Octofitter.ref(ctx, obs.target))
    @test Octofitter.simulate(obs, ctx).Δpmra_dr3 == 0.0
    @test isfinite(lnl(sys, nt))
end

@testset "channels= restricts the datum, not the cost" begin
    # The "overhead" of G23HObs is not in its channel set. `transit_priorities`
    # marginalizes the epoch selection over the *whole* forecast pool and is
    # declared regardless of which channels survive, so a two-channel
    # observation carries one sampled dimension per forecast transit.
    Random.seed!(1234)
    sampled = with_logger(COMP_QUIET) do
        G23HObs(; target=:B, blends=(), ref=Barycentre,
            gaia_id=COMP_CAT.gaia_source_id, catalog=COMP_CAT,
            forecast_table=COMP_FORECAST,
            channels=(:ra_dr3, :dec_dr3), ueva_mode=:none,
            include_rv=false, include_iad=false,
            frame_shift=false, freeze_epochs=false)
    end
    @test sampled.table.kind == [:ra_dr3, :dec_dr3]
    sys = comp_system(sampled; a=750.0, mass_A=1.0)
    model = with_logger(COMP_QUIET) do
        Octofitter.LogDensityModel(sys; verbosity=0)
    end
    # Two channels, and still exactly one sampled dimension per pooled
    # transit (every body variable in `comp_system` is fixed, so
    # `transit_priorities` is the entire model dimension).
    @test model.D == length(sampled.gaia_table.epoch)
end
