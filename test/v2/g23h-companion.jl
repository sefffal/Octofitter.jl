# G23HObs pointed at a *resolved stellar companion's own* Gaia source.
#
# Will's ask on PR #140 was to be able to "attach a full set of Gaia
# observations to a particular companion without the overhead of G23HObs
# etc." — the case `feature/gaia-dr3-companion` (PR #104) built a separate
# `GaiaDR3CompanionObs` type for. On the v9 surface that needs no new type:
# `host=` names whichever body the catalog source is centred on, so a
# resolved secondary's own DR3 row is `G23HObs(host=B, …)` restricted to the
# DR3 channels. See `design/gaia-dr3-companion-v2-decision.md`.
#
# This file is the evidence for that claim, and for the two traps in it.
# Everything runs offline against the `fixtures/` G23H catalog row and scan
# forecast, with `hip_id` set to NaN so the source looks like what it is
# meant to represent: a faint secondary with a Gaia entry and no Hipparcos
# one.

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
pair those fluxes must not reach the observation, which is the point of the
`fluxratio` argument below.
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
function comp_obs(; companions, seed=1234)
    Random.seed!(seed)
    with_logger(COMP_QUIET) do
        G23HObs(; host=:B, companions, ref=Barycentre,
            gaia_id=COMP_CAT.gaia_source_id, catalog=COMP_CAT,
            forecast_table=COMP_FORECAST,
            channels=(:ra_dr3, :dec_dr3), ueva_mode=:none,
            include_rv=false, include_iad=false,
            frame_shift=false, freeze_epochs=true)
    end
end

"""
`fluxratio` is the observation's blending vector, indexed by `companions=`.
`[0.0]` is a *resolved* pair — A contributes no light to B's source;
`nothing` falls through to the bodies' own `flux_G`, which for this
configuration makes the bright primary dominate the faint secondary's
photocentre.
"""
function comp_system(obs; a, mass_A, fluxratio=[0.0])
    A, B = comp_bodies(; a, mass_A)
    vars = isnothing(fluxratio) ?
        (@variables begin
            plx = $(COMP_CAT.parallax); ra = $(COMP_CAT.ra); dec = $(COMP_CAT.dec)
            pmra = $(COMP_CAT.pmra_dr3); pmdec = $(COMP_CAT.pmdec_dr3)
            rv = 0.0; ref_epoch = 57388.5
        end) :
        (@variables begin
            plx = $(COMP_CAT.parallax); ra = $(COMP_CAT.ra); dec = $(COMP_CAT.dec)
            pmra = $(COMP_CAT.pmra_dr3); pmdec = $(COMP_CAT.pmdec_dr3)
            rv = 0.0; ref_epoch = 57388.5
            fluxratio = $fluxratio
        end)
    with_logger(COMP_QUIET) do
        System(name="src", bodies=[A, B], observations=[obs], verbosity=0, variables=vars)
    end
end

comp_lnlike(sys) = (θ = Octofitter.make_arr2nt(sys)(Float64[]);
                    Octofitter.make_ln_like(sys, θ)(sys, θ))

@testset "companions=() silently drops the orbital signal" begin
    # THE TRAP. `host=B, companions=()` reads as "this catalog source is just
    # body B", and constructs happily — but `_g23h_simulate!` gates the whole
    # sky-path computation on `any_active = any(active)`, where `active` runs
    # over `companions=` only. With none declared it takes the "every
    # per-transit perturbation is exactly zero" branch, which is true when the
    # host *is* the body the barycentre sits on and false here.
    #
    # If this test ever fails, the fast path has been fixed — read
    # `design/gaia-dr3-companion-v2-decision.md` §"Gap 1" and delete this
    # testset rather than adjusting it.
    obs = comp_obs(companions=())
    lls = [comp_lnlike(comp_system(obs; a, mass_A=1.0)) for a in (10.0, 750.0, 3000.0)]
    @test all(==(lls[1]), lls)
end

@testset "companions=(A,) makes the source move with the orbit" begin
    obs = comp_obs(companions=(:A,))
    lls = [comp_lnlike(comp_system(obs; a, mass_A=1.0)) for a in (10.0, 750.0, 3000.0)]
    @test allunique(lls)
    # Wider orbit ⇒ smaller reflex ⇒ closer to the (unperturbed) catalog row.
    @test issorted(lls)

    # The blending spelling is load-bearing: with the bodies' own fluxes the
    # bright primary dominates B's photocentre, and the likelihood is a
    # different function. Neither is a default the user can ignore.
    blended = [comp_lnlike(comp_system(obs; a, mass_A=1.0, fluxratio=nothing))
               for a in (10.0, 750.0, 3000.0)]
    @test all(abs.(blended .- lls) .> 1.0)
end

@testset "recovers an injected orbit and host mass" begin
    a_true, m_true = 750.0, 1.0
    obs = comp_obs(companions=(:A,))
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

@testset "channels= restricts the datum, not the cost" begin
    # The "overhead" of G23HObs is not in its channel set. `transit_priorities`
    # marginalizes the epoch selection over the *whole* forecast pool and is
    # declared regardless of which channels survive, so a two-channel
    # observation carries one sampled dimension per forecast transit.
    Random.seed!(1234)
    sampled = with_logger(COMP_QUIET) do
        G23HObs(; host=:B, companions=(:A,), ref=Barycentre,
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
