# The standalone Hipparcos IAD likelihood.
#
# Runs offline against the same checked-in fixtures `v2/g23h.jl` uses — a real
# Hipparcos IAD table and five-parameter solution for HIP 1475, extracted by
# `fixtures/gen-g23h-fixtures.jl`. The loader itself is gated in `v2/g23h.jl`
# ("the Hipparcos IAD loader reproduces v1"); what is gated here is the
# likelihood built on top of it.

using Octofitter: Octofitter, HipparcosIADObs, G23HObs
using Octofitter.TypedTables: Table, FlexTable
using Octofitter.CSV
using PlanetOrbits
using StaticArrays, Distributions, LinearAlgebra, Random, Statistics, Test
using FiniteDiff

# Self-contained fixture access, so this file does not depend on `v2/g23h.jl`
# having been included first.
const HIP_FIX = joinpath(@__DIR__, "fixtures")

_hip_fixture(name) = Table(CSV.File(joinpath(HIP_FIX, name)))

function _hip_onerow(t)
    cols = Tuple(Octofitter.Tables.columnnames(t))
    return NamedTuple{cols}(map(c -> getproperty(t, c)[1], cols))
end

const HIP_IAD = (; table=_hip_fixture("g23h-hip-iad.csv"),
    hip_sol=_hip_onerow(_hip_fixture("g23h-hip-sol.csv")))
const HIP_SOL = HIP_IAD.hip_sol
const HIP_PLX_CATALOG = HIP_SOL.plx
const HIP_CAT = _hip_onerow(_hip_fixture("g23h-catalog-row.csv"))
const HIP_EVAL = let t = _hip_fixture("g23h-eval-point.csv")
    Dict(Symbol(n) => v for (n, v) in zip(t.name, t.value))
end
const HIP_NPOOL = let t = _hip_fixture("g23h-v1-reference.csv")
    Int(Dict(Symbol(n) => v for (n, v) in zip(t.name, t.value))[:n_pool])
end

"""
A one-companion Hipparcos-only model. `f_hp` is the companion's Hp-band
contrast against the host; `use_body_fluxes` decides whether it is spelled as
the bodies' own `flux_Hp` variables or as the observation's `fluxratio_hip`
vector, which must give the same numbers.
"""
function hip_model(; mass=0.03, f_hp=0.05, a=4.7, plx=HIP_SOL.plx,
                   use_body_fluxes=true, recalibrate=false, jitter=0.01)
    hostvars = @variables begin
        mass = $(0.486 - mass)
    end
    compvars = @variables begin
        mass = $mass
        M = 0.486
        a = $a; e = 0.23; i = 0.94; ω = 1.31; Ω = 2.44; tp = 52000.0
    end
    if use_body_fluxes
        hostvars = vcat(hostvars, @variables begin
            flux_Hp = 1.0
        end)
        compvars = vcat(compvars, @variables begin
            flux_Hp = $f_hp
        end)
    end
    host = Octofitter.Body(name="A", variables=hostvars)
    b = Octofitter.Body(name="b", about=host, variables=compvars)

    obsvars = @variables begin
        hip_iad_jitter = $jitter
        iad_Δra = 0.0
        iad_Δdec = 0.0
        iad_Δplx = 0.0
        iad_pmra = $(HIP_SOL.pm_ra)
        iad_pmdec = $(HIP_SOL.pm_de)
    end
    if !use_body_fluxes
        obsvars = vcat(obsvars, @variables begin
            fluxratio_hip = ($f_hp,)
        end)
    end
    obs = HipparcosIADObs(; host, companions=(b,), iad=HIP_IAD,
        recalibrate, variables=obsvars)
    sys = Octofitter.System(name="hiponly", bodies=[host, b], observations=[obs],
        variables=(@variables begin
            plx = $plx
        end), observing_geometry=false)
    return (; sys, obs, host, b)
end

function hip_context(sys, obs, θ=Float64[])
    nt = Octofitter.make_arr2nt(sys)(θ)
    lnl = Octofitter.make_ln_like(sys, nt)
    posys = lnl.build(nt)
    eps, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, eps; method=sys.method,
        observing_geometry=sys.observing_geometry,
        barycentric_lighttime=sys.barycentric_lighttime)
    θ_obs = getproperty(nt.observations,
        Symbol(Octofitter.normalizename(Octofitter.likelihoodname(obs))))
    return (; nt, lnl, ctx=Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs]))
end

@testset "construction" begin
    m = hip_model()
    obs = m.obs
    @test length(Octofitter.epochs(obs)) == length(HIP_IAD.table.epoch)
    @test Octofitter.epochs(obs) == collect(Float64, HIP_IAD.table.epoch)
    @test Octofitter._refnames.(Octofitter.refspecs(obs)) == ((:A,), (:b,), ())
    @test occursin("A (blended with b) vs Barycentre",
        sprint(show, MIME"text/plain"(), obs))
    @test obs.hip_sol.plx == HIP_SOL.plx

    # Either the data or an id, not neither.
    @test_throws r"needs either `hip_id`" HipparcosIADObs(host=:A)

    # A typo'd body name is caught at model-build time.
    @test_throws r"references :nope" Octofitter.System(
        name="bad", bodies=[m.host, m.b],
        observations=[HipparcosIADObs(; host=:nope, companions=(:b,), iad=HIP_IAD,
            variables=(obs.priors, obs.derived))],
        variables=@variables begin plx = 278.0 end)

    # The default variable set is the frame-offset block plus the jitter, and
    # `iad_Δplx` is fixed: the system's own `plx` is the anchor, so leaving it
    # free would turn the one thing these data measure best into a nuisance.
    d = HipparcosIADObs(; host=:A, companions=(:b,), iad=HIP_IAD)
    for k in (:hip_iad_jitter, :iad_Δra, :iad_Δdec, :iad_Δpmra, :iad_Δpmdec)
        @test k ∈ keys(d.priors.priors)
    end
    @test :iad_Δplx ∈ keys(d.derived.variables)
    @test :iad_Δplx ∉ keys(d.priors.priors)
end

@testset "Hp flux ratios default to the bodies' own flux_Hp" begin
    # The refs + `flux_Hp` replacement for v1's mandatory positional
    # `fluxratio_hip ~ Product([...])`. Both spellings, same arithmetic.
    a = hip_model(use_body_fluxes=true)
    b = hip_model(use_body_fluxes=false)
    @test :fluxratio_hip ∉ keys(a.obs.derived.variables)
    @test :fluxratio_hip ∈ keys(b.obs.derived.variables)
    ca = hip_context(a.sys, a.obs); cb = hip_context(b.sys, b.obs)
    sa = Octofitter.simulate(a.obs, ca.ctx); sb = Octofitter.simulate(b.obs, cb.ctx)
    @test sa.resid == sb.resid
    @test sa.σ_inflation == sb.σ_inflation
    @test ca.lnl(a.sys, ca.nt) == cb.lnl(b.sys, cb.nt)
    # …and the companion's light really is doing something
    @test maximum(sa.σ_inflation) > 1.0
end

@testset "the sky-path model, against G23H's :iad_hip channel" begin
    # With no companion perturbation the two must agree exactly: G23H's
    # abscissa channel is the same five-parameter frame offset projected on
    # the same scans, and its catalog-bias correction is zero when there is
    # nothing to bias. Anchor the standalone's parallax on the Hipparcos
    # catalog value (G23H's choice) by giving the system that `plx`, and
    # recalibrate its table the way G23H does unconditionally.
    m = hip_model(mass=0.0, f_hp=0.0, plx=HIP_SOL.plx, recalibrate=true)
    c = hip_context(m.sys, m.obs)
    s = Octofitter.simulate(m.obs, c.ctx)

    cat = HIP_CAT
    gobs = G23HObs(; host=m.host, companions=(m.b,),
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=_hip_fixture("g23h-forecast.csv"), hipparcos=HIP_IAD,
        dr2_transits_catalog=_hip_fixture("g23h-dr2-transits.csv"),
        variables=@variables begin
            σ_AL = $(HIP_EVAL[:sigma_AL])
            σ_att = $(HIP_EVAL[:sigma_att])
            σ_calib = $(HIP_EVAL[:sigma_calib])
            fluxratio = (0.0,)
            fluxratio_hip = (0.0,)
            hip_iad_jitter = 0.01
            iad_Δra = 0.0
            iad_Δdec = 0.0
            iad_Δplx = 0.0
            iad_pmra = $(HIP_SOL.pm_ra)
            iad_pmdec = $(HIP_SOL.pm_de)
            σ_rv_per_transit = $(HIP_EVAL[:sigma_rv_per_transit])
            pmra = $(cat.pmra_dr3)
            pmdec = $(cat.pmdec_dr3)
            transit_priorities = $(collect(range(-1.0, 1.0, length=HIP_NPOOL)))
            transits = $(collect(1:Int(HIP_EVAL[:n_transits])))
            transits_dr2 = $(collect(1:Int(HIP_EVAL[:n_transits_dr2])))
            transits_rv = $(collect(1:Int(HIP_EVAL[:n_transits_rv])))
        end)
    gsys = Octofitter.System(name="g23hhip", bodies=[m.host, m.b], observations=[gobs],
        variables=(@variables begin
            plx = $(HIP_SOL.plx)
        end), observing_geometry=false)
    gc = hip_context(gsys, gobs)
    gs = Octofitter.simulate(gobs, gc.ctx)
    @test abs.(s.resid) == gs.iad_resid
    @test all(iszero, s.Δα_mas)
    @test all(≈(1.0), s.σ_inflation)
end

@testset "the grating response, and its wide-separation limit" begin
    # A close luminous companion modulates the abscissa and inflates σ …
    close = hip_model(a=4.7, f_hp=0.2)
    c = hip_context(close.sys, close.obs)
    s = Octofitter.simulate(close.obs, c.ctx)
    @test maximum(s.σ_inflation) > 1.0
    @test maximum(abs, s.Δα_mas) > 0

    # … and beyond the resolution scale it tapers away entirely, leaving the
    # host's reflex about the barycentre and nothing else. Same statement as
    # the G23H wide-separation test, made through the standalone type.
    wide = hip_model(a=4000.0, f_hp=0.3)
    cw = hip_context(wide.sys, wide.obs)
    sw = Octofitter.simulate(wide.obs, cw.ctx)
    hostref = Octofitter.ref(cw.ctx, wide.obs.host)
    reference = Octofitter.ref(cw.ctx, wide.obs.ref)
    n = length(wide.obs.table.epoch)
    host_along = [raoff(Octofitter.solutionat(cw.ctx, i), hostref, reference) *
                  wide.obs.table.cosϕ[i] +
                  decoff(Octofitter.solutionat(cw.ctx, i), hostref, reference) *
                  wide.obs.table.sinϕ[i] for i in 1:n]
    @test sw.Δα_mas ≈ host_along .* wide.obs.table.cosϕ rtol = 1e-12
    @test all(≈(1.0), sw.σ_inflation)

    # A massless companion contributes neither reflex nor light, whatever its
    # declared flux says.
    dark = hip_model(mass=0.0, f_hp=0.3)
    cd = hip_context(dark.sys, dark.obs)
    sd = Octofitter.simulate(dark.obs, cd.ctx)
    @test all(iszero, sd.Δα_mas)
    @test all(≈(1.0), sd.σ_inflation)
end

@testset "the likelihood, subsets and simulation" begin
    m = hip_model()
    c = hip_context(m.sys, m.obs)
    ll = c.lnl(m.sys, c.nt)
    @test isfinite(ll)

    # Rejected scans (`sres ≤ 0` in the catalog) contribute nothing.
    @test any(m.obs.table.reject)
    n_used = count(!, m.obs.table.reject)
    keep = findall(!, collect(m.obs.table.reject))
    sub = Octofitter.likeobj_from_epoch_subset(m.obs, keep)
    @test length(sub.table.epoch) == n_used
    subsys = Octofitter.System(name="hiponly", bodies=[m.host, m.b], observations=[sub],
        variables=(m.sys.priors, m.sys.derived), observing_geometry=false)
    ntsub = Octofitter.make_arr2nt(subsys)(Float64[])
    @test Octofitter.make_ln_like(subsys, ntsub)(subsys, ntsub) ≈ ll rtol = 1e-12

    # Noiseless generation puts the data on the model, so every residual is
    # zero and the regenerated likelihood beats the one it was drawn from.
    gen = Octofitter.generate_from_params(m.obs, c.ctx; add_noise=false)
    gsys = Octofitter.System(name="hiponly", bodies=[m.host, m.b], observations=[gen],
        variables=(m.sys.priors, m.sys.derived), observing_geometry=false)
    gc = hip_context(gsys, gen)
    @test maximum(abs, Octofitter.simulate(gen, gc.ctx).resid) < 1e-9
    @test gc.lnl(gsys, gc.nt) > ll

    # With noise it still evaluates, and the residuals are the right size.
    Random.seed!(7)
    gn = Octofitter.generate_from_params(m.obs, c.ctx; add_noise=true)
    gnsys = Octofitter.System(name="hiponly", bodies=[m.host, m.b], observations=[gn],
        variables=(m.sys.priors, m.sys.derived), observing_geometry=false)
    gnc = hip_context(gnsys, gn)
    r = Octofitter.simulate(gn, gnc.ctx).resid
    σ = gn.table.sres_renorm
    @test isfinite(gnc.lnl(gnsys, gnc.nt))
    @test 0.3 < std(r[.!collect(gn.table.reject)] ./ σ[.!collect(gn.table.reject)]) < 3.0
end

@testset "the likelihood is differentiable" begin
    # A Hipparcos-only fit: the parallax, the companion's mass and its
    # semi-major axis, with the frame's position and proper motion as
    # nuisances. This is the use case the type exists for.
    host = Octofitter.Body(name="A", variables=@variables begin
        mass = 0.486 - system.m_b
        flux_Hp = 1.0
    end)
    b = Octofitter.Body(name="b", about=host, variables=@variables begin
        mass = system.m_b
        M = 0.486
        flux_Hp = system.f_b
        a ~ LogUniform(1.0, 30.0)
        e = 0.23; i = 0.94; ω = 1.31; Ω = 2.44; tp = 52000.0
    end)
    obs = HipparcosIADObs(; host, companions=(b,), iad=HIP_IAD,
        variables=@variables begin
            hip_iad_jitter ~ LogUniform(0.01, 10.0)
            iad_Δra ~ Uniform(-50, 50)
            iad_Δdec ~ Uniform(-50, 50)
            iad_Δplx = 0.0
            iad_Δpmra ~ Uniform(-50, 50)
            iad_Δpmdec ~ Uniform(-50, 50)
            iad_pmra = $(HIP_SOL.pm_ra) + iad_Δpmra
            iad_pmdec = $(HIP_SOL.pm_de) + iad_Δpmdec
        end)
    sys = Octofitter.System(name="hipgrad", bodies=[host, b], observations=[obs],
        variables=(@variables begin
            m_b ~ LogUniform(1e-4, 0.2)
            f_b ~ LogUniform(1e-4, 0.5)
            # A global, not `$(HIP_SOL.plx)`: `@variables` only quasiquotes the
            # right-hand side of `=` lines, so `$` in a `~` line is a syntax
            # error. The escaped expression sees this file's scope anyway.
            plx ~ truncated(Normal(HIP_PLX_CATALOG, 1.0), lower=1.0)
        end), observing_geometry=false)

    model = Octofitter.LogDensityModel(sys; verbosity=0)

    # A fixed, physically sensible point rather than a prior draw — the
    # gradient check below needs a well-scaled step, and a random draw of
    # `iad_Δra`/`iad_Δpmra` puts the model tens of mas off its own data.
    # Order is the flat parameter vector's: system, then bodies, then
    # observations. Pinned by the `arr2nt` assertions so it cannot drift.
    θ = [0.03, 0.05, HIP_SOL.plx, 4.7, 0.5, 0.0, 0.0, 0.0, 0.0]
    nt = model.arr2nt(θ)
    @test nt.m_b == 0.03 && nt.plx == HIP_SOL.plx
    @test nt.bodies.b.a == 4.7
    @test nt.observations.Hipparcos_IAD.hip_iad_jitter == 0.5

    θt = model.link(θ)
    v, g = model.∇ℓπcallback(θt)
    @test isfinite(v)
    @test all(isfinite, g)

    # `relstep=1e-8`, not FiniteDiff's default. 145 transits at ~2 mas each
    # constrain the parallax to ~1e-4 of its prior scale, so the log-density
    # varies on that scale in the linked coordinate; the default central step
    # (∛eps ≈ 6e-6) is *inside* the curvature and mis-estimates the two stiff
    # components (plx, a) by ~1%. That is finite differencing failing, not AD:
    # the disagreement shrinks monotonically to 2e-8 as the step is reduced to
    # 1e-8, then grows again from roundoff. Do not "fix" this by loosening the
    # tolerance — a loose tolerance here would pass a genuinely wrong gradient.
    gfd = FiniteDiff.finite_difference_gradient(
        model.ℓπcallback, θt, Val{:central}; relstep=1e-8, absstep=1e-8)
    @test norm(g .- gfd) / norm(g) < 1e-6

    # The parallax is genuinely constrained: `iad_Δplx = 0` means the system's
    # own `plx` sets the abscissa's parallax signature.
    i_plx = findfirst(==(:plx), collect(keys(model.system.priors.priors)))
    @test !isnothing(i_plx)
end
