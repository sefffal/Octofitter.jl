# G23H — the joint Gaia/Hipparcos catalog observation.
#
# Everything here runs offline against the checked-in fixtures in
# `fixtures/`, which were extracted from a real G23H catalog row (HIP 1475 /
# Gaia DR3 385334230892516480), a real cached Gaia scan forecast, and the real
# Hipparcos IAD, by `fixtures/gen-g23h-fixtures.jl`. That script also produced
# the Octofitter v1 reference values the port is gated against; v1 is deleted
# on this branch, so they are recorded numbers in the same spirit as the
# other v1 regression fixtures.

using Octofitter: Octofitter, G23HObs
using Octofitter.TypedTables: Table, FlexTable
using Octofitter.CSV
using PlanetOrbits
using StaticArrays, Distributions, LinearAlgebra, Random, Test
using Octofitter.DataDeps

const G23H_FIX = joinpath(@__DIR__, "fixtures")

_g23h_fixture(name) = Table(CSV.File(joinpath(G23H_FIX, name)))

function g23h_catalog_row()
    t = _g23h_fixture("g23h-catalog-row.csv")
    cols = Tuple(Octofitter.Tables.columnnames(t))
    return NamedTuple{cols}(map(c -> getproperty(t, c)[1], cols))
end

g23h_forecast() = _g23h_fixture("g23h-forecast.csv")

function g23h_hipparcos()
    tab = _g23h_fixture("g23h-hip-iad.csv")
    s = _g23h_fixture("g23h-hip-sol.csv")
    hip_sol = NamedTuple{Tuple(Octofitter.Tables.columnnames(s))}(
        map(c -> getproperty(s, c)[1], Tuple(Octofitter.Tables.columnnames(s))))
    return (; table=tab, hip_sol)
end

g23h_dr2_sidecar() = _g23h_fixture("g23h-dr2-transits.csv")

function g23h_named_values(name)
    t = _g23h_fixture(name)
    return Dict(Symbol(n) => v for (n, v) in zip(t.name, t.value))
end

const V1_REF = g23h_named_values("g23h-v1-reference.csv")
const EVAL = g23h_named_values("g23h-eval-point.csv")

# The epoch selections the reference was evaluated at (see the generator).
const N_POOL = Int(V1_REF[:n_pool])
const TRANSITS = collect(1:Int(EVAL[:n_transits]))
const TRANSITS_DR2 = collect(1:Int(EVAL[:n_transits_dr2]))
const TRANSITS_RV = collect(1:Int(EVAL[:n_transits_rv]))
const PRIORITIES = collect(range(-1.0, 1.0, length=N_POOL))

"""
Build the v2 model that reproduces the v1 reference configuration: one
luminous companion, a parallax-only frame (v1's `Visual{KepOrbit}`, i.e. the
`absolute_orbits=false` branch), and `observing_geometry=false`, which is
v1's sky geometry exactly.

`fluxratios`/`masses` are per-companion tuples, so the same builder covers
the multi-companion divergence measurement.
"""
function g23h_model(; n_companions=1,
                    masses=(EVAL[:M_comp],),
                    fluxratios=(EVAL[:fluxratio_G],),
                    fluxratios_hip=(EVAL[:fluxratio_Hp],),
                    elements=((a=EVAL[:a], e=EVAL[:e], i=EVAL[:i],
                               ω=EVAL[:omega], Ω=EVAL[:Omega], tp=EVAL[:tp]),),
                    observing_geometry=false,
                    ueva_mode=:RUWE, include_rv=true)
    cat = g23h_catalog_row()
    M_tot = EVAL[:M_tot]
    host = Octofitter.Body(name="A", variables=@variables begin
        mass = $(M_tot - sum(masses))
    end)
    comps = map(1:n_companions) do k
        el = elements[k]
        nm = ("b", "c", "d")[k]
        Octofitter.Body(name=nm, about=host, variables=@variables begin
            mass = $(masses[k])
            # v1 drove Kepler's third law from the total system mass for every
            # companion; `M=` is the compatibility override for exactly that.
            M = $M_tot
            a = $(el.a); e = $(el.e); i = $(el.i)
            ω = $(el.ω); Ω = $(el.Ω); tp = $(el.tp)
        end)
    end

    obs = G23HObs(;
        host=host, companions=Tuple(comps),
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(),
        ueva_mode, include_rv,
        variables=@variables begin
            σ_AL = $(EVAL[:sigma_AL])
            σ_att = $(EVAL[:sigma_att])
            σ_calib = $(EVAL[:sigma_calib])
            fluxratio = $(Tuple(fluxratios))
            fluxratio_hip = $(Tuple(fluxratios_hip))
            hip_iad_jitter = $(EVAL[:hip_iad_jitter])
            iad_Δra = $(EVAL[:iad_dra])
            iad_Δdec = $(EVAL[:iad_ddec])
            iad_Δplx = $(EVAL[:iad_dplx])
            iad_pmra = $(EVAL[:iad_pmra])
            iad_pmdec = $(EVAL[:iad_pmdec])
            σ_rv_per_transit = $(EVAL[:sigma_rv_per_transit])
            # v2 reserves the system-level `pmra`/`pmdec` names for the
            # absolute frame, and a *partial* frame is rejected at build time,
            # so a parallax-only model declares the reference point's proper
            # motion on the observation. v1 read it from the system block; the
            # values are identical.
            pmra = $(cat.pmra_dr3)
            pmdec = $(cat.pmdec_dr3)
            transit_priorities = $PRIORITIES
            transits = $TRANSITS
            transits_dr2 = $TRANSITS_DR2
            transits_rv = $TRANSITS_RV
        end)

    sys = Octofitter.System(name="g23h", bodies=[host, comps...], observations=[obs],
        variables=(@variables begin
            plx = $(cat.parallax)
        end), observing_geometry=observing_geometry)
    return (; sys, obs, host, comps)
end

# One solved context, for calling `simulate` directly.
function g23h_context(sys, obs, θ=Float64[])
    nt = Octofitter.make_arr2nt(sys)(θ)
    lnl = Octofitter.make_ln_like(sys, nt)
    posys = lnl.build(nt)
    eps, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, eps; method=sys.method,
        observing_geometry=sys.observing_geometry,
        barycentric_lighttime=sys.barycentric_lighttime)
    θ_obs = getproperty(nt.observations, Symbol(Octofitter.normalizename(Octofitter.likelihoodname(obs))))
    return (; nt, lnl, ctx=Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs]))
end

@testset "construction and the epoch layout" begin
    m = g23h_model()
    obs = m.obs
    @test length(obs.gaia_table.epoch) == N_POOL
    @test obs.n_hip == Int(V1_REF[:n_hip])
    # hip transits + the Gaia pool + six catalog reference epochs
    @test length(Octofitter.epochs(obs)) == obs.n_hip + obs.n_gaia + 6
    @test obs.table.kind == [:iad_hip, :ra_hip, :dec_hip, :ra_hg, :dec_hg,
        :ra_dr2, :dec_dr2, :ra_dr32, :dec_dr32, :ra_dr3, :dec_dr3,
        :ueva_dr3, :rv_dr3]
    # `host=`/`companions=` are declared, so they show up for name validation
    @test Octofitter._refnames.(Octofitter.refspecs(obs)) == ((:A,), (:b,), ())
    # A typo'd body name is caught at model-build time, not in the sampler
    @test_throws r"references :nope" Octofitter.System(
        name="bad", bodies=[m.host, m.comps[1]],
        observations=[G23HObs(; host=:nope, companions=(:b,),
            gaia_id=g23h_catalog_row().gaia_source_id, catalog=g23h_catalog_row(),
            forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
            dr2_transits_catalog=g23h_dr2_sidecar(), variables=(obs.priors, obs.derived))],
        variables=@variables begin plx = 278.0 end)
end

@testset "single luminous companion: v1 agreement to roundoff" begin
    # The legacy per-companion coefficient (−m + f·m_host)/(M(1+f)) is
    # algebraically identical to the exact flux-weighted photocentre offset
    # when there is exactly one luminous companion, so v1 and v2 do the same
    # arithmetic and must agree to roundoff — on every channel, not just the
    # total.
    m = g23h_model()
    c = g23h_context(m.sys, m.obs)
    sim = Octofitter.simulate(m.obs, c.ctx)

    pairs = [
        (:mu_h_ra, sim.μ_h[1]), (:mu_h_dec, sim.μ_h[2]),
        (:mu_hg_ra, sim.μ_hg[1]), (:mu_hg_dec, sim.μ_hg[2]),
        (:mu_dr2_ra, sim.μ_dr2[1]), (:mu_dr2_dec, sim.μ_dr2[2]),
        (:mu_dr32_ra, sim.μ_dr32[1]), (:mu_dr32_dec, sim.μ_dr32[2]),
        (:mu_dr3_ra, sim.μ_dr3[1]), (:mu_dr3_dec, sim.μ_dr3[2]),
        (:UEVA_model, sim.UEVA_model), (:UEVA_unc, sim.UEVA_unc),
        (:mu_1_3, sim.μ_1_3), (:sample_variance, sim.sample_variance),
        (:deflation_factor_dr3, sim.deflation_factor_dr3),
        (:hip_bias_pm_sq, sim.hip_bias_pm_sq),
        (:delta_alpha_dr3, sim.Δα_dr3), (:delta_delta_dr3, sim.Δδ_dr3),
        (:delta_pmra_dr3, sim.Δpmra_dr3), (:delta_pmdec_dr3, sim.Δpmdec_dr3),
    ]
    # Measured worst relative residual over these channels: 7.5e-16, and the
    # total log-likelihood is bit-identical. The gate is set an order of
    # magnitude above that, not at "close enough".
    for (k, got) in pairs
        want = V1_REF[k]
        @test abs(got - want) <= 1e-14 * max(abs(want), 1e-12)
    end
    @test c.lnl(m.sys, c.nt) ≈ V1_REF[:ll] rtol = 1e-14
end

# ── the Hipparcos instrument response, against the legacy math ────────────
#
# Transcribed from `src/legacy/gaia-utils-fitting.jl`
# (`_simulate_skypath_hippacentre_combined!`), with the per-companion orbit
# solutions replaced by the per-epoch separations they produced. That is the
# whole content of the routine: given the same separations, host reflex,
# fluxes and scan angles it must produce the same Δν_B and f_σ.
function legacy_hippacentre!(Δα, Δδ, σ_inflation, cosϕ, sinϕ,
                             host_ra, host_dec, sep_ra, sep_dec, f, s)
    N = length(f)
    inv_2π_over_s = s / (2π)
    two_π_over_s = (2π) / s
    inv_res_mas2 = 1 / (1000 * Octofitter.HIPPARCOS_RESOLUTION_ARCSEC)^2
    for i in eachindex(cosϕ)
        Re = 1.0; Im = 0.0; f_total = 0.0
        host_along = host_ra[i] * cosϕ[i] + host_dec[i] * sinϕ[i]
        for k in 1:N
            ra_p = sep_ra[k][i]
            dec_p = sep_dec[k][i]
            ρ_pk = ra_p * cosϕ[i] + dec_p * sinϕ[i]
            α_k = exp(-(ra_p^2 + dec_p^2) * inv_res_mas2)
            ζ_k = two_π_over_s * ρ_pk
            f_k = f[k] * α_k
            sin_ζk, cos_ζk = isfinite(ζ_k) ? sincos(ζ_k) : (NaN, NaN)
            Re += f_k * cos_ζk
            Im += f_k * sin_ζk
            f_total += f_k
        end
        Δν_B = inv_2π_over_s * atan(Im, Re) + host_along
        Δα[i] += Δν_B * cosϕ[i]
        Δδ[i] += Δν_B * sinϕ[i]
        σ_inflation[i] *= (1 + f_total) / sqrt(Re^2 + Im^2)
    end
    return
end

@testset "the Hippacentre response matches the legacy math" begin
    for (nc, f) in ((1, (0.037,)), (2, (0.037, 0.011)), (3, (0.05, 0.02, 0.004)))
        m = g23h_model(n_companions=nc,
            masses=ntuple(k -> (0.03, 0.012, 0.004)[k], nc),
            fluxratios=ntuple(k -> f[k], nc),
            fluxratios_hip=ntuple(k -> f[k], nc),
            elements=((a=4.7, e=0.23, i=0.94, ω=1.31, Ω=2.44, tp=52000.0),
                (a=1.9, e=0.05, i=1.7, ω=0.4, Ω=5.1, tp=53100.0),
                (a=12.0, e=0.4, i=0.3, ω=3.0, Ω=1.0, tp=50000.0)))
        c = g23h_context(m.sys, m.obs)
        obs = m.obs
        n = obs.n_hip
        hostref = Octofitter.ref(c.ctx, obs.host)
        reference = Octofitter.ref(c.ctx, obs.ref)
        comps = map(x -> Octofitter.ref(c.ctx, x), obs.companions)

        host_ra = [raoff(Octofitter.solutionat(c.ctx, i), hostref, reference) for i in 1:n]
        host_dec = [decoff(Octofitter.solutionat(c.ctx, i), hostref, reference) for i in 1:n]
        sep_ra = [[raoff(Octofitter.solutionat(c.ctx, i), ck, hostref) for i in 1:n] for ck in comps]
        sep_dec = [[decoff(Octofitter.solutionat(c.ctx, i), ck, hostref) for i in 1:n] for ck in comps]

        Δα_l = zeros(n); Δδ_l = zeros(n); σ_l = ones(n)
        legacy_hippacentre!(Δα_l, Δδ_l, σ_l, obs.hip_table.cosϕ, obs.hip_table.sinϕ,
            host_ra, host_dec, sep_ra, sep_dec, collect(f),
            Octofitter.HIPPARCOS_GRID_STEP_ARCSEC)

        Δα_n = zeros(n); Δδ_n = zeros(n); σ_n = ones(n)
        Octofitter._hippacentre!(Δα_n, Δδ_n, σ_n, c.ctx, 1:n,
            obs.hip_table.cosϕ, obs.hip_table.sinϕ, hostref, reference,
            comps, ntuple(k -> f[k], nc), ntuple(_ -> true, nc),
            Octofitter.HIPPARCOS_GRID_STEP_ARCSEC)

        @test Δα_n ≈ Δα_l rtol = 1e-14
        @test Δδ_n ≈ Δδ_l rtol = 1e-14
        @test σ_n ≈ σ_l rtol = 1e-14
        # It is a real inflation, not a no-op that the tolerance hides
        @test maximum(σ_n) > 1.0
    end

    # Wide-separation limit: with every companion beyond the resolution scale,
    # the modulation tapers away and only the host reflex remains.
    m = g23h_model(masses=(0.03,), fluxratios=(0.3,), fluxratios_hip=(0.3,),
        elements=((a=4000.0, e=0.0, i=0.0, ω=0.0, Ω=0.0, tp=52000.0),))
    c = g23h_context(m.sys, m.obs)
    obs = m.obs
    n = obs.n_hip
    hostref = Octofitter.ref(c.ctx, obs.host)
    reference = Octofitter.ref(c.ctx, obs.ref)
    comps = map(x -> Octofitter.ref(c.ctx, x), obs.companions)
    Δα = zeros(n); Δδ = zeros(n); σ = ones(n)
    Octofitter._hippacentre!(Δα, Δδ, σ, c.ctx, 1:n, obs.hip_table.cosϕ, obs.hip_table.sinϕ,
        hostref, reference, comps, (0.3,), (true,), Octofitter.HIPPARCOS_GRID_STEP_ARCSEC)
    host_only = [raoff(Octofitter.solutionat(c.ctx, i), hostref, reference) * obs.hip_table.cosϕ[i] +
                 decoff(Octofitter.solutionat(c.ctx, i), hostref, reference) * obs.hip_table.sinϕ[i]
                 for i in 1:n]
    @test Δα ≈ host_only .* obs.hip_table.cosϕ rtol = 1e-12
    @test all(≈(1.0), σ)
end

@testset "multiple luminous companions: the documented divergence" begin
    # Legacy normalized each companion's photocentre term by its own (1 + f_k)
    # and superposed; the exact photocentre normalizes by (1 + Σ f_k). The
    # difference is real and is a correctness fix, so it is measured here
    # rather than reproduced.
    m = g23h_model(n_companions=2, masses=(0.03, 0.012),
        fluxratios=(0.037, 0.011), fluxratios_hip=(0.037, 0.011),
        elements=((a=4.7, e=0.23, i=0.94, ω=1.31, Ω=2.44, tp=52000.0),
            (a=1.9, e=0.05, i=1.7, ω=0.4, Ω=5.1, tp=53100.0)))
    c = g23h_context(m.sys, m.obs)
    obs = m.obs
    hostref = Octofitter.ref(c.ctx, obs.host)
    reference = Octofitter.ref(c.ctx, obs.ref)
    comps = map(x -> Octofitter.ref(c.ctx, x), obs.companions)
    M_tot = EVAL[:M_tot]
    masses = (0.03, 0.012)
    f = (0.037, 0.011)

    n = obs.n_gaia
    exact_ra = zeros(n); exact_dec = zeros(n)
    legacy_ra = zeros(n); legacy_dec = zeros(n)
    photo = PlanetOrbits.photocentre(SVector(1.0, f[1], f[2]))
    for j in 1:n
        sol = Octofitter.solutionat(c.ctx, obs.n_hip + j)
        exact_ra[j] = raoff(sol, photo, reference)
        exact_dec[j] = decoff(sol, photo, reference)
        # v1: Σ_k raoff(sol_k) · (−m_k + f_k·(M − m_k)) / (M(1 + f_k))
        for k in 1:2
            coeff = (-masses[k] + f[k] * (M_tot - masses[k])) / (M_tot * (1 + f[k]))
            legacy_ra[j] += raoff(sol, comps[k], hostref) * coeff
            legacy_dec[j] += decoff(sol, comps[k], hostref) * coeff
        end
    end
    Δ = maximum(hypot.(exact_ra .- legacy_ra, exact_dec .- legacy_dec))
    scale = maximum(hypot.(exact_ra, exact_dec))

    # …and what that does to the quantity actually compared against the
    # catalog: the five-parameter refit's proper motion.
    A = obs.A_prepared_5_dr3
    cϕ = obs.gaia_table.cosϕ
    sϕ = obs.gaia_table.sinϕ
    pe = Octofitter.fit_5param_prepared(A, cϕ, sϕ, exact_ra, exact_dec).parameters
    pl = Octofitter.fit_5param_prepared(A, cϕ, sϕ, legacy_ra, legacy_dec).parameters
    Δpm = hypot(pe[3] - pl[3], pe[4] - pl[4])
    @info("multi-companion photocentre divergence (2 luminous companions)",
        max_skypath_mas = Δ, max_signal_mas = scale, relative = Δ / scale,
        delta_pm_mas_per_yr = Δpm,
        catalog_sigma_pm_dr3 = g23h_catalog_row().pmra_dr3_error)
    # It is not zero (the two really differ) and it is not absurd (both are
    # the same physics to first order in f).
    @test Δ > 1e-4
    @test Δ / scale < 1.0
    @test Δpm > 0
end

@testset "the flux-ratio contract" begin
    T = Float64
    # correct length: accepted, in `companions=` order
    @test Octofitter._g23h_fluxratios((; fluxratio=(0.1, 0.2)), :fluxratio, (true, true), T) == (0.1, 0.2)
    @test Octofitter._g23h_fluxratios((; fluxratio=SVector(0.1, 0.2)), :fluxratio, (true, true), T) == (0.1, 0.2)
    @test Octofitter._g23h_fluxratios((; fluxratio=[0.1, 0.2]), :fluxratio, (true, true), T) == (0.1, 0.2)
    # an inactive companion is dark, whatever the variable says
    @test Octofitter._g23h_fluxratios((; fluxratio=(0.1, 0.2)), :fluxratio, (true, false), T) == (0.1, 0.0)
    # absent entirely: all dark
    @test Octofitter._g23h_fluxratios((;), :fluxratio, (true, true), T) == (0.0, 0.0)
    # a scalar is unambiguous only for one companion
    @test Octofitter._g23h_fluxratios((; fluxratio=0.3), :fluxratio, (true,), T) == (0.3,)
    @test_throws r"declares 2 companions" Octofitter._g23h_fluxratios(
        (; fluxratio=0.3), :fluxratio, (true, true), T)
    @test_throws r"received 3 flux ratios" Octofitter._g23h_fluxratios(
        (; fluxratio=(0.1, 0.2, 0.3)), :fluxratio, (true, true), T)
end

@testset "the all-inactive fast path" begin
    # Every companion massless: the Gaia perturbations vanish, the Hipparcos
    # LSQ collapses to the cached projection of the constant residuals, and
    # the RV sample variance is zero. The likelihood must still be finite —
    # this is ~half the prior volume in a variable-companion-count model.
    m = g23h_model(masses=(0.0,))
    c = g23h_context(m.sys, m.obs)
    sim = Octofitter.simulate(m.obs, c.ctx)
    @test sim.hip_bias_pm_sq == 0.0
    @test sim.sample_variance == 0.0
    @test all(iszero, sim.Δα_mas_hip)
    @test all(≈(1.0), sim.σ_inflation_hip)
    @test sim.Δpmra_dr3 == 0.0
    ll = c.lnl(m.sys, c.nt)
    @test isfinite(ll)

    # …and it agrees with a companion whose mass and flux are merely tiny,
    # which is the check that the fast path is a shortcut and not a
    # different model.
    m2 = g23h_model(masses=(1e-14,), fluxratios=(0.0,), fluxratios_hip=(0.0,))
    c2 = g23h_context(m2.sys, m2.obs)
    @test c2.lnl(m2.sys, c2.nt) ≈ ll rtol = 1e-6
end

@testset "kind-based epoch subsets" begin
    m = g23h_model()
    obs = m.obs
    kinds = collect(obs.table.kind)
    # Drop the Hipparcos abscissa channel, as the production driver does
    keep = findall(k -> k !== :iad_hip, kinds)
    sub = Octofitter.likeobj_from_epoch_subset(obs, keep)
    @test :iad_hip ∉ sub.table.kind
    # …and its six nuisance parameters go with it
    for k in (:hip_iad_jitter, :iad_Δra, :iad_Δdec, :iad_Δplx, :iad_Δpmra, :iad_Δpmdec)
        @test k ∉ keys(sub.priors.priors)
        @test k ∉ keys(sub.derived.variables)
    end
    @test :iad_pmra ∉ keys(sub.derived.variables)

    # Drop every Hipparcos channel: the catalog distributions go with them
    keep2 = findall(k -> k ∉ (:iad_hip, :ra_hip, :dec_hip, :ra_hg, :dec_hg), kinds)
    sub2 = Octofitter.likeobj_from_epoch_subset(obs, keep2)
    @test isnothing(sub2.catalog.dist_hip)
    @test isnothing(sub2.catalog.dist_hg)

    # A Gaia-only subset still evaluates
    sys2 = Octofitter.System(name="g23hsub", bodies=[m.host, m.comps...],
        observations=[sub2], variables=(@variables begin
            plx = $(g23h_catalog_row().parallax)
        end), observing_geometry=false)
    nt2 = Octofitter.make_arr2nt(sys2)(Float64[])
    @test isfinite(Octofitter.make_ln_like(sys2, nt2)(sys2, nt2))
end

@testset "ueva_mode = :none" begin
    m = g23h_model(ueva_mode=:none)
    @test :ueva_dr3 ∉ m.obs.table.kind
    c = g23h_context(m.sys, m.obs)
    sim = Octofitter.simulate(m.obs, c.ctx)
    # No UEVA datum means no UEVA-driven deflation of the DR3 block.
    @test sim.deflation_factor_dr3 == 1.0
    @test isfinite(c.lnl(m.sys, c.nt))
    @test_throws r"must be :RUWE, :EAN or :none" g23h_model(ueva_mode=:bogus)
end

@testset "generate_from_params round-trips" begin
    m = g23h_model()
    c = g23h_context(m.sys, m.obs)
    obs2 = Octofitter.generate_from_params(m.obs, c.ctx; add_noise=false)
    @test obs2 isa G23HObs
    @test obs2.table.kind == m.obs.table.kind
    sys2 = Octofitter.System(name="g23h", bodies=[m.host, m.comps...], observations=[obs2],
        variables=(@variables begin
            plx = $(g23h_catalog_row().parallax)
        end), observing_geometry=false)
    nt2 = Octofitter.make_arr2nt(sys2)(Float64[])
    c2 = g23h_context(sys2, obs2)
    sim2 = Octofitter.simulate(obs2, c2.ctx)
    # Noiseless generation puts the catalog on the model, so the catalog block
    # of the likelihood is evaluated at zero residual: every modelled proper
    # motion equals the regenerated catalog value it is compared against.
    cat2 = obs2.catalog
    @test sim2.μ_dr3[1] ≈ cat2.pmra_dr3 rtol = 1e-10
    @test sim2.μ_dr3[2] ≈ cat2.pmdec_dr3 rtol = 1e-10
    @test sim2.μ_dr2[1] ≈ cat2.pmra_dr2 rtol = 1e-10
    @test sim2.μ_h[1] ≈ cat2.pmra_hip rtol = 1e-10
    @test sim2.μ_hg[1] ≈ cat2.pmra_hg rtol = 1e-10
    # …and the regenerated likelihood beats the original one it was drawn from
    @test c2.lnl(sys2, nt2) > c.lnl(m.sys, c.nt)
end

@testset "the hot path is Bumper-only apart from the two least-squares solves" begin
    # `ln_like` carves every per-sample buffer out of the context's arena. The
    # one exception is the DR2/DR3 five-parameter refits: their design matrix
    # is a *sampled* row subset, so it cannot be pre-factorized, and
    # `fit_5param_prepared` falls through to `Base.\`, whose LAPACK QR
    # allocates a workspace. Measure that reference here and assert the rest
    # of the likelihood adds nothing to it — the gate that catches a new
    # allocation, which is what this is for.
    #
    # (The Hipparcos side has a fixed design matrix and per-epoch σ, so it
    # uses the cached pseudo-inverse and allocates nothing.)
    m1 = g23h_model()
    c1 = g23h_context(m1.sys, m1.obs)
    m3 = g23h_model(n_companions=3, masses=(0.03, 0.012, 0.004),
        fluxratios=(0.037, 0.011, 0.004), fluxratios_hip=(0.021, 0.008, 0.002),
        elements=((a=4.7, e=0.23, i=0.94, ω=1.31, Ω=2.44, tp=52000.0),
            (a=1.9, e=0.05, i=1.7, ω=0.4, Ω=5.1, tp=53100.0),
            (a=12.0, e=0.4, i=0.3, ω=3.0, Ω=1.0, tp=50000.0)))
    c3 = g23h_context(m3.sys, m3.obs)

    obs = m1.obs
    n3 = length(TRANSITS)
    n2 = length(TRANSITS_DR2)
    A3 = obs.A_prepared_5_dr3[TRANSITS, :]
    A2 = obs.A_prepared_5_dr2[TRANSITS_DR2, :]
    cϕ3 = obs.gaia_table.cosϕ[TRANSITS]; sϕ3 = obs.gaia_table.sinϕ[TRANSITS]
    cϕ2 = obs.gaia_table.cosϕ[TRANSITS_DR2]; sϕ2 = obs.gaia_table.sinϕ[TRANSITS_DR2]
    z3 = zeros(n3); z2 = zeros(n2)
    function lsq_reference()
        Octofitter.fit_5param_prepared(A3, cϕ3, sϕ3, z3, z3, 0.0, 0.3; include_chi2=Val(true))
        Octofitter.fit_5param_prepared(A2, cϕ2, sϕ2, z2, z2)
        return nothing
    end
    lsq_reference()
    ref = @allocated lsq_reference()

    f1() = Octofitter.ln_like(m1.obs, c1.ctx)
    f3() = Octofitter.ln_like(m3.obs, c3.ctx)
    f1(); f3()
    a1 = @allocated f1()
    a3 = @allocated f3()
    @info "G23H ln_like allocations" lsq_workspace = ref one_companion = a1 three_companions = a3 n_hip = m1.obs.n_hip n_gaia = m1.obs.n_gaia
    if Base.JLOptions().check_bounds != 1
        # `Pkg.test` forces --check-bounds=yes, whose throw branches spill
        # stack buffers to the heap; the gate is meaningful only at default
        # bounds (§12).
        @test a1 <= ref + 1024
        # Nothing scales with the number of companions.
        @test a3 <= a1 + 1024
    end
end

@testset "the likelihood is differentiable" begin
    # The whole point of the port is that this runs inside a sampler.
    #
    # NB `include_rv=false`. `Distributions.NoncentralChisq`'s log-density
    # does not accept ForwardDiff `Dual`s, so with the :rv_dr3 channel active
    # the *whole* log-density is -Inf under AD while the primal path stays
    # finite — see the next testset, and the port notes. That is inherited
    # behaviour, not something this port introduced.
    cat = g23h_catalog_row()
    host = Octofitter.Body(name="A", variables=@variables begin
        mass = system.M_tot - system.m_b
    end)
    b = Octofitter.Body(name="b", about=host, variables=@variables begin
        mass = system.m_b
        M = system.M_tot
        a ~ LogUniform(1.0, 30.0)
        e = 0.23
        i = 0.94
        ω = 1.31
        Ω = 2.44
        tp = 52000.0
    end)
    obs = G23HObs(; host=host, companions=(b,),
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(), include_rv=false,
        variables=@variables begin
            σ_AL ~ truncated(Normal(cat.sig_AL, cat.sig_AL_sigma), lower=1e-4, upper=10.0)
            σ_att = $(EVAL[:sigma_att])
            σ_calib = $(EVAL[:sigma_calib])
            fluxratio = (system.f_b,)
            fluxratio_hip = (system.f_b,)
            hip_iad_jitter = $(EVAL[:hip_iad_jitter])
            iad_Δra = 0.0
            iad_Δdec = 0.0
            iad_Δplx = 0.0
            iad_pmra = $(EVAL[:iad_pmra])
            iad_pmdec = $(EVAL[:iad_pmdec])
            σ_rv_per_transit = $(EVAL[:sigma_rv_per_transit])
            pmra = $(cat.pmra_dr3)
            pmdec = $(cat.pmdec_dr3)
            transit_priorities = $PRIORITIES
            transits = $TRANSITS
            transits_dr2 = $TRANSITS_DR2
            transits_rv = $TRANSITS_RV
        end)
    sys = Octofitter.System(name="g23hgrad", bodies=[host, b], observations=[obs],
        variables=(@variables begin
            M_tot = 0.486
            m_b ~ LogUniform(1e-4, 0.2)
            f_b ~ LogUniform(1e-4, 0.5)
            plx = $(cat.parallax)
        end), observing_geometry=false)
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θt = model.link(Octofitter.sample_priors(Random.Xoshiro(11), sys))
    v, g = model.∇ℓπcallback(θt)
    @test isfinite(v)
    @test all(isfinite, g)
    gfd = FiniteDiff.finite_difference_gradient(model.ℓπcallback, θt)
    @test norm(g .- gfd) / norm(g) < 1e-4
end

@testset "the Gaia RV channel is not differentiable (inherited)" begin
    # Recorded deliberately: `Distributions.NoncentralChisq` cannot evaluate
    # its log-density at a `ForwardDiff.Dual`, so the :rv_dr3 term falls into
    # the rejection branch and takes the whole log-density with it. The
    # primal path is unaffected, which is why a derivative-free sampler (the
    # production driver uses Pigeons) never sees this. If this test starts
    # failing because the value came out finite, `NoncentralChisq` grew AD
    # support and the warning in `g23h.jl` can go.
    m = g23h_model()
    c = g23h_context(m.sys, m.obs)
    @test :rv_dr3 ∈ m.obs.table.kind
    @test isfinite(c.lnl(m.sys, c.nt))
    f(x) = Octofitter.ln_like(m.obs, Octofitter.ObsContext(
        c.ctx.θ_system, merge(c.ctx.θ_obs, (; σ_rv_per_transit=x[1])),
        c.ctx.system, c.ctx.traj, c.ctx.epoch_index, c.ctx.buf))
    @test isfinite(f([EVAL[:sigma_rv_per_transit]]))
    g = (@test_logs (:warn, r"Gaia RV variability channel") match_mode = :any Octofitter.ForwardDiff.gradient(f, [EVAL[:sigma_rv_per_transit]]))
    @test !all(isfinite, g) || iszero(g)
end

@testset "end to end: a lumcomp-shaped model samples" begin
    # A miniature of the production model this observation exists for: a host
    # plus two companions, a *sampled* resolved-flag latent per companion at
    # system level gating that companion out of the Gaia photocentre, and the
    # resulting effective flux ratios read by the observation.
    #
    # The point is the wiring, not the science: `f_b` depends on `b.a`, which
    # makes it a *deferred* system variable (it mentions a body), and an
    # observation reading a deferred system variable is exactly the case the
    # three-tier design turns on — the blending state never round-trips
    # through a body's own `flux_G`, which would be the cycle codegen rejects.
    cat = g23h_catalog_row()
    host = Octofitter.Body(name="A", variables=@variables begin
        mass = system.M_tot - system.m_b - system.m_c
    end)
    b = Octofitter.Body(name="b", about=host, variables=@variables begin
        mass = system.m_b
        M = system.M_tot
        a ~ LogUniform(0.5, 20.0)
        e = 0.1
        i ~ Sine()
        ω = 1.0
        Ω = 0.5
        tp = 52000.0
    end)
    c = Octofitter.Body(name="c", about=host, variables=@variables begin
        mass = system.m_c
        M = system.M_tot
        a ~ LogUniform(0.5, 20.0)
        e = 0.05
        i ~ Sine()
        ω = 2.0
        Ω = 3.0
        tp = 51000.0
    end)

    obs = G23HObs(; host=host, companions=(b, c),
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(),
        # The RV channel is not differentiable; see the testset above.
        include_rv=false,
        variables=@variables begin
            σ_AL = $(EVAL[:sigma_AL])
            σ_att = $(EVAL[:sigma_att])
            σ_calib = $(EVAL[:sigma_calib])
            # Tier 2: the observation's own effective flux ratios, one per
            # declared companion, read from (deferred) system variables.
            fluxratio = (system.f_b, system.f_c)
            fluxratio_hip = (system.fh_b, system.fh_c)
            hip_iad_jitter ~ LogUniform(0.01, 10.0)
            iad_Δra = 0.0
            iad_Δdec = 0.0
            iad_Δplx = 0.0
            iad_pmra = $(EVAL[:iad_pmra])
            iad_pmdec = $(EVAL[:iad_pmdec])
            pmra ~ Uniform(cat.pmra_dr3 - 20, cat.pmra_dr3 + 20)
            pmdec ~ Uniform(cat.pmdec_dr3 - 20, cat.pmdec_dr3 + 20)
            transit_priorities = $PRIORITIES
            transits = $TRANSITS
            transits_dr2 = $TRANSITS_DR2
        end)

    sys = Octofitter.System(name="lumcomp", bodies=[host, b, c], observations=[obs],
        variables=(@variables begin
            M_tot = 0.486
            m_b ~ LogUniform(1e-3, 0.1)
            m_c ~ LogUniform(1e-3, 0.1)
            plx = $(cat.parallax)
            # Sampled resolved-flag gating, as continuous latents. `b.a` makes
            # these deferred lines.
            u_res_b ~ Uniform(0, 1)
            u_res_c ~ Uniform(0, 1)
            resolved_b = u_res_b < min(b.a / 40, 0.9)
            resolved_c = u_res_c < min(c.a / 40, 0.9)
            f_b = resolved_b ? 0.0 : 4 * m_b
            f_c = resolved_c ? 0.0 : 4 * m_c
            fh_b = 3 * m_b
            fh_c = 3 * m_c
        end), observing_geometry=false)

    @test :resolved_b ∈ sys.deferred
    @test :f_b ∈ sys.deferred

    model = Octofitter.LogDensityModel(sys; verbosity=0)
    @test model.D == 11
    θ = Octofitter.sample_priors(Random.Xoshiro(5), sys)
    nt = model.arr2nt(θ)
    @test length(nt.observations.G23H.fluxratio) == 2
    lnl = Octofitter.make_ln_like(sys, nt)
    @test isfinite(lnl(sys, nt))

    # generate_from_params round-trips through a full model rebuild
    c_ = g23h_context(sys, obs, θ)
    obs_gen = Octofitter.generate_from_params(obs, c_.ctx; add_noise=true)
    sys_gen = Octofitter.System(name="lumcomp", bodies=[host, b, c], observations=[obs_gen],
        variables=(sys.priors, sys.derived), observing_geometry=false)
    @test isfinite(Octofitter.make_ln_like(sys_gen, nt)(sys_gen, nt))

    # …and a short chain completes. Expect divergences: `resolved_b` is a step
    # function of a sampled latent, so the log-density is genuinely
    # discontinuous and HMC has nothing to follow across the jump. That is a
    # property of the model, not of the likelihood — the production driver
    # samples this shape with Pigeons. The assertion is that the machinery
    # runs end to end and produces finite values.
    chain = octofit(Random.Xoshiro(2), model; adaptation=60, iterations=60, verbosity=0)
    @test size(chain, 1) == 60
    @test haskey(chain, :b_a)
    @test haskey(chain, :G23H_pmra)
    @test all(isfinite, chain[:b_a])
end

@testset "the Hipparcos IAD loader reproduces v1" begin
    # Opt-in: the reconstruction needs the `Hipparcos_IAD` (~330 MB) and
    # `DE440_Ephemeris` (~130 MB) data dependencies, so it runs only where
    # they are already present — never as a download in CI. The fixture it
    # compares against is v1's own output for the same star, so this is the
    # gate on the loader itself; every other test here feeds the fixture
    # straight in and therefore cannot see a loader regression.
    #
    # Measured 2026-08-02: bit-identical on all 20 columns and the
    # five-parameter header.
    have_data = try
        isdir(joinpath(DataDeps.@datadep_str("Hipparcos_IAD"), "ResRec_JavaTool_2014")) &&
            isfile(joinpath(DataDeps.@datadep_str("DE440_Ephemeris"), "de440.bsp"))
    catch
        false
    end
    if !have_data
        @info "skipping the Hipparcos IAD loader gate: the data dependencies are not present"
    else
        ref = _g23h_fixture("g23h-hip-iad.csv")
        got = Octofitter.hipparcos_iad(hip_id=1475)
        @test length(got.table.epoch) == length(ref.epoch)
        for c in Octofitter.Tables.columnnames(ref)
            @test collect(getproperty(got.table, c)) == collect(getproperty(ref, c))
        end
        sref = _g23h_fixture("g23h-hip-sol.csv")
        for c in Octofitter.Tables.columnnames(sref)
            @test getproperty(got.hip_sol, c) == getproperty(sref, c)[1]
        end
    end
end
