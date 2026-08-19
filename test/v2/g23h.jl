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

`body_fluxes=true` moves the flux ratios out of the observation and onto the
bodies as `flux_G`/`flux_Hp` variables (the host at 1.0, so every other body's
flux *is* its contrast ratio). The two spellings must give the same numbers;
that equivalence is what the "flux ratios default to the bodies' own fluxes"
testset checks.
"""
function g23h_model(; n_companions=1,
                    masses=(EVAL[:M_comp],),
                    fluxratios=(EVAL[:fluxratio_G],),
                    fluxratios_hip=(EVAL[:fluxratio_Hp],),
                    elements=((a=EVAL[:a], e=EVAL[:e], i=EVAL[:i],
                               ω=EVAL[:omega], Ω=EVAL[:Omega], tp=EVAL[:tp]),),
                    observing_geometry=false,
                    body_fluxes=false, obs_ratios=nothing, channels=nothing,
                    body_fluxratios=nothing, body_fluxratios_hip=nothing,
                    ueva_mode=:RUWE, include_rv=true)
    # `body_fluxes` puts the ratios on the bodies; `obs_ratios` puts them on the
    # observation. Setting both is the override case (an explicit vector wins).
    obs_ratios = isnothing(obs_ratios) ? !body_fluxes : obs_ratios
    body_fluxratios = isnothing(body_fluxratios) ? fluxratios : body_fluxratios
    body_fluxratios_hip = isnothing(body_fluxratios_hip) ? fluxratios_hip : body_fluxratios_hip
    cat = g23h_catalog_row()
    M_tot = EVAL[:M_tot]
    hostvars = @variables begin
        mass = $(M_tot - sum(masses))
    end
    if body_fluxes
        hostvars = vcat(hostvars, @variables begin
            flux_G = 1.0
            flux_Hp = 1.0
        end)
    end
    host = Octofitter.Body(name="A", variables=hostvars)
    comps = map(1:n_companions) do k
        el = elements[k]
        nm = ("b", "c", "d")[k]
        cvars = @variables begin
            mass = $(masses[k])
            # v1 drove Kepler's third law from the total system mass for every
            # companion; `M=` is the compatibility override for exactly that.
            M = $M_tot
            a = $(el.a); e = $(el.e); i = $(el.i)
            ω = $(el.ω); Ω = $(el.Ω); tp = $(el.tp)
        end
        if body_fluxes
            cvars = vcat(cvars, @variables begin
                flux_G = $(body_fluxratios[k])
                flux_Hp = $(body_fluxratios_hip[k])
            end)
        end
        Octofitter.Body(name=nm, about=host, variables=cvars)
    end

    obsvars = @variables begin
        σ_AL = $(EVAL[:sigma_AL])
        σ_att = $(EVAL[:sigma_att])
        σ_calib = $(EVAL[:sigma_calib])
    end
    if obs_ratios
        obsvars = vcat(obsvars, @variables begin
            fluxratio = $(Tuple(fluxratios))
            fluxratio_hip = $(Tuple(fluxratios_hip))
        end)
    end
    obs = G23HObs(;
        host=host, companions=Tuple(comps),
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(),
        ueva_mode, include_rv, channels,
        variables=vcat(obsvars, @variables begin
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
        end))

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

# Override the three UEVA σ nuisances inside this observation's own namespace,
# leaving the rest of the parameter tree untouched. Used to show they are inert
# under `ueva_mode = :none`.
function _g23h_with_sigmas(nt, obs, v)
    key = Symbol(Octofitter.normalizename(Octofitter.likelihoodname(obs)))
    θ_obs = merge(getproperty(nt.observations, key), (; σ_AL=v, σ_att=v, σ_calib=v))
    return merge(nt, (; observations=merge(nt.observations, NamedTuple{(key,)}((θ_obs,)))))
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

@testset "the coupled block's layout is derived, not hardcoded" begin
    # The block used to be assembled with literal index arithmetic — `Val(11)`,
    # `Σ_full[5:6, 9:10]`, `mask[11]` — in five places that had to be edited in
    # step. It is now derived from `_g23h_blocks`. This pins the derivation to
    # the layout that arithmetic encoded, so an edit to the block list that
    # silently reorders or resizes a channel is caught here rather than as a
    # wrong likelihood.
    @test Octofitter._g23h_channel_kinds == (
        :ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr2, :dec_dr2,
        :ra_dr32, :dec_dr32, :ra_dr3, :dec_dr3, :ueva_dr3)
    @test Octofitter._g23h_nchan == 11
    @test Octofitter._g23h_block_offset == (0, 2, 4, 6, 8, 10)
    @test Octofitter._g23h_offsetof(:dr2) == 4    # was the literal `5:6`
    @test Octofitter._g23h_offsetof(:dr3) == 8    # was the literal `9:10`
    # Offsets and widths agree: each sub-block starts where the previous ended.
    let o = 0
        for b in Octofitter._g23h_blocks
            @test Octofitter._g23h_offsetof(b.name) == o
            o += length(b.kinds)
        end
        @test o == Octofitter._g23h_nchan
    end
    # Every coupled channel is reachable through `channels=`, and the two
    # separate (uncorrelated) terms are deliberately not in the block.
    @test all(k -> k ∈ Octofitter._g23h_all_kinds, Octofitter._g23h_channel_kinds)
    @test :iad_hip ∉ Octofitter._g23h_channel_kinds
    @test :rv_dr3 ∉ Octofitter._g23h_channel_kinds
end

@testset "a Gaia-only source constructs (no Hipparcos row at all)" begin
    # v8 guarded each Hipparcos epoch conversion individually; folding all
    # eleven into one `_e(y) = years2mjd(y)` helper during the v9 port dropped
    # the guard, and `years2mjd(NaN)` throws `InexactError: Int64(NaN)` from
    # `Dates.Date` — before the Hipparcos rows are dropped. That made G23HObs
    # unconstructible for *any* source with no Hipparcos entry, which is every
    # faint secondary in a wide pair (ups And B, and the reason multi-source
    # work needs this).
    #
    # `hip_id=NaN` alone does not reproduce it: the fixture row is a real
    # Hipparcos star, so its epoch columns stay finite and never reach the
    # throwing conversion. A genuinely Gaia-only catalog row has both.
    cat = merge(g23h_catalog_row(), (; hip_id=NaN,
        epoch_ra_hip=NaN, epoch_dec_hip=NaN, epoch_ra_hg=NaN, epoch_dec_hg=NaN))
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.29
        flux_G = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 0.05
        flux_G = 0.03
        a = 5.0; e = 0.2; i = 1.0; ω = 1.3; Ω = 2.4; tp = 52000.0
    end)
    obs = G23HObs(; host=A, companions=(b,), ref=Barycentre,
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=g23h_forecast(), hipparcos=nothing,
        dr2_transits_catalog=g23h_dr2_sidecar(),
        variables=@variables begin
            σ_AL = $(EVAL[:sigma_AL])
            σ_att = $(EVAL[:sigma_att])
            σ_calib = $(EVAL[:sigma_calib])
            σ_rv_per_transit = $(EVAL[:sigma_rv_per_transit])
            transit_priorities = $PRIORITIES
            transits = $TRANSITS
            transits_dr2 = $TRANSITS_DR2
            transits_rv = $TRANSITS_RV
        end)

    # Every Hipparcos and Hipparcos–Gaia channel is gone, the Gaia ones remain.
    @test Tuple(obs.table.kind) ==
          (:ra_dr2, :dec_dr2, :ra_dr32, :dec_dr32, :ra_dr3, :dec_dr3, :ueva_dr3, :rv_dr3)
    # No NaN epoch survives into a *kept* row — the guard drops the value, and
    # the row it belonged to is dropped too.
    @test all(isfinite, obs.table.epoch)
    @test all(isfinite, obs.table.start_epoch)
    @test all(isfinite, obs.table.stop_epoch)
    # `_g23h_restrict` reads "no Hipparcos channels" as "no Hipparcos", so the
    # simulator's `isnothing(dist_hip)` branch is the one that runs.
    @test isnothing(obs.catalog.dist_hip)
    @test isnothing(obs.catalog.dist_hg)
    @test obs.n_hip == 0

    # ...and it evaluates, which is the part §2.1 stopped short of auditing.
    sys = Octofitter.System(name="gaiaonly", bodies=[A, b], observations=[obs],
        variables=(@variables begin
            plx = $(cat.parallax)
            ra = $(cat.ra)
            dec = $(cat.dec)
            pmra = $(cat.pmra_dr3)
            pmdec = $(cat.pmdec_dr3)
            rv = 0.0
            ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
        end), observing_geometry=false)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    lnl = Octofitter.make_ln_like(sys, nt)
    @test isfinite(lnl(sys, nt))

    # The coupled block really is the Gaia-only subset, assembled at the
    # offsets the *full* layout defines — the case the fixed-width assembly
    # could only reach by masking rows it had already built from NaNs.
    ctx = g23h_context(sys, obs).ctx
    mom = Octofitter._g23h_catalog_moments(obs, ctx)
    @test mom.labels == [:ra_dr2, :dec_dr2, :ra_dr32, :dec_dr32, :ra_dr3, :dec_dr3, :ueva_dr3]
    @test all(isfinite, mom.sigma)
    @test all(isfinite, mom.Σ)
    # …and the mask bookkeeping agrees with the layout: seven of the eleven
    # coupled channels, at the offsets the full block defines, with the four
    # Hipparcos rows absent rather than present-and-zero.
    @test length(mom.labels) == length(mom.sigma) == size(mom.Σ, 1) == 7
    @test !any(k -> k in (:ra_hip, :dec_hip, :ra_hg, :dec_hg), mom.labels)

    # §11: `generate_from_params` has its own `has_hip` branch, never
    # exercised for a Gaia-only source. It must produce finite synthetic data
    # for every remaining channel, and — the actual hazard — data the
    # likelihood then scores at its maximum, i.e. the model and the generator
    # must agree about which corrections apply.
    gen = Octofitter.generate_from_params(obs, ctx; add_noise=false)
    @test Tuple(gen.table.kind) == Tuple(obs.table.kind)
    # `:ueva_dr3` and `:rv_dr3` carry NaN in `pm`/`σ_pm` by design — they are a
    # scatter statistic and an RV-variability statistic, read from the catalog
    # row rather than from these columns — so the assertion is that the
    # generated table has NaN in exactly the places the original does, not
    # that it is finite everywhere.
    @test map(isfinite, gen.table.pm) == map(isfinite, obs.table.pm)
    @test map(isfinite, gen.table.σ_pm) == map(isfinite, obs.table.σ_pm)
    @test all(isfinite, gen.table.pm[map(k -> k ∉ (:ueva_dr3, :rv_dr3), gen.table.kind)])
    sys_gen = Octofitter.System(name="gaiaonly-gen", bodies=[A, b], observations=[gen],
        variables=(sys.priors, sys.derived), observing_geometry=false)
    nt_gen = Octofitter.make_arr2nt(sys_gen)(Float64[])
    lnl_gen = Octofitter.make_ln_like(sys_gen, nt_gen)
    # Noiseless data generated at θ is scored at θ: every pull is zero, so the
    # log-likelihood is the normalization alone and strictly beats the value
    # at the same model against the *real* catalog row.
    @test isfinite(lnl_gen(sys_gen, nt_gen))
    @test lnl_gen(sys_gen, nt_gen) > lnl(sys, nt)
end

@testset "the perspective-acceleration guard is finiteness, not has_hip" begin
    # §11. `μ_hg` used to pick up `cat.nonlinear_dpm*` under
    # `absolute && !isnothing(dist_hip)`. The two rows anyone had checked
    # happened to carry finite nonlinear terms exactly where they had
    # Hipparcos, so a row with a Hipparcos entry and null nonlinear terms
    # would have added NaN to the model proper motion and NaN'd the whole log
    # density, with nothing in the message to say why.
    cat = g23h_catalog_row()
    @test Octofitter._g23h_nonlinear(cat, true, true) ==
          (SVector(cat.nonlinear_dpmra, cat.nonlinear_dpmdec),
           2 * SVector(cat.nonlinear_dpmra, cat.nonlinear_dpmdec))
    # Not applied without an absolute frame, nor without Hipparcos…
    @test all(iszero, Octofitter._g23h_nonlinear(cat, false, true))
    @test all(iszero, Octofitter._g23h_nonlinear(cat, true, false))
    # …and not applied — rather than propagated — when the terms are null.
    bad = merge(cat, (; nonlinear_dpmra=NaN, nonlinear_dpmdec=NaN))
    @test all(iszero, Octofitter._g23h_nonlinear(bad, true, true))
    @test all(isfinite, Octofitter._g23h_nonlinear(bad, true, true)[1])

    # A row like that is rejected at construction, where it can be named,
    # rather than silently losing a correction the model is entitled to.
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 0.01
        a = 5.0; e = 0.2; i = 1.0; ω = 1.3; Ω = 2.4; tp = 52000.0
    end)
    mk(c; kw...) = G23HObs(; host=A, companions=(b,), gaia_id=c.gaia_source_id,
        catalog=c, forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(), kw...)
    @test_throws r"nonlinear_dpmra" mk(bad)
    # …unless the channels that would use it are not in the model at all.
    @test mk(bad; channels=(:ra_dr3, :dec_dr3, :ueva_dr3)) isa G23HObs
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
        (:mu_1_3, sim.μ_1_3),
        (:deflation_factor_dr3, sim.deflation_factor_dr3),
        (:hip_bias_pm_sq, sim.hip_bias_pm_sq),
        (:delta_alpha_dr3, sim.Δα_dr3), (:delta_delta_dr3, sim.Δδ_dr3),
        (:delta_pmra_dr3, sim.Δpmra_dr3), (:delta_pmdec_dr3, sim.Δpmdec_dr3),
    ]
    # Measured worst relative residual over these channels: 7.5e-16. The gate
    # is set an order of magnitude above that, not at "close enough".
    for (k, got) in pairs
        want = V1_REF[k]
        @test abs(got - want) <= 1e-14 * max(abs(want), 1e-12)
    end

    # The two RV-derived quantities are the exception, and deliberately so.
    # The Gaia RVS channel is a *spectroscopic* measurement, so v2 predicts it
    # with `radvel`, which now carries the Einstein (second-order Doppler plus
    # gravitational-redshift) term that v1 had no notion of. The term varies
    # over the orbit, so it moves the epoch-RV sample variance — by 4.5e-6
    # relative here — and the total log-likelihood with it, by 2.9e-11. Both
    # gates are an order of magnitude above the measured departure: this is a
    # bounded, understood difference from v1, not a free tolerance.
    @test sim.sample_variance ≈ V1_REF[:sample_variance] rtol = 1e-5
    @test !isapprox(sim.sample_variance, V1_REF[:sample_variance]; rtol=1e-14)
    @test c.lnl(m.sys, c.nt) ≈ V1_REF[:ll] rtol = 1e-9
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

# Two-companion elements reused by several of the testsets below.
const TWO_EL = ((a=4.7, e=0.23, i=0.94, ω=1.31, Ω=2.44, tp=52000.0),
                (a=1.9, e=0.05, i=1.7, ω=0.4, Ω=5.1, tp=53100.0))

@testset "the flux-ratio contract" begin
    T = Float64
    # A solved system to resolve fluxes against. Its bodies carry `flux_G` and
    # `flux_Hp` (host 1.0), so the body-flux tier has something to find.
    m = g23h_model(n_companions=2, masses=(0.03, 0.012),
        fluxratios=(0.037, 0.011), fluxratios_hip=(0.021, 0.008),
        elements=TWO_EL, body_fluxes=true)
    c = g23h_context(m.sys, m.obs)
    sys = c.ctx.system
    hostidx = Octofitter.ref(c.ctx, m.obs.host).idx
    cidx = map(x -> Octofitter.ref(c.ctx, x).idx, m.obs.companions)
    f(θ_obs, θ_sys, key, band, active) =
        Octofitter._g23h_fluxratios(θ_obs, θ_sys, key, Val(band), sys, hostidx, cidx, active, T)

    # Tier 1 — the observation's own vector, in `companions=` order.
    @test f((; fluxratio=(0.1, 0.2)), (;), :fluxratio, :G, (true, true)) == (0.1, 0.2)
    @test f((; fluxratio=SVector(0.1, 0.2)), (;), :fluxratio, :G, (true, true)) == (0.1, 0.2)
    @test f((; fluxratio=[0.1, 0.2]), (;), :fluxratio, :G, (true, true)) == (0.1, 0.2)
    # an inactive companion is dark, whatever the variable says
    @test f((; fluxratio=(0.1, 0.2)), (;), :fluxratio, :G, (true, false)) == (0.1, 0.0)
    # a scalar is unambiguous only for one companion
    @test Octofitter._g23h_fluxratios((; fluxratio=0.3), (;), :fluxratio, Val(:G),
        sys, hostidx, (cidx[1],), (true,), T) == (0.3,)
    @test_throws r"declares 2 companions" f((; fluxratio=0.3), (;), :fluxratio, :G, (true, true))
    @test_throws r"received 3 flux ratios" f((; fluxratio=(0.1, 0.2, 0.3)), (;),
        :fluxratio, :G, (true, true))

    # Tier 2 — a system-level vector of the same name, which is what the old
    # default-variables block forwarded.
    @test f((;), (; fluxratio=(0.3, 0.4)), :fluxratio, :G, (true, true)) == (0.3, 0.4)
    # …and tier 1 wins over it.
    @test f((; fluxratio=(0.1, 0.2)), (; fluxratio=(0.3, 0.4)), :fluxratio, :G, (true, true)) == (0.1, 0.2)

    # Tier 3 — the bodies' own fluxes, as companion/host.
    @test f((;), (;), :fluxratio, :G, (true, true)) == (0.037, 0.011)
    @test f((;), (;), :fluxratio_hip, :Hp, (true, true)) == (0.021, 0.008)
    @test f((;), (;), :fluxratio, :G, (true, false)) == (0.037, 0.0)
    # A band no body declares, while other bands exist, is a name mismatch —
    # `flux_H` for `flux_Hp` silently modelling a dark companion is exactly the
    # failure mode this refuses.
    @test_throws r"declare" f((;), (;), :fluxratio, :H, (true, true))

    # A model with no photometry at all: every companion dark, as v1 was when
    # the flux-ratio vector was omitted.
    m0 = g23h_model(n_companions=2, masses=(0.03, 0.012), elements=TWO_EL,
        fluxratios=(0.0, 0.0), fluxratios_hip=(0.0, 0.0))
    c0 = g23h_context(m0.sys, m0.obs)
    @test Octofitter._g23h_fluxratios((;), (;), :fluxratio, Val(:G), c0.ctx.system,
        1, (2, 3), (true, true), T) == (0.0, 0.0)
end

@testset "flux ratios default to the bodies' own fluxes" begin
    # The §3.5 refinement: the positional vectors stop being mandatory. A model
    # that declares `flux_G`/`flux_Hp` on its bodies must produce *exactly* the
    # numbers the equivalent explicit vectors produce — same arithmetic, one
    # fewer thing to keep in the right order.
    kw = (; n_companions=2, masses=(0.03, 0.012), fluxratios=(0.037, 0.011),
        fluxratios_hip=(0.021, 0.008), elements=TWO_EL)
    m_vec = g23h_model(; kw...)
    m_flx = g23h_model(; kw..., body_fluxes=true)
    # The observation really does declare no flux variables in the second case
    @test :fluxratio ∈ keys(m_vec.obs.derived.variables)
    @test :fluxratio ∉ keys(m_flx.obs.derived.variables)
    @test :fluxratio_hip ∉ keys(m_flx.obs.derived.variables)

    c_vec = g23h_context(m_vec.sys, m_vec.obs)
    c_flx = g23h_context(m_flx.sys, m_flx.obs)
    s_vec = Octofitter.simulate(m_vec.obs, c_vec.ctx)
    s_flx = Octofitter.simulate(m_flx.obs, c_flx.ctx)
    for k in (:μ_h, :μ_hg, :μ_dr2, :μ_dr32, :μ_dr3)
        @test getproperty(s_flx, k) == getproperty(s_vec, k)
    end
    @test s_flx.UEVA_model == s_vec.UEVA_model
    @test s_flx.σ_inflation_hip == s_vec.σ_inflation_hip
    @test s_flx.Δpmra_dr3 == s_vec.Δpmra_dr3
    @test c_flx.lnl(m_flx.sys, c_flx.nt) == c_vec.lnl(m_vec.sys, c_vec.nt)
    # …and it is not a trivial agreement between two all-dark models
    @test maximum(s_flx.σ_inflation_hip) > 1.0

    # An explicit vector still overrides the bodies' fluxes, which is what
    # per-draw resolved-flag gating needs: the bodies keep their `flux_G` /
    # `flux_Hp`, and this draw's vector gates both companions dark.
    m_gate = g23h_model(; kw..., body_fluxes=true, obs_ratios=true,
        fluxratios=(0.0, 0.0), fluxratios_hip=(0.0, 0.0),
        body_fluxratios=(0.037, 0.011), body_fluxratios_hip=(0.021, 0.008))
    c_gate = g23h_context(m_gate.sys, m_gate.obs)
    s_gate = Octofitter.simulate(m_gate.obs, c_gate.ctx)
    @test all(≈(1.0), s_gate.σ_inflation_hip)
    @test s_gate.Δpmra_dr3 != s_flx.Δpmra_dr3
    # …and it matches the model that never declared body fluxes at all
    m_none = g23h_model(; kw..., fluxratios=(0.0, 0.0), fluxratios_hip=(0.0, 0.0))
    c_none = g23h_context(m_none.sys, m_none.obs)
    @test Octofitter.simulate(m_none.obs, c_none.ctx).Δpmra_dr3 == s_gate.Δpmra_dr3
end

@testset "the channels= keyword is the subset path" begin
    full = g23h_model().obs
    chans = (:ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr3, :dec_dr3)
    ctor = g23h_model(channels=chans).obs
    # Constructing with `channels=` must land on the same place as subsetting
    # after the fact — that is the whole reason it filters `table.kind`.
    post = Octofitter.likeobj_from_epoch_subset(full, findall(k -> k ∈ chans, collect(full.table.kind)))
    @test collect(ctor.table.kind) == collect(chans)
    @test collect(post.table.kind) == collect(chans)
    @test Octofitter.Tables.columntable(ctor.table) == Octofitter.Tables.columntable(post.table)
    # …including the two consequences: the abscissa nuisance block is gone,
    # and the Hipparcos catalog distributions survive (Hipparcos channels stay)
    for k in (:hip_iad_jitter, :iad_Δra, :iad_Δdec, :iad_Δplx, :iad_Δpmra, :iad_Δpmdec)
        @test k ∉ keys(ctor.priors.priors)
        @test k ∉ keys(ctor.derived.variables)
    end
    @test !isnothing(ctor.catalog.dist_hip)
    @test !isnothing(ctor.catalog.dist_hg)

    # Dropping every Hipparcos channel drops the distributions, via the same
    # helper the post-hoc subset uses.
    gaia_only = g23h_model(channels=(:ra_dr3, :dec_dr3, :ra_dr2, :dec_dr2)).obs
    @test isnothing(gaia_only.catalog.dist_hip)
    @test isnothing(gaia_only.catalog.dist_hg)

    # A dropped :rv_dr3 must not leave σ_rv_per_transit sampled but unused: the
    # defaults are built against the restricted table.
    defaults = G23HObs(; host=:A, companions=(:b,),
        gaia_id=g23h_catalog_row().gaia_source_id, catalog=g23h_catalog_row(),
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(), channels=chans)
    @test :σ_rv_per_transit ∉ keys(defaults.priors.priors)
    @test :transits_rv ∉ keys(defaults.derived.variables)
    @test :hip_iad_jitter ∉ keys(defaults.priors.priors)

    @test_throws r"is not a G23HObs channel" g23h_model(channels=(:ra_dr3, :nope))
    @test_throws r"selects no channels" g23h_model(channels=())
    @test_throws r"selected none of this source" g23h_model(channels=(:rv_dr3,), include_rv=false)
end

@testset "HGCAObs is a G23HObs over the HGCA's channels" begin
    cat = g23h_catalog_row()
    m = g23h_model()
    obs = Octofitter.HGCAObs(; host=m.host, companions=(m.comps[1],),
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(),
        variables=@variables begin
            σ_AL = $(EVAL[:sigma_AL])
            σ_att = $(EVAL[:sigma_att])
            σ_calib = $(EVAL[:sigma_calib])
            fluxratio = $((EVAL[:fluxratio_G],))
            fluxratio_hip = $((EVAL[:fluxratio_Hp],))
            # The Hipparcos abscissa nuisances are declared and must be
            # stripped anyway: `channels=` drops the `:iad_hip` row.
            hip_iad_jitter = $(EVAL[:hip_iad_jitter])
            iad_Δra = 0.0
            iad_Δdec = 0.0
            iad_Δplx = 0.0
            iad_pmra = $(EVAL[:iad_pmra])
            iad_pmdec = $(EVAL[:iad_pmdec])
            pmra = $(cat.pmra_dr3)
            pmdec = $(cat.pmdec_dr3)
            transit_priorities = $PRIORITIES
            transits = $TRANSITS
            transits_dr2 = $TRANSITS_DR2
        end)
    @test obs isa G23HObs
    for k in (:hip_iad_jitter, :iad_Δra, :iad_Δdec, :iad_Δplx, :iad_pmra, :iad_pmdec)
        @test k ∉ keys(obs.derived.variables)
    end
    # Hipparcos, Hipparcos-Gaia and DR3 — no DR2, no DR32, no UEVA, no RV, no
    # abscissae.
    @test collect(obs.table.kind) == [:ra_hip, :dec_hip, :ra_hg, :dec_hg, :ra_dr3, :dec_dr3]
    # :none, not false: the datum AND the UEVA-driven DR3 deflation both go.
    @test obs.ueva_mode === :none
    @test obs.include_iad == false
    @test Octofitter.likelihoodname(obs) == "HGCA"

    sys = Octofitter.System(name="hgca", bodies=[m.host, m.comps...], observations=[obs],
        variables=(@variables begin
            plx = $(cat.parallax)
        end), observing_geometry=false)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    @test isfinite(Octofitter.make_ln_like(sys, nt)(sys, nt))

    c = g23h_context(sys, obs)
    sim = Octofitter.simulate(obs, c.ctx)
    @test sim.deflation_factor_dr3 == 1.0     # :none ⇒ no deflation
    @test all(isfinite, sim.μ_hg)

    # The retired types name themselves rather than failing as undefined.
    @test_throws r"no longer exists" Octofitter.HGCAInstantaneousObs()
    @test_throws r"no longer exists" Octofitter.GaiaCatalogFitObs()
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
    ruwe = g23h_model(ueva_mode=:RUWE)
    @test :ueva_dr3 ∈ ruwe.obs.table.kind
    @test :ueva_dr3 ∉ m.obs.table.kind
    @test length(m.obs.table.kind) == length(ruwe.obs.table.kind) - 1
    c = g23h_context(m.sys, m.obs)
    sim = Octofitter.simulate(m.obs, c.ctx)
    # No UEVA datum means no UEVA-driven deflation of the DR3 block.
    @test sim.deflation_factor_dr3 == 1.0
    @test isfinite(c.lnl(m.sys, c.nt))
    @test_throws r"must be :RUWE, :EAN or :none" g23h_model(ueva_mode=:bogus)

    # The three σ nuisances are *fixed*, not sampled, under `:none` — in the
    # **default** variable set, which is what a user who does not pass
    # `variables=` gets. (`g23h_model` always passes one, so it cannot show
    # this.) They reach the likelihood only through the UEVA channel, which
    # `:none` switches off, so sampling them would add three unconstrained
    # dimensions; worse, it makes the observation unconstructable for the very
    # sources `:none` exists to serve — the brightest stars, whose
    # sig_AL/sig_att_radec/sig_cal calibration is NaN, and
    # `truncated(Normal(NaN, NaN), …)` cannot be built. Ported from main's
    # test_g23h_ueva_none.jl, which was written against the v1 API.
    cat = g23h_catalog_row()
    defaults(mode) = G23HObs(; host=m.host, companions=Tuple(m.comps),
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(),
        ueva_mode=mode, include_rv=false)
    d_none, d_ruwe = defaults(:none), defaults(:RUWE)
    for σ in (:σ_AL, :σ_att, :σ_calib)
        @test !haskey(d_none.priors.priors, σ)      # fixed…
        @test haskey(d_none.derived.variables, σ)
        @test haskey(d_ruwe.priors.priors, σ)       # …and control: :RUWE samples them
    end

    # …and they are genuinely inert: perturbing them cannot move the likelihood
    # under `:none`, while `:RUWE` responds. This is the assertion that would
    # catch a future channel quietly starting to read σ_formal.
    ll_base = c.lnl(m.sys, c.nt)
    for probe in (0.05, 0.5, 2.0)
        nt_p = _g23h_with_sigmas(c.nt, m.obs, probe)
        @test c.lnl(m.sys, nt_p) == ll_base              # bit-identical
    end
    c_r = g23h_context(ruwe.sys, ruwe.obs)
    ll_r = c_r.lnl(ruwe.sys, c_r.nt)
    @test c_r.lnl(ruwe.sys, _g23h_with_sigmas(c_r.nt, ruwe.obs, 2.0)) != ll_r
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

# ── several catalog sources in one system ─────────────────────────────────
#
# Two Gaia/Hipparcos sources 2.8″ apart, each with its own inner companion,
# modelled off ONE trajectory against ONE shared frame. Nothing here is new
# machinery — `refspecs` has always been variadic and the epoch union has
# always been shared — so these tests exist to pin that it stays true.
function g23h_two_source_model(; a_wide=12.0, observing_geometry=false, frame_shift=false)
    cat = g23h_catalog_row()
    Aa = Octofitter.Body(name="Aa", variables=@variables begin
        mass = 0.40
        flux_G = 1.0
        flux_Hp = 1.0
    end)
    Ab = Octofitter.Body(name="Ab", about=Aa, variables=@variables begin
        mass = 0.05
        flux_G = 0.04
        flux_Hp = 0.02
        a = 1.0; e = 0.1; i = 0.9; ω = 1.2; Ω = 2.0; tp = 52000.0
    end)
    Ba = Octofitter.Body(name="Ba", about=Aa, variables=@variables begin
        mass = 0.30
        flux_G = 0.5
        flux_Hp = 0.4
        a = $a_wide; e = 0.2; i = 1.1; ω = 0.3; Ω = 1.0; tp = 40000.0
    end)
    Bb = Octofitter.Body(name="Bb", about=Ba, variables=@variables begin
        mass = 0.03
        flux_G = 0.01
        flux_Hp = 0.005
        a = 1.5; e = 0.05; i = 0.6; ω = 2.2; Ω = 4.0; tp = 51000.0
    end)

    # Both sources declare exactly the same variables, and NONE of them is a
    # per-source offset on a Gaia channel: the two sources share the system's
    # frame, and that shared frame is what makes the wide orbit identifiable.
    srcvars() = @variables begin
        σ_AL = $(EVAL[:sigma_AL])
        σ_att = $(EVAL[:sigma_att])
        σ_calib = $(EVAL[:sigma_calib])
        hip_iad_jitter = $(EVAL[:hip_iad_jitter])
        iad_Δra = 0.0
        iad_Δdec = 0.0
        iad_Δplx = 0.0
        iad_pmra = $(EVAL[:iad_pmra])
        iad_pmdec = $(EVAL[:iad_pmdec])
        σ_rv_per_transit = $(EVAL[:sigma_rv_per_transit])
        transit_priorities = $PRIORITIES
        transits = $TRANSITS
        transits_dr2 = $TRANSITS_DR2
        transits_rv = $TRANSITS_RV
    end
    # `frame_shift=false` is mandatory here and the model errors without it —
    # see the guard's own testset. The shift redefines the system's `pmra` to
    # mean one observation's host body, and two observations cannot both do
    # that to the same parameter.
    src(host, comps, nm) = G23HObs(;
        host, companions=comps, ref=Barycentre, name=nm, frame_shift,
        gaia_id=cat.gaia_source_id, catalog=cat,
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(), variables=srcvars())
    obsA = src(Aa, (Ab,), "srcA")
    obsB = src(Ba, (Bb,), "srcB")

    sys = Octofitter.System(name="ABquad", bodies=[Aa, Ab, Ba, Bb],
        observations=[obsA, obsB],
        variables=(@variables begin
            plx = $(cat.parallax)
            ra = $(cat.ra)
            dec = $(cat.dec)
            pmra = $(cat.pmra_dr3)
            pmdec = $(cat.pmdec_dr3)
            rv = 0.0
            ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
        end), observing_geometry=observing_geometry)
    return (; sys, obsA, obsB, Aa, Ab, Ba, Bb)
end

@testset "several catalog sources in one system" begin
    m = g23h_two_source_model()
    @test Octofitter._refnames.(Octofitter.refspecs(m.obsA)) == ((:Aa,), (:Ab,), ())
    @test Octofitter._refnames.(Octofitter.refspecs(m.obsB)) == ((:Ba,), (:Bb,), ())
    @test sprint(show, MIME"text/plain"(), m.obsB) |>
          s -> occursin("Ba (blended with Bb) vs Barycentre", s)

    # One trajectory: the two sources' epochs are solved as a single union.
    eps, maps = Octofitter.epoch_plan(m.sys)
    @test eps == sort(unique(vcat(Octofitter.epochs(m.obsA), Octofitter.epochs(m.obsB))))
    @test length(maps[m.obsA]) == length(Octofitter.epochs(m.obsA))
    @test length(maps[m.obsB]) == length(Octofitter.epochs(m.obsB))

    nt = Octofitter.make_arr2nt(m.sys)(Float64[])
    lnl = Octofitter.make_ln_like(m.sys, nt)
    @test isfinite(lnl(m.sys, nt))

    # One shared frame: neither source declares its own reference proper
    # motion, and both read the system's. Per-source Gaia offsets are exactly
    # what would dilute the constraint that binds the two pairs together, so
    # their absence is the assertion.
    for nm in (:srcA, :srcB)
        θ_obs = getproperty(nt.observations, nm)
        @test !hasproperty(θ_obs, :pmra)
        @test !hasproperty(θ_obs, :pmdec)
        @test Octofitter._g23h_pm(nt, θ_obs, Float64) ==
              (Float64(nt.pmra), Float64(nt.pmdec))
    end

    # …and the wide orbit really is constrained by both sources at once: move
    # it and *both* sources' modelled proper motions move, without either
    # source knowing anything about the other.
    cA = g23h_context(m.sys, m.obsA); cB = g23h_context(m.sys, m.obsB)
    sA = Octofitter.simulate(m.obsA, cA.ctx); sB = Octofitter.simulate(m.obsB, cB.ctx)
    m2 = g23h_two_source_model(a_wide=18.0)
    c2A = g23h_context(m2.sys, m2.obsA); c2B = g23h_context(m2.sys, m2.obsB)
    s2A = Octofitter.simulate(m2.obsA, c2A.ctx); s2B = Octofitter.simulate(m2.obsB, c2B.ctx)
    # `Δpmra_dr3` is the perturbation the DR3 refit sees. With `frame_shift`
    # off it also reaches `μ_dr3`, which is the property the next testset is
    # about; here all that is claimed is that moving the wide orbit moves both
    # sources.
    @test sA.Δpmra_dr3 != s2A.Δpmra_dr3
    @test sB.Δpmra_dr3 != s2B.Δpmra_dr3
    @test sA.μ_hg != s2A.μ_hg
    @test sB.μ_hg != s2B.μ_hg
    # The two sources are genuinely different sources, not two copies
    @test sA.Δpmra_dr3 != sB.Δpmra_dr3
    @test sA.μ_hg != sB.μ_hg
end

@testset "frame_shift: what it does, and why two sources cannot share it" begin
    # § What it does. The shift subtracts the model's own DR3-epoch reflex of
    # the host from *every* channel, so `μ_dr3` becomes exactly the frame
    # `pmra`/`pmdec` and every other channel is translated by the same
    # constant. Both halves of that are asserted, because the first alone
    # would also hold if the shift were applied to DR3 only.
    on = g23h_model()
    off = g23h_model()
    # `g23h_model` is parallax-only, so build the shifted/unshifted pair by
    # hand off the same catalog row rather than through the frame.
    shifted = Octofitter.simulate(on.obs, g23h_context(on.sys, on.obs).ctx)
    obs_off = G23HObs(; host=off.host, companions=(off.comps[1],), frame_shift=false,
        gaia_id=g23h_catalog_row().gaia_source_id, catalog=g23h_catalog_row(),
        forecast_table=g23h_forecast(), hipparcos=g23h_hipparcos(),
        dr2_transits_catalog=g23h_dr2_sidecar(),
        variables=(off.obs.priors, off.obs.derived))
    sys_off = Octofitter.System(name="g23h", bodies=[off.host, off.comps...],
        observations=[obs_off], variables=(@variables begin
            plx = $(g23h_catalog_row().parallax)
        end), observing_geometry=false)
    plain = Octofitter.simulate(obs_off, g23h_context(sys_off, obs_off).ctx)

    Δ = SVector(plain.Δpmra_dr3, plain.Δpmdec_dr3)
    @test shifted.Δpmra_dr3 ≈ plain.Δpmra_dr3      # the shift itself is unchanged
    @test shifted.μ_dr3 ≈ plain.μ_dr3 .- Δ
    @test shifted.μ_dr2 ≈ plain.μ_dr2 .- Δ
    @test shifted.μ_hg ≈ plain.μ_hg .- Δ
    @test shifted.μ_h ≈ plain.μ_h .- Δ
    @test shifted.μ_dr32 ≈ plain.μ_dr32 .- Δ
    # It is a real difference, not a no-op on this configuration.
    @test !(shifted.μ_dr3 ≈ plain.μ_dr3)
    # …and the UEVA channel, which is a scatter statistic rather than a
    # velocity, is untouched by a common-mode translation.
    @test shifted.UEVA_model ≈ plain.UEVA_model

    # § Why two sources cannot. Each observation would subtract its own Δpm
    # from the one shared `pmra`, so both DR3 channels would predict the same
    # number. Caught at model-build time, by name, rather than as a quietly
    # wrong likelihood.
    err = try
        g23h_two_source_model(frame_shift=true); nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err)
    @test occursin("srcA", err) && occursin("srcB", err)
    @test occursin("frame_shift", err) && occursin("AnchoredFrame", err)

    # One source with the shift on is untouched — this is still the default,
    # and `g23h_model` above exercises it throughout the file.
    @test g23h_model().obs.frame_shift
    # …and so is a *pair* that shares no frame: with no absolute frame each
    # observation carries its own reference proper motion, so nothing is
    # redefined out from under anyone.
    @test Octofitter.check_siblings(on.obs, (on.obs, obs_off),
        (; name=:x, framevars=Symbol[], sysnames=Symbol[])) === nothing
end

@testset "a shared system-level fluxratio names the observations that clash" begin
    # §10. Tier 2 of `_g23h_fluxratios` is keyed by the bare name, so a
    # system-level `fluxratio` is read by *every* G23HObs, each validating it
    # against its own companion count. Two observations with different counts
    # can never both be right, and the runtime failure ("received 3 flux
    # ratios but the observation declares 1 companion") names neither of them
    # nor says why they are competing.
    cat = g23h_catalog_row()
    Aa = Octofitter.Body(name="Aa", variables=@variables begin mass = 0.4 end)
    Ab = Octofitter.Body(name="Ab", about=Aa, variables=@variables begin
        mass = 0.05; a = 1.0; e = 0.1; i = 0.9; ω = 1.2; Ω = 2.0; tp = 52000.0
    end)
    Ac = Octofitter.Body(name="Ac", about=Aa, variables=@variables begin
        mass = 0.02; a = 3.0; e = 0.1; i = 0.9; ω = 1.2; Ω = 2.0; tp = 52000.0
    end)
    Ba = Octofitter.Body(name="Ba", about=Aa, variables=@variables begin
        mass = 0.3; a = 12.0; e = 0.2; i = 1.1; ω = 0.3; Ω = 1.0; tp = 40000.0
    end)
    srcvars() = @variables begin
        σ_AL = $(EVAL[:sigma_AL]); σ_att = $(EVAL[:sigma_att])
        σ_calib = $(EVAL[:sigma_calib])
        hip_iad_jitter = $(EVAL[:hip_iad_jitter])
        iad_Δra = 0.0; iad_Δdec = 0.0; iad_Δplx = 0.0
        iad_pmra = $(EVAL[:iad_pmra]); iad_pmdec = $(EVAL[:iad_pmdec])
        σ_rv_per_transit = $(EVAL[:sigma_rv_per_transit])
        transit_priorities = $PRIORITIES
        transits = $TRANSITS; transits_dr2 = $TRANSITS_DR2; transits_rv = $TRANSITS_RV
    end
    src(host, comps, nm; vars=srcvars()) = G23HObs(;
        host, companions=comps, name=nm, frame_shift=false,
        gaia_id=cat.gaia_source_id, catalog=cat, forecast_table=g23h_forecast(),
        hipparcos=g23h_hipparcos(), dr2_transits_catalog=g23h_dr2_sidecar(),
        variables=vars)
    build(obsA, obsB; extra=()) = Octofitter.System(
        name="clash", bodies=[Aa, Ab, Ac, Ba], observations=[obsA, obsB],
        variables=(@variables begin
            plx = $(cat.parallax); ra = $(cat.ra); dec = $(cat.dec)
            pmra = $(cat.pmra_dr3); pmdec = $(cat.pmdec_dr3); rv = 0.0
            ref_epoch = $(Octofitter.meta_gaia_DR3.ref_epoch_mjd)
            fluxratio = (0.04, 0.01)
        end), observing_geometry=false, verbosity=0)

    twoc = src(Aa, (Ab, Ac), "srcA")
    onec = src(Ba, (), "srcB")
    err = try
        build(twoc, onec); nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err)
    @test occursin("srcA", err) && occursin("srcB", err)
    @test occursin("fluxratio", err)
    # The counts are in the message, which is the part the runtime error had.
    @test occursin("2 companion", err) && occursin("0 companion", err)

    # Equal companion counts are fine — one vector, one meaning.
    @test build(src(Aa, (Ab,), "srcA"), src(Ba, (Ab,), "srcB")) isa Octofitter.System
    # …and so is a clash where one observation declares its own `fluxratio`,
    # because tier 1 wins and it never reaches the system-level vector.
    own = vcat(srcvars(), @variables begin fluxratio = (0.5,) end)
    @test build(src(Aa, (Ab,), "srcA"; vars=own), onec) isa Octofitter.System
end

@testset "differential parallax is per body, already" begin
    # `_observe_pass!` writes `cart2angle[k, j] = rad2mas / (d_au + pz)` — per
    # epoch AND per body, using that body's own retarded line-of-sight depth.
    # So each component's AU→mas scale is its own, and nothing extra is needed
    # for a wide pair whose components sit at measurably different distances.
    # (The parallax *factors* are per-source by a different route: they arrive
    # as data columns from each source's own forecast/IAD.)
    m = g23h_two_source_model(observing_geometry=true)
    c = g23h_context(m.sys, m.obsA)
    traj = c.ctx.traj
    c2a = traj.cart2angle
    @test size(c2a, 2) == 4
    k = firstindex(traj.epochs)
    for j in 1:4
        @test c2a[k, j] ≈ PlanetOrbits.rad2mas / (traj.d_au[k] + traj.z[k, j])
    end
    # It is a real per-body spread, not four copies of one number.
    @test !all(≈(c2a[k, 1]), @view c2a[k, :])

    # …and it collapses to one shared scale when the geometry is skipped,
    # which is what `observing_geometry=false` means.
    mflat = g23h_two_source_model(observing_geometry=false)
    cflat = g23h_context(mflat.sys, mflat.obsA)
    c2aflat = cflat.ctx.traj.cart2angle
    @test all(≈(c2aflat[k, 1]), @view c2aflat[k, :])
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

@testset "catalog convention and the correction impact" begin
    # G23H does NOT declare `reduced_lighttime_free`, though its inputs are
    # catalog reductions. Its dominant channels (Hipparcos−Gaia, DR3−DR2) are
    # differences BETWEEN reduction windows, and the path between two windows is
    # not governed by either one's parameterization — so the flag is measured
    # like any other correction rather than pinned. See the comment beside
    # `G23HObs`'s correction-flag section, and `reduced_lighttime_free`'s
    # docstring for the distinction the trait actually draws.
    @test !Octofitter.reduced_lighttime_free(G23HObs)
    m = g23h_model()
    d = only(filter(x -> x.flag === :barycentric_lighttime,
                    m.sys.corrections.decisions))
    # Decided by the impact test, not by a convention pin: whichever way it
    # resolves, it got there by taking draws.
    @test d.source !== :data
    @test d.ndraws > 0

    # ...and G23H can now answer the impact question for the *other* flag,
    # instead of vetoing every model it appears in with "cannot say": the five
    # PM channel pairs, against the tightest catalog σ.
    @test Octofitter.has_correction_impact(G23HObs)
    sys, obs = m.sys, m.obs
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    lnl = Octofitter.make_ln_like(sys, nt)
    posys = lnl.build(nt)
    eps, maps = Octofitter.epoch_plan(sys)
    θ_obs = getproperty(nt.observations,
        Symbol(Octofitter.normalizename(Octofitter.likelihoodname(obs))))
    mkctx(og) = Octofitter.ObsContext(nt, θ_obs, posys,
        orbitsolve(posys, eps; method=sys.method, observing_geometry=og,
                   barycentric_lighttime=false), maps[obs])
    r = Octofitter.correction_impact(obs, mkctx(true), mkctx(false))
    @test r.n == 10                       # five (ra, dec) channel pairs
    @test isfinite(r.delta) && r.delta >= 0
    @test isfinite(r.sigma) && r.sigma > 0
    # Identical contexts must measure an identical model: zero impact.
    r0 = Octofitter.correction_impact(obs, mkctx(true), mkctx(true))
    @test r0.delta == 0
end
