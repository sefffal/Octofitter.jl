#=
Under `ueva_mode = :none`, σ_AL / σ_att / σ_calib must not reach the likelihood.

They are declared only to keep `σ_formal = √(σ_att² + σ_AL²)` finite for the
DR3 5-parameter refit, where a *scalar* weight cancels out of the least-squares
solution. Every other consumer (the UEVA datum, its change-of-variables
Jacobian, its covariance entry, the DR3 covariance deflation) is switched off by
:none. This test pins that down by evaluating the SAME compiled ln_like closure
at the same parameter vector, changing only those three values inside θ_obs:
under :none the result must be bit-identical, under :RUWE it must move.

Standalone (like test_g23h_simulation.jl) — needs the G23H catalog DataDep.
=#

using Octofitter, Distributions, Random, Test

ref_epoch = Octofitter.meta_gaia_DR3.ref_epoch_mjd

function build_probe_system(absastrom, name)
    planet = Planet(name="b", basis=Visual{KepOrbit}, observations=[],
        variables=@variables begin
            a    ~ LogUniform(0.1, 100)
            e    ~ Uniform(0, 0.99)
            ω    ~ Uniform(0, 2pi)
            i    ~ Sine()
            Ω    ~ Uniform(0, 2pi)
            θ    ~ Uniform(0, 2pi)
            tp   = θ_at_epoch_to_tperi(θ, $ref_epoch; M=system.M, e, a, i, ω, Ω)
            mass ~ LogUniform(0.01, 1000)
        end)
    plx0 = Float64(absastrom.catalog.parallax)
    eplx = ismissing(absastrom.catalog.parallax_error) ? 0.01 * plx0 :
           Float64(absastrom.catalog.parallax_error)
    pmra0 = Float64(absastrom.catalog.pmra_dr3)
    pmdec0 = Float64(absastrom.catalog.pmdec_dr3)
    return System(name=name, companions=[planet], observations=[absastrom],
        variables=@variables begin
            M   ~ truncated(Normal(0.71, 0.04), lower=0.1)
            plx ~ truncated(Normal(plx0, eplx), lower=max(0.1, plx0 - 10*eplx))
            pmra  ~ Uniform(pmra0 - 100, pmra0 + 100)
            pmdec ~ Uniform(pmdec0 - 100, pmdec0 + 100)
            dec = $(absastrom.catalog.dec)
            ra  = $(absastrom.catalog.ra)
            ref_epoch = $ref_epoch
        end)
end

"θnt with only the three σ nuisances replaced."
function with_sigmas(θnt, σ_AL, σ_att, σ_calib)
    g = θnt.observations.G23H
    g2 = merge(g, (; σ_AL=σ_AL, σ_att=σ_att, σ_calib=σ_calib))
    return merge(θnt, (; observations=merge(θnt.observations, (; G23H=g2))))
end

# Spans ~200x, well beyond any plausible calibration.
const PROBES = ((0.132, 0.0779, 0.0795), (0.01, 0.01, 0.01), (2.0, 1.0, 1.0))

@testset "G23H ueva_mode=:none — sigma nuisances are inert" begin
    for mode in (:none, :RUWE)
        obs = G23HObs(; hip_id=18512, freeze_epochs=true, include_rv=false,
                        ueva_mode=mode)
        sys = build_probe_system(obs, "probe_$mode")
        Random.seed!(20260728)
        θnt = Octofitter.make_arr2nt(sys)(Octofitter.sample_priors(sys))
        lnlike = Octofitter.make_ln_like(sys, θnt)
        ll_base = lnlike(sys, θnt)
        @test isfinite(ll_base)

        @testset "$mode / σ=$p" for p in PROBES
            ll = lnlike(sys, with_sigmas(θnt, p...))
            if mode == :none
                @test ll == ll_base          # bit-identical
            else
                @test ll != ll_base          # control: :RUWE does respond
            end
        end
    end
end
