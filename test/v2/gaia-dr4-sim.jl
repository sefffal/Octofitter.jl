# Simulated Gaia DR4 epoch astrometry: the measured per-source noise budget
# (`g23h_scan_uncertainty`) and the transit-level table built from a scan
# forecast (`gaia_dr4_transit_template`).
#
# Everything here is offline. The uncertainty model reads the same real G23H
# catalog row fixture `g23h.jl` uses (HIP 1475 / Gaia DR3 385334230892516480),
# and the scan geometry is the same cached GOST forecast fixture — which is the
# point of the shared `forecast_table` contract: one table drives both a
# `G23HObs` and a simulated DR4 observation.

using Octofitter
using Octofitter.TypedTables: Table
using Octofitter.CSV
using PlanetOrbits
using Distributions, Random, Statistics, Test

const DR4SIM_FIX = joinpath(@__DIR__, "fixtures")
_dr4sim_fixture(name) = Table(CSV.File(joinpath(DR4SIM_FIX, name)))

function dr4sim_catalog_row()
    t = _dr4sim_fixture("g23h-catalog-row.csv")
    cols = Tuple(Octofitter.Tables.columnnames(t))
    return NamedTuple{cols}(map(c -> getproperty(t, c)[1], cols))
end

const DR4SIM_CAT = dr4sim_catalog_row()
const DR4SIM_FORECAST = _dr4sim_fixture("g23h-forecast.csv")

@testset "g23h_scan_uncertainty" begin
    σ = g23h_scan_uncertainty(gaia_id=DR4SIM_CAT.gaia_source_id, catalog=DR4SIM_CAT)

    # The three catalog terms are passed through untouched...
    @test σ.σ_AL == DR4SIM_CAT.sig_AL
    @test σ.σ_att == DR4SIM_CAT.sig_att_radec
    @test σ.σ_calib == DR4SIM_CAT.sig_cal
    @test σ.gaia_id == DR4SIM_CAT.gaia_source_id
    @test σ.phot_g_mean_mag == DR4SIM_CAT.phot_g_mean_mag_dr3

    # ...and combined the way `G23HObs` combines them. `σ_formal` is the
    # per-CCD quantity G23H forms in `_g23h_simulate!`:
    @test σ.σ_formal ≈ sqrt(DR4SIM_CAT.sig_att_radec^2 + DR4SIM_CAT.sig_AL^2)

    # `n_ccd` is G23H's own `N_AL = N / N_FoV`, and the per-transit formal
    # uncertainty is `σ_formal / √N_AL` — the effective per-transit variance in
    # G23H's DR3 channel, where a transit-level χ² computed with `σ_formal` is
    # multiplied by `N_AL` before it is compared to the catalog.
    n_ccd = DR4SIM_CAT.astrometric_n_good_obs_al_dr3 /
            DR4SIM_CAT.astrometric_matched_transits_dr3
    @test σ.n_ccd ≈ n_ccd
    @test 5 < σ.n_ccd < 12                      # Gaia sees SM + AF1..AF9
    @test σ.σ_transit_formal ≈ σ.σ_formal / sqrt(n_ccd)

    # The calibration term is per transit and does *not* average down with the
    # CCD observations, so it dominates the true per-transit scatter.
    @test σ.σ_transit_true ≈ hypot(σ.σ_transit_formal, DR4SIM_CAT.sig_cal)
    @test σ.σ_transit_true > σ.σ_transit_formal

    # Recorded values, so a change in the arithmetic has to be deliberate.
    @test σ.σ_formal ≈ 0.10702 atol = 1e-5
    @test σ.σ_transit_formal ≈ 0.03630 atol = 1e-5
    @test σ.σ_transit_true ≈ 0.15193 atol = 1e-5

    # `hip_id` resolves the same row.
    @test g23h_scan_uncertainty(hip_id=Int(DR4SIM_CAT.hip_id),
        catalog=DR4SIM_CAT).σ_formal == σ.σ_formal

    @test_throws Exception g23h_scan_uncertainty(catalog=DR4SIM_CAT)
    # A source outside the calibration's footprint carries NaN, and must say so
    # rather than propagating a NaN error bar into a simulated table.
    @test_throws Exception g23h_scan_uncertainty(gaia_id=DR4SIM_CAT.gaia_source_id,
        catalog=merge(DR4SIM_CAT, (; sig_cal=NaN)))
    @test_throws Exception g23h_scan_uncertainty(gaia_id=DR4SIM_CAT.gaia_source_id,
        catalog=(; sig_AL=0.05, sig_att_radec=0.05, sig_cal=0.1))  # no transit counts
end

@testset "gaia_dr4_transit_template" begin
    σ_al = 0.05
    tbl = gaia_dr4_transit_template(; σ_al, forecast_table=DR4SIM_FORECAST)
    n = length(DR4SIM_FORECAST.epoch)

    @test length(tbl.epoch) == n
    @test issorted(tbl.epoch)
    @test issubset(Octofitter.dr4_cols, Octofitter.Tables.columnnames(tbl))
    @test tbl.epoch == sort(collect(DR4SIM_FORECAST.epoch))
    @test all(iszero, tbl.centroid_pos_al)
    @test all(==(σ_al), tbl.centroid_pos_error_al)
    @test all(iszero, tbl.outlier_flag)

    # Degrees, not radians: the unit `GaiaDR4AstromObs` ingests.
    perm = sortperm(collect(DR4SIM_FORECAST.epoch))
    @test tbl.scan_pos_angle ≈ rad2deg.(collect(DR4SIM_FORECAST.scanAngle_rad)[perm])
    @test maximum(abs, tbl.scan_pos_angle) > 2π   # would be false in radians
    @test tbl.parallax_factor_al ≈
          collect(DR4SIM_FORECAST.parallaxFactorAlongScan)[perm]

    # A GOST-shaped table (Julian days, `scanAngle_rad_`) gives the same thing.
    gost_shaped = Table(
        ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_=mjd2jd.(collect(DR4SIM_FORECAST.epoch)),
        scanAngle_rad_=collect(DR4SIM_FORECAST.scanAngle_rad),
        parallaxFactorAlongScan=collect(DR4SIM_FORECAST.parallaxFactorAlongScan),
    )
    tbl2 = gaia_dr4_transit_template(; σ_al, forecast_table=gost_shaped)
    @test tbl2.epoch ≈ tbl.epoch
    @test tbl2.scan_pos_angle ≈ tbl.scan_pos_angle
    @test tbl2.parallax_factor_al ≈ tbl.parallax_factor_al

    # Per-transit σ vectors are taken in the caller's row order and permuted
    # with the epochs, not silently mismatched.
    σvec = collect(range(0.03, 0.09, length=n))
    tblv = gaia_dr4_transit_template(; σ_al=σvec, forecast_table=DR4SIM_FORECAST)
    @test tblv.centroid_pos_error_al == σvec[perm]

    @test_throws Exception gaia_dr4_transit_template(; σ_al=σvec[1:end-1],
        forecast_table=DR4SIM_FORECAST)
    @test_throws Exception gaia_dr4_transit_template(; σ_al=-0.05,
        forecast_table=DR4SIM_FORECAST)
    # Neither a forecast nor anything to fetch one with: fail before the network.
    @test_throws Exception gaia_dr4_transit_template(; σ_al)
    @test_throws Exception gaia_dr4_transit_template(; σ_al,
        forecast_table=Table(epoch=[1.0, 2.0], scanAngle_rad=[0.1, 0.2]))
end

# A simulated table has to be fittable: it must construct a `GaiaDR4AstromObs`,
# give a finite log-likelihood, and round-trip through `generate_from_params`.
@testset "a simulated DR4 observation" begin
    σ = g23h_scan_uncertainty(gaia_id=DR4SIM_CAT.gaia_source_id, catalog=DR4SIM_CAT)
    tbl = gaia_dr4_transit_template(; σ_al=σ.σ_transit_true,
        forecast_table=DR4SIM_FORECAST)
    mkobs(t) = GaiaDR4AstromObs(t; target=Photocentre, ref=Barycentre, name="GaiaDR4",
        variables=@variables begin
            astrometric_jitter = 0.0
            ra_offset_mas = 0.4
            dec_offset_mas = -0.2
            pmra = 12.0
            pmdec = -8.0
            ref_epoch = 57388.5
        end)

    # The scan angles round-trip through the degree convention.
    obs = mkobs(tbl)
    @test obs.sinψ ≈ sin.(deg2rad.(tbl.scan_pos_angle))
    @test obs.cosψ ≈ cos.(deg2rad.(tbl.scan_pos_angle))

    function mksys(o)
        A = Octofitter.Body(name="A", variables=@variables begin
            mass = 0.4
            flux = 1.0
        end)
        b = Octofitter.Body(name="b", about=A, variables=@variables begin
            mass = 0.01
            flux = 0.0
            a = 1.5
            e = 0.1
            ω = 0.3
            i = 0.8
            Ω = 0.6
            θ = 1.1
            epoch = 57388.5
        end)
        return Octofitter.System(name="sim", bodies=[A, b], observations=[o],
            variables=@variables begin
                plx = 20.0
            end)
    end

    sys = mksys(obs)
    lnlike(s) = Octofitter.make_ln_like(s)(s, Octofitter.make_arr2nt(s)(Float64[]))
    @test isfinite(lnlike(sys))

    # ...and the same table goes into a real sampled model. (`sys` above fixes
    # every variable, which is what makes the analytic check below exact, but a
    # model with no free variable cannot be compiled.)
    let A = Octofitter.Body(name="A", variables=@variables begin
            mass = 0.4
            flux = 1.0
        end),
        b = Octofitter.Body(name="b", about=A, variables=@variables begin
            mass ~ LogUniform(0.001, 0.1)
            flux = 0.0
            a ~ LogUniform(0.5, 5.0)
            e = 0.1; ω = 0.3; i = 0.8; Ω = 0.6; θ = 1.1
            epoch = 57388.5
        end)

        sysfree = Octofitter.System(name="simfree", bodies=[A, b],
            observations=[mkobs(tbl)], variables=@variables begin
                plx ~ truncated(Normal(20.0, 1.0), lower=1.0)
            end)
        model = Octofitter.LogDensityModel(sysfree)
        @test model.D == 3
        @test isfinite(model.ℓπcallback(model.link(model.sample_priors(Random.Xoshiro(7)))))
    end

    θ = Octofitter.drawfrompriors(sys)

    # Noiseless: the data become exactly the modelled abscissae, so the
    # residual term vanishes and the log-likelihood is the normalization alone.
    clean = Octofitter.generate_from_params(sys, θ; add_noise=false)
    n = length(tbl.epoch)
    expected = -n * log(2π) / 2 - sum(log.(tbl.centroid_pos_error_al .^ 2)) / 2
    @test lnlike(clean) ≈ expected rtol = 1e-10
    @test !all(iszero, clean.observations[1].table.centroid_pos_al)

    # Noised: scattered by the table's own error column, so the χ² per transit
    # is ~1 and the log-likelihood is finite but lower.
    Random.seed!(1234)
    noisy = Octofitter.generate_from_params(sys, θ; add_noise=true)
    resid = noisy.observations[1].table.centroid_pos_al .-
            clean.observations[1].table.centroid_pos_al
    @test isfinite(lnlike(noisy))
    @test lnlike(noisy) < lnlike(clean)
    @test 0.5 < std(resid) / σ.σ_transit_true < 1.6
    # The uncertainties, epochs and scan geometry survive the round trip.
    @test noisy.observations[1].table.centroid_pos_error_al == tbl.centroid_pos_error_al
    @test noisy.observations[1].table.scan_pos_angle == tbl.scan_pos_angle
end
