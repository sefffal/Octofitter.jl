# Integration test for GaiaDR3CompanionObs: inject a synthetic Gaia DR3
# companion solution and recover the dark central mass with MCMC.
using Octofitter
using Distributions
using PlanetOrbits
using Test
using Random
Random.seed!(42)

const _mjup2msol = Octofitter.mjup2msol
const _REF = Octofitter.meta_gaia_DR3.ref_epoch_mjd
const _M_central_true = 1.2
const _M_comp_msol = 0.3
const _M_comp_mjup = _M_comp_msol / _mjup2msol
const _ra0, _dec0 = 150.0, -30.0

# --- generate synthetic data at truth ---
gen_planet = Planet(
    name="B", basis=AbsoluteVisual{KepOrbit},
    observations=[GaiaDR3CompanionObs(;
        ra=_ra0, dec=_dec0, parallax=50.0, pmra=80.0, pmdec=-40.0,
        ra_error=0.04, dec_error=0.04, parallax_error=0.03,
        pmra_error=0.04, pmdec_error=0.04, inflation_factor=1.37)],
    variables=@variables begin
        mass = $_M_comp_mjup
        a = 20.0; e = 0.3; i = 0.6; ω = 0.4; Ω = 1.0
        M = system.M
        tp = 50000.0
    end)
gen_sys = System(
    name="darkmass_gen", companions=[gen_planet], observations=[],
    variables=@variables begin
        M_central = $_M_central_true
        M = M_central + $_M_comp_msol
        plx = 50.0; pmra = 80.0; pmdec = -40.0; rv = 10_000.0
        ra = $_ra0; dec = $_dec0; ref_epoch = $_REF
    end)
sim_sys = Octofitter.generate_from_params(gen_sys; add_noise=true)
sim_obs = sim_sys.planets[1].observations[1]

# --- inference model: broad prior on the dark mass ---
fit_planet = Planet(
    name="B", basis=AbsoluteVisual{KepOrbit}, observations=[sim_obs],
    variables=@variables begin
        mass = $_M_comp_mjup
        a ~ truncated(Normal(20.0, 0.3), lower=1)
        e ~ truncated(Normal(0.3, 0.02), lower=0, upper=0.99)
        i ~ truncated(Normal(0.6, 0.02), lower=0, upper=pi)
        ω ~ Normal(0.4, 0.02)
        Ω ~ Normal(1.0, 0.02)
        M = system.M
        tp = 50000.0
    end)
fit_sys = System(
    name="darkmass_fit", companions=[fit_planet], observations=[],
    variables=@variables begin
        M_central ~ LogUniform(0.3, 5.0)
        M = M_central + $_M_comp_msol
        plx ~ Normal(50.0, 0.05)
        pmra ~ Normal(80.0, 0.2)
        pmdec ~ Normal(-40.0, 0.2)
        rv = 10_000.0
        ra = $_ra0; dec = $_dec0; ref_epoch = $_REF
    end)
fit_model = Octofitter.LogDensityModel(fit_sys)

initialize!(fit_model)
chain = octofit(fit_model; adaptation=400, iterations=1500, verbosity=0)

Mc = vec(chain[:M_central])
q = sort(Mc)
lo = q[max(1, round(Int, 0.16*length(q)))]
hi = q[round(Int, 0.84*length(q))]
mn = sum(Mc)/length(Mc)
@info "GaiaDR3CompanionObs recovery" recovered=mn ci68=(lo, hi) truth=_M_central_true

@testset "dark central mass recovery" begin
    @test lo < _M_central_true < hi
    @test abs(mn - _M_central_true) < 0.15
end
