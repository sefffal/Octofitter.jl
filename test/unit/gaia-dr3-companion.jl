# Unit tests for GaiaDR3CompanionObs — fast, deterministic, no network / MCMC.
using Octofitter
using Distributions
using PlanetOrbits
using Test

const _mjup2msol = Octofitter.mjup2msol
const _REF = Octofitter.meta_gaia_DR3.ref_epoch_mjd

# Representative system: dark central object + resolved bound stellar companion.
const _sys = (; plx=50.0, ra=150.0, dec=-30.0, pmra=80.0, pmdec=-40.0, rv=10_000.0, ref_epoch=_REF)
const _M_central = 1.2
const _planet = (; mass=0.3/_mjup2msol, M=_M_central + 0.3, a=20.0, e=0.3, i=0.6, ω=0.4, Ω=1.0, tp=50000.0)
const _orbit = AbsoluteVisual{KepOrbit}(; merge(_sys, _planet)...)

function _ctx(obs, sysv, planetv)
    θ_planet = planetv
    θ_system = merge(sysv, (; M=planetv.M, planets=(; B=θ_planet)))
    orbit = AbsoluteVisual{KepOrbit}(; merge(sysv, planetv)...)
    sols = [orbitsolve(orbit, e) for e in obs.table.epoch]
    Octofitter.PlanetObservationContext(θ_system, θ_planet, (;), (orbit,), (sols,), 1, 0)
end

_obs(; kw...) = GaiaDR3CompanionObs(;
    ra=_sys.ra, dec=_sys.dec, parallax=_sys.plx, pmra=_sys.pmra, pmdec=_sys.pmdec,
    ra_error=0.05, dec_error=0.05, parallax_error=0.03, pmra_error=0.05, pmdec_error=0.05,
    inflation_factor=1.0, kw...)

@testset "construction & inflation" begin
    o = GaiaDR3CompanionObs(; ra=_sys.ra, dec=_sys.dec, parallax=_sys.plx,
        pmra=_sys.pmra, pmdec=_sys.pmdec, ra_error=0.04, dec_error=0.04,
        parallax_error=0.03, pmra_error=0.05, pmdec_error=0.06, inflation_factor=1.37)
    @test o.σ_vec[3] ≈ 1.37*0.03 rtol=1e-10
    @test o.σ_vec[4] ≈ 1.37*0.05 rtol=1e-10
    @test length(o.table.epoch) == 3
    @test o.table.epoch[2] == _REF
    # Missing parameters must error in the user-specified path.
    @test_throws ArgumentError GaiaDR3CompanionObs(; ra=1.0, dec=2.0)
end

@testset "companion offset physics" begin
    sol = orbitsolve(_orbit, _REF)
    M_comp = _planet.mass*_mjup2msol
    M_host = totalmass(_orbit) - M_comp
    Δα_bary = (sol.compensated.ra2 - _sys.ra)*3.6e6*cosd(_sys.dec)
    Δδ_bary = (sol.compensated.dec2 - _sys.dec)*3.6e6
    Δα, Δδ, plx = Octofitter._dr3_companion_offset(sol, M_comp, _sys.ra, _sys.dec)
    @test (Δα - Δα_bary) ≈ (M_host/totalmass(_orbit))*raoff(sol)  rtol=1e-8
    @test (Δδ - Δδ_bary) ≈ (M_host/totalmass(_orbit))*decoff(sol) rtol=1e-8
    @test plx ≈ sol.compensated.parallax2
    # Massless companion sits at the full separation from the barycentre.
    Δα0, Δδ0, _ = Octofitter._dr3_companion_offset(sol, 0.0, _sys.ra, _sys.dec)
    @test (Δα0 - Δα_bary) ≈ raoff(sol)  rtol=1e-8
end

@testset "proper-motion finite difference" begin
    obs = _obs()
    ctx = _ctx(obs, _sys, _planet)
    sim = Octofitter.simulate(obs, ctx.θ_system, ctx.θ_planet, ctx.θ_obs, ctx.orbits,
        ctx.orbit_solutions, ctx.i_planet, ctx.orbit_solutions_i_epoch_start)
    M_comp = _planet.mass*_mjup2msol
    dt = 5.0
    a1 = Octofitter._dr3_companion_offset(orbitsolve(_orbit, _REF-dt), M_comp, _sys.ra, _sys.dec)
    a2 = Octofitter._dr3_companion_offset(orbitsolve(_orbit, _REF+dt), M_comp, _sys.ra, _sys.dec)
    @test sim.pmra  ≈ (a2[1]-a1[1])/((2dt)/365.25) rtol=2e-3
    @test sim.pmdec ≈ (a2[2]-a1[2])/((2dt)/365.25) rtol=2e-3
end

@testset "AbsoluteVisual orbit required" begin
    obs = _obs()
    vis = Visual{KepOrbit}(; plx=_sys.plx, M=_planet.M, a=_planet.a, e=_planet.e,
        i=_planet.i, ω=_planet.ω, Ω=_planet.Ω, tp=_planet.tp)
    sols = [orbitsolve(vis, e) for e in obs.table.epoch]
    θ_planet = _planet
    θ_system = merge(_sys, (; M=_planet.M, planets=(; B=θ_planet)))
    ctx = Octofitter.PlanetObservationContext(θ_system, θ_planet, (;), (vis,), (sols,), 1, 0)
    @test_throws Exception Octofitter.ln_like(obs, ctx)
end

@testset "profile likelihood peaks at injected dark mass" begin
    # Generate noiseless data at truth.
    obs0 = _obs()
    ctx0 = _ctx(obs0, _sys, _planet)
    sim = Octofitter.simulate(obs0, ctx0.θ_system, ctx0.θ_planet, ctx0.θ_obs, ctx0.orbits,
        ctx0.orbit_solutions, ctx0.i_planet, ctx0.orbit_solutions_i_epoch_start)
    obs = _obs(ra=_sys.ra+sim.Δα_mas/3.6e6/cosd(_sys.dec), dec=_sys.dec+sim.Δδ_mas/3.6e6,
        parallax=sim.parallax, pmra=sim.pmra, pmdec=sim.pmdec)
    Mgrid = range(0.8, 1.6, length=17)
    lls = map(Mgrid) do Mc
        pv = merge(_planet, (; M = Mc + _planet.mass*_mjup2msol))
        Octofitter.ln_like(obs, _ctx(obs, _sys, pv))
    end
    @test abs(Mgrid[argmax(lls)] - _M_central) < 0.06
end

@testset "optional Gaia RV term" begin
    obs = GaiaDR3CompanionObs(; ra=_sys.ra, dec=_sys.dec, parallax=_sys.plx,
        pmra=_sys.pmra, pmdec=_sys.pmdec, ra_error=0.05, dec_error=0.05,
        parallax_error=0.03, pmra_error=0.05, pmdec_error=0.05, inflation_factor=1.0,
        include_rv=true, radial_velocity=11.0, radial_velocity_error=0.5)
    ctx = _ctx(obs, _sys, _planet)
    sim = Octofitter.simulate(obs, ctx.θ_system, ctx.θ_planet, ctx.θ_obs, ctx.orbits,
        ctx.orbit_solutions, ctx.i_planet, ctx.orbit_solutions_i_epoch_start)
    sol0 = ctx.orbit_solutions[1][2]
    M_comp = _planet.mass*_mjup2msol
    M_host = totalmass(_orbit) - M_comp
    rv_expected = (sol0.compensated.rv2 + (M_host/totalmass(_orbit))*radvel(sol0.sol))/1000
    @test sim.rv_kms ≈ rv_expected rtol=1e-10
    # include_rv=true requires RV data
    @test_throws ArgumentError GaiaDR3CompanionObs(; ra=_sys.ra, dec=_sys.dec,
        parallax=_sys.plx, pmra=_sys.pmra, pmdec=_sys.pmdec, ra_error=0.05, dec_error=0.05,
        parallax_error=0.03, pmra_error=0.05, pmdec_error=0.05, include_rv=true)
end
