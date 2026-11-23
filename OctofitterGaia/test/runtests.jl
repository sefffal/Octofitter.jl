# OctofitterGaia — the emitter layer.
#
# These tests run on two small fixtures cut from the full products (see
# fixtures/make_fixtures.jl): a few MIST isochrones around solar age and
# metallicity, and the main-sequence part of a spectral library. They cover
# what the layer promises: that the tables are reproduced at their nodes,
# that everything in between is continuous and differentiable, that photon
# rates are on Gaia's scale, and that the pieces compose into an Octofitter
# model.

using Test
using OctofitterGaia
using Octofitter
using Octofitter: Table
using Distributions
using ForwardDiff
using Random

const OG = OctofitterGaia
const FIXTURES = joinpath(@__DIR__, "fixtures")

evol = MISTGrid(joinpath(FIXTURES, "mist_fixture.h5"))
lib = SpectralLibrary(joinpath(FIXTURES, "sed_fixture.h5"))

central(f, x, i; h=1e-6) = (f(x .+ h .* (eachindex(x) .== i)) .- f(x .- h .* (eachindex(x) .== i))) ./ 2h

@testset "OctofitterGaia" begin

    @testset "Catmull–Rom stencil" begin
        nodes = [-1.0, -0.4, 0.0, 0.25, 0.5, 1.5]          # deliberately non-uniform
        for x in (-1.0, -0.9, -0.4, -0.1, 0.1, 0.3, 0.9, 1.5)
            idx, w, _ = OG._stencil(nodes, x)
            @test sum(w) ≈ 1
            @test sum(w .* nodes[collect(idx)]) ≈ x atol = 1e-12    # reproduces linear functions
        end
        # interpolating: on a node, that node alone
        idx, w, t = OG._stencil(nodes, 0.25)
        @test sum(w .* (nodes[collect(idx)] .== 0.25)) ≈ 1
        @test t == 0 || t == 1
        # exact for quadratics away from the ends, where the tangents are centred
        idx, w, t = OG._stencil(nodes, 0.1)
        @test t ≈ 0.4
        @test sum(w .* nodes[collect(idx)] .^ 2) ≈ 0.1^2 atol = 1e-12
        @test OG._stencil(nodes, 1.6) === nothing
        @test OG._stencil(nodes, -1.1) === nothing
    end

    @testset "MIST evolution" begin
        # the stored isochrones are reproduced on the nodes
        worst = 0.0
        for j in eachindex(evol.feh), k in eachindex(evol.log_age)
            for e in evol.eep_lo[k, j]+1:9:evol.eep_hi[k, j]-1
                (evol.log_mass[e+1, k, j] > evol.log_mass[e, k, j] > evol.log_mass[e-1, k, j]) || continue
                s = evol(exp10(evol.log_mass[e, k, j]), evol.log_age[k], evol.feh[j])
                worst = max(worst, abs(log10(s.teff) - evol.log_teff[e, k, j]),
                    abs(s.logL - evol.log_L[e, k, j]), abs(s.logg - evol.logg[e, k, j]))
            end
        end
        @test worst < 1e-12

        # a solar-mass star at 4.5 Gyr. MIST's own model, not the Sun: it is
        # 9 % more luminous than the Sun at that age.
        s = evol(1.0, 9.65, 0.0)
        @test 5800 < s.teff < 5900
        @test 1.0 < s.radius < 1.04
        @test s.mass_now ≈ 1.0 atol = 1e-3
        # log g, radius and present mass are one consistent set
        @test s.logg ≈ OG.LOGG_SUN + log10(s.mass_now / s.radius^2)
        @test s.logL ≈ 2 * log10(s.radius) + 4 * log10(s.teff / OG.TEFF_SUN)

        # a subgiant, past the terminal-age main sequence (EEP 454)
        s = evol(1.2, 9.75, 0.0)
        @test s.eep > 454
        @test s.logg < 3.6
        @test s.radius > 3

        # outside the grid, and dead stars
        @test isnan(evol(3.0, 9.7, 0.0).teff)       # long dead at 5 Gyr
        @test isnan(evol(0.09, 9.7, 0.0).teff)      # below the 0.1 M☉ floor
        @test isnan(evol(1.0, 9.9, 0.0).teff)       # age outside the fixture
        @test isnan(evol(1.0, 9.7, 0.4).teff)       # [Fe/H] outside the fixture
        @test isnan(evol(-1.0, 9.7, 0.0).teff)

        f(v) = (s = evol(v[1], v[2], v[3]); [s.teff, s.logg, s.radius, s.mass_now, s.mh])

        # continuous across age and metallicity nodes, below and above the
        # intermediate-age main sequence, where the correspondence between
        # isochrones changes
        for m in (0.3, 0.7, 1.0, 1.1, 1.2)
            for (x, i) in (([m, 9.65, 0.07], 2), ([m, 9.63, 0.0], 3))
                lo = f(x .- 1e-9 .* (1:3 .== i))
                hi = f(x .+ 1e-9 .* (1:3 .== i))
                @test all(isapprox.(lo, hi; rtol=1e-6))
            end
        end

        # … and with a continuous gradient in age and metallicity there
        for m in (0.3, 0.7)
            for (x, i) in (([m, 9.65, 0.07], 2), ([m, 9.63, 0.0], 3))
                below = ForwardDiff.jacobian(f, x .- 1e-7 .* (1:3 .== i))
                above = ForwardDiff.jacobian(f, x .+ 1e-7 .* (1:3 .== i))
                @test all(isapprox.(below[:, 2:3], above[:, 2:3]; rtol=1e-3, atol=1e-6))
            end
        end

        # continuous in mass from the main sequence through the IAMS to the
        # giant branch: halving the mass step halves the largest step in
        # radius, which a jump anywhere would not allow
        largest_step(n) = maximum(abs, diff([evol(m, 9.71, 0.06).radius for m in range(0.95, 1.25; length=n)]))
        @test largest_step(8001) < 0.6 * largest_step(4001) < 0.36 * largest_step(2001)
        eeps = [evol(m, 9.71, 0.06).eep for m in range(0.95, 1.25; length=2001)]
        @test !any(isnan, eeps)
        @test issorted(eeps)
        @test minimum(eeps) < OG.MIST_EEP_IAMS && maximum(eeps) > 454
        # the bottom of the grid is reachable, at every age and metallicity
        @test evol(0.1, 9.71, 0.06).teff > 2500
        @test evol(0.1, 9.65, 0.0).teff > 2500

        # derivatives: ForwardDiff against central differences, below and above the IAMS
        for x in ([0.6, 9.62, -0.1], [1.05, 9.72, 0.12], [1.2, 9.74, 0.05])
            J = ForwardDiff.jacobian(f, x)
            Jfd = hcat((central(f, x, i) for i in 1:3)...)
            @test all(isfinite, J)
            @test isapprox(J, Jfd; rtol=2e-4, atol=1e-5)
        end
    end

    @testset "passbands" begin
        pb = gaia_passbands()
        @test pb.G.zp ≈ 25.6874 atol = 1e-4
        @test pb.BP.zp ≈ 25.3385 atol = 1e-4
        @test pb.RP.zp ≈ 24.7479 atol = 1e-4
        # The published AB zero points follow from the curves and the pupil
        # area alone: the rate of a flat f_ν source of AB magnitude 0.
        f = SEDFunctionals(collect(320.0:0.5:1100.0))
        λ = (f.edges[1:end-1] .+ f.edges[2:end]) ./ 2 .* 1e-7          # cm
        fλ = exp10(-48.6 / 2.5) .* OG.C_LIGHT ./ λ .^ 2
        @test 2.5 * log10(photon_rate(f, :G, fλ)) ≈ 25.8010 atol = 1e-3
        @test 2.5 * log10(photon_rate(f, :BP, fλ)) ≈ 25.3540 atol = 1e-3
        @test 2.5 * log10(photon_rate(f, :RP, fλ)) ≈ 25.1040 atol = 1e-3
        # G_RVS: the band-averaged flux density of Vega is 2219 Jy (SVO), which
        # must come out near magnitude 0 with the published zero point
        @test OG.magnitude(pb.RVS, photon_rate(f, :RVS, fλ .* 2219.05e-23 ./ exp10(-48.6 / 2.5))) ≈ 0 atol = 0.01
        # the n = 0 chromatic row is the G rate
        @test f.chroma[1, :] == f.G
        # coarse bins lose nothing for a smooth spectrum
        fc = SEDFunctionals(collect(320.0:5.0:1100.0))
        λc = (fc.edges[1:end-1] .+ fc.edges[2:end]) ./ 2 .* 1e-7
        @test photon_rate(fc, :G, exp10(-48.6 / 2.5) .* OG.C_LIGHT ./ λc .^ 2) ≈ photon_rate(f, :G, fλ) rtol = 1e-4
    end

    @testset "spectral library" begin
        n = OG.nbins(lib)
        sed = zeros(n)
        # nodes are reproduced
        for (i, j, k) in ((11, 2, 2), (1, 1, 1), (length(lib.teff), length(lib.logg), length(lib.mh)))
            surface_sed!(sed, lib, lib.teff[i], lib.logg[j], lib.mh[k])
            @test sed ≈ exp10.(Float64.(lib.log_flux[:, i, j, k])) rtol = 1e-6
        end
        # outside the grid
        @test all(isnan, surface_sed!(sed, lib, 2900.0, 4.5, 0.0))
        @test all(isnan, surface_sed!(sed, lib, 5000.0, 3.9, 0.0))
        @test all(isnan, surface_sed!(sed, lib, 5000.0, 4.5, 0.6))
        @test_throws DimensionMismatch surface_sed!(zeros(n + 1), lib, 5000.0, 4.5, 0.0)

        # a missing model makes the cells it touches invalid, and only those
        holed = SpectralLibrary(lib.name, lib.teff, lib.logg, lib.mh, lib.edges, lib.log_flux,
            copy(lib.present), lib.path)
        holed.present[11, 2, 2] = false
        @test all(isnan, surface_sed!(sed, holed, lib.teff[11] + 10, lib.logg[2] + 0.1, lib.mh[2] + 0.1))
        @test !any(isnan, surface_sed!(sed, holed, lib.teff[13] + 10, lib.logg[2] + 0.1, lib.mh[2] + 0.1))

        # more flux from a hotter star, in every bin
        cool = copy(surface_sed!(sed, lib, 4000.0, 4.6, 0.0))
        hot = copy(surface_sed!(sed, lib, 4100.0, 4.6, 0.0))
        @test all(hot .> cool)

        g(v) = (T = eltype(v); out = zeros(T, n); surface_sed!(out, lib, v[1], v[2], v[3]); log.(out[[5, 60, 150]]))
        # ForwardDiff against central differences
        x = [4321.0, 4.63, 0.13]
        @test isapprox(ForwardDiff.jacobian(g, x), hcat((central(g, x, i; h=1e-4) for i in 1:3)...); rtol=1e-5, atol=1e-9)
        # continuous first derivatives across a node of each axis
        for (x, i, h) in (([4300.0, 4.63, 0.13], 1, 1e-6), ([4321.0, 4.5, 0.13], 2, 1e-9), ([4321.0, 4.63, 0.0], 3, 1e-9))
            below = ForwardDiff.jacobian(g, x .- h .* (1:3 .== i))
            above = ForwardDiff.jacobian(g, x .+ h .* (1:3 .== i))
            @test isapprox(below, above; rtol=1e-5, atol=1e-9)
        end
    end

    @testset "emitters" begin
        star = Star(evol, lib)
        f = SEDFunctionals(lib)
        n = OG.nbins(star)
        @test collect(scene_points(star, (;), (;), 57000.0)) == [(0.0, 0.0, 1.0)]

        body = (; teff=5772.0, logg=OG.LOGG_SUN, mh=0.0, radius=1.0)
        sys = (; plx=100.0, av=0.0)
        sun = sed!(zeros(n), star, body, sys)

        # A star with the Sun's parameters at 10 pc. The literature value is
        # M_G = 4.67 and BP - RP = 0.82 (Casagrande & VandenBerg 2018).
        MG = OG.magnitude(f.passbands.G, photon_rate(f, :G, sun))
        BPRP = OG.magnitude(f.passbands.BP, photon_rate(f, :BP, sun)) -
               OG.magnitude(f.passbands.RP, photon_rate(f, :RP, sun))
        @test MG ≈ 4.67 atol = 0.04
        @test BPRP ≈ 0.82 atol = 0.03

        # dilution
        @test sed!(zeros(n), star, body, (; plx=50.0, av=0.0)) ≈ sun ./ 4
        @test sed!(zeros(n), star, merge(body, (; radius=2.0)), sys) ≈ sun .* 4
        # extinction: A_G / A_0 for a solar-type star is about 0.84
        red = sed!(zeros(n), star, body, (; plx=100.0, av=1.0))
        @test all(red .< sun)
        @test 0.78 < -2.5 * log10(photon_rate(f, :G, red) / photon_rate(f, :G, sun)) < 0.90
        @test red[1] / sun[1] < red[end] / sun[end]            # reddening, not grey

        # BandRate is the same number
        for band in (:G, :BP, :RP, :RVS)
            r = BandRate(star, band; functionals=f)
            @test r(body.teff, body.logg, body.mh, body.radius, sys.plx, sys.av) ≈ photon_rate(f, band, sun)
        end
        rate_G = BandRate(star, :G)
        h(v) = rate_G(v...)
        x = [4321.0, 4.63, 0.13, 0.7, 40.0, 0.1]
        @test isapprox(ForwardDiff.gradient(h, x), [central(h, x, i; h=1e-5)[1] for i in 1:6]; rtol=1e-5)
        @test (@allocated rate_G(4321.0, 4.63, 0.13, 0.7, 40.0, 0.1)) == 0

        # chromatic moments: Legendre moments lie in [-1, 1], a redder star has
        # more of its photons at low wavenumber, and the moments of a blend
        # are the rate-weighted means of its components' — which is the
        # reason the image model carries moments and not a colour
        cool = sed!(zeros(n), star, (; teff=3400.0, logg=4.9, mh=0.0, radius=0.35), sys)
        c_sun, c_cool, c_blend = zeros(6), zeros(6), zeros(6)
        r_sun = chromatic_moments!(c_sun, f, sun)
        r_cool = chromatic_moments!(c_cool, f, cool)
        r_blend = chromatic_moments!(c_blend, f, sun .+ cool)
        @test r_sun ≈ photon_rate(f, :G, sun)
        @test all(abs.(c_sun) .<= 1) && all(abs.(c_cool) .<= 1)
        @test c_cool[1] < c_sun[1]
        @test r_blend ≈ r_sun + r_cool
        @test c_blend ≈ (r_sun .* c_sun .+ r_cool .* c_cool) ./ (r_sun + r_cool)
    end

    # Magnitudes integrated from the library against MIST's own, along the
    # main sequence of an isochrone. The two use different model atmospheres
    # and the comparison is not expected to be exact. What PHOENIX-ACES gives,
    # and what the bounds below record, is G fainter than MIST's by 0.04 to
    # 0.06 mag at 4500–6200 K, falling to zero at 3500–4000 K, with BP - RP
    # within 0.02 mag throughout. (Below 3500 K, not tested here, the library
    # is up to 0.05 mag brighter in G and 0.2 mag bluer.) A change of flux
    # scale, units or passband would break this.
    @testset "library magnitudes against MIST" begin
        star = Star(evol, lib)
        f = SEDFunctionals(lib)
        mist = OG.mist_magnitudes(evol)
        j = findfirst(==(0.0), evol.feh)
        k = findfirst(≈(9.7), evol.log_age)
        sed = zeros(OG.nbins(lib))
        warm = (; dG=Float64[], dcol=Float64[])        # 4500–6200 K
        cool = (; dG=Float64[], dcol=Float64[])        # 3500–4500 K
        for e in evol.eep_lo[k, j]:evol.eep_hi[k, j]
            e - 1 < 454 || break               # main sequence only
            teff = exp10(evol.log_teff[e, k, j])
            radius = exp10(evol.log_L[e, k, j] / 2) * (OG.TEFF_SUN / teff)^2
            OG.observed_sed!(sed, star, teff, evol.logg[e, k, j], evol.feh_surf[e, k, j], radius, 100.0, 0.0)
            isnan(sed[1]) && continue
            G, BP, RP = (OG.magnitude(f.passbands[b], photon_rate(f, b, sed)) for b in (:G, :BP, :RP))
            bin = 4500 <= teff <= 6200 ? warm : 3500 <= teff < 4500 ? cool : continue
            push!(bin.dG, G - mist.MG[e, k, j])
            push!(bin.dcol, (BP - RP) - (mist.MBP[e, k, j] - mist.MRP[e, k, j]))
        end
        @test length(warm.dG) > 50 && length(cool.dG) > 30
        @test all(x -> 0.03 < x < 0.07, warm.dG)
        @test all(x -> -0.01 < x < 0.045, cool.dG)
        @test all(x -> abs(x) < 0.03, warm.dcol)
        @test all(x -> abs(x) < 0.03, cool.dcol)
    end

    @testset "in an Octofitter model" begin
        star = Star(evol, lib)
        rate_G = BandRate(star, :G)
        A = Body(name="A", variables=@variables begin
            mass ~ truncated(Normal(0.8, 0.05), lower=0.6, upper=1.0)
            star = $evol(mass, system.logage, system.feh)
            teff = star.teff
            logg = star.logg
            mh = star.mh
            radius = star.radius
            flux_G = $rate_G(teff, logg, mh, radius, system.plx, system.av)
        end)
        b = Body(name="b", about=A, variables=@variables begin
            mass ~ Uniform(0.4, 0.6)
            star = $evol(mass, system.logage, system.feh)
            teff = star.teff
            logg = star.logg
            mh = star.mh
            radius = star.radius
            flux_G = $rate_G(teff, logg, mh, radius, system.plx, system.av)
            a ~ LogUniform(1.0, 10.0)
            e ~ Uniform(0.0, 0.5)
            i ~ Sine()
            ω ~ Uniform(0, 2pi)
            Ω ~ Uniform(0, 2pi)
            θ ~ Uniform(0, 2pi)
            epoch = 57000.0
        end)
        phot = PhotometryObs(Table(band=[:G], phot=[9.0e6], σ_phot=[9.0e4]); target=A, band=:G, name="Gphot")
        sys = System(name="emitters", bodies=[A, b], observations=[phot], variables=@variables begin
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            logage ~ Uniform(9.6, 9.75)
            feh ~ truncated(Normal(0.0, 0.05), lower=-0.2, upper=0.2)
            av = 0.0
        end)
        model = Octofitter.LogDensityModel(sys; verbosity=0)
        θ = Octofitter.sample_priors(Random.Xoshiro(2), sys)
        nt = model.arr2nt(θ)
        # the derived variables are the emitter's inputs, and the flux is what
        # the emitter itself computes from them
        @test nt.bodies.A.teff == evol(nt.bodies.A.mass, nt.logage, nt.feh).teff
        @test nt.bodies.A.flux_G > nt.bodies.b.flux_G > 0
        f = SEDFunctionals(lib)
        sed = sed!(zeros(OG.nbins(star)), star, nt.bodies.b, nt)
        @test photon_rate(f, :G, sed) ≈ nt.bodies.b.flux_G
        # the log density and its gradient
        θt = model.link(θ)
        @test isfinite(model.ℓπcallback(θt))
        v, g = model.∇ℓπcallback(θt)
        @test isfinite(v) && all(isfinite, g)
        gfd = [central(x -> [model.ℓπcallback(x)], θt, i; h=1e-6)[1] for i in eachindex(θt)]
        @test isapprox(g, gfd; rtol=1e-3, atol=1e-4)
    end
end
