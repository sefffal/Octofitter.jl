# OctofitterGaia: the emitter layer, the RVS likelihood and the image basis.
#
# These tests run on small fixtures (see fixtures/): a few MIST isochrones
# around solar age and metallicity, the main-sequence part of a spectral
# library, four BHAC15 isochrones, and spectra and images simulated by
# GaiaTwin. They cover what the package promises: that the tables are
# reproduced at their nodes, that everything in between is continuous and
# differentiable, that photon rates are on Gaia's scale, that the likelihood
# of the RVS epoch spectra is the Gaussian it claims to be and recovers what
# GaiaTwin put in, that the image model's mean is GaiaTwin's, and that the
# pieces compose into an Octofitter model.

using Test
using OctofitterGaia
using Octofitter
using Octofitter: Table
using Distributions
using ForwardDiff
using Random
using LinearAlgebra
import Bumper
using Bumper: @no_escape, @alloc

const OG = OctofitterGaia
const FIXTURES = joinpath(@__DIR__, "fixtures")

evol = MISTGrid(joinpath(FIXTURES, "mist_fixture.h5"))
lib = SpectralLibrary(joinpath(FIXTURES, "sed_fixture.h5"))
rvs = RVSLibrary(joinpath(FIXTURES, "rvs_fixture.h5"))

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

    @testset "RVS library" begin
        n, nc = OG.nwave(rvs), OG.ncomponents(rvs)
        λ = OG.wavelengths(rvs)
        @test length(λ) == n
        @test all(diff(log.(λ)) .≈ rvs.dln_wave)
        @test rvs.dln_wave * OG.C_LIGHT / 1e5 ≈ 2.0             # km/s per sample
        # the eigenspectra are orthonormal over the full window; over the
        # fixture's part of it their Gram matrix is at least positive definite
        @test all(>(0), [sum(abs2, rvs.basis[c, :]) for c in 1:nc])

        # weights are reproduced on the nodes
        w = zeros(nc)
        rvs_weights!(w, rvs, rvs.teff[11], rvs.logg[2], rvs.mh[2])
        @test w ≈ rvs.weights[:, 11, 2, 2]
        @test all(isnan, rvs_weights!(w, rvs, 2900.0, 4.5, 0.0))
        @test all(isnan, rvs_weights!(w, rvs, 5000.0, 4.5, 0.6))
        @test_throws DimensionMismatch rvs_weights!(zeros(nc + 1), rvs, 5000.0, 4.5, 0.0)

        # a K dwarf: the deepest feature of the window is the middle line of
        # the Ca II triplet, at 854.44 nm in vacuum
        spec = rvs_surface_spectrum!(zeros(n), rvs, 4800.0, 4.6, 0.0)
        @test all(>(0), spec)
        @test λ[argmin(spec)] ≈ 854.44 atol = 0.03
        @test minimum(spec) < 0.6 * maximum(spec)
        # the other line in the fixture's range, at 850.04 nm
        near = findall(x -> 849.5 < x < 850.5, λ)
        @test λ[near[argmin(spec[near])]] ≈ 850.04 atol = 0.03

        # the two node products are on one flux scale: the RVS spectrum
        # averaged over a bin of the low-resolution SED is that bin's value.
        # The fixture keeps 12 components, good to a few per cent.
        sed = surface_sed!(zeros(OG.nbins(lib)), lib, 4800.0, 4.6, 0.0)
        b = findfirst(==(850.0), lib.edges)
        inbin = findall(x -> lib.edges[b] <= x < lib.edges[b+1], λ)
        @test length(inbin) > 800
        @test sum(spec[inbin] .* λ[inbin]) / sum(λ[inbin]) ≈ sed[b] rtol = 0.03

        # ForwardDiff against central differences, and continuous first
        # derivatives across a node of each axis
        g(v) = (T = eltype(v); out = zeros(T, n); rvs_surface_spectrum!(out, rvs, v[1], v[2], v[3]); log.(out[[200, 700, 1100]]))
        x = [4321.0, 4.63, 0.13]
        @test isapprox(ForwardDiff.jacobian(g, x), hcat((central(g, x, i; h=1e-4) for i in 1:3)...); rtol=1e-5, atol=1e-9)
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

    @testset "banded Cholesky" begin
        rng = Random.Xoshiro(3)
        n, b = 40, 3
        M = [abs(i - j) <= b ? randn(rng) : 0.0 for i in 1:n, j in 1:n]
        A = M * M' + 5I                               # bandwidth 2b
        F = OG.BandedCholesky(A, 2b)
        L = cholesky(Symmetric(A)).L
        @test F.logdet ≈ logdet(A)
        x = randn(rng, n)
        @test OG.whiten!(copy(x), F) ≈ L \ x
        @test OG.colour!(copy(x), F) ≈ L * x
        # dual numbers go through the substitution
        d = ForwardDiff.derivative(t -> sum(abs2, OG.whiten!(x .* t, F)), 1.0)
        @test d ≈ 2 * dot(x, A \ x)
        @test_throws ArgumentError OG.BandedCholesky(A - 100I, 2b)
        # the small dense factorisation used for the marginalised modes
        S = Matrix(Symmetric(A[1:6, 1:6]))
        C = copy(S)
        @test OG._cholesky_small!(C) ≈ logdet(S)
        @test OG._forward_small!(x[1:6], C) ≈ cholesky(S).L \ x[1:6]
        @test OG._cholesky_small!(-S) == -Inf
    end

    @testset "RVS components" begin
        grid = LogWaveGrid(rvs)
        @test length(grid) == OG.nwave(rvs)
        h = OG.velocity_step(grid)
        @test h ≈ 2.0 atol = 1e-9
        raw = rvs_surface_spectrum!(zeros(length(grid)), rvs, 5777.0, 4.44, 0.0)
        # below half a sample nothing happens; flux is conserved; lines get shallower
        @test rotational_broaden!(similar(raw), raw, 0.9, grid) == raw
        b20 = rotational_broaden!(similar(raw), raw, 20.0, grid)
        inner = 30:length(raw)-30
        @test sum(b20[inner]) ≈ sum(raw[inner]) rtol = 1e-3
        @test minimum(b20) > minimum(raw)
        # the kernel is that of a limb-darkened sphere: a delta function comes
        # out with the profile's variance, vsini² ((1 - ε)/4 + 2ε/15) / (1 - ε/3) for ε = 0.6
        δ = zeros(length(grid))
        δ[200] = 1.0
        k = rotational_broaden!(similar(δ), δ, 30.0, grid)
        @test sum(k) ≈ 1
        v = ((1:length(k)) .- 200) .* h
        @test sum(k .* v) ≈ 0 atol = 1e-9
        @test sum(k .* v .^ 2) ≈ 30.0^2 * (0.4 / 4 + 1.2 / 15) / (1 - 0.6 / 3) rtol = 3e-3
        # continuous, with a continuous derivative, in vsini
        s(vs) = rotational_broaden!(zeros(typeof(vs), length(raw)), raw, vs, grid)[300]
        for vs in (1.0, 3.0, 5.0, 12.0)          # 1, 3, 5: where the kernel gains samples
            @test s(vs - 1e-7) ≈ s(vs + 1e-7) rtol = 1e-9
            @test ForwardDiff.derivative(s, vs - 1e-6) ≈ ForwardDiff.derivative(s, vs + 1e-6) rtol = 1e-3 atol = 1e-3 * abs(raw[300])
        end
        @test all(isfinite, ForwardDiff.derivative(vs -> rotational_broaden!(zeros(typeof(vs), length(raw)), raw, vs, grid), 3.0))
        # derivatives in velocity: a shift by δv changes the sampled spectrum by dflux δv
        c = rvs_component!(similar(raw), similar(raw), similar(raw), raw, 6.0, grid)
        @test c.level ≈ sum(c.flux) / length(c.flux)
        i = 300
        @test c.dflux[i] ≈ -(c.flux[i+1] - c.flux[i-1]) / 2h
        @test c.d2flux[i] ≈ (c.flux[i+1] - 2c.flux[i] + c.flux[i-1]) / h^2
        # physics helpers
        @test gravitational_redshift(4.438, 1.0) ≈ 0.636 atol = 2e-3
        inst = RVSInstrument()
        @test rvs_velocity(inst, 10.0, 200.0) ≈ 10.0 + OG.RVS_KAPPA * 200
        @test abs(rvs_velocity(inst, 0.0, 200.0)) > 26       # 200 mas along scan is more than a resolution element
        @test window_fraction(inst, 0.0) ≈ 1 atol = 1e-6
        @test window_fraction(inst, inst.ac_halfwidth) ≈ 0.5 atol = 1e-6
        @test window_fraction(inst, 300.0) > window_fraction(inst, 800.0) > window_fraction(inst, 1200.0)
    end

    @testset "RVS epoch spectra" begin
        grid = LogWaveGrid(rvs)
        inst = RVSInstrument()
        wave = collect(849.0:0.025:855.0)
        function component(teff, logg, mh, vsini, scale=1.0)
            T = promote_type(typeof(teff), typeof(logg), typeof(mh), typeof(vsini), typeof(scale), Float64)
            raw = rvs_surface_spectrum!(Vector{T}(undef, length(grid)), rvs, teff, logg, mh) .* scale
            rvs_component!(similar(raw), similar(raw), similar(raw), raw, vsini, grid; inst.limb_darkening)
        end
        cA, cB = component(5800.0, 4.4, 0.0, 5.0), component(4000.0, 4.7, 0.0, 3.0, 0.3)
        nt = 6
        vA = [10.0, -20.0, 35.0, 0.0, 5.0, -12.0]
        vB = [-40.0, 60.0, -80.0, 20.0, -5.0, 44.0]
        σ = fill(0.01, length(wave), nt)
        blank = RVSEpochSpectra(ones(length(wave), nt), σ, grid; instrument=inst, wavelength=wave)
        rng = Random.Xoshiro(11)

        # the model spectrum: normalised, and a single star shifts as a whole
        m0 = rvs_transit_model(blank, 1, (cA,), (0.0,), (1.0,))
        @test 0.9 < sum(m0) / length(m0) < 1.1
        @test rvs_transit_model(blank, 1, (cA,), (0.0,), (7.0,)) ≈ m0           # one emitter: its weight cancels
        dv = 0.1
        m1 = rvs_transit_model(blank, 1, (cA,), (dv,), (1.0,))
        B = zeros(length(wave), OG.nmodes(inst))
        OG._rvs_design!(B, blank, 1, (cA,), (0.0,), (1.0,))
        @test maximum(abs, m1 .- m0 .- dv .* B[:, end-1]) < 0.02 * maximum(abs, m1 .- m0)     # the zero-point mode is ∂m/∂v
        # the width mode is what a broader line-spread function does
        @test maximum(abs, B[:, end]) > 0.05
        # a blend is the weighted sum of its parts
        mA = rvs_transit_model(blank, 1, (cA,), (vA[1],), (1.0,))
        mB = rvs_transit_model(blank, 1, (cB,), (vB[1],), (1.0,))
        mAB = rvs_transit_model(blank, 1, (cA, cB), (vA[1], vB[1]), (1.0, 0.9))
        fB = 0.9 * cB.level / (cA.level + 0.9 * cB.level)
        @test mAB ≈ (1 - fB) .* mA .+ fB .* mB
        # velocities that take the bins off the library grid are impossible, not an error
        @test first(rvs_transit_loglike(RVSEpochSpectra(ones(length(wave), 1), σ[:, 1:1], grid; instrument=inst, wavelength=wave, p_bad=0.0),
            1, (cA,), (1500.0,), (1.0,))) == -Inf

        # simulate, with gaps, and compare with the dense Gaussian the model defines
        flux = reduce(hcat, [rvs_transit_simulate(rng, blank, t, (cA, cB), (vA[t], vB[t]), (1.0, 0.9)) for t in 1:nt])
        flux[40:47, 2] .= NaN
        data = RVSEpochSpectra(flux, σ, grid; instrument=inst, wavelength=wave)
        @test OG.ntransits(data) == nt
        @test length(data.transits[2].bins) == length(wave) - 8
        for t in (1, 2)
            tr = data.transits[t]
            n, J = length(tr.bins), OG.nmodes(inst)
            Bt = zeros(n, J)
            OG._rvs_design!(Bt, data, t, (cA, cB), (vA[t], vB[t]), (1.0, 0.9))
            K, _ = OG._rvs_covariance(tr.bins, view(σ, :, t), view(data.mismatch, :, t), wave, inst)
            Λ = Diagonal([1 / OG._mode_precision(inst, j) for j in 1:J])
            C = Symmetric(K, :L) + Bt * Λ * Bt'
            r = flux[tr.bins, t] .- Bt[:, 1]
            dense = -(dot(r, C \ r) + logdet(C) + n * log(2π)) / 2
            ll, ll_good = rvs_transit_loglike(data, t, (cA, cB), (vA[t], vB[t]), (1.0, 0.9))
            @test ll_good ≈ dense rtol = 1e-9
            @test ll ≈ log((1 - 1e-3) * exp(ll_good - dense) + 1e-3 * exp(tr.ll_bad - dense)) + dense
            @test tr.ll_bad < ll_good - 100          # a good spectrum is not mistaken for a bad one
        end
        @test rvs_transit_simulate(rng, data, 3, (cA, cB), (vA[3], vB[3]), (1.0, 0.9); add_noise=false) ==
              rvs_transit_model(data, 3, (cA, cB), (vA[3], vB[3]), (1.0, 0.9))
        # a spectrum unlike the model falls to the bad state instead of -Inf
        junk = copy(flux)
        junk[:, 1] .= 1 .+ 0.5 .* randn(rng, length(wave))
        jd = RVSEpochSpectra(junk, σ, grid; instrument=inst, wavelength=wave)
        llj, llj_good = rvs_transit_loglike(jd, 1, (cA, cB), (vA[1], vB[1]), (1.0, 0.9))
        @test llj ≈ log(1e-3) + jd.transits[1].ll_bad
        @test llj_good < llj - 1000
        # subsetting by transit
        sub = data[[2, 5]]
        @test first(rvs_transit_loglike(sub, 1, (cA, cB), (vA[2], vB[2]), (1.0, 0.9))) ==
              first(rvs_transit_loglike(data, 2, (cA, cB), (vA[2], vB[2]), (1.0, 0.9)))
        @test_throws ArgumentError RVSEpochSpectra(fill(NaN, length(wave), 1), σ[:, 1:1], grid; wavelength=wave)

        # the total as a function of every kind of parameter: gradient, and the truth is near the top
        total(x) = begin
            a, b = component(x[1], x[2], x[3], x[4]), component(x[5], x[6], x[7], x[8], x[9])
            sum(first(rvs_transit_loglike(data, t, (a, b), (vA[t] + x[10], vB[t] + x[11]), (one(x[12]), x[12]))) for t in 1:nt)
        end
        x0 = [5800.0, 4.4, 0.0, 5.0, 4000.0, 4.7, 0.0, 3.0, 0.3, 0.0, 0.0, 0.9]
        g = ForwardDiff.gradient(total, x0)
        @test all(isfinite, g)
        gfd = map(eachindex(x0)) do i
            e = zeros(length(x0))
            e[i] = 1e-5 * max(abs(x0[i]), 1.0)
            (total(x0 .+ e) - total(x0 .- e)) / 2e[i]
        end
        @test isapprox(g, gfd; rtol=2e-3, atol=2e-3)
        H = ForwardDiff.hessian(x -> total(vcat(x0[1:9], x, x0[12])), [0.0, 0.0])
        step = -(H \ g[10:11])
        σv = sqrt.(diag(inv(-H)))
        @test all(abs.(step) .< 4 .* σv)
        @test σv[1] < 0.5 && σv[2] < 3.0
        @test total(setindex!(copy(x0), 3.0, 10)) < total(x0) - 20
        @test (@allocated rvs_transit_loglike(data, 1, (cA, cB), (vA[1], vB[1]), (1.0, 0.9))) == 0
    end

    # GaiaTwin simulates a double-lined pair: `composite_pixels` is its
    # noise-free blend on detector pixels, `dr4_epoch_spectrum` the published
    # format (three CCDs at random sub-pixel phases, noise, linear
    # interpolation onto the 961 bins, average, normalisation). The fixture
    # holds GaiaTwin's own library spectra, so the instrument is the same and
    # what is tested is the block's arithmetic and its noise model.
    @testset "RVS epoch spectra against GaiaTwin" begin
        fx = OG.h5open(joinpath(FIXTURES, "gaiatwin_rvs_fixture.h5")) do f
            (; (Symbol(k) => read(f[k]) for k in keys(f))...)
        end
        grid = LogWaveGrid(fx.ln_wave0, fx.dln_wave, size(fx.spectra, 1))
        nt = size(fx.velocity, 2)
        c = ntuple(2) do k
            raw = Float64.(fx.spectra[:, k])
            rvs_component!(similar(raw), similar(raw), similar(raw), raw, 0.0, grid)
        end
        v(t, dv=(0.0, 0.0)) = (fx.velocity[1, t] + dv[1], fx.velocity[2, t] + dv[2])
        w(t, f=1.0) = (fx.electrons[1, t], f * fx.electrons[2, t])

        # the mean, on GaiaTwin's detector pixels: equal up to the overall scale.
        # (GaiaTwin interpolates its library linearly, the block cubically.)
        pix = RVSEpochSpectra(ones(length(fx.pixel_wave), nt), fill(0.01, length(fx.pixel_wave), nt), grid; wavelength=fx.pixel_wave)
        for t in 1:nt
            ratio = fx.composite[:, t] ./ rvs_transit_model(pix, t, c, v(t), w(t))
            @test maximum(abs, ratio ./ (sum(ratio) / length(ratio)) .- 1) < 1e-4
        end

        # the published format. GaiaTwin's spectra have no template mismatch.
        inst = RVSInstrument(mismatch_amplitude=1e-4)
        data = RVSEpochSpectra(Float64.(fx.flux), Float64.(fx.flux_error), grid; instrument=inst)
        @test data.wavelength == fx.wave
        ll(x) = sum(first(rvs_transit_loglike(data, t, c, v(t, (x[1], x[2])), w(t, x[3]))) for t in 1:nt)
        x = [0.0, 0.0, 1.0]
        H = ForwardDiff.hessian(ll, x)
        for _ in 1:4
            H = ForwardDiff.hessian(ll, x)
            x -= H \ ForwardDiff.gradient(ll, x)
        end
        σx = sqrt.(diag(inv(-H)))
        # both velocities and the flux ratio are recovered within their errors
        # (over 40 realisations of this scene the pulls have mean 0.2, -0.1,
        # 0.4 and standard deviation 0.6, 1.0, 0.95)
        @test all(abs.(x .- [0, 0, 1]) .< 3.5 .* σx)
        @test σx[1] < 0.25 && σx[2] < 0.7 && σx[3] < 0.03
        # the normalisation and the smoothing by interpolation land in the modes
        tr = data.transits[1]
        J = OG.nmodes(inst)
        B = zeros(length(tr.bins), J)
        OG._rvs_design!(B, data, 1, c, v(1), w(1))
        Z = copy(B)
        foreach(j -> OG.whiten!(view(Z, :, j), tr.chol), 1:J)
        P = Diagonal([OG._mode_precision(inst, j) for j in 1:J])
        ĉ = (P + Z'Z) \ (Z' * tr.white + P * [1; zeros(J - 1)])
        @test 0.9 < ĉ[1] < 1.0                   # the brightest 30 % of bins sit above the mean
        @test 0.03 < ĉ[end] < 0.15               # interpolation broadens the line-spread function
        # residuals: unit variance, neighbours correlated as the block assumes
        r = (fx.flux[tr.bins, 1] .- B * ĉ) ./ fx.flux_error[tr.bins, 1]
        @test 0.9 < sqrt(sum(abs2, r) / length(r)) < 1.1
        ρ = sum(r[1:end-1] .* r[2:end]) / sum(abs2, r)
        @test abs(ρ - inst.interp_corr) < 0.1
        @test 0.9 < sum(abs2, tr.white .- Z * ĉ) / length(r) < 1.1
        # ignoring the correlation claims more than the data hold
        d0 = RVSEpochSpectra(Float64.(fx.flux), Float64.(fx.flux_error), grid; instrument=RVSInstrument(mismatch_amplitude=1e-4, interp_corr=0.0))
        H0 = ForwardDiff.hessian(y -> sum(first(rvs_transit_loglike(d0, t, c, v(t, (y[1], y[2])), w(t, y[3]))) for t in 1:nt), x)
        @test all(sqrt.(diag(inv(-H0))) .< σx)
    end

    @testset "library mismatch" begin
        wave = collect(849.0:0.025:855.0)
        same = library_mismatch(rvs, rvs, 5777.0, 4.44, 0.0; wavelength=wave)
        @test maximum(abs, same.residual) < 1e-5
        @test abs(same.velocity) < 1e-3
        @test all(==(0.003), same.amplitude)
        # a second library that differs by a feature 5 % deep at 852 nm
        λ = OG.wavelengths(rvs)
        other = RVSLibrary(rvs.name, rvs.teff, rvs.logg, rvs.mh, rvs.ln_wave0, rvs.dln_wave,
            rvs.mean .+ log10.(1 .- 0.05 .* exp.(-0.5 .* ((λ .- 852.0) ./ 0.05) .^ 2)),
            rvs.basis, rvs.weights, rvs.present, rvs.rms_wave, rvs.resolving_power, rvs.path)
        diff = library_mismatch(rvs, other, 5777.0, 4.44, 0.0; wavelength=wave)
        i = argmax(abs.(diff.residual))
        @test abs(wave[i] - 852.0) < 0.03
        @test 0.03 < abs(diff.residual[i]) < 0.055
        @test argmax(diff.amplitude) in i-3:i+3
        @test diff.amplitude[1] < 0.01 < diff.amplitude[i]
        @test abs(diff.velocity) < 0.5
        @test_throws ArgumentError library_mismatch(rvs, other, 9000.0, 4.44, 0.0; wavelength=wave)
        # the profile as the mismatch of a data set: one number per bin, used in every transit
        grid = LogWaveGrid(rvs)
        σ = fill(0.01, length(wave), 2)
        d1 = RVSEpochSpectra(ones(length(wave), 2), σ, grid; wavelength=wave, mismatch=diff.amplitude)
        d2 = RVSEpochSpectra(ones(length(wave), 2), σ, grid; wavelength=wave, mismatch=repeat(diff.amplitude, 1, 2))
        @test d1.mismatch == d2.mismatch
        @test d1.transits[2].chol.logdet == d2.transits[2].chol.logdet
        flat = RVSEpochSpectra(ones(length(wave), 2), σ, grid; wavelength=wave, mismatch=0.003)
        @test d1.transits[1].chol.logdet > flat.transits[1].chol.logdet
        @test_throws DimensionMismatch RVSEpochSpectra(ones(length(wave), 2), σ, grid; wavelength=wave, mismatch=[0.01, 0.02])
    end

    @testset "seamed libraries" begin
        # the same library, 0.02 dex brighter, joined to itself over 5000–5600 K
        hot = SpectralLibrary(lib.name * " (hot)", lib.teff, lib.logg, lib.mh, lib.edges, lib.log_flux .+ 0.02f0, lib.present, lib.path)
        seam = SeamedLibrary(lib, hot; teff=(5000, 5600))
        n = OG.nbins(lib)
        @test OG.nbins(seam) == n
        f(l, T) = surface_sed!(zeros(typeof(T), n), l, T, 4.5, 0.0)
        @test f(seam, 4900.0) == f(lib, 4900.0)
        @test f(seam, 5000.0) == f(lib, 5000.0)
        @test f(seam, 5600.0) == f(hot, 5600.0)
        @test f(seam, 6000.0) == f(hot, 6000.0)
        @test f(seam, 5300.0) ≈ f(lib, 5300.0) .* exp10(0.01) rtol = 1e-5      # half way: the geometric mean
        @test all(f(lib, 5150.0) .< f(seam, 5150.0) .< f(hot, 5150.0))
        for T in (5000.0, 5600.0)                                              # value and slope are continuous at both ends
            @test f(seam, T - 1e-6) ≈ f(seam, T + 1e-6) rtol = 1e-7
            g(x) = log(sum(f(seam, x)))
            @test ForwardDiff.derivative(g, T - 1e-6) ≈ ForwardDiff.derivative(g, T + 1e-6) rtol = 1e-4
        end
        @test all(isnan, f(seam, 5300.0 + 0 * NaN))
        @test_throws ArgumentError SeamedLibrary(lib, hot; teff=(6000, 7000))   # beyond the cool library
        @test_throws ArgumentError SeamedLibrary(lib, hot; teff=(5600, 5000))
        # a star and its photometry take either kind
        star, plain = Star(evol, seam), Star(evol, lib)
        @test SEDFunctionals(seam).G == SEDFunctionals(lib).G
        rate, rate0 = BandRate(star, :G), BandRate(plain, :G)
        @test rate(4800.0, 4.5, 0.0, 0.8, 50.0, 0.0) == rate0(4800.0, 4.5, 0.0, 0.8, 50.0, 0.0)
        @test rate(5900.0, 4.5, 0.0, 0.8, 50.0, 0.0) ≈ rate0(5900.0, 4.5, 0.0, 0.8, 50.0, 0.0) * exp10(0.02) rtol = 1e-5
        # the RVS window
        hotrvs = RVSLibrary(rvs.name * " (hot)", rvs.teff, rvs.logg, rvs.mh, rvs.ln_wave0, rvs.dln_wave, rvs.mean .+ 0.02,
            rvs.basis, rvs.weights, rvs.present, rvs.rms_wave, rvs.resolving_power, rvs.path)
        rseam = SeamedLibrary(rvs, hotrvs; teff=(5000, 5600))
        @test LogWaveGrid(rseam) == LogWaveGrid(rvs)
        r(l, T) = rvs_surface_spectrum!(zeros(length(LogWaveGrid(rvs))), l, T, 4.5, 0.0)
        @test r(rseam, 4800.0) == r(rvs, 4800.0)
        # into an array that lives in a Bumper arena, as the likelihood's scratch does
        @no_escape Bumper.default_buffer() begin
            tmp = @alloc(Float64, length(LogWaveGrid(rvs)))
            rvs_surface_spectrum!(tmp, rseam, 5300.0, 4.5, 0.0)
            @test tmp ≈ r(rvs, 5300.0) .* exp10(0.01)
            surf = @alloc(Float64, n)
            surface_sed!(surf, seam, 5300.0, 4.5, 0.0)
            @test surf == f(seam, 5300.0)
            nothing
        end
        @test r(rseam, 5300.0) ≈ r(rvs, 5300.0) .* exp10(0.01)
        @test r(rseam, 5800.0) == r(hotrvs, 5800.0)
        both = Star(evol, seam, rseam)
        @test rvs_grid(both) == LogWaveGrid(rvs)
        @test_throws ArgumentError SeamedLibrary(rvs, RVSLibrary(rvs.name, rvs.teff, rvs.logg, rvs.mh, rvs.ln_wave0 + 1e-6, rvs.dln_wave,
            rvs.mean, rvs.basis, rvs.weights, rvs.present, rvs.rms_wave, rvs.resolving_power, rvs.path))
    end

    @testset "BHAC15 and the mass seam" begin
        low = BHAC15Grid(joinpath(FIXTURES, "bhac15_fixture.txt"))
        # the table at its nodes: 0.1 M☉ at 5 Gyr is 2811 K, 0.124 R☉, log g 5.25
        s = low(0.1, log10(5e9))
        @test s.teff ≈ 2811 atol = 0.01
        @test s.radius ≈ 0.124 atol = 1e-6
        @test s.logg ≈ 5.25 atol = 0.01
        @test low(0.2, log10(4e9)).teff ≈ 3262 atol = 0.01
        @test 2811 < low(0.105, 9.68).teff < 2916                   # between nodes
        @test isnan(low(0.05, log10(5e9)).teff)                      # the 5 Gyr isochrone starts at 0.07 M☉
        @test isnan(low(0.1, 9.3).teff) && isnan(low(0.1, 10.0).teff)
        @test isnan(low(2.0, 9.7).radius) && isnan(low(-1.0, 9.7).radius)
        @test all(isfinite, ForwardDiff.gradient(x -> low(x[1], x[2]).radius, [0.09, 9.65]))

        seam = SeamedEvolution(low, evol)
        age = 9.66
        # the stellar model, untouched, from the top of the overlap up
        @test seam(0.5, age, 0.1) == evol(0.5, age, 0.1)
        @test seam(0.15, age, 0.1) == evol(0.15, age, 0.1)
        # MIST alone stops at 0.1 M☉; the seamed model goes on
        @test isnan(evol(0.09, age, 0.0).teff)
        x = seam(0.09, age, 0.0)
        @test 2500 < x.teff < 2800 && 0.10 < x.radius < 0.12 && 5.2 < x.logg < 5.4
        @test x.mass_now ≈ 0.09 rtol = 0.02
        @test x.logL ≈ 2log10(x.radius) + 4log10(x.teff / OG.TEFF_SUN)
        @test isnan(x.eep) && !isnan(seam(0.12, age, 0.0).eep)
        @test isnan(seam(0.04, age, 0.0).teff)                       # below BHAC15 at this age
        @test isnan(seam(0.09, 9.9, 0.0).teff)                       # outside the MIST fixture's ages
        # continuous, with a continuous slope, at both ends of the overlap and at every metallicity
        for feh in (-0.2, 0.0, 0.2), m in (0.10, 0.15), q in (:teff, :radius, :logg)
            f(mass) = log(getproperty(seam(mass, age, feh), q))
            @test f(m * (1 - 1e-7)) ≈ f(m * (1 + 1e-7)) atol = 1e-6
            @test ForwardDiff.derivative(f, m * (1 - 1e-6)) ≈ ForwardDiff.derivative(f, m * (1 + 1e-6)) rtol = 1e-3 atol = 1e-3
        end
        # in the overlap the two models are close at solar metallicity, so the blend stays between them
        for m in (0.11, 0.13)
            a, b, c = low(m, age).teff, evol(m, age, 0.0).teff, seam(m, age, 0.0).teff
            @test abs(a - b) < 60 && min(a, b) - 30 < c < max(a, b) + 30
        end
        # monotonic in mass through the seam
        ms = 0.075:0.0025:0.2
        @test issorted([seam(m, age, 0.0).teff for m in ms])
        @test issorted([seam(m, age, 0.0).radius for m in ms])
        # below the join the metallicity dependence is that of the star at its top
        ratio(m) = seam(m, age, -0.2).teff / seam(m, age, 0.0).teff
        @test ratio(0.09) ≈ ratio(0.15) rtol = 1e-10
        @test ratio(0.15) > 1.01
        @test seam(0.09, age, -0.2).mh == evol(0.15, age, -0.2).mh
        g = ForwardDiff.gradient(v -> seam(v[1], v[2], v[3]).teff, [0.12, age, 0.05])
        @test all(isfinite, g) && g[1] > 0 && g[3] < 0
        @test_throws ArgumentError SeamedEvolution(low, evol; mass=(0.2, 0.1))
        # as the evolution model of a star
        star = Star(seam, lib)
        rate = BandRate(star, :G)
        lo, hi = seam(0.125, age, 0.0), seam(0.2, age, 0.0)
        @test 0 < rate(lo.teff, lo.logg, lo.mh, lo.radius, 50.0, 0.0) < rate(hi.teff, hi.logg, hi.mh, hi.radius, 50.0, 0.0)
    end

    # GaiaTwin's wavefront optics are the stand-in instrument. The fixture holds
    # the chromatic basis of one device built from its monochromatic images,
    # and its own polychromatic images of two stars in an 18 × 12 window.
    @testset "image basis against GaiaTwin" begin
        path = joinpath(FIXTURES, "gaiatwin_image_fixture.h5")
        basis = ImageBasis(path, "fov1_ccd5")
        @test OG.nterms(basis) == 7
        @test_throws ArgumentError ImageBasis(path, "fov3_ccd1")
        # the instrument and the spectra expand over the same range of wavenumber
        f = SEDFunctionals(lib)
        @test collect(basis.wavenumber_range) ≈ collect(f.wavenumber_range)
        @test collect(wavenumber_range(f.passbands.G)) == collect(f.wavenumber_range)
        ref = OG.h5open(g -> (; (Symbol(k) => read(g["reference"][k]) for k in keys(g["reference"]))...), path)
        ν_lo, ν_hi = basis.wavenumber_range
        x = 2 .* (1 ./ ref.lam .- ν_lo) ./ (ν_hi - ν_lo) .- 1
        moments(k, degree=6) = [n <= degree ? sum(ref.weights[:, k] .* OG._legendre.(n, x)) : 0.0 for n in 1:6]
        image(k, c) = [image_fraction(basis, c, ref.w[i, j] - ref.position[1, k], ref.z[i, j] - ref.position[2, k]) for i in 1:18, j in 1:12]
        centroid(a) = sum(a .* ref.w) / sum(a)
        for k in 1:2
            truth = ref.image[:, :, k]
            img = image(k, moments(k))
            # (2e-4 and 7e-4 of the peak with the product's 4 × 12 mas table; the fixture's is 12 × 30)
            @test maximum(abs, img .- truth) < 2e-3 * maximum(truth)
            @test sum(img) ≈ sum(truth) rtol = 1e-3
            @test abs(centroid(img) - centroid(truth)) < 0.005                 # mas
            # one colour is not enough, and none is far off
            @test abs(centroid(image(k, moments(k, 1))) - centroid(truth)) > 0.015
            @test abs(centroid(image(k, moments(k, 0))) - centroid(truth)) > 0.1
        end
        # the package's own spectra give nearly the moments GaiaTwin's do (5 nm against 30 nm bins)
        sed = surface_sed!(zeros(OG.nbins(lib)), lib, 5800.0, 4.5, 0.0)
        c = zeros(6)
        chromatic_moments!(c, f, sed)
        @test isapprox(c, moments(1); atol=0.02)
        # the terms, and the fraction as their sum
        t = image_terms!(zeros(7), basis, 37.0, -112.0)
        @test image_fraction(basis, c, 37.0, -112.0) ≈ t[1] + dot(c, t[2:end])
        @test t[1] > 0.02 && all(abs.(t[2:end]) .< t[1])
        # zero beyond the table, in either direction
        @test image_fraction(basis, c, 700.0, 0.0) == 0 && image_fraction(basis, c, 0.0, -1400.0) == 0
        @test all(iszero, image_terms!(ones(7), basis, -700.0, 0.0))
        @test_throws DimensionMismatch image_fraction(basis, zeros(3), 0.0, 0.0)
        # continuous with continuous derivatives across the table's nodes (x = 36 mas, y = 60 mas are nodes)
        g(v) = image_fraction(basis, c, v[1], v[2])
        for v in ([36.0, 17.0], [5.0, 60.0])
            e = [1e-7, 1e-7]
            @test g(v .- e) ≈ g(v .+ e) rtol = 1e-6
            @test ForwardDiff.gradient(g, v .- e) ≈ ForwardDiff.gradient(g, v .+ e) rtol = 1e-3
        end
        v = [21.3, -77.7]
        gfd = [(g(v .+ h) - g(v .- h)) / 2e-4 for h in ([1e-4, 0.0], [0.0, 1e-4])]
        @test ForwardDiff.gradient(g, v) ≈ gfd rtol = 1e-5
        @test all(isfinite, ForwardDiff.gradient(m -> image_fraction(basis, m, 21.3, -77.7), c))
        # a window: background plus each source's electrons times its own image
        srcs = ntuple(k -> (; w=ref.position[1, k], z=ref.position[2, k], electrons=(4.0e5, 3.0e4)[k], moments=moments(k)), 2)
        μ = window_mean!(zeros(18, 12), basis, ref.w, ref.z, 12.5, srcs)
        @test μ ≈ 12.5 .+ 4.0e5 .* image(1, moments(1)) .+ 3.0e4 .* image(2, moments(2))
        @test maximum(abs, μ .- (12.5 .+ 4.0e5 .* ref.image[:, :, 1] .+ 3.0e4 .* ref.image[:, :, 2])) < 2e-3 * 4.0e5 * maximum(ref.image[:, :, 1])
        @test window_mean!(zeros(18, 12), basis, ref.w, ref.z, 0.0, (srcs[1],)) ≈ 4.0e5 .* image(1, moments(1))
        # positions and fluxes as dual numbers
        d = ForwardDiff.gradient(p -> sum(abs2, window_mean!(zeros(eltype(p), 18, 12), basis, ref.w, ref.z, p[1],
                ((; w=p[2], z=p[3], electrons=p[4], moments=moments(1)),))), [12.5, 0.0, 0.0, 4.0e5])
        @test all(isfinite, d) && all(!iszero, d)
        @test_throws DimensionMismatch window_mean!(zeros(18, 11), basis, ref.w, ref.z, 0.0, srcs)
    end

    # The RVS block on the model's absolute frame: positions as `GaiaDR4AstromObs`
    # places them, velocities from the frame's own radial velocity.
    @testset "GaiaDR4RawObs" begin
        RA0, DEC0, REF = 209.506326888, 31.695499700, 57936.875
        star = Star(evol, lib, rvs)
        A = Body(name="A", variables=@variables begin
            mass ~ truncated(Normal(0.9, 0.03), lower=0.8, upper=1.0)
            star = $evol(mass, system.logage, system.feh)
            teff = star.teff
            logg = star.logg
            mh = star.mh
            radius = star.radius
            vsini ~ Uniform(0, 15)
            drv ~ Normal(0, 300)
        end)
        b = Body(name="b", about=A, variables=@variables begin
            mass ~ Uniform(0.6, 0.7)
            star = $evol(mass, system.logage, system.feh)
            teff = star.teff
            logg = star.logg
            mh = star.mh
            radius = star.radius
            vsini ~ Uniform(0, 15)
            drv ~ Normal(0, 300)
            a ~ LogUniform(0.5, 2.0)
            e ~ Uniform(0.0, 0.4)
            i ~ Sine()
            ω ~ Uniform(0, 2pi)
            Ω ~ Uniform(0, 2pi)
            θ ~ Uniform(0, 2pi)
            epoch = 57000.0
        end)
        rng = Random.Xoshiro(5)
        nt = 10
        epoch = sort(56900 .+ 2000 .* rand(rng, nt))
        ψ = 360 .* rand(rng, nt)
        f_al, f_ac = 1.4 .* rand(rng, nt) .- 0.7, 1.4 .* rand(rng, nt) .- 0.7
        track = catalogue_track(epoch, ψ, f_al, f_ac; parallax=40.0, pmra=120.0, pmdec=-80.0, ref_epoch=REF)
        wave = collect(849.0:0.025:855.0)
        table(; ref_al=track.al, flux=[ones(length(wave)) for _ in 1:nt]) = Table(; epoch, scan_pos_angle=ψ,
            parallax_factor_al=f_al, parallax_factor_ac=f_ac, ref_pos_al=ref_al, ref_pos_ac=track.ac,
            flux, flux_error=[fill(0.01, length(wave)) for _ in 1:nt])
        raw(; kw...) = GaiaDR4RawObs(; rvs=table(; kw...), emitters=(A=star, b=star), tangent_point=(RA0, DEC0), wavelength=wave)
        frame = @variables begin
            plx ~ truncated(Normal(40.0, 0.05), lower=1.0)
            Δra_mas ~ Normal(0, 5)
            Δdec_mas ~ Normal(0, 5)
            ra = $RA0 + Δra_mas / (3.6e6 * cosd($DEC0))
            dec = $DEC0 + Δdec_mas / 3.6e6
            pmra ~ Normal(120.0, 1.0)
            pmdec ~ Normal(-80.0, 1.0)
            rv ~ Normal(15_000, 3_000)
            ref_epoch = $REF
            logage ~ Uniform(9.6, 9.75)
            feh ~ truncated(Normal(0.0, 0.05), lower=-0.2, upper=0.2)
            av = 0.0
        end
        astrom = GaiaDR4AstromObs(Table(; epoch, scan_pos_angle=ψ, parallax_factor_al=f_al, centroid_pos_al=zeros(nt),
                centroid_pos_error_al=fill(0.1, nt)); target=A, ref=Barycentre, tangent_point=(RA0, DEC0), name="astrom")
        system(obs...) = System(name="rawpair", bodies=[A, b], observations=collect(obs), variables=frame)
        obs = raw()
        sys = system(obs, astrom)
        @test Octofitter.refspecs(obs) == (Octofitter.refspec(:A), Octofitter.refspec(:b), Barycentre)
        Random.seed!(3)
        θ = Octofitter.drawfrompriors(sys)
        function contexts(sys, θ)
            unique_ep, maps = Octofitter.epoch_plan(sys)
            posys = Octofitter._build_po_system(sys, θ)
            traj = Octofitter.orbitsolve(posys, unique_ep; Octofitter._solvekw(sys)...)
            return map(o -> Octofitter._obsctx(sys, θ, (;), posys, traj, maps[o]), sys.observations)
        end
        ctx, ctx_astrom = contexts(sys, θ)
        sim = Octofitter.simulate(obs, ctx)
        @test size(sim.flux) == (length(wave), nt) && size(sim.velocity) == (2, nt)

        # positions: the same local plane as the DR4 astrometry of the same body
        refA, refb, bary = Octofitter.ref(ctx, Octofitter.refspec(:A)), Octofitter.ref(ctx, Octofitter.refspec(:b)), Octofitter.ref(ctx, Barycentre)
        al = Octofitter.simulate(astrom, ctx_astrom).along_scan
        pos(body, t) = transit_position(Octofitter.solutionat(ctx, t), body, bary, obs.tangent_basis, obs.sinψ[t], obs.cosψ[t], f_al[t], f_ac[t])
        @test [pos(refA, t)[1] for t in 1:nt] ≈ al rtol = 1e-12
        # the pipeline's track is the frame's to within what the prior allows the two to differ by
        @test maximum(abs, al .- track.al) < 60
        # velocities: frame + orbit + redshift + template zero point + κ × along-scan offset
        inst = obs.rvs.instrument
        for t in (1, nt), (k, body, θb) in ((1, refA, θ.bodies.A), (2, refb, θ.bodies.b))
            sol = Octofitter.solutionat(ctx, t)
            w, z = pos(body, t)
            v = (Octofitter.frame_rv(sol) + Octofitter.radvel(sol, body, bary)) / 1000 +
                gravitational_redshift(θb.logg, θb.radius) + θb.drv / 1000 + inst.kappa * (w - track.al[t])
            @test sim.velocity[k, t] ≈ v rtol = 1e-12
            @test sim.weight[k, t] ≈ window_fraction(inst, z - track.ac[t])
        end
        @test abs(sim.velocity[1, 1] - θ.rv / 1000) < 40                # the systemic velocity is the frame's
        @test all(x -> 0.99 < x <= 1, sim.weight)                        # a close pair is whole in the window
        # moving the assumed position along scan is a velocity shift for every emitter
        shifted = Octofitter.simulate(raw(ref_al=track.al .+ 100), contexts(system(raw(ref_al=track.al .+ 100), astrom), θ)[1])
        @test shifted.velocity ≈ sim.velocity .- inst.kappa * 100
        @test abs(inst.kappa * 100) > 14
        # an emitter the pipeline subtracted is not in the spectrum
        mem = ones(2, nt)
        mem[2, 3] = 0
        deb = GaiaDR4RawObs(; rvs=table(), emitters=(A=star, b=star), tangent_point=(RA0, DEC0), wavelength=wave, membership=mem)
        simd = Octofitter.simulate(deb, contexts(system(deb, astrom), θ)[1])
        @test simd.weight[2, 3] == 0 && simd.weight[2, 4] == sim.weight[2, 4]
        @test simd.flux[:, 4] == sim.flux[:, 4] && simd.flux[:, 3] != sim.flux[:, 3]
        single = GaiaDR4RawObs(; rvs=table(), emitters=(A=star,), tangent_point=(RA0, DEC0), wavelength=wave)
        @test Octofitter.simulate(single, contexts(system(single, astrom), θ)[1]).flux[:, 3] ≈ simd.flux[:, 3]

        # data simulated by the model are most probable at the parameters that made them
        Random.seed!(11)
        truth = Octofitter.generate_from_params(sys, θ; add_noise=true)
        tobs = truth.observations[1]
        @test tobs isa GaiaDR4RawObs && tobs.table.flux[1] != obs.table.flux[1]
        @test tobs.rvs.p_bad == obs.rvs.p_bad && tobs.membership == obs.membership
        ll = Octofitter.make_ln_like(truth)
        ℓ0 = ll(truth, θ)
        @test isfinite(ℓ0)
        # (over 40 such data sets the pull of `rv` has mean -0.14 and standard deviation 1.19, with σ = 0.24 km/s)
        @test ll(truth, merge(θ, (; rv=θ.rv + 3000))) < ℓ0 - 20
        @test ll(truth, merge(θ, (; rv=θ.rv - 3000))) < ℓ0 - 20
        @test ll(truth, merge(θ, (; bodies=merge(θ.bodies, (; A=merge(θ.bodies.A, (; drv=θ.bodies.A.drv + 2000))))))) < ℓ0 - 10
        # transits are independent given the parameters
        tctx = contexts(truth, θ)[1]
        part(inds) = Octofitter.ln_like(Octofitter.likeobj_from_epoch_subset(tobs, inds),
            contexts(system(Octofitter.likeobj_from_epoch_subset(tobs, inds), astrom), θ)[1])
        @test part(1:4) + part(5:nt) ≈ Octofitter.ln_like(tobs, tctx)
        # the log density and its gradient through the whole model
        model = Octofitter.LogDensityModel(truth; verbosity=0)
        θv = Octofitter.sample_priors(Random.Xoshiro(7), truth)
        θt = model.link(θv)
        v, g = model.∇ℓπcallback(θt)
        @test isfinite(v) && all(isfinite, g)
        gfd = [central(x -> [model.ℓπcallback(x)], θt, i; h=1e-6)[1] for i in eachindex(θt)]
        @test isapprox(g, gfd; rtol=2e-3, atol=1e-3 * maximum(abs, g))

        # what the model refuses
        @test_throws ErrorException System(name="noframe", bodies=[A, b], observations=[raw()], variables=@variables begin
            plx ~ Uniform(30, 50)
            logage ~ Uniform(9.6, 9.75)
            feh = 0.0
            av = 0.0
        end)
        @test_throws ErrorException GaiaDR4RawObs(; rvs=table(), emitters=(A=star, b=star), wavelength=wave)        # no tangent point
        @test_throws ErrorException GaiaDR4RawObs(; rvs=table()[nt:-1:1], emitters=(A=star,), tangent_point=(RA0, DEC0), wavelength=wave)
        @test_throws ArgumentError GaiaDR4RawObs(; rvs=table(), emitters=(A=star, b=star), tangent_point=(RA0, DEC0), wavelength=wave, membership=zeros(2, nt))
        @test_throws MethodError GaiaDR4RawObs(; rvs=table(), emitters=(A=Star(evol, lib),), tangent_point=(RA0, DEC0), wavelength=wave)
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
