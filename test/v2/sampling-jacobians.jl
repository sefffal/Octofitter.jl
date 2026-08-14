# Sampling-coordinate Jacobians.
#
# `logjac_cartesian_to_campbell` is the potential term that converts a flat
# prior in Cartesian relative phase space into a flat prior in Campbell
# elements. Two things need gating and they are different in kind:
#
#   1. the closed form is the actual Jacobian of the actual map — measured
#      against PlanetOrbits' own elements→state path, not re-derived on paper;
#   2. the distributional consequence the docstring claims (thermal e,
#      isotropic cos i, and the Mtot^{3/2} tilt) really follows — and that
#      reweighting by the term removes it.
#
# (1) is the sharp test; (2) is the one that would catch a correct formula
# documented with the wrong sign or the wrong direction.

using Octofitter, PlanetOrbits, Distributions, LinearAlgebra, Random, Test
using Octofitter: logjac_cartesian_to_campbell
import ForwardDiff

const JAC_EPOCH = 57388.0

# Raw dynamical state of the secondary relative to the primary: no
# `observe_pass!`, so these are the propagator's states at the common emission
# epoch, which is what the canonical (Delaunay) variables refer to. Going
# through `orbitsolve` instead would retard each body to its own emission time
# and measure the Jacobian of a *different*, v/c-perturbed map.
function _jac_state(p, Mtot)
    a, e, cosi, M0, ω, Ω = p
    i = acos(clamp(cosi, -one(cosi), one(cosi)))
    mA = Mtot * 0.85
    A = PlanetOrbits.Body(mass=mA, name=:A)
    B = PlanetOrbits.Body(mass=Mtot - mA, name=:B)
    o = PlanetOrbits.Orbit(B, about=A; a, e, i, ω, Ω, M0, epoch=JAC_EPOCH)
    sys = PlanetOrbits.System((A, B), (o,))
    traj = Trajectory(sys, [JAC_EPOCH])
    PlanetOrbits.frame_pass!(traj, sys.frame)
    PlanetOrbits.propagate!(traj, sys, KeplerianApprox())
    sol = traj[1]
    return [posx(sol, :B, :A), posy(sol, :B, :A), posz(sol, :B, :A),
            velx(sol, :B, :A), vely(sol, :B, :A), velz(sol, :B, :A)]
end

@testset "logjac_cartesian_to_campbell" begin

    @testset "the closed form is the map's actual Jacobian" begin
        # |∂(x, v)/∂(a, e, cos i, M, ω, Ω)| = ½ · μ^{3/2} · √a · e
        # — the product of the diagonal of a triangular (L, G, H) transform,
        # as derived in the docstring.
        closed(a, e, Mtot) =
            0.5 * (PlanetOrbits.GM_sun_au3_julianyr2 * Mtot)^1.5 * sqrt(a) * e
        worst = 0.0
        for a in (0.8, 5.0, 750.0), e in (0.05, 0.31, 0.77),
            cosi in (-0.6, 0.2, 0.9), Mtot in (0.5, 1.5156, 3.0),
            (M0, ω, Ω) in ((0.7, 1.3, 2.9), (4.1, 5.5, 0.4))

            J = ForwardDiff.jacobian(q -> _jac_state(q, Mtot), [a, e, cosi, M0, ω, Ω])
            worst = max(worst, abs(abs(det(J)) / closed(a, e, Mtot) - 1))
        end
        @test worst < 1e-9

        # ...and the potential term is minus its log, up to a constant. Test
        # *differences*, which is all an additive log-density term can mean —
        # and which is exactly what would catch a sign error.
        for (p, q) in (((5.0, 0.31, 1.5), (12.0, 0.77, 2.0)),
                       ((0.8, 0.05, 0.5), (750.0, 0.6, 3.0)))
            Δjac = log(closed(q...)) - log(closed(p...))
            Δterm = logjac_cartesian_to_campbell(q...) -
                    logjac_cartesian_to_campbell(p...)
            @test Δterm ≈ -Δjac rtol = 1e-12
        end
    end

    # Weighted-CDF maximum deviation — a KS statistic for importance weights.
    function _wks(xs, ws, cdf)
        p = sortperm(xs)
        acc = 0.0
        dev = 0.0
        for k in p
            lo = acc
            acc += ws[k]
            dev = max(dev, abs(cdf(xs[k]) - lo), abs(cdf(xs[k]) - acc))
        end
        return dev
    end

    # The distributional claim, tested over a **product region in element
    # space** and weighted by the **numerically measured** Jacobian.
    #
    # Both choices are deliberate. Sampling a box in (x…vz) instead — the
    # obvious thing — does not work, and the reason is worth recording because
    # it looks like a tolerance problem and is not: a Cartesian box is not a
    # product region in element space. A velocity bound truncates high-speed
    # phases, which cuts small-r states harder at large `a` (so the marginal in
    # `a` is √a times the region's own a-dependent volume) and cuts high-e
    # orbits near periapsis (so `e` is not thermal either). Measured: the √a
    # marginal was off by 20%, and `e` still failed at the 0.1% level after
    # conditioning on a thin a-shell. The truncation is real, it is a property
    # of the region rather than of the density, and the docstring says so.
    #
    # Using ForwardDiff for the weights rather than the closed form keeps this
    # test independent of the formula it is checking: it asserts that the
    # *actual* map's Jacobian produces the *claimed* marginals.
    @testset "flat-Cartesian implies thermal e, isotropic cos i, p(a) ∝ √a" begin
        rng = Random.Xoshiro(20260812)
        Mtot = 1.5156
        alo, ahi = 2.0, 12.0
        n = 6000
        as = alo .+ (ahi - alo) .* rand(rng, n)
        es = 0.9 .* rand(rng, n)
        cosis = 2 .* rand(rng, n) .- 1
        M0s = 2π .* rand(rng, n)
        ωs = 2π .* rand(rng, n)
        Ωs = 2π .* rand(rng, n)

        # Uniform draws over the element box, importance-weighted by |det J| —
        # which turns them into draws from a flat density in Cartesian phase
        # space, restricted to this box's image.
        w = map(1:n) do k
            J = ForwardDiff.jacobian(q -> _jac_state(q, Mtot),
                [as[k], es[k], cosis[k], M0s[k], ωs[k], Ωs[k]])
            abs(det(J))
        end
        w ./= sum(w)
        ess = 1 / sum(abs2, w)
        @test ess > 2000
        tol = 2.5 / sqrt(ess)

        # p(e) ∝ e on [0, b]   ⇒ F(e) = e²/b²          (thermal eccentricity)
        # p(a) ∝ √a on [l, u]  ⇒ F(a) = (a^{3/2}−l^{3/2})/(u^{3/2}−l^{3/2})
        # p(cos i) flat        ⇒ F(x) = (x + 1)/2      (isotropic orientation)
        @test _wks(es, w, x -> (x / 0.9)^2) < tol
        @test _wks(as, w, x -> (x^1.5 - alo^1.5) / (ahi^1.5 - alo^1.5)) < tol
        @test _wks(cosis, w, x -> (x + 1) / 2) < tol
        # ω and Ω pick up nothing: they are canonical angles, so the Jacobian
        # is independent of them.
        @test _wks(ωs, w, x -> x / 2π) < tol
        @test _wks(Ωs, w, x -> x / 2π) < tol

        # ...and the weighting is discriminating: e and a are emphatically not
        # flat under it, which is the whole point of the correction below.
        # Measured separations are ~6.7σ (e) and ~2.3σ (a); `a` is the weaker
        # of the two because √a varies by only 2.4× across this box where e
        # varies without bound at its lower end.
        @test _wks(es, w, x -> x / 0.9) > 4tol
        @test _wks(as, w, x -> (x - alo) / (ahi - alo)) > 2tol

        # Adding the potential term cancels the tilt exactly: the combined
        # weight must be flat in all three, recovering the uniform element box
        # the draws came from.
        w2 = w .* exp.([logjac_cartesian_to_campbell(as[k], es[k], Mtot) for k in 1:n])
        w2 ./= sum(w2)
        ess2 = 1 / sum(abs2, w2)
        @test ess2 > 0.99n          # a perfect cancellation leaves equal weights
        tol2 = 2.5 / sqrt(ess2)
        @test _wks(es, w2, x -> x / 0.9) < tol2
        @test _wks(as, w2, x -> (x - alo) / (ahi - alo)) < tol2
        @test _wks(cosis, w2, x -> (x + 1) / 2) < tol2
    end

    @testset "the mass tilt is the documented one" begin
        # The Mtot^{3/2} factor is the surprise the docstring calls out in
        # bold, so it gets its own assertion rather than riding along inside
        # the determinant test.
        for (m1, m2) in ((0.5, 2.0), (1.0, 1.5156))
            Δ = logjac_cartesian_to_campbell(5.0, 0.3, m2) -
                logjac_cartesian_to_campbell(5.0, 0.3, m1)
            @test Δ ≈ -1.5 * (log(m2) - log(m1)) rtol = 1e-12
        end
        # a and e enter with the exponents the derivation gives.
        @test logjac_cartesian_to_campbell(20.0, 0.3, 1.0) -
              logjac_cartesian_to_campbell(5.0, 0.3, 1.0) ≈ -0.5 * log(4) rtol = 1e-12
        @test logjac_cartesian_to_campbell(5.0, 0.6, 1.0) -
              logjac_cartesian_to_campbell(5.0, 0.3, 1.0) ≈ -log(2) rtol = 1e-12
    end

    @testset "usable as an `LL +=` term in a model" begin
        # The whole point is that it is a potential term in a `@variables`
        # block, so gate that it actually compiles and moves the log-density.
        # The primary's mass reaches the companion's block through a *system*
        # variable: a body block sees `system.*`, never a sibling, so
        # `system.bodies.A.mass` is not available there.
        primary() = Octofitter.Body(name="A", variables=@variables begin
            mass = system.M_A
        end)
        wrap(A, B) = Octofitter.System(name="wide", bodies=[A, B], observations=[],
            variables=@variables begin
                plx = 74.19
                M_A = 1.29
            end)

        A1 = primary()
        plain = wrap(A1, Octofitter.Body(name="B", about=A1, variables=@variables begin
            mass = 0.2256
            x ~ Uniform(-40, 40)
            y ~ Uniform(-40, 40)
            z ~ Uniform(-40, 40)
            vx = 0.35
            vy = 1.1
            vz = -0.2
            epoch = $JAC_EPOCH
        end))

        A2 = primary()
        tilted = wrap(A2, Octofitter.Body(name="B", about=A2, variables=@variables begin
            mass = 0.2256
            x ~ Uniform(-40, 40)
            y ~ Uniform(-40, 40)
            z ~ Uniform(-40, 40)
            vx = 0.35
            vy = 1.1
            vz = -0.2
            epoch = $JAC_EPOCH
            LL += logjac_cartesian_to_campbell(x, y, z, vx, vy, vz,
                                               mass + system.M_A)
        end))

        θ = [0.21, -0.34, 0.12]
        ld1 = Octofitter.LogDensityModel(plain; verbosity=0)
        ld2 = Octofitter.LogDensityModel(tilted; verbosity=0)
        @test ld1.D == 3 && ld2.D == 3
        lp1, lp2 = ld1.ℓπcallback(θ), ld2.ℓπcallback(θ)
        @test isfinite(lp1) && isfinite(lp2)

        # The difference is exactly the closed form at this draw — computed
        # here from the state the model actually built, so this checks the
        # plumbing (units, mass argument, sign) and not just that "it moved".
        st = ld1.invlink(θ)
        nt = Octofitter.make_arr2nt(plain)(st)
        b = nt.bodies.B
        expect = logjac_cartesian_to_campbell(b.x, b.y, b.z, b.vx, b.vy, b.vz,
                                              b.mass + nt.bodies.A.mass)
        @test lp2 - lp1 ≈ expect rtol = 1e-10
    end
end
