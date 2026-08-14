# The interim system exposed to deferred lines, and `AnchoredFrame` on top of
# it.
#
# Two things are being pinned here, and they are different in kind. The
# codegen half (`system_interim`) is a scoping and cost question: who can see
# it, what it is built from, and that a model which never mentions it pays
# nothing. The `AnchoredFrame` half is a *physics* question: that the frame
# the model ends up with is the one whose anchor body reproduces the sampled
# catalogue solution.

using Octofitter
using Octofitter: INTERIM_SYSTEM_VAR, INTERIM_PARALLAX_VAR
# Two testsets build `Derived` blocks by hand rather than with `@variables`,
# because `$` in a `@variables` block is *value* interpolation: `$INTERIM_SYSTEM_VAR`
# would splice the Symbol in as a captured constant and the expression would
# never mention `system_interim` at all.
using Octofitter.OrderedCollections: OrderedDict
using PlanetOrbits
using StaticArrays, Distributions, LinearAlgebra, Random, Test
using ForwardDiff

const REF_EP = 57388.0     # Gaia DR3's reference epoch, MJD

# A host plus one companion at `a` AU. Deliberately not a fixture: everything
# here is about the relationship between two systems built from the same
# bodies, so the bodies only have to be nontrivial.
function anchored_bodies(; a=5.0, mass_b=0.02, mass_A=1.29)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = $mass_A
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = $mass_b
        a = $a; e = 0.35; i = 1.1; ω = 0.7; Ω = 2.3; tp = 51000.0
    end)
    return (A, b)
end

const CATALOG = (ra=24.19928, dec=41.40546, plx=74.19,
                 pmra=-172.57, pmdec=-381.32, rv=-26900.0)

function anchored_system(; a=5.0, mass_b=0.02, kwargs...)
    A, b = anchored_bodies(; a, mass_b)
    sys = Octofitter.System(name="anchored", bodies=[A, b], verbosity=0,
        variables=AnchoredFrame(A; ref_epoch=REF_EP, kwargs..., variables=@variables begin
            ra_A = $(CATALOG.ra)
            dec_A = $(CATALOG.dec)
            plx_A = $(CATALOG.plx)
            pmra_A = $(CATALOG.pmra)
            pmdec_A = $(CATALOG.pmdec)
            rv_A = $(CATALOG.rv)
        end))
    return (; sys, A, b)
end

# θ and the built PlanetOrbits.System for a model with no free parameters.
function anchored_solved(sys)
    θ = Octofitter.make_arr2nt(sys)(Float64[])
    return (; θ, posys=Octofitter.make_ln_like(sys, θ).build(θ))
end

@testset "the interim costs nothing when nothing reads it" begin
    A, b = anchored_bodies()
    plain = Octofitter.System(name="plain", bodies=[A, b], verbosity=0,
        variables=@variables begin plx = 74.19 end)
    @test !Octofitter._uses_interim(plain)
    src = string(Base.remove_linenums!(Octofitter.arr2nt_expr(plain)))
    @test !occursin("_interim", src)

    # …and appears exactly once when something does. Once, not once per line
    # that reads it: six anchored quantities share one build and one solve.
    m = anchored_system()
    @test Octofitter._uses_interim(m.sys)
    src2 = string(Base.remove_linenums!(Octofitter.arr2nt_expr(m.sys)))
    @test count("_interim = let T", src2) == 1
    # Once, not once per line that reads it: six anchored quantities share one
    # build and one call.
    @test count("anchored_frame(", src2) == 1
end

@testset "the interim is built from the same declarations as the final system" begin
    # This is the reason the interim lives in codegen rather than being
    # documented as a user-space pattern: one generator emits both, so they
    # cannot drift. Asserted by comparing the interim's own masses and
    # elements against the final system's, which are built from θ.
    m = anchored_system(a=7.5, mass_b=0.031)
    sysvars = (Octofitter.Priors(), Octofitter.Derived(OrderedDict{Symbol,Any}(
        :plx_interim => 74.19,
        :snapshot => :($INTERIM_SYSTEM_VAR),
        :plx => 74.19), (), ()))
    A, b = anchored_bodies(a=7.5, mass_b=0.031)
    snap = Octofitter.System(name="snap", bodies=[A, b], variables=sysvars, verbosity=0)
    @test :snapshot in snap.deferred          # inferred purely from the mention
    interim = Octofitter.make_arr2nt(snap)(Float64[]).snapshot
    @test interim isa PlanetOrbits.System

    final = anchored_solved(m.sys).posys
    @test PlanetOrbits._names(interim) == PlanetOrbits._names(final) == (:A, :b)
    @test collect(interim.masses) == collect(final.masses)
    for f in (:a, :e, :i, :ω, :Ω, :tp, :M)
        @test getfield(interim.rows[1], f) == getfield(final.rows[1], f)
    end
    @test interim.Ainv == final.Ainv
    # …and the one thing that *is* different is the frame, which is the point.
    @test interim.frame isa PlanetOrbits.Parallax
    @test final.frame isa PlanetOrbits.AbsoluteFrame
    @test interim.frame.plx == 74.19
end

@testset "the anchor body reproduces the sampled catalogue solution" begin
    # The acceptance criterion, and the reason `anchored_frame` runs two
    # passes. Solving the *final* system and adding the anchor's reflex back
    # must return the sampled catalogue values.
    for (a, mass_b) in ((5.0, 0.02), (750.0, 0.2256))
        m = anchored_system(; a, mass_b)
        (; θ, posys) = anchored_solved(m.sys)
        sol = orbitsolve(posys, REF_EP)
        bary = PlanetOrbits.barycentre(posys)
        @test PlanetOrbits.frame_pmra(sol) + PlanetOrbits.pmra(sol, :A, bary) ≈
              CATALOG.pmra atol = 1e-6
        @test PlanetOrbits.frame_pmdec(sol) + PlanetOrbits.pmdec(sol, :A, bary) ≈
              CATALOG.pmdec atol = 1e-6
        @test PlanetOrbits.frame_rv(sol) + PlanetOrbits.radvel(sol, :A, bary) ≈
              CATALOG.rv atol = 1e-3         # m/s, against an rv of 2.7e4 of them
        @test PlanetOrbits.frame_ra(sol) ≈ θ.ra atol = 1e-12

        # The corrections themselves are real at both separations — the point
        # of the exercise is that they are large next to a catalogue σ
        # (σ(μ_DR3) ~ 0.04 mas/yr), so the closure above is not trivially
        # "nothing happened".
        @test abs(θ.pmra - CATALOG.pmra) > 1e-3
        @test θ.plx != θ.plx_A
    end
end

@testset "the second pass is what closes it, and by the documented mechanism" begin
    # Both terms the `Parallax`-level interim cannot see, isolated and
    # compared against their closed forms. Written this way rather than as a
    # tolerance because a tolerance would pass just as happily if the second
    # pass silently stopped running.
    m = anchored_system(a=750.0, mass_b=0.2256)
    (; θ, posys) = anchored_solved(m.sys)
    interim = PlanetOrbits.reframe(posys; plx=θ.plx_interim)
    one_pass = Octofitter._anchor_map(interim, :A, REF_EP, ntuple(_ -> true, 6),
        (; ra=CATALOG.ra, dec=CATALOG.dec, plx=CATALOG.plx,
           pmra=CATALOG.pmra, pmdec=CATALOG.pmdec, rv=CATALOG.rv))

    # § Proper motion: the interim's distance is constant, so it misses the
    # angular offset shrinking at ρ⋅(ṙ/d) as the system recedes.
    sol = orbitsolve(posys, REF_EP)
    bary = PlanetOrbits.barycentre(posys)
    d_au = 1000 * PlanetOrbits.pc2au / θ.plx
    rdot = CATALOG.rv / PlanetOrbits.au2m * PlanetOrbits.year2sec_julian   # AU/yr
    # `one_pass − θ` is (2nd-pass offset − 1st-pass offset), and the final
    # system's reflex is the interim's *minus* ρ⋅(ṙ/d) — hence the sign.
    predicted = PlanetOrbits.raoff(sol, :A, bary) * rdot / d_au
    @test one_pass.pmra - θ.pmra ≈ -predicted rtol = 0.01
    @test abs(predicted) > 1e-3          # 10% of a DR3 σ: not a rounding term

    # § Radial velocity: ½|v_frame|²/c, absent from an interim with no space
    # motion. Read off the Einstein terms themselves (radvel − velz), because
    # the *total* one-pass rv error also carries the transverse coupling and
    # is 2.05 m/s rather than the 2.40 this accounts for.
    interim_sol = orbitsolve(interim, REF_EP)
    au_yr2ms = PlanetOrbits.au2m / PlanetOrbits.year2sec_julian
    ein(s, sys) = PlanetOrbits.radvel(s, :A, PlanetOrbits.barycentre(sys)) -
                  PlanetOrbits.velz(s, :A, PlanetOrbits.barycentre(sys)) * au_yr2ms
    v_tan = hypot(θ.pmra, θ.pmdec) / θ.plx * 4.740470446                   # km/s
    v_tot = hypot(v_tan * 1e3, CATALOG.rv)                                 # m/s
    @test ein(sol, posys) - ein(interim_sol, interim) ≈
          v_tot^2 / 2 / PlanetOrbits.c_light_ms rtol = 0.02
    @test abs(one_pass.rv - θ.rv) > 1.0                                    # m/s

    # And the residual left *after* the second pass is the second-order
    # feedback of the frame on its own offsets — four orders down from the
    # first-pass error, i.e. 10⁻⁵ of a DR3 σ.
    @test abs(PlanetOrbits.frame_pmra(sol) + PlanetOrbits.pmra(sol, :A, bary) -
              CATALOG.pmra) < 1e-6

    # `refine=false` is the escape hatch, and it reproduces the one-pass map
    # exactly — which is what makes the Jacobian below exactly unit.
    coarse = anchored_frame(interim, :A, REF_EP; refine=false,
        ra=CATALOG.ra, dec=CATALOG.dec, plx=CATALOG.plx,
        pmra=CATALOG.pmra, pmdec=CATALOG.pmdec, rv=CATALOG.rv)
    @test coarse.pmra == one_pass.pmra
end

@testset "the wide pair: two sources on one frame differ by the model's motion" begin
    # The structural claim behind multi-source astrometry. With the frame
    # anchored on A, the *absolute* proper motions the model predicts for A
    # and for a wide companion B differ by exactly the model's relative
    # proper motion — which is what a per-observation frame shift deletes.
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.29 end)
    B = Octofitter.Body(name="B", about=A, variables=@variables begin
        mass = 0.2256
        a = 750.0; e = 0.3; i = 1.2; ω = 0.4; Ω = 1.9; tp = 30000.0
    end)
    sys = Octofitter.System(name="wide", bodies=[A, B], verbosity=0,
        variables=AnchoredFrame(A; ref_epoch=REF_EP, variables=@variables begin
            ra_A = $(CATALOG.ra)
            dec_A = $(CATALOG.dec)
            plx_A = $(CATALOG.plx)
            pmra_A = $(CATALOG.pmra)
            pmdec_A = $(CATALOG.pmdec)
            rv_A = $(CATALOG.rv)
        end))
    (; θ, posys) = anchored_solved(sys)
    sol = orbitsolve(posys, REF_EP)
    bary = PlanetOrbits.barycentre(posys)
    pm_abs(nm) = (PlanetOrbits.frame_pmra(sol) + PlanetOrbits.pmra(sol, nm, bary),
                  PlanetOrbits.frame_pmdec(sol) + PlanetOrbits.pmdec(sol, nm, bary))
    pmA = pm_abs(:A); pmB = pm_abs(:B)
    rel = (PlanetOrbits.pmra(sol, :A, :B), PlanetOrbits.pmdec(sol, :A, :B))
    @test pmA[1] - pmB[1] ≈ rel[1] rtol = 1e-12
    @test pmA[2] - pmB[2] ≈ rel[2] rtol = 1e-12
    # And it is a difference a catalogue could see: σ(μ_DR3) is ~0.04 mas/yr,
    # so anything at this scale is many σ, which is exactly the signal the
    # shared frame exists to capture.
    @test hypot(rel...) > 0.1
    # A's own value is still the anchored one, so the pair is pinned to the
    # anchor and B floats relative to it — not the other way round.
    @test pmA[1] ≈ CATALOG.pmra atol = 1e-6
end

@testset "ForwardDiff runs through the whole anchored path" begin
    A = Octofitter.Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.29, 0.02), lower=0)
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass ~ LogUniform(1e-4, 1e-1)
        a ~ Uniform(1.0, 10.0)
        e = 0.35; i = 1.1; ω = 0.7; Ω = 2.3; tp = 51000.0
    end)
    sys = Octofitter.System(name="ad", bodies=[A, b], verbosity=0,
        variables=AnchoredFrame(A; ref_epoch=REF_EP, variables=@variables begin
            ra_A = $(CATALOG.ra)
            dec_A = $(CATALOG.dec)
            plx_A ~ Normal(CATALOG.plx, 0.2)
            pmra_A ~ Normal(CATALOG.pmra, 0.05)
            pmdec_A = $(CATALOG.pmdec)
            rv_A = $(CATALOG.rv)
        end))
    arr2nt = Octofitter.make_arr2nt(sys)
    x0 = collect(Octofitter.make_prior_sampler(sys)(Random.Xoshiro(3)))
    # The barycentric frame quantities are differentiable functions of every
    # sampled parameter, including the *body* ones — which is the whole
    # coupling the interim introduces.
    g = ForwardDiff.gradient(x -> arr2nt(x).pmra, x0)
    @test all(isfinite, g)
    @test !all(iszero, g)
    # ∂pmra/∂pmra_A is exactly 1: the correction subtracted from it is built
    # from the body variables and `plx_A`, never from `pmra_A` itself.
    names = Octofitter.list_parameter_names(sys)
    @test g[findfirst(==(:pmra_A), names)] ≈ 1.0 atol = 1e-9
    # …and it does depend on the companion's mass, which is the coupling.
    @test g[findfirst(==(:b_mass), names)] != 0.0

    # Through the log density too, not only through arr2nt.
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    @test all(isfinite, ForwardDiff.gradient(model.ℓπcallback, model.link(x0)))
end

@testset "the anchored → barycentric map is triangular with unit diagonal" begin
    # The Jacobian claim in `AnchoredFrame`'s docstring, measured — and the
    # price the refinement pass charges for it, measured too. The one-pass map
    # is exactly triangular with unit diagonal (each correction is built from
    # the body variables and the anchored parallax and nothing else); the
    # refined map couples every output weakly to every input.
    A, b = anchored_bodies(a=750.0, mass_b=0.2256)
    QS = (:ra, :dec, :plx, :pmra, :pmdec, :rv)
    x0 = Float64[CATALOG.ra, CATALOG.dec, CATALOG.plx,
                 CATALOG.pmra, CATALOG.pmdec, CATALOG.rv]

    function jac(refine)
        sys = Octofitter.System(name="jac", bodies=[A, b], verbosity=0,
            variables=AnchoredFrame(A; ref_epoch=REF_EP, refine, variables=@variables begin
                ra_A ~ Normal(CATALOG.ra, 1e-5)
                dec_A ~ Normal(CATALOG.dec, 1e-5)
                plx_A ~ Normal(CATALOG.plx, 0.2)
                pmra_A ~ Normal(CATALOG.pmra, 0.05)
                pmdec_A ~ Normal(CATALOG.pmdec, 0.04)
                rv_A ~ Normal(CATALOG.rv, 500.0)
            end))
        arr2nt = Octofitter.make_arr2nt(sys)
        names = Octofitter.list_parameter_names(sys)
        @test names == [:ra_A, :dec_A, :plx_A, :pmra_A, :pmdec_A, :rv_A]
        J = ForwardDiff.jacobian(x -> [getproperty(arr2nt(x), q) for q in QS], x0)
        return (J, arr2nt(x0))
    end

    # § One pass: exact.
    J1, θ1 = jac(false)
    for (k, q) in enumerate(QS)
        q === :plx && continue
        @test J1[k, k] == 1.0
    end
    @test J1[3, 3] ≈ (θ1.plx / θ1.plx_A)^2 rtol = 1e-10
    # Triangular: parallax feeds everything (it sets the AU→mas scale) and
    # nothing feeds parallax but itself.
    @test all(iszero, J1[3, [1, 2, 4, 5, 6]])
    @test det(J1) == prod(diag(J1))
    # …and the parallax entry is measurably not 1 at this separation, so the
    # check above is discriminating rather than a rounding coincidence. `dz`
    # is the *anchor's* offset from the barycentre — 750 au × m_b/M — not the
    # pair's separation, which is why it is 3e-5 and not 7e-3.
    @test abs(J1[3, 3] - 1) > 1e-5

    # § Two passes: the same map to within the feedback, which is what the
    # docstring quotes. Nothing here is exact any more, and nothing here is
    # large enough to matter either.
    J2, θ2 = jac(true)
    for (k, q) in enumerate(QS)
        q === :plx && continue
        @test J2[k, k] ≈ 1.0 atol = 2e-4
    end
    # Not a no-op: refinement does move the diagonal, though `ra` and `dec`
    # move by less than a Float64 ulp at 1.0 and so stay exactly 1.
    @test any(k -> J2[k, k] != J1[k, k], 1:6)
    @test J2[3, 3] ≈ J1[3, 3] rtol = 1e-6
    @test maximum(abs, J2[3, [1, 2, 4, 5, 6]]) < 1e-6
    @test det(J2) ≈ prod(diag(J2)) rtol = 1e-8
    # The whole volume factor stays within 1e-4 of the one-pass value, i.e.
    # the prior a user declares on the anchored variables is transported
    # essentially unchanged either way.
    @test det(J2) ≈ det(J1) rtol = 1e-4

    # For a planet-mass companion an AU out the parallax entry is 1 to 1e-6,
    # which is the other number the docstring quotes.
    let m = anchored_system(a=1.0, mass_b=0.001)
        θ = anchored_solved(m.sys).θ
        @test (θ.plx / θ.plx_A)^2 ≈ 1 atol = 1e-6
        @test (θ.plx / θ.plx_A)^2 != 1
    end
end

@testset "per-quantity anchoring, and the frame-free interim" begin
    # Anchoring proper motion but not parallax: `plx` is declared directly, so
    # nothing derives it, and `plx_interim` is not emitted — the interim falls
    # back to the non-deferred `plx`, which is exactly what it should use.
    A, b = anchored_bodies()
    sys = Octofitter.System(name="mixed", bodies=[A, b], verbosity=0,
        variables=AnchoredFrame(A; ref_epoch=REF_EP, plx=false, ra=false, dec=false,
            variables=@variables begin
                ra = $(CATALOG.ra)
                dec = $(CATALOG.dec)
                plx = $(CATALOG.plx)
                pmra_A = $(CATALOG.pmra)
                pmdec_A = $(CATALOG.pmdec)
                rv_A = $(CATALOG.rv)
            end))
    @test Octofitter._interim_plx(sys) === :plx
    @test !(INTERIM_PARALLAX_VAR in keys(sys.derived.variables))
    @test !(:plx in sys.deferred) && :pmra in sys.deferred
    θ = Octofitter.make_arr2nt(sys)(Float64[])
    @test θ.plx == CATALOG.plx           # untouched: not anchored
    @test θ.pmra != CATALOG.pmra         # anchored

    # An un-anchored quantity is passed through the map untouched — it still
    # goes *in*, because the second pass needs a complete frame to reframe
    # onto, but its `correct` slot is false and it comes back unchanged.
    sys2 = Octofitter.System(name="physical", bodies=[A, b], verbosity=0,
        variables=AnchoredFrame(A; ref_epoch=REF_EP, ra=false, dec=false,
            pmra=false, pmdec=false, variables=@variables begin
                ra = $(CATALOG.ra)
                dec = $(CATALOG.dec)
                pmra = $(CATALOG.pmra)
                pmdec = $(CATALOG.pmdec)
                plx_A = $(CATALOG.plx)
                rv_A = $(CATALOG.rv)
            end))
    @test occursin("(false, false, true, false, false, true)",
        string(sys2.derived.variables[:_frame]))
    θ2 = Octofitter.make_arr2nt(sys2)(Float64[])
    @test θ2.rv != CATALOG.rv
    @test θ2.pmra == CATALOG.pmra           # not anchored: untouched
    @test isfinite(θ2.plx)

    # Renaming the sampled variable.
    sys3 = Octofitter.System(name="renamed", bodies=[A, b], verbosity=0,
        variables=AnchoredFrame(A; ref_epoch=REF_EP, pmra=:pmra_A_dr3,
            pmdec=:pmdec_A_dr3, variables=@variables begin
                ra_A = $(CATALOG.ra); dec_A = $(CATALOG.dec); plx_A = $(CATALOG.plx)
                pmra_A_dr3 = $(CATALOG.pmra)
                pmdec_A_dr3 = $(CATALOG.pmdec)
                rv_A = $(CATALOG.rv)
            end))
    @test occursin("pmra_A_dr3", string(sys3.derived.variables[:_frame]))
    @test sys3.derived.variables[:pmra] == :(_frame.pmra)
end

@testset "the interim is deferred-scope only, and says so" begin
    A, b = anchored_bodies()
    # Note these blocks are built by hand rather than with `@variables`: `$` in
    # a `@variables` block is *value* interpolation, so `$INTERIM_SYSTEM_VAR`
    # would splice the Symbol as a captured constant and the expression would
    # never mention `system_interim` at all.
    blk(pairs...) = (Octofitter.Priors(),
        Octofitter.Derived(OrderedDict{Symbol,Any}(pairs...), (), ()))

    # A body block cannot see it: it is built *from* the bodies.
    badbody = Octofitter.Body(name="c", about=A, variables=blk(
        :mass => 0.02, :a => :(system_interim === nothing ? 5.0 : 5.0),
        :e => 0.1, :i => 1.0, :ω => 1.0, :Ω => 1.0, :tp => 50000.0))
    @test_throws r"only exists in \*deferred system lines\*" Octofitter.System(
        name="x", bodies=[A, badbody], verbosity=0,
        variables=@variables begin plx = 74.19 end)

    # Neither can an observation block: observations are evaluated against the
    # final system, so the interim is never the object they want.
    obs = Octofitter.RelAstromObs(
        (epoch=[57000.0], ra=[100.0], dec=[100.0], σ_ra=[1.0], σ_dec=[1.0], cor=[0.0]);
        target=b, ref=A, name="astrom",
        variables=blk(:jitter => :(system_interim === nothing ? 0.0 : 0.0)))
    @test_throws r"only exists in deferred system lines" Octofitter.System(
        name="y", bodies=[A, b], observations=[obs], verbosity=0,
        variables=@variables begin plx = 74.19 end)

    # And `plx_interim` cannot itself depend on the interim, because it is
    # what the interim is built at.
    @test_throws r"must be known before the bodies" Octofitter.System(
        name="z", bodies=[A, b], verbosity=0,
        variables=blk(
            :plx_interim => :(b.a > 0 ? 74.19 : 74.19),
            :plx => :(anchor_offsets(system_interim, :A, 57388.0).dz + 74.19)))
end

@testset "AnchoredFrame rejects the mistakes it can see" begin
    A, b = anchored_bodies()
    mk(; kwargs...) = AnchoredFrame(A; ref_epoch=REF_EP, kwargs...)
    full() = @variables begin
        ra_A = $(CATALOG.ra); dec_A = $(CATALOG.dec); plx_A = $(CATALOG.plx)
        pmra_A = $(CATALOG.pmra); pmdec_A = $(CATALOG.pmdec); rv_A = $(CATALOG.rv)
    end
    @test mk(variables=full()) isa Tuple

    # An un-anchored quantity with nothing defining it.
    @test_throws r"`rv` is not anchored" mk(rv=false, variables=@variables begin
        ra_A = $(CATALOG.ra); dec_A = $(CATALOG.dec); plx_A = $(CATALOG.plx)
        pmra_A = $(CATALOG.pmra); pmdec_A = $(CATALOG.pmdec)
    end)
    # A missing sampled variable, named, with the three ways out.
    @test_throws r"defines no `pmra_A`" mk(variables=@variables begin
        ra_A = $(CATALOG.ra); dec_A = $(CATALOG.dec); plx_A = $(CATALOG.plx)
        pmdec_A = $(CATALOG.pmdec); rv_A = $(CATALOG.rv)
    end)
    # Defining the barycentric quantity *and* anchoring it.
    @test_throws r"directly \*and\* anchors it" mk(variables=vcat(full(),
        @variables begin pmra = $(CATALOG.pmra) end))
    # Nothing anchored at all.
    @test_throws r"every quantity was passed as `false`" mk(
        ra=false, dec=false, plx=false, pmra=false, pmdec=false, rv=false,
        variables=full())
    # A keyword that is neither a name nor a switch.
    @test_throws r"takes `true`" mk(pmra=1.0, variables=full())
    # An anchor that is not a body.
    @test_throws r"must be a `Body` model node" AnchoredFrame(3.0; ref_epoch=REF_EP,
        variables=full())
end

@testset "anchor_offsets and barycentre_parallax on their own" begin
    # The two helpers are public, so pin their contracts directly rather than
    # only through the codegen path.
    A, b = anchored_bodies(a=750.0, mass_b=0.2256)
    poA = PlanetOrbits.Body(mass=1.29, name=:A)
    pob = PlanetOrbits.Body(mass=0.2256, name=:b)
    posys = PlanetOrbits.System((poA, pob),
        (PlanetOrbits.Orbit(pob; about=poA, a=750.0, e=0.35, i=1.1,
                            ω=0.7, Ω=2.3, tp=51000.0),); plx=CATALOG.plx)
    off = anchor_offsets(posys, :A, REF_EP)
    sol = orbitsolve(posys, REF_EP)
    bary = PlanetOrbits.barycentre(posys)
    @test off.ra_cosdec ≈ PlanetOrbits.raoff(sol, :A, bary)
    @test off.pmra ≈ PlanetOrbits.pmra(sol, :A, bary)
    @test off.rv ≈ PlanetOrbits.radvel(sol, :A, bary)
    @test off.dz ≈ PlanetOrbits.posz(sol, :A, bary)
    # A `Body` model node and its name are the same anchor.
    @test anchor_offsets(posys, A, REF_EP).pmra == off.pmra
    # The angular fields need a parallax, and the map says so by name rather
    # than by MethodError.
    nofr = PlanetOrbits.reframe(posys)
    @test_throws r"has no parallax" anchor_offsets(nofr, :A, REF_EP)
    @test_throws r"needs the interim system to carry a parallax" anchored_frame(
        nofr, :A, REF_EP; ra=CATALOG.ra, dec=CATALOG.dec, plx=CATALOG.plx,
        pmra=CATALOG.pmra, pmdec=CATALOG.pmdec, rv=CATALOG.rv)

    # `anchored_frame` with everything switched off is the identity — nothing
    # is corrected, so nothing changes, two passes or not.
    idn = anchored_frame(posys, :A, REF_EP, ntuple(_ -> false, 6);
        ra=CATALOG.ra, dec=CATALOG.dec, plx=CATALOG.plx,
        pmra=CATALOG.pmra, pmdec=CATALOG.pmdec, rv=CATALOG.rv)
    @test all(q -> getproperty(idn, q) == getproperty(CATALOG, q),
              (:ra, :dec, :plx, :pmra, :pmdec, :rv))

    # `barycentre_parallax` inverts the distance, exactly.
    plx = barycentre_parallax(CATALOG.plx, off.dz)
    d_anchor = 1000 * PlanetOrbits.pc2au / CATALOG.plx
    @test 1000 * PlanetOrbits.pc2au / plx ≈ d_anchor - off.dz rtol = 1e-14
    @test barycentre_parallax(CATALOG.plx, 0.0) == CATALOG.plx
    # A body *behind* the barycentre gives a smaller parallax, and vice versa.
    @test barycentre_parallax(CATALOG.plx, 1000.0) > CATALOG.plx
    @test barycentre_parallax(CATALOG.plx, -1000.0) < CATALOG.plx
end
