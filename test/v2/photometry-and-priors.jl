# Photometry and the configuration priors — `design/observation-types-migration.md`
# §3.4 (PhotometryObs), §3.7 (the dynamical priors and O'Neil's Jacobian).
#
# The priors are the interesting half to test: they carry no data, so the only
# thing worth asserting is that the number they add to the log density is the
# number v1 added. Each v1 formula is transcribed here from
# `src/legacy/likelihoods/` and compared against the port on a fixed
# configuration — exactly, where the arithmetic is unchanged.

# --- shared scaffolding -------------------------------------------------------

# A two-companion system whose orbital elements are all fixed, so every model
# below evaluates at a single known point with an empty parameter vector.
function pp_system(; obs, a_b=3.0, e_b=0.15, a_c=9.0, e_c=0.05,
                   m_b=2mjup, m_c=4mjup, flux_H_b=15.2, jacobi=false)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
        flux_H = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = $m_b
        flux_H = $flux_H_b
        a = $a_b; e = $e_b; i = 0.9; ω = 1.2; Ω = 2.0; tp = 55000.0
    end)
    c = Octofitter.Body(name="c", about=(jacobi ? (A, b) : A), variables=@variables begin
        mass = $m_c
        a = $a_c; e = $e_c; i = 1.1; ω = 0.3; Ω = 1.0; tp = 56000.0
    end)
    return Octofitter.System(name="pp", bodies=[A, b, c], observations=obs,
        variables=@variables begin plx = 25.0 end)
end

# The evaluation context one observation of `sys` would be handed, at the
# model's single parameter point.
function pp_ctx(sys, obs)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep)
    nm = Symbol(Octofitter.normalizename(Octofitter.likelihoodname(obs)))
    θ_obs = hasproperty(nt.observations, nm) ? getproperty(nt.observations, nm) : (;)
    return Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs])
end

pp_lnlike(sys) = Octofitter.make_ln_like(sys)(sys, Octofitter.make_arr2nt(sys)(Float64[]))

const PP_PHOT = Table(phot=[15.0, 14.8, 15.6], σ_phot=[3.0, 0.5, 1.2])

# ---------------------------------------------------------------------------
# PhotometryObs
# ---------------------------------------------------------------------------

@testset "PhotometryObs: construction and interface" begin
    b = Octofitter.Body(name="b", variables=@variables begin mass = 1.0 end)
    obs = PhotometryObs(PP_PHOT; target=b, band=:H, name="NIRC2")

    @test Octofitter.refspecs(obs) === (Octofitter.BodyRefSpec{:b}(),)
    @test Octofitter._photband(obs) === :H
    @test Octofitter._photvar(obs) === :flux_H
    # It is an observation, not a prior term: that is the whole reason the type
    # survives the refactor (§3.4).
    @test Octofitter._isprior(obs) == false
    # It touches no epoch, so it contributes nothing to the trajectory solve
    # even though a user's table may carry an :epoch column for bookkeeping.
    @test Octofitter.epochs(obs) == Float64[]
    @test Octofitter.epochs(PhotometryObs(
        Table(epoch=[57000.0, 57100.0], phot=[1.0, 2.0], σ_phot=[1.0, 1.0]);
        target=b, band=:H, name="withepoch")) == Float64[]

    # Default band is the plain `flux` variable, matching `_flux_band`.
    plain = PhotometryObs(PP_PHOT; target=b, name="plain")
    @test Octofitter._photvar(plain) === :flux
    @test Octofitter._photband(plain) === :default

    @test occursin("b.flux_H", sprint(show, MIME"text/plain"(), obs))

    @test_throws r"Expected columns" PhotometryObs(
        Table(phot=[1.0]); target=b, name="x")
    # A photocentre is a point on the sky, not a brightness.
    @test_throws r"must be a `Body`" PhotometryObs(
        PP_PHOT; target=Photocentre(:H), name="x")
    @test_throws r"must be a `Body`" PhotometryObs(
        PP_PHOT; target=Barycentre, name="x")
end

@testset "PhotometryObs: reproduces v1's log-likelihood bit-for-bit" begin
    obs = PhotometryObs(PP_PHOT; target=:b, band=:H, name="NIRC2")
    sys = pp_system(obs=[obs], flux_H_b=15.2)
    ctx = pp_ctx(sys, obs)

    # v1's `ln_like(::PhotometryObs, ::PlanetObservationContext)`, transcribed
    # from src/legacy/likelihoods/photometry.jl. The only difference is where
    # the flux comes from: v1 read `θ_obs.flux`, declared on the observation;
    # v2 reads the *body's* `flux_H`, which is the flux/band unification.
    function v1_phot_ll(flux, tab)
        ll = 0.0
        for i in eachindex(tab.phot)
            resid = flux - tab.phot[i]
            σ² = tab.σ_phot[i]^2
            ll += -(1 / 2) * resid^2 / σ² - log(sqrt(2π * σ²))
        end
        return ll
    end

    @test Octofitter.ln_like(obs, ctx) === v1_phot_ll(15.2, PP_PHOT)
    # …and the same number reaches the model's total log density.
    @test pp_lnlike(sys) ≈ v1_phot_ll(15.2, PP_PHOT) rtol = 1e-14

    # The flux really is the body's variable, not a constant.
    sys2 = pp_system(obs=[obs], flux_H_b=14.0)
    @test pp_lnlike(sys2) ≈ v1_phot_ll(14.0, PP_PHOT) rtol = 1e-14
end

@testset "PhotometryObs: subsetting, simulation, and a missing flux variable" begin
    obs = PhotometryObs(PP_PHOT; target=:b, band=:H, name="NIRC2")
    sys = pp_system(obs=[obs], flux_H_b=15.2)
    ctx = pp_ctx(sys, obs)

    # Cross-validation needs this; a `~` line in the body's block could not
    # provide it, which is §3.4's argument for keeping the type.
    sub = Octofitter.likeobj_from_epoch_subset(obs, 2:3)
    @test sub isa PhotometryObs
    @test vec(sub.table.phot) == PP_PHOT.phot[2:3]
    @test Octofitter._photvar(sub) === :flux_H
    @test Octofitter.likelihoodname(sub) == "NIRC2"

    @test Octofitter.simulate(obs, ctx).phot_model == 15.2

    # …and so does simulation.
    gen = Octofitter.generate_from_params(obs, ctx; add_noise=false)
    @test gen isa PhotometryObs
    @test all(==(15.2), gen.table.phot)
    @test vec(gen.table.σ_phot) == PP_PHOT.σ_phot
    noisy = Octofitter.generate_from_params(obs, ctx; add_noise=true)
    @test noisy.table.phot != gen.table.phot

    # `PlanetOrbits.fluxes` would report a body that declares no flux in the
    # band as zero, silently comparing 0.0 against the data forever. Reading the
    # variable turns that into a message naming what is missing.
    badobs = PhotometryObs(PP_PHOT; target=:b, band=:Ks, name="bad")
    badsys = pp_system(obs=[badobs])
    @test_throws r"declares no `flux_Ks`" Octofitter.ln_like(badobs, pp_ctx(badsys, badobs))
end

@testset "PhotometryObs is differentiable in the body's flux" begin
    obs = PhotometryObs(PP_PHOT; target=:b, band=:H, name="NIRC2")
    sys = pp_system(obs=[obs], flux_H_b=15.2)
    ctx = pp_ctx(sys, obs)

    # A Dual flux against a plain-`Float64` solved system: the mixed case, and
    # the one the port has to survive, since a sampled flux touches θ but not
    # the trajectory. (See the `@test_broken` below for why it has to be built
    # by hand rather than by sampling `flux_H`.)
    FD = Octofitter.ForwardDiff
    dflux = FD.Dual{Nothing}(15.2, 1.0)
    dbodies = merge(ctx.θ_system.bodies,
        (; b=merge(ctx.θ_system.bodies.b, (; flux_H=dflux))))
    dctx = Octofitter.ObsContext(merge(ctx.θ_system, (; bodies=dbodies)),
        ctx.θ_obs, ctx.system, ctx.traj, ctx.epoch_index)
    ll = Octofitter.ln_like(obs, dctx)
    @test FD.value(ll) === Octofitter.ln_like(obs, ctx)
    # d/df Σ −(f − pᵢ)²/2σᵢ² = −Σ (f − pᵢ)/σᵢ²
    @test FD.partials(ll)[1] ≈
          -sum((15.2 .- PP_PHOT.phot) ./ PP_PHOT.σ_phot .^ 2) rtol = 1e-12

    # A non-finite flux (an atmosphere-model interpolation off the end of its
    # grid) rejects the sample rather than poisoning the log density.
    nanbodies = merge(ctx.θ_system.bodies,
        (; b=merge(ctx.θ_system.bodies.b, (; flux_H=NaN))))
    nanctx = Octofitter.ObsContext(merge(ctx.θ_system, (; bodies=nanbodies)),
        ctx.θ_obs, ctx.system, ctx.traj, ctx.epoch_index)
    @test Octofitter.ln_like(obs, nanctx) == -Inf
end

@testset "sampling a body flux, with nothing else free" begin
    # `PlanetOrbits.System` used to promote its scalar type over the masses, the
    # orbital elements and the frame but *not* over the bodies' fluxes, so
    # `_collect_fluxes` narrowed a Dual back to Float64, `T(::Dual)` threw, and
    # `make_ln_like`'s guard turned the whole log density into `-Inf` with a
    # zero gradient. Fixed upstream (`_flux_scalar_type` now joins that
    # `promote_type` call), which is what makes a sampled contrast usable.
    #
    # This was always the degenerate case — it bit only when the flux was the
    # *sole* sampled quantity entering the promotion. Sample any mass, element
    # or parallax alongside it, as every real model does, and the parameter
    # vector was already Dual. That path is covered in `v2/api-smoke.jl`; this
    # testset is the narrow one that used to fail.
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
        flux_H = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 2mjup
        flux_H ~ Uniform(0.0, 40.0)
        a = 3.0; e = 0.1; i = 0.5; ω = 0.2; Ω = 0.3; tp = 55000.0
    end)
    obs = PhotometryObs(PP_PHOT; target=b, band=:H, name="NIRC2")
    sys = Octofitter.System(name="phot", bodies=[A, b], observations=[obs],
        variables=@variables begin plx = 25.0 end)
    arr2nt = Octofitter.make_arr2nt(sys)
    lnl = Octofitter.make_ln_like(sys)
    @test isfinite(lnl(sys, arr2nt([15.2])))          # the Float64 path is fine
    g = Octofitter.ForwardDiff.derivative(x -> lnl(sys, arr2nt([x])), 15.2)
    @test isfinite(g)
    # d/df Σ −(f − pᵢ)²/2σᵢ² = −Σ (f − pᵢ)/σᵢ², plus the Uniform prior's zero
    # contribution, so the gradient is the photometric residual sum.
    @test g ≈ -sum((15.2 .- PP_PHOT.phot) ./ PP_PHOT.σ_phot .^ 2) rtol = 1e-10
end

# ---------------------------------------------------------------------------
# OrbitOrderPrior
# ---------------------------------------------------------------------------

@testset "OrbitOrderPrior" begin
    prior = OrbitOrderPrior(:b, :c)
    @test Octofitter._isprior(prior)
    @test Octofitter.refspecs(prior) ===
          (Octofitter.BodyRefSpec{:b}(), Octofitter.BodyRefSpec{:c}())
    @test Octofitter.epochs(prior) == Float64[]
    @test Octofitter.likeobj_from_epoch_subset(prior, 1:1) === prior
    @test occursin("b < c", sprint(show, MIME"text/plain"(), prior))
    @test PlanetOrderPrior === OrbitOrderPrior      # v1 name still resolves

    @test_throws r"at least two bodies" OrbitOrderPrior(:b)
    @test_throws r"must be a `Body`" OrbitOrderPrior(:b, Barycentre)

    ordered = pp_system(obs=[prior], a_b=3.0, a_c=9.0)
    swapped = pp_system(obs=[prior], a_b=9.0, a_c=3.0)
    @test Octofitter.ln_like(prior, pp_ctx(ordered, prior)) == 0.0
    @test Octofitter.ln_like(prior, pp_ctx(swapped, prior)) == -Inf
    # Equal semi-major axes are rejected, as in v1 (`>=`).
    equal = pp_system(obs=[prior], a_b=5.0, a_c=5.0)
    @test Octofitter.ln_like(prior, pp_ctx(equal, prior)) == -Inf

    # Reversing the argument order reverses the assertion.
    rev = OrbitOrderPrior(:c, :b)
    revsys = pp_system(obs=[rev], a_b=9.0, a_c=3.0)
    @test Octofitter.ln_like(rev, pp_ctx(revsys, rev)) == 0.0

    # The named body's *placing row* is what is read, so the same physical
    # system gives the same answer under either hierarchy convention.
    jac = pp_system(obs=[prior], a_b=3.0, a_c=9.0, jacobi=true)
    @test Octofitter.ln_like(prior, pp_ctx(jac, prior)) == 0.0

    # A root body has no orbit of its own.
    rootprior = OrbitOrderPrior(:A, :b)
    @test_throws r"not placed by any orbit" Octofitter.ln_like(
        rootprior, pp_ctx(pp_system(obs=[rootprior]), rootprior))
end

# ---------------------------------------------------------------------------
# NonCrossingPrior / LimitClosestApproachAUPrior
# ---------------------------------------------------------------------------

# v1's `ln_like(::LimitClosestApproachAUPrior, …)`, transcribed from
# src/legacy/likelihoods/prior-non-crossing.jl. `orbits` is a vector of (a, e).
function v1_closest_approach(orbits, hard, soft)
    length(orbits) == 1 && return 0.0
    ll = 0.0
    srt = sort(collect(orbits), by=first)
    for (oa, ob) in zip(srt[1:end-1], srt[2:end])
        sep_farthest_inner = oa[1] * (1 + oa[2])
        sep_closest_outer = ob[1] * (1 - ob[2])
        closest_approach = sep_closest_outer - sep_farthest_inner
        if closest_approach <= hard
            return -Inf
        end
        if closest_approach < soft
            ll -= 1 / (closest_approach - soft)^2
        end
    end
    return ll
end

@testset "NonCrossingPrior / LimitClosestApproachAUPrior matches v1" begin
    prior = NonCrossingPrior()
    @test Octofitter._isprior(prior)
    @test Octofitter.refspecs(prior) === ()
    @test occursin("every orbit", sprint(show, MIME"text/plain"(), prior))
    @test prior isa LimitClosestApproachAUPrior
    @test prior.hard_closest_approach_au == 0.0
    @test prior.soft_closest_approach_au == 0.0

    # Well-separated: both give exactly zero.
    wide = ((3.0, 0.15), (9.0, 0.05))
    sysw = pp_system(obs=[prior], a_b=3.0, e_b=0.15, a_c=9.0, e_c=0.05)
    @test Octofitter.ln_like(prior, pp_ctx(sysw, prior)) ==
          v1_closest_approach(collect(wide), 0.0, 0.0)

    # Crossing: both give -Inf. (b's apoapsis 3·1.6 = 4.8 > c's periapsis
    # 5.0·0.4 = 2.0.)
    cross = ((3.0, 0.6), (5.0, 0.6))
    sysx = pp_system(obs=[prior], a_b=3.0, e_b=0.6, a_c=5.0, e_c=0.6)
    @test Octofitter.ln_like(prior, pp_ctx(sysx, prior)) == -Inf
    @test v1_closest_approach(collect(cross), 0.0, 0.0) == -Inf

    # The soft repulsive term, bit-for-bit — this is the only branch with any
    # arithmetic in it, so it is the one worth pinning exactly.
    soft = LimitClosestApproachAUPrior(2.0)
    @test soft.hard_closest_approach_au == 0.0
    @test soft.soft_closest_approach_au == 2.0
    near = ((3.0, 0.1), (4.5, 0.1))
    sysn = pp_system(obs=[soft], a_b=3.0, e_b=0.1, a_c=4.5, e_c=0.1)
    v1 = v1_closest_approach(collect(near), 0.0, 2.0)
    @test v1 < 0 && isfinite(v1)
    @test Octofitter.ln_like(soft, pp_ctx(sysn, soft)) === v1

    # A non-zero hard floor rejects what the soft term would only penalise.
    hard = LimitClosestApproachAUPrior(1.0, 2.0)
    sysh = pp_system(obs=[hard], a_b=3.0, e_b=0.1, a_c=4.5, e_c=0.1)
    @test Octofitter.ln_like(hard, pp_ctx(sysh, hard)) == -Inf
    @test v1_closest_approach(collect(near), 1.0, 2.0) == -Inf

    # The sort is by semi-major axis, not by declaration order: swapping which
    # body is inner must not change the answer.
    swapped = pp_system(obs=[soft], a_b=4.5, e_b=0.1, a_c=3.0, e_c=0.1)
    @test Octofitter.ln_like(soft, pp_ctx(swapped, soft)) === v1
end

@testset "LimitClosestApproachAUPrior: bodies= and no allocations" begin
    # `bodies=` restricts the prior to the rows placing those bodies; with no
    # list every row is included, which is what v1 always did.
    restricted = NonCrossingPrior(bodies=(:b, :c))
    sysx = pp_system(obs=[restricted], a_b=3.0, e_b=0.6, a_c=5.0, e_c=0.6)
    @test Octofitter.ln_like(restricted, pp_ctx(sysx, restricted)) == -Inf
    @test occursin("orbits of b, c", sprint(show, MIME"text/plain"(), restricted))
    # One orbit named: no adjacent pair, so nothing to say.
    single = NonCrossingPrior(bodies=(:b,))
    @test Octofitter.ln_like(single, pp_ctx(pp_system(obs=[single]), single)) == 0.0
    @test_throws r"must be a `Body`" NonCrossingPrior(bodies=(Barycentre,))

    # v1 built a `sort(collect(orbits))` Vector on *every* evaluation and
    # carried a `# TODO: would be nice to make this non-allocating`. The
    # fixed-width walk that replaced it allocates nothing.
    prior = NonCrossingPrior()
    sys = pp_system(obs=[prior])
    ctx = pp_ctx(sys, prior)
    Octofitter.ln_like(prior, ctx)         # compile
    @test (@allocated Octofitter.ln_like(prior, ctx)) == 0
end

# ---------------------------------------------------------------------------
# HillStabilityPrior
# ---------------------------------------------------------------------------

@testset "HillStabilityPrior" begin
    prior = HillStabilityPrior()
    @test Octofitter._isprior(prior)
    @test Octofitter.epochs(prior) == Float64[]

    # The Gladman criterion, written out independently of the port.
    function hill_ok(a_in, a_out, m_in, m_out, M★)
        R_H = a_out * cbrt((m_in + m_out) / (3M★))
        return (a_out - a_in) > 2 * sqrt(3) * R_H
    end

    m_b, m_c = 2mjup, 4mjup
    tight = pp_system(obs=[prior], a_b=3.0, a_c=3.2, m_b=m_b, m_c=m_c)
    loose = pp_system(obs=[prior], a_b=3.0, a_c=30.0, m_b=m_b, m_c=m_c)
    @test !hill_ok(3.0, 3.2, m_b, m_c, 1.0)
    @test hill_ok(3.0, 30.0, m_b, m_c, 1.0)
    @test Octofitter.ln_like(prior, pp_ctx(tight, prior)) == -Inf
    @test Octofitter.ln_like(prior, pp_ctx(loose, prior)) == 0.0

    # M★ is the mass interior to the outer row *excluding* the inner row's
    # bodies, so the astrocentric and Jacobi spellings of the same physical
    # system agree. v1's `θ_planet.M` did not have this property.
    jac = pp_system(obs=[prior], a_b=3.0, a_c=30.0, m_b=m_b, m_c=m_c, jacobi=true)
    @test Octofitter.ln_like(prior, pp_ctx(jac, prior)) == 0.0
    jact = pp_system(obs=[prior], a_b=3.0, a_c=3.2, m_b=m_b, m_c=m_c, jacobi=true)
    @test Octofitter.ln_like(prior, pp_ctx(jact, prior)) == -Inf

    # Heavier companions push the Hill radius out, so a separation that passes
    # at low mass fails at high mass. This is exactly what v1 could not see:
    # its `θ_planet_b = θ_system.planets[idx_a]` line overwrote the outer
    # planet's parameters with the inner one's, so `m_out` never entered.
    a_out = 8.0
    @test hill_ok(3.0, a_out, m_b, 4mjup, 1.0)
    @test !hill_ok(3.0, a_out, m_b, 400mjup, 1.0)
    light = pp_system(obs=[prior], a_b=3.0, a_c=a_out, m_b=m_b, m_c=4mjup)
    heavy = pp_system(obs=[prior], a_b=3.0, a_c=a_out, m_b=m_b, m_c=400mjup)
    @test Octofitter.ln_like(prior, pp_ctx(light, prior)) == 0.0
    @test Octofitter.ln_like(prior, pp_ctx(heavy, prior)) == -Inf

    ctx = pp_ctx(loose, prior)
    Octofitter.ln_like(prior, ctx)
    @test (@allocated Octofitter.ln_like(prior, ctx)) == 0
end

# ---------------------------------------------------------------------------
# ObsPriorONeil2019
# ---------------------------------------------------------------------------

const PP_EPOCHS = [50000.0, 50100.0, 50250.0, 50600.0]

function pp_astrom(; target=:b, ref=:A, name="astrom")
    return RelAstromObs((epoch=PP_EPOCHS,
            ra=[100.0, 110.0, 118.0, 130.0], dec=[50.0, 55.0, 61.0, 70.0],
            σ_ra=fill(1.0, 4), σ_dec=fill(1.0, 4)); target, ref, name)
end

# The O'Neil et al. 2019 Jacobian, written from the paper's expression rather
# than from the port: mean/eccentric anomalies from an independent Newton solve
# of Kepler's equation, and the mean anomaly built as 2π(t − tp)/P.
function v1_oneil_term(P_days, e, tp, epochs)
    jac = 0.0
    for t in epochs
        # v1 fed an *unwrapped* `n(t − tp)` to `kepler_solver`, whose Markley
        # branch wraps it to [−π, π] before solving, and then rebuilt the mean
        # anomaly from the solved E. The Jacobian is not periodic in M, so the
        # wrap is part of the definition, not an implementation detail.
        MA = rem2pi(2π * (t - tp) / P_days, RoundNearest)
        E = MA
        for _ in 1:100
            E -= (E - e * sin(E) - MA) / (1 - e * cos(E))
        end
        M = E - e * sin(E)
        jac += abs(3M * (e + cos(E)) + 2 * (-2 + e^2 + e * cos(E)) * sin(E))
    end
    jac *= cbrt(P_days / 365.25) / sqrt(1 - e^2)
    return 2log(jac)
end

@testset "ObsPriorONeil2019: interface" begin
    astrom = pp_astrom()
    prior = ObsPriorONeil2019(astrom)
    @test prior.wrapped_like === astrom
    @test prior.orbit === (Octofitter.BodyRefSpec{:b}(),)      # inherited target
    @test Octofitter.likelihoodname(prior) == "obspri_astrom"
    # It wraps data, so it is a likelihood, not a prior term — v1 had this too.
    @test Octofitter._isprior(prior) == false
    @test Octofitter.epochs(prior) == PP_EPOCHS
    @test ObsPriorAstromONeil2019 === ObsPriorONeil2019          # v1 name

    sub = Octofitter.likeobj_from_epoch_subset(prior, 1:2)
    @test sub isa ObsPriorONeil2019
    @test sub.wrapped_like.table.epoch == PP_EPOCHS[1:2]
    @test sub.orbit === prior.orbit

    @test_throws r"must be a `Body`" ObsPriorONeil2019(astrom; orbit=Barycentre)
end

@testset "ObsPriorONeil2019: adds exactly the paper's Jacobian" begin
    astrom = pp_astrom()
    prior = ObsPriorONeil2019(astrom)
    sysw = pp_system(obs=[astrom])       # the wrapped likelihood alone
    sysp = pp_system(obs=[prior])        # …and wrapped in the prior

    posys = Octofitter.make_ln_like(sysw).build(Octofitter.make_arr2nt(sysw)(Float64[]))
    P_days = PlanetOrbits.period(posys, 1)
    expected = v1_oneil_term(P_days, 0.15, 55000.0, PP_EPOCHS)
    @test isfinite(expected)

    # The wrapper's contribution over the wrapped likelihood is the Jacobian.
    @test pp_lnlike(sysp) - pp_lnlike(sysw) ≈ expected rtol = 1e-9

    # The Jacobian is a property of the orbit and the epochs, so wrapping a
    # *radial-velocity* likelihood at the same epochs adds the same term. This
    # is the collapse of v1's separate `prior-observable-rv.jl` dispatch.
    rvs = RadialVelocityObs((epoch=PP_EPOCHS, rv=[10.0, -5.0, 3.0, 8.0],
        σ_rv=fill(2.0, 4)); target=:A, ref=Barycentre, name="rv")
    rvprior = ObsPriorONeil2019(rvs; orbit=:b)
    sysr = pp_system(obs=[rvs])
    sysrp = pp_system(obs=[rvprior])
    @test pp_lnlike(sysrp) - pp_lnlike(sysr) ≈ expected rtol = 1e-9

    # Two orbits sum, as v1's system-attached method did over every planet.
    two = ObsPriorONeil2019(astrom; orbit=(:b, :c))
    systwo = pp_system(obs=[two])
    expected_c = v1_oneil_term(PlanetOrbits.period(posys, 2), 0.05, 56000.0, PP_EPOCHS)
    @test pp_lnlike(systwo) - pp_lnlike(sysw) ≈ expected + expected_c rtol = 1e-9

    # A stellar-reflex RV has no orbit of its own to inherit, and says so.
    bad = ObsPriorONeil2019(rvs)
    @test bad.orbit === (Octofitter.BodyRefSpec{:A}(),)
    @test_throws r"not placed by any orbit" Octofitter.ln_like(
        bad, pp_ctx(pp_system(obs=[bad]), bad))
end

@testset "ObsPriorONeil2019 is finite and differentiable" begin
    astrom = pp_astrom()
    prior = ObsPriorONeil2019(astrom)
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 2mjup
        P_yr ~ Uniform(0.1, 50.0)              # years — the prior assumes this
        a = cbrt(system.M_tot * P_yr^2)
        e ~ Uniform(0.0, 0.9)
        i = 0.9; ω = 1.2; Ω = 2.0; tp = 55000.0
    end)
    sys = Octofitter.System(name="oneil", bodies=[A, b], observations=[prior],
        variables=@variables begin
            plx = 25.0
            M_tot = 1.0
        end)
    arr2nt = Octofitter.make_arr2nt(sys)
    lnl = Octofitter.make_ln_like(sys)
    f = x -> lnl(sys, arr2nt(x))
    x0 = [4.0, 0.3]
    @test isfinite(f(x0))
    g = Octofitter.ForwardDiff.gradient(f, x0)
    @test all(isfinite, g)
    @test g ≈ FiniteDiff.finite_difference_gradient(f, x0) rtol = 1e-4
end

@testset "all five terms compose in one system" begin
    astrom = pp_astrom()
    obs = [ObsPriorONeil2019(astrom),
           PhotometryObs(PP_PHOT; target=:b, band=:H, name="NIRC2"),
           NonCrossingPrior(),
           OrbitOrderPrior(:b, :c),
           HillStabilityPrior()]
    sys = pp_system(obs=obs)
    ll = pp_lnlike(sys)
    @test isfinite(ll)

    # Each term is type-stable and allocation-free on the hot path — the
    # priors because the sort was replaced by a fixed-width walk, the O'Neil
    # wrapper because the epoch scan was replaced by `solutionat`.
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep)
    for o in obs
        nm = Symbol(Octofitter.normalizename(Octofitter.likelihoodname(o)))
        θ_obs = hasproperty(nt.observations, nm) ? getproperty(nt.observations, nm) : (;)
        ctx = Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[o])
        @test only(Base.return_types(Octofitter.ln_like, (typeof(o), typeof(ctx)))) === Float64
        Octofitter.ln_like(o, ctx)
        @test (@allocated Octofitter.ln_like(o, ctx)) == 0
        @test !isempty(sprint(show, MIME"text/plain"(), o))
    end
end
