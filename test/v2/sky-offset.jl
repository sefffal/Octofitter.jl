# The shared observation front-ends: `resolverefs`, `sky_offset!`/`sky_offset`
# and the structured `show` description. Seven likelihood ports build on
# these, so the conventions are pinned here rather than in each of them.
#
# `reference_system`, `model_system` and `EPOCHS_A` come from v2/likelihoods.jl.

# Build an `ObsContext` for `obs` inside `sys`, at the parameters' point values.
function _ctx_for(sys, obs)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep)
    sym = Symbol(Octofitter.likelihoodname(obs))
    # An observation that declares no variables gets no entry at all.
    θ_obs = hasproperty(nt.observations, sym) ? getproperty(nt.observations, sym) : (;)
    return Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs])
end

function _astrom_obs(; name="astrom", variables=(Octofitter.Priors(), Octofitter.Derived()))
    n = length(EPOCHS_A)
    return RelAstromObs((epoch=EPOCHS_A, ra=fill(1.0, n), dec=fill(1.0, n),
        σ_ra=fill(1.0, n), σ_dec=fill(1.0, n)); target=:c, ref=:b, name, variables)
end

@testset "resolverefs" begin
    posys = reference_system()
    obs = _astrom_obs()
    ctx = _ctx_for(model_system(obs=[obs]), obs)

    specs = (Octofitter.refspec(:b), Barycentre, Photocentre)
    refs = Octofitter.resolverefs(ctx, specs)
    @test refs isa Tuple{PlanetOrbits.BodyRef,PlanetOrbits.WeightedPoint,
                         PlanetOrbits.WeightedPoint}
    @test refs === map(s -> Octofitter.ref(ctx, s), specs)

    # The whole reason it exists: a tuple stays concrete and allocation-free
    # even when the resolved types differ, where a comprehension would widen.
    @test isconcretetype(typeof(refs))
    f(ctx, specs) = Octofitter.resolverefs(ctx, specs)
    f(ctx, specs)
    @test (@allocated f(ctx, specs)) == 0

    # A vector of specs is rejected rather than silently deoptimized.
    @test_throws r"takes a `Tuple`" Octofitter.resolverefs(ctx, collect(specs))
end

@testset "sky_offset: the identity calibration is bit-exact" begin
    obs = _astrom_obs()
    ctx = _ctx_for(model_system(obs=[obs]), obs)
    t, r = Octofitter.ref(ctx, obs.target), Octofitter.ref(ctx, obs.ref)
    n = length(EPOCHS_A)

    Δα, Δδ = zeros(n), zeros(n)
    Octofitter.sky_offset!(Δα, Δδ, ctx, obs.target, obs.ref)
    for i in 1:n
        sol = Octofitter.solutionat(ctx, i)
        # `===` not `≈`: with platescale 1 and northangle 0 the rotate-and-
        # scale block must not perturb the last bit, or every port that routes
        # through it stops reproducing v1 exactly.
        @test Δα[i] === raoff(sol, t, r)
        @test Δδ[i] === decoff(sol, t, r)
    end

    # Specs and already-resolved refs are interchangeable at the call site.
    Δα2, Δδ2 = zeros(n), zeros(n)
    Octofitter.sky_offset!(Δα2, Δδ2, ctx, t, r)
    @test Δα2 == Δα && Δδ2 == Δδ

    # The scalar variant agrees with the vectorized one.
    sol = Octofitter.solutionat(ctx, 3)
    @test Octofitter.sky_offset(sol, t, r) === (Δα[3], Δδ[3])

    # Buffers are overwritten, not accumulated into (unlike
    # `accumulate_offsets!`, whose callers pre-load a linear motion).
    fill!(Δα2, 99.0); fill!(Δδ2, 99.0)
    Octofitter.sky_offset!(Δα2, Δδ2, ctx, obs.target, obs.ref)
    @test Δα2 == Δα && Δδ2 == Δδ

    @test_throws r"indexed by table row" Octofitter.sky_offset!(
        zeros(n - 1), zeros(n - 1), ctx, obs.target, obs.ref)
end

@testset "sky_offset: platescale and northangle convention" begin
    obs = _astrom_obs()
    ctx = _ctx_for(model_system(obs=[obs]), obs)
    t, r = Octofitter.ref(ctx, obs.target), Octofitter.ref(ctx, obs.ref)
    sol = Octofitter.solutionat(ctx, 4)

    ra0, dec0 = Octofitter.sky_offset(sol, t, r)
    sep0, pa0 = hypot(ra0, dec0), atan(ra0, dec0)   # PA: north through east

    # The helper carries the model onto the detector, so it applies the
    # *inverse* of the correction a user applies to reported data:
    #   sep_true = sep_reported × platescale,  pa_true = pa_reported + northangle
    ps, na = 1.07, 0.13
    ra1, dec1 = Octofitter.sky_offset(sol, t, r; platescale=ps, northangle=na)
    @test hypot(ra1, dec1) ≈ sep0 / ps
    @test mod2pi(atan(ra1, dec1) - (pa0 - na) + π) - π ≈ 0 atol = 1e-12

    # …which is the same rotation `RelAstromObs` applies to its data, in the
    # opposite direction: rotating the data by +northangle and rotating the
    # model by −northangle land on the same relative geometry.
    @test ra1 * ps ≈ sep0 * sin(pa0 - na)
    @test dec1 * ps ≈ sep0 * cos(pa0 - na)

    # Non-allocating and ForwardDiff-clean in both arguments.
    g = Octofitter.ForwardDiff.derivative(
        p -> first(Octofitter.sky_offset(sol, t, r; platescale=p, northangle=na)), ps)
    @test isfinite(g) && !iszero(g)
    n = length(EPOCHS_A)
    Δα, Δδ = zeros(n), zeros(n)
    h(Δα, Δδ, ctx, obs) = Octofitter.sky_offset!(Δα, Δδ, ctx, obs.target, obs.ref;
                                                 platescale=1.07, northangle=0.13)
    h(Δα, Δδ, ctx, obs)
    @test (@allocated h(Δα, Δδ, ctx, obs)) == 0
end

@testset "sky_calibration defaults to the identity" begin
    plain = _astrom_obs()
    ctx = _ctx_for(model_system(obs=[plain]), plain)
    @test Octofitter.sky_calibration(ctx) === (1.0, 0.0)

    cal = _astrom_obs(name="cal", variables=@variables begin
        platescale = 1.02
        northangle = 0.05
    end)
    ctxc = _ctx_for(model_system(obs=[cal]), cal)
    @test Octofitter.sky_calibration(ctxc) === (1.02, 0.05)
end

@testset "show describes targets, then the reference" begin
    obs = _astrom_obs()
    s = sprint(show, MIME"text/plain"(), obs)
    @test occursin("c vs b", s)
    @test Octofitter._refdesc(obs) == "c vs b"

    # The variadic reading: everything but the last spec is a target. v1's
    # flat `join(…, " vs ")` printed `A vs Ab vs Barycentre`, which reads as a
    # chain of differences.
    @test Octofitter._refdesc((Octofitter.refspec(:b), Octofitter.refspec(:c)),
                              Octofitter.refspec(:A)) == "(b, c) vs A"
    @test Octofitter._refdesc_default((Octofitter.refspec(:b),
                                       Octofitter.refspec(:c), Barycentre)) ==
          "(b, c) vs Barycentre"
    @test Octofitter._refdesc_default(()) == ""
end
