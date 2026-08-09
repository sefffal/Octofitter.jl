# Performance and differentiability gates on the hot path.

function gradient_testmodel()
    A = Octofitter.Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.1), lower=0.2)
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass_jup ~ LogUniform(0.1, 100)
        mass = mass_jup * mjup
        a ~ LogUniform(1, 30)
        e ~ Uniform(0, 0.9)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        # Parametrize tp as an offset from a base epoch: a raw MJD ~5.7e4
        # gives FiniteDiff a relative step of ~0.4 d against a short period,
        # which looks like an AD bug and is not one.
        dtp ~ Uniform(-500, 500)
        tp = 56000.0 + dtp
    end)
    epochs = collect(range(56000.0, 59000.0, length=14))
    astrom = RelAstromObs((epoch=epochs,
            ra=100 .* sin.(epochs ./ 500), dec=100 .* cos.(epochs ./ 500),
            σ_ra=fill(3.0, length(epochs)), σ_dec=fill(3.0, length(epochs)));
        target=b, ref=A, name="astrom")
    rv = RadialVelocityObs((epoch=epochs, rv=10 .* sin.(epochs ./ 700),
            σ_rv=fill(2.0, length(epochs)));
        target=A, ref=Barycentre, name="rv",
        variables=@variables begin
            offset ~ Normal(0, 50)
            jitter ~ LogUniform(0.5, 20)
        end)
    return Octofitter.System(name="grad", bodies=[A, b], observations=[astrom, rv],
        variables=@variables begin plx ~ Uniform(5, 60) end)
end

@testset "ForwardDiff agrees with finite differences" begin
    sys = gradient_testmodel()
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θt = model.link(Octofitter.sample_priors(Random.Xoshiro(3), sys))
    _, g = model.∇ℓπcallback(θt)
    gfd = FiniteDiff.finite_difference_gradient(model.ℓπcallback, θt)
    # Compare against the gradient *norm*: individual components can be
    # near-zero, where a relative comparison means nothing.
    @test norm(g .- gfd) / norm(g) < 1e-5
end

@testset "the hot path allocates nothing" begin
    sys = gradient_testmodel()
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θ = Octofitter.sample_priors(Random.Xoshiro(4), sys)
    θt = model.link(θ)
    lnlike = Octofitter.make_ln_like(sys)
    nt = model.arr2nt(θ)

    model.arr2nt(θ); lnlike(sys, nt); model.ℓπcallback(θt)
    # `broken` below 1.12, measured rather than guessed. The escape analysis
    # that folds these away landed in 1.12, not 1.11: on 1.11 *every* absolute
    # allocation gate in the suite still allocates, and by the same amounts
    # 1.10 did — 752 B here, 352 B at the subset model, 128 B for `arr2nt`,
    # 64/48/32/16 B for the sky-offset and prior gates. Marked rather than
    # skipped, so the summary keeps reporting them and an unexpected pass
    # fails loudly if a 1.11 patch starts folding them.
    @test (@allocated model.arr2nt(θ)) == 0 broken=(VERSION < v"1.12")
    @test (@allocated lnlike(sys, nt)) == 0 broken=(VERSION < v"1.12")
    @test (@allocated model.ℓπcallback(θt)) == 0 broken=(VERSION < v"1.12")
end

@testset "type stability through the whole path" begin
    sys = gradient_testmodel()
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θ = Octofitter.sample_priors(Random.Xoshiro(5), sys)
    nt = model.arr2nt(θ)
    lnlike = Octofitter.make_ln_like(sys)
    @test isconcretetype(Core.Compiler.return_type(model.arr2nt, Tuple{typeof(θ)}))
    @test Core.Compiler.return_type(lnlike, Tuple{typeof(sys),typeof(nt)}) === Float64
    @test Core.Compiler.return_type(model.ℓπcallback, Tuple{Vector{Float64}}) === Float64
    # The PlanetOrbits system built per sample must be concretely typed, or
    # every observable read downstream pays for a dynamic dispatch.
    @test isconcretetype(typeof(lnlike.build(nt)))
end

@testset "one solve, not one per body" begin
    # The trajectory covers every epoch in the model exactly once, so adding a
    # companion adds columns rather than another pass over the epochs.
    sys = gradient_testmodel()
    ep, maps = Octofitter.epoch_plan(sys)
    @test length(ep) == 14                      # both observations share epochs
    @test all(m -> m == 1:14, values(maps))
end

@testset "the scratch arena is sized to the model" begin
    # Bumper's default buffer keeps only its first 1 MB slab across
    # `@no_escape` boundaries and frees the rest, so a model whose per-sample
    # scratch overflows one slab re-faults the excess every single evaluation
    # (measured at 1.8-2.4x on the gradient past ~1000 epochs). `_slab_size`
    # prices the trajectory at model-build time so the arena covers it in one
    # slab. See `_slab_size` in src/model/codegen.jl.

    # A small model must not pay for an arena it does not need: `nothing`
    # means "keep Bumper's default buffer", so no second arena is held.
    small = gradient_testmodel()
    @test Octofitter._slab_size(small, first(Octofitter.epoch_plan(small))) === nothing

    # A model with enough epochs to overflow the default slab gets its own.
    function bigmodel(nep)
        A = Octofitter.Body(name="A", variables=@variables begin
            mass ~ Uniform(0.5, 2.0)
        end)
        b = Octofitter.Body(name="b", about=A, variables=@variables begin
            mass = 0.001
            a ~ LogUniform(0.1, 100)
            e ~ Uniform(0, 0.9)
            ω ~ Uniform(0, 2pi)
            i ~ Sine()
            Ω ~ Uniform(0, 2pi)
            θ ~ Uniform(0, 2pi)
            epoch = 57388.0
        end)
        ep = collect(range(55000.0, 60000.0, length=nep))
        rv = RadialVelocityObs((epoch=ep, rv=10 .* sin.(ep ./ 700), σ_rv=fill(2.0, nep));
            target=A, ref=Barycentre, name="rv",
            variables=@variables begin
                jitter ~ LogUniform(0.5, 20)
            end)
        return Octofitter.System(name="big", bodies=[A, b], observations=[rv],
            variables=@variables begin plx ~ Uniform(5, 60) end)
    end

    big = bigmodel(2000)
    slab = Octofitter._slab_size(big, first(Octofitter.epoch_plan(big)))
    @test slab isa Int
    D = length(Octofitter._list_priors(big))
    # It has to actually cover one Dual{D} trajectory — being a hair too small
    # costs the whole penalty, which is the failure mode this guards.
    posys = Octofitter._example_posys(big)
    need = PlanetOrbits.trajectory_storage(
        Octofitter.ForwardDiff.Dual{Nothing,Float64,D}, posys, first(Octofitter.epoch_plan(big)))
    @test slab > need
    @test ispow2(slab)

    # And the buffer it hands back is that size, task-locally, and stable
    # across calls (a fresh buffer per evaluation would defeat the point).
    buf1 = Octofitter._scratch_buffer(Val(slab))
    buf2 = Octofitter._scratch_buffer(Val(slab))
    @test buf1 === buf2
    @test buf1 isa Octofitter.Bumper.SlabBuffer{slab}
    @test Octofitter._scratch_buffer(Val(nothing)) === Octofitter.Bumper.default_buffer()

    # The sizing is an optimization and must never change the answer.
    m = Octofitter.LogDensityModel(big; verbosity=0)
    θt = m.link(Octofitter.sample_priors(Random.Xoshiro(11), big))
    lp = m.ℓπcallback(θt)
    @test isfinite(lp)
    @test (@allocated m.ℓπcallback(θt)) == 0 broken=(VERSION < v"1.12")
    v, g = m.∇ℓπcallback(θt)
    @test v ≈ lp
    @test all(isfinite, g)
end

@testset "ForwardDiff chunk width" begin
    # The rule, checked against the widths that actually measured fastest on
    # the N-companion models (see the table above `_chunk_heuristic` in
    # src/logdensitymodel.jl). Both sets of numbers are hardware-specific, so
    # each is asserted against the constants that produced it rather than
    # against whatever host happens to run CI.
    let single = Octofitter.CHUNK_SINGLE_MAX[], width = Octofitter.CHUNK_WIDTH_MAX[]
        try
            # AVX2 measurements: single==width==24.
            Octofitter.CHUNK_SINGLE_MAX[] = 24
            Octofitter.CHUNK_WIDTH_MAX[]  = 24
            @test Octofitter._chunk_heuristic(23) == 23   # 23 beats 12/16/20
            @test Octofitter._chunk_heuristic(30) == 15   # 15 beats 10/12/18/20/24/30
            @test Octofitter._chunk_heuristic(37) == 19   # 19 beats 13/16/22/25/30/37

            # AVX-512 (Zen 5) measurements: the single chunk keeps winning to
            # 59, then falls off a 5x cliff at 60, after which narrow splits
            # win. See the second table in src/logdensitymodel.jl.
            Octofitter.CHUNK_SINGLE_MAX[] = 59
            Octofitter.CHUNK_WIDTH_MAX[]  = 24
            @test Octofitter._chunk_heuristic(27) == 27   # 27 beats 14/16/20/24
            @test Octofitter._chunk_heuristic(37) == 37   # 37 beats 16/19/24/32
            @test Octofitter._chunk_heuristic(47) == 47   # 47 beats 16/24/32
            @test Octofitter._chunk_heuristic(57) == 57   # 57 beats 20/29/40
            @test Octofitter._chunk_heuristic(59) == 59   # last width before the cliff
            @test Octofitter._chunk_heuristic(60) != 60   # 60 is past it (5.1x)
            @test Octofitter._chunk_heuristic(67) == 23   # 23 beats 34/45/67
        finally
            Octofitter.CHUNK_SINGLE_MAX[] = single
            Octofitter.CHUNK_WIDTH_MAX[]  = width
        end
    end

    # Single chunk right up to the bound, never past it, and always balanced.
    for D in 1:200
        c = Octofitter._chunk_heuristic(D)
        @test 1 <= c <= D
        @test c <= max(D, Octofitter.CHUNK_SINGLE_MAX[], Octofitter.CHUNK_WIDTH_MAX[])
        n = cld(D, c)
        if n > 1
            # Split: balanced, and dropping to one fewer chunk would exceed the
            # width cap. (A single chunk is exempt -- that branch is governed by
            # CHUNK_SINGLE_MAX, which may exceed CHUNK_WIDTH_MAX.)
            @test c <= Octofitter.CHUNK_WIDTH_MAX[]
            @test cld(D, n - 1) > Octofitter.CHUNK_WIDTH_MAX[]
        end
    end
    @test Octofitter._chunk_heuristic(Octofitter.CHUNK_SINGLE_MAX[]) == Octofitter.CHUNK_SINGLE_MAX[]

    # The measured probe is opt-in and must agree with the default's answer.
    sys = gradient_testmodel()
    D = length(Octofitter._list_priors(sys))
    m_default = Octofitter.LogDensityModel(sys; verbosity=0)
    m_probed = Octofitter.LogDensityModel(sys; verbosity=0, chunk_sizes=[2, 4, D])
    θt = m_default.link(Octofitter.sample_priors(Random.Xoshiro(12), sys))
    v1, g1 = m_default.∇ℓπcallback(θt)
    v2, g2 = m_probed.∇ℓπcallback(θt)
    @test v1 ≈ v2
    @test g1 ≈ g2

    # An explicit `autodiff` still wins over both.
    m_explicit = Octofitter.LogDensityModel(sys; verbosity=0,
        autodiff=Octofitter.AutoForwardDiff(chunksize=3))
    v3, g3 = m_explicit.∇ℓπcallback(θt)
    @test v3 ≈ v1
    @test g3 ≈ g1

    # Nonsense candidates fall back rather than failing model construction.
    m_bad = (@test_logs (:warn,) match_mode=:any Octofitter.LogDensityModel(
        sys; verbosity=0, chunk_sizes=[0, D + 50]))
    @test m_bad.∇ℓπcallback(θt)[1] ≈ v1
end

@testset "threaded trajectory solve is bit-identical to serial" begin
    # `orbitsolve!(...; threads=n)` splits the epoch union into contiguous
    # chunks solved on concurrent tasks. Every pass is epoch-local and each
    # chunk writes a disjoint range of the same storage, so this is not an
    # approximately-equal comparison: the likelihood and every gradient
    # component must come out bit-for-bit identical, or the chunking touched
    # something it must not. Meaningful on the multithreaded CI job; with one
    # thread the solve runs serial and the comparison is trivially true.
    A = Octofitter.Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.1), lower=0.2)
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 0.0
        a ~ LogUniform(1, 30)
        e ~ Uniform(0, 0.9)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        dtp ~ Uniform(-500, 500)
        tp = 56000.0 + dtp
    end)
    # Enough epochs that the solve actually chunks (2 × MIN_EPOCHS_PER_TASK).
    epochs = collect(range(56000.0, 59000.0, length=2*PlanetOrbits.MIN_EPOCHS_PER_TASK + 100))
    rv = RadialVelocityObs((epoch=epochs, rv=10 .* sin.(epochs ./ 700),
            σ_rv=fill(2.0, length(epochs)));
        target=A, ref=Barycentre, name="rv",
        variables=@variables begin
            offset ~ Normal(0, 50)
            jitter ~ LogUniform(0.5, 20)
        end)
    sys = Octofitter.System(name="threadgate", bodies=[A, b], observations=[rv],
        variables=@variables begin plx ~ Uniform(5, 60) end)
    model = Octofitter.LogDensityModel(sys; verbosity=0)

    was = Octofitter._kepsolve_use_threads[]
    try
        rng = Random.Xoshiro(4)
        for _ in 1:5
            θt = model.link(Octofitter.sample_priors(rng, sys))
            Octofitter._kepsolve_use_threads[] = false
            lp_s = model.ℓπcallback(θt)
            v_s, g_s = model.∇ℓπcallback(θt)
            g_s = copy(g_s)
            Octofitter._kepsolve_use_threads[] = true
            lp_t = model.ℓπcallback(θt)
            v_t, g_t = model.∇ℓπcallback(θt)
            @test lp_s === lp_t
            @test v_s === v_t
            @test all(g_s .=== g_t)
        end
    finally
        Octofitter._kepsolve_use_threads[] = was
    end
end
