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
    @test (@allocated model.arr2nt(θ)) == 0
    @test (@allocated lnlike(sys, nt)) == 0
    @test (@allocated model.ℓπcallback(θt)) == 0
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
