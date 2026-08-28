# The Pigeons extension.
#
# `octofit_pigeons` is exported from core with no methods; they come from
# `ext/OctofitterPigeonsExt.jl`, which loads when Pigeons does. That extension
# was left unported through most of the v2 work and nothing noticed until the
# documentation authors tried to run the tutorials — so this file exists mainly
# to make "the function has no methods" a test failure.
#
# The runs here are deliberately tiny (2^4 scans, 4 chains, no variational
# reference). What is being checked is that the plumbing between Octofitter's
# `LogDensityModel` and Pigeons' `initialization` / `sample_iid!` /
# `default_reference` / `Chains` hooks still lines up, not that tempering
# converges.

using Test
using Octofitter
using Octofitter: Octofitter, Body, System
using Octofitter.TypedTables: Table
using Distributions
using PlanetOrbits
using Random
using Pigeons

const PG_ASTROM = Table(
    epoch = [56000.0, 56400.0, 56800.0, 57200.0, 57600.0, 58000.0],
    ra    = [100.0, 115.0, 120.0, 118.0, 105.0, 85.0],
    dec   = [50.0, 35.0, 15.0, -5.0, -25.0, -40.0],
    σ_ra  = fill(5.0, 6), σ_dec = fill(5.0, 6),
)

function pg_system()
    A = Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.1)
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass = 0.0
        a ~ Uniform(1.0, 20.0)
        e ~ Uniform(0.0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 56000.0
    end)
    obs = RelAstromObs(PG_ASTROM; target=b, ref=A, name="astrom")
    return System(name="pigtest", bodies=[A, b], observations=[obs],
        variables=@variables begin
            plx ~ truncated(Normal(25.0, 0.1), lower=1.0)
        end)
end

@testset "octofit_pigeons has methods once Pigeons is loaded" begin
    # The regression, stated plainly.
    @test !isempty(methods(octofit_pigeons))
    @test any(m -> m.module === Base.get_extension(Octofitter, :OctofitterPigeonsExt),
              methods(octofit_pigeons))
end

@testset "a tempered fit runs and returns a usable chain" begin
    sys = pg_system()
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    chain, pt = octofit_pigeons(model; n_rounds=4, n_chains=4,
        n_chains_variational=0, variational=nothing)

    @test chain isa Octofitter.MCMCChains.Chains
    @test size(chain, 1) == 2^4
    # Model parameters and the sampler's own diagnostics both come through, and
    # the diagnostics are marked internal (v1 dropped that distinction on
    # reload; the Pigeons path sets it explicitly).
    for k in (:plx, :A_mass, :b_a, :b_e, :b_i)
        @test k in keys(chain)
    end
    @test :pigeons_logpotential in keys(chain)
    @test :loglike in keys(chain)
    @test all(isfinite, chain[:loglike])
    @test all(isfinite, chain[:logprior])
    # `Chains(model, pt)` recomputes the likelihood rather than trusting the
    # tempered log potential, so logpost == loglike + logprior exactly.
    @test collect(chain[:logpost]) ≈ collect(chain[:loglike]) .+ collect(chain[:logprior])

    @test isfinite(Pigeons.stepping_stone(pt))
    @test chain.info.sampler == "pigeons"
end

@testset "cores= runs the same fit in worker processes" begin
    # The `cores` path serializes the whole `LogDensityModel` into freshly
    # launched worker processes (Pigeons `ChildProcess`, local MPI) and reads
    # the chain back from the round checkpoints. What is being checked is the
    # full round trip — serialize, launch, deserialize (RuntimeGeneratedFunctions
    # included), sample, load, convert — on every OS in CI, notably Windows,
    # where the MPI runtime is a different vendor (MicrosoftMPI_jll).
    sys = pg_system()
    model = Octofitter.LogDensityModel(sys; verbosity=0)

    @test_throws ArgumentError octofit_pigeons(model; n_rounds=3, cores=0)
    @test_throws ArgumentError octofit_pigeons(model; n_rounds=3, cores=2, threads_per_process=0)
    @test_throws ArgumentError octofit_pigeons(model; n_rounds=3, cores=2, threads_per_process=4)

    # ChildProcess checkpoints under `results/` in the working directory;
    # run from a scratch directory so the repo stays clean.
    chain = cd(mktempdir()) do
        chain, pt = octofit_pigeons(model; n_rounds=4, n_chains=4,
            n_chains_variational=0, variational=nothing, cores=2)
        chain
    end
    @test chain isa Octofitter.MCMCChains.Chains
    @test size(chain, 1) == 2^4
    for k in (:plx, :A_mass, :b_a, :b_e, :b_i)
        @test k in keys(chain)
    end
    @test all(isfinite, chain[:loglike])
    @test chain.info.sampler == "pigeons"
end

@testset "the reference model is sampled IID when it can be" begin
    sys = pg_system()

    # `prior_only_model` replaces each observation with a `BlankLikelihood`,
    # which v2 marks `_isprior = true` even though it contributes exactly zero.
    # Counting it would send every reference down the explorer path instead of
    # drawing from the priors — v1 got this right by accident, having never
    # defined `_isprior` for `BlankLikelihood` at all.
    ext = Base.get_extension(Octofitter, :OctofitterPigeonsExt)
    @test !isnothing(ext)

    bare = Octofitter.LogDensityModel(
        Octofitter.prior_only_model(sys, exclude_all=true); verbosity=0)
    @test !ext._has_non_sampleable_priors(bare)

    # …while the default `prior_only_model` keeps the `UnitLengthPrior` behind
    # each `UniformCircular`, which genuinely reshapes the prior.
    kept = Octofitter.LogDensityModel(Octofitter.prior_only_model(sys); verbosity=0)
    @test ext._has_non_sampleable_priors(kept)

    # And the full model, which also carries a real likelihood.
    @test ext._has_non_sampleable_priors(Octofitter.LogDensityModel(sys; verbosity=0))

    # The documented log-evidence recipe: `stepping_stone` on the bare
    # prior-only model is the reference normalization `log_Z0`.
    _, pt0 = octofit_pigeons(bare; n_rounds=3, n_chains=3,
        n_chains_variational=0, variational=nothing)
    @test isfinite(Pigeons.stepping_stone(pt0))
end
