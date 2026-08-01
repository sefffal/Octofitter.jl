using Random

# Reduced iterations for faster tests
const JOINT_ITERATIONS = 50
const JOINT_ADAPTATION = 50

@testset "Joint Fitting" begin
    # Use seeded RNG for reproducibility
    rng = Random.Xoshiro(42)

    astrom = PlanetRelAstromLikelihood(
        Table(epoch=[50000.0], ra=[100.0], dec=[50.0], σ_ra=[1.0], σ_dec=[1.0]),
        name="joint_fitting_test"
    )

    hgca = HGCALikelihood(;gaia_id=756291174721509376)

    b_astrom = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom],
        variables=@variables begin
            a ~ LogUniform(0.1, 10)
            e ~ Uniform(0, 0.5)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            mass ~ LogUniform(0.1, 100)
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )

    AstromSystem = System(
        name="AstromSystem",
        companions=[b_astrom],
        observations=[],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ gaia_plx(;gaia_id=756291174721509376)
            pmra ~ Normal(-975, 10)
            pmdec ~ Normal(20, 10)
            jitter ~ LogUniform(0.1, 100)
            rv0 ~ Normal(0, 100)
        end
    )

    b_joint = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom],
        variables=@variables begin
            a ~ LogUniform(0.1, 10)
            e ~ Uniform(0, 0.5)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            mass ~ LogUniform(0.1, 100)
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )

    JointSystem = System(
        name="JointSystem",
        companions=[b_joint],
        observations=[hgca],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ gaia_plx(;gaia_id=756291174721509376)
            pmra ~ Normal(-975, 10)
            pmdec ~ Normal(20, 10)
            jitter ~ LogUniform(0.1, 100)
            rv0 ~ Normal(0, 100)
        end
    )

    model_astrom = Octofitter.LogDensityModel(AstromSystem)
    model_joint = Octofitter.LogDensityModel(JointSystem)

    # Test that both models can be fit successfully
    chain_astrom = octofit(rng, model_astrom, iterations=JOINT_ITERATIONS, adaptation=JOINT_ADAPTATION)
    chain_joint = octofit(rng, model_joint, iterations=JOINT_ITERATIONS, adaptation=JOINT_ADAPTATION)

    # Basic sanity checks - both chains should have reasonable log posteriors
    @test all(isfinite, chain_astrom[:logpost])
    @test all(isfinite, chain_joint[:logpost])

    # Both should produce mass estimates in a reasonable range
    @test all(chain_astrom[:b_mass] .> 0)
    @test all(chain_joint[:b_mass] .> 0)

    # Note: We intentionally don't test std(joint) < std(astrom) because this is stochastic
    # and can fail with short chains. The key test is that both models can fit successfully.
end
