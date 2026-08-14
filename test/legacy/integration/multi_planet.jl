using Random

# Reduced iterations for faster tests
const MULTI_PLANET_ITERATIONS = 50
const MULTI_PLANET_ADAPTATION = 50

@testset "Multi-Planet Systems" begin
    # Create test data for two planets
    astrom_b = PlanetRelAstromLikelihood(
        Table(epoch=[50000.0, 50100.0], ra=[100.0, 110.0], dec=[50.0, 55.0], σ_ra=[1.0, 1.0], σ_dec=[1.0, 1.0]),
        name="astrom_b"
    )
    astrom_c = PlanetRelAstromLikelihood(
        Table(epoch=[50000.0, 50100.0], ra=[-200.0, -210.0], dec=[-100.0, -110.0], σ_ra=[1.0, 1.0], σ_dec=[1.0, 1.0]),
        name="astrom_c"
    )

    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom_b],
        variables=@variables begin
            a ~ LogUniform(0.1, 10)
            e ~ Uniform(0, 0.5)
            i = system.i
            ω ~ UniformCircular()
            Ω = system.Ω
            P = 2*system.P_nominal * P_mul
            P_mul ~ Normal(1, 0.1)
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )

    c = Planet(
        name="c",
        basis=Visual{KepOrbit},
        observations=[astrom_c],
        variables=@variables begin
            a ~ LogUniform(0.1, 10)
            e ~ Uniform(0, 0.5)
            i = system.i
            ω ~ UniformCircular()
            Ω = system.Ω
            P = system.P_nominal * P_mul
            P_mul ~ truncated(Normal(1, 0.1), lower=0.1)
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )

    TwoPlanetSystem = System(
        name="TwoPlanetSystem",
        companions=[b, c],
        observations=[],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            i ~ Sine()
            Ω ~ UniformCircular()
            P_nominal ~ LogUniform(50, 300)
        end
    )

    model = Octofitter.LogDensityModel(TwoPlanetSystem)
    chain = octofit(model, iterations=MULTI_PLANET_ITERATIONS, adaptation=MULTI_PLANET_ADAPTATION)

    # Test that period ratio is approximately preserved (with tolerance for short chains)
    period_ratio = chain[:b_P][:] ./ chain[:c_P][:]
    @test mean(period_ratio) ≈ 2.0 rtol=0.5  # Relaxed tolerance for short chains

    # Test coplanarity is preserved
    @test all(chain[:b_i][:] .== chain[:c_i][:])
    @test all(chain[:b_Ω][:] .== chain[:c_Ω][:])
end
