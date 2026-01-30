using Random

# Reduced iterations for faster tests
const CV_ITERATIONS = 50
const CV_ADAPTATION = 50

@testset "Cross Validation" begin
    epochs = collect(range(50000, 51000, length=10))
    astrom1 = PlanetRelAstromLikelihood(
        Table(
            epoch=epochs,
            ra=[100.0+t/100 for t in epochs],
            dec=[50.0+t/200 for t in epochs],
            σ_ra=fill(1.0, length(epochs)),
            σ_dec=fill(1.0, length(epochs))
        ),
        name="astrom1"
    )
    astrom2 = PlanetRelAstromLikelihood(
        Table(
            epoch=epochs,
            ra=[100.0+t/100 for t in epochs],
            dec=[50.0+t/200 for t in epochs],
            σ_ra=fill(1.0, length(epochs)),
            σ_dec=fill(1.0, length(epochs))
        ),
        name="astrom2"
    )

    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom1, astrom2],
        variables=@variables begin
            a  = 1.0
            e  = 0.0
            i  = 1.0
            ω  = 1.0
            Ω  = 1.0
            tp = 50000
        end
    )

    Sys = System(
        name="Sys",
        companions=[b],
        observations=[],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx = 100.
        end
    )

    model = Octofitter.LogDensityModel(Sys)
    chain = octofit(model, iterations=CV_ITERATIONS, adaptation=CV_ADAPTATION)

    @testset "Pointwise likelihood" begin
        like_mat, cv_epochs = Octofitter.pointwise_like(model, chain)
        @test size(like_mat, 2) == 20  # One column per epoch
        @test length(cv_epochs) == 20
        @test all(isfinite, like_mat)
    end

    @testset "K-fold systems generation" begin
        kfold_systems = Octofitter.generate_kfold_systems(model.system)
        @test length(kfold_systems) == 2  # One per dataset
    end

    @testset "Per-epoch system generation" begin
        per_epoch_systems, cv_epochs = Octofitter.generate_system_per_epoch(model.system)
        @test length(per_epoch_systems) == 20
        @test length(cv_epochs) == 20
    end
end
