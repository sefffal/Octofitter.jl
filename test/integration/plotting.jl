using Random
using PairPlots
using CairoMakie

# Reduced iterations for faster tests
const PLOT_TEST_ITERATIONS = 50
const PLOT_TEST_ADAPTATION = 50

@testset "Plotting" begin
    astrom_like = PlanetRelAstromLikelihood(
        Table(
            epoch = [50000, 50120, 50240, 50360, 50480, 50600, 50720, 50840],
            ra = [-494.4, -495.0, -493.7, -490.4, -485.2, -478.1, -469.1, -458.3],
            dec = [-76.7, -44.9, -12.9, 19.1, 51.0, 82.8, 114.3, 145.3],
            σ_ra = [12.6, 10.4, 9.9, 8.7, 8.0, 6.9, 5.8, 4.2],
            σ_dec = [12.6, 10.4, 9.9, 8.7, 8.0, 6.9, 5.8, 4.2],
            cor = [0.2, 0.5, 0.1, -0.8, 0.3, -0.0, 0.1, -0.2]
        ),
        name="plotting_test"
    )

    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[astrom_like],
        variables=@variables begin
            a ~ LogUniform(0.1, 100)
            e ~ Uniform(0, 0.99)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
            mass = 10
        end
    )

    gaia_id = 756291174721509376
    hgca = HGCALikelihood(;gaia_id=gaia_id)

    TestSystem = System(
        name="TestSystem",
        companions=[b],
        observations=[hgca],
        variables=@variables begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            pmra ~ Normal(-137, 10)
            pmdec ~ Normal(2, 10)
        end
    )

    model = Octofitter.LogDensityModel(TestSystem)
    Random.seed!(0)
    Octofitter.default_initializer!(model, nruns=1)
    chain = octofit(model, iterations=PLOT_TEST_ITERATIONS, adaptation=PLOT_TEST_ADAPTATION)

    @testset "octoplot" begin
        fig = octoplot(model, chain)
        @test fig isa Makie.Figure

        # Test individual plot options
        fig_astrom = octoplot(
            model, chain,
            N = 1,
            show_astrom=true,
            show_physical_orbit=false,
            show_astrom_time=false,
            show_pma=false,
            show_mass=false,
            show_rv=false,
            show_relative_rv=false,
            show_hipparcos=false,
        )
        @test fig_astrom isa Makie.Figure

        fig_phys_orb = octoplot(
            model, chain,
            N = 1,
            show_astrom=false,
            show_physical_orbit=true,
            show_astrom_time=false,
            show_pma=false,
            show_mass=false,
            show_rv=false,
            show_relative_rv=false,
            show_hipparcos=false,
        )
        @test fig_phys_orb isa Makie.Figure
    end

    @testset "octocorner" begin
        fig_corner = octocorner(model, chain, small=true)
        @test fig_corner isa Makie.Figure
    end

    # Clean up generated files
    rm("$(model.system.name)-plot-grid.png", force=true)
end
