# Additional imports needed for this test file
# (base imports come from runtests.jl)
using Pigeons
using CairoMakie
using MCMCChains
using Octofitter: TypedTables

@testset "Gaia DR4 Integration" begin

    # Load real test data
    test_data_path = joinpath(@__DIR__, "..", "data", "gaia_dr4_test.csv")
    df = CSV.read(test_data_path, DataFrame)

    # Create model for all tests
    ref_epoch_mjd = 57936.375

    priors, derived = @variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)
        ra_offset_mas ~ Normal(0, 10000)
        dec_offset_mas ~ Normal(0, 10000)
        pmra ~ Uniform(-1000, 1000)
        pmdec ~ Uniform(-1000, 1000)
        plx = system.plx
        ref_epoch = $ref_epoch_mjd
    end

    table = TypedTables.Table(df)
    if hasproperty(table, :obs_time_tcb) && !hasproperty(table, :epoch)
        table = TypedTables.Table(table; epoch=Octofitter.jd2mjd.(table.obs_time_tcb))
    end
    xyz = TypedTables.Table(Octofitter.geocentre_position_query.(table.epoch))
    table = TypedTables.Table(table; xyz)

    mock_gaia_sol = (
        ra = 180.0,
        dec = 45.0,
        pmra = 5.0,
        pmdec = -10.0,
        parallax = 20.0,
    )

    gaia_obs = Octofitter.GaiaDR4AstromObs{typeof(table), typeof(mock_gaia_sol)}(
        table,
        123456789,
        mock_gaia_sol,
        priors,
        derived,
        "GaiaDR4"
    )

    orbit_ref_epoch = mean(gaia_obs.table.epoch)

    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[],
        variables=@variables begin
            a ~ LogUniform(0.1, 50)
            e ~ Uniform(0, 0.9)
            ω ~ Uniform(0, 2π)
            i ~ Sine()
            Ω ~ Uniform(0, 2π)
            θ ~ Uniform(0, 2π)
            tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
            mass ~ LogUniform(0.1, 500)
        end
    )

    sys = System(
        name="gaia_dr4_integration_test",
        companions=[b],
        observations=[gaia_obs],
        variables=@variables begin
            M = 1.0
            plx ~ Uniform(1, 100)
        end
    )

    model = Octofitter.LogDensityModel(sys, verbosity=0)

    @testset "Initialization and short MCMC" begin
        # Test initialization
        init_chain = initialize!(model)
        @test init_chain isa MCMCChains.Chains
        @test size(init_chain, 1) > 0

        # Run a very short MCMC chain
        chain, pt = octofit_pigeons(model, n_rounds=4)
        @test chain isa MCMCChains.Chains
        @test size(chain, 1) > 0
    end

    @testset "gaiastarplot" begin
        init_chain = initialize!(model)

        # Test gaiastarplot with default arguments
        fig = Octofitter.gaiastarplot(model, init_chain)
        @test fig isa Makie.Figure

        # Test gaiastarplot with explicit sample index
        fig2 = Octofitter.gaiastarplot(model, init_chain, 1)
        @test fig2 isa Makie.Figure

        # Clean up generated files
        rm("gaia_dr4_integration_test-gaiastarplot.png", force=true)
    end

    @testset "skytrackplot" begin
        init_chain = initialize!(model)

        # Test skytrackplot with default arguments
        fig = Octofitter.skytrackplot(model, init_chain, 1)
        @test fig isa Makie.Figure

        # Test skytrackplot with keplerian_mult
        fig2 = Octofitter.skytrackplot(model, init_chain, 1, keplerian_mult=5.0)
        @test fig2 isa Makie.Figure

        # Clean up generated files
        rm("gaia_dr4_integration_test-skytrack.png", force=true)
    end

    @testset "gaiastarplot with MCMC chain" begin
        # Run short MCMC
        chain, pt = octofit_pigeons(model, n_rounds=4)

        # Test with MCMC chain (uses MAP by default)
        fig = Octofitter.gaiastarplot(model, chain)
        @test fig isa Makie.Figure

        # Test with specific sample from chain
        idx = rand(1:size(chain, 1))
        fig2 = Octofitter.gaiastarplot(model, chain, idx)
        @test fig2 isa Makie.Figure

        # Clean up
        rm("gaia_dr4_integration_test-gaiastarplot.png", force=true)
    end

    @testset "skytrackplot with MCMC chain" begin
        chain, pt = octofit_pigeons(model, n_rounds=4)

        idx = rand(1:size(chain, 1))
        fig = Octofitter.skytrackplot(model, chain, idx)
        @test fig isa Makie.Figure

        # Clean up
        rm("gaia_dr4_integration_test-skytrack.png", force=true)
    end

    @testset "mcmcchain2result with Gaia DR4" begin
        chain, pt = octofit_pigeons(model, n_rounds=4)

        # Test single sample extraction
        sample_idx = 1
        result = Octofitter.mcmcchain2result(model, chain, sample_idx)
        @test result isa NamedTuple
        @test haskey(result, :plx)
        @test haskey(result, :planets)
        @test haskey(result.planets, :b)
        @test haskey(result, :observations)
        @test haskey(result.observations, :GaiaDR4)
    end
end
