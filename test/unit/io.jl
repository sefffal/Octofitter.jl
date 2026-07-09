using MCMCChains: Chains

@testset "Data I/O" begin
    @testset "FITS Chain basic I/O" begin
        # Test handling of unicode characters in parameter names
        test_chain = Chains(
            randn(100, 5, 1),  # 100 samples, 5 parameters, 1 chain
            [:a, :e, :ω, :Ω, :θ]  # Unicode parameter names
        )

        # Save and reload with FITS
        fname = tempname() * ".fits"
        Octofitter.savechain(fname, test_chain)
        loaded_chain = Octofitter.loadchain(fname)

        @test keys(test_chain) == keys(loaded_chain)
        for name in keys(test_chain)
            @test test_chain[name] ≈ loaded_chain[name]
        end
    end

    @testset "CSV Astrometry Data I/O" begin
        # Create a test CSV file
        test_data = """
        epoch,ra,dec,σ_ra,σ_dec,cor
        50000.0,100.0,50.0,1.0,1.0,0.0
        50100.0,110.0,55.0,1.0,1.0,0.1
        50200.0,120.0,60.0,1.0,1.0,-0.1
        """

        fname = tempname() * ".csv"
        open(fname, "w") do io
            write(io, test_data)
        end

        # Load using CSV.read
        csv_table = CSV.read(fname, Table)
        astrom_data = PlanetRelAstromLikelihood(
            csv_table,
            name="csv_test"
        )

        @test length(astrom_data.table) == 3
        @test astrom_data.table.epoch[1] ≈ 50000.0
        @test astrom_data.table.ra[2] ≈ 110.0
        @test astrom_data.table.dec[3] ≈ 60.0
        @test astrom_data.table.cor[2] ≈ 0.1
    end
end
