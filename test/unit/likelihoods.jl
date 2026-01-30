@testset "Likelihood Objects" begin
    @testset "PlanetRelAstromLikelihood" begin
        # Test RA/Dec format
        data_radec = PlanetRelAstromLikelihood(
            Table(epoch=[5000.0, 5100.0], ra=[100.0, 110.0], dec=[50.0, 55.0], σ_ra=[1.0, 1.0], σ_dec=[1.0, 1.0]),
            name="test_radec"
        )
        @test data_radec isa PlanetRelAstromLikelihood
        @test length(data_radec.table) == 2
        @test hasproperty(data_radec.table, :ra)

        # Test sep/PA format
        data_seppa = PlanetRelAstromLikelihood(
            Table(epoch=[5000.0, 5100.0], sep=[100.0, 110.0], pa=[1.0, 1.1], σ_sep=[1.0, 1.0], σ_pa=[0.1, 0.1]),
            name="test_seppa"
        )
        @test data_seppa isa PlanetRelAstromLikelihood
        @test length(data_seppa.table) == 2
        @test hasproperty(data_seppa.table, :sep)

        # Test that invalid column combinations throw errors
        @test_throws Exception PlanetRelAstromLikelihood(
            Table(epoch=[5000.0], ra=[100.0], pa=[1.0], σ_ra=[1.0], σ_pa=[0.1]),
            name="test_invalid"
        )

        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(data_seppa, 1:1)
        @test length(subset.table) == 1
    end

    @testset "PhotometryLikelihood" begin
        phot = PhotometryLikelihood(
            Table(band=[:Z, :J], phot=[15.0, 14.0], σ_phot=[0.1, 0.2]),
            name="test_phot"
        )
        @test phot isa PhotometryLikelihood
        @test length(phot.table) == 2
        @test all([:band, :phot, :σ_phot] .∈ Ref(propertynames(phot.table)))

        subset = Octofitter.likeobj_from_epoch_subset(phot, 1:1)
        @test length(subset.table) == 1
    end

    @testset "HGCALikelihood" begin
        gaia_id = 756291174721509376
        hgca = HGCALikelihood(;gaia_id=gaia_id)
        @test hgca isa HGCALikelihood
    end
end
