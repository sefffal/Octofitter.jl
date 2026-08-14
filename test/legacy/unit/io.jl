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

    @testset "Chain reconstruction preserves vector-valued derived variables" begin
        # Regression test for #115. `mcmcchain2result` reconstructs the
        # per-sample named tuple from a saved chain. `result2mcmcchain`
        # flattens and saves AbstractArray-valued derived variables (e.g.
        # G23HObs's `transits_dr2 :: Vector{Int}`), but the reconstruction
        # side previously only mapped Number/Tuple values back, silently
        # dropping vectors. That made e.g. `octoplot` crash on a reloaded /
        # sampled chain (`G23HObs requires a transits_dr2 observation
        # variable`). The two directions must stay symmetric.
        import Random

        planet = Planet(
            name="b",
            basis=Visual{KepOrbit},
            observations=Octofitter.AbstractLikelihood[],
            variables=@variables begin
                a ~ Uniform(1, 10)
                e ~ Uniform(0, 0.5)
                ω ~ Uniform(0, 2π)
                i ~ Sine()
                Ω ~ Uniform(0, 2π)
                θ ~ Uniform(0, 2π)
                tp = θ_at_epoch_to_tperi(θ, 57388.5; M=system.M, e, a, i, ω, Ω)
                mass ~ Uniform(1, 10)
                pl_vec = [a, 2a, 3a]          # planet-level vector-valued derived var
            end
        )
        sys = System(
            name="sys_rt",
            companions=[planet],
            observations=Octofitter.AbstractLikelihood[],
            variables=@variables begin
                M ~ Uniform(0.5, 2)
                plx ~ Uniform(10, 100)
                sys_vec = [M, plx]            # system-level vector-valued derived var
                sys_tup = (M, plx)            # tuple control: always round-tripped
            end
        )
        model = Octofitter.LogDensityModel(sys, verbosity=0)

        nts = [model.arr2nt(model.sample_priors(Random.default_rng())) for _ in 1:4]
        chain = Octofitter.result2mcmcchain(nts)
        recon = Octofitter.mcmcchain2result(model, chain)

        r, n = recon[1], nts[1]
        # System-level vector-valued derived var must survive reconstruction.
        @test hasproperty(r, :sys_vec)
        @test collect(Float64, r.sys_vec) ≈ collect(Float64, n.sys_vec)
        # Tuple control must keep working.
        @test hasproperty(r, :sys_tup)
        # Planet-level vector-valued derived var must survive too.
        @test hasproperty(r.planets.b, :pl_vec)
        @test collect(Float64, r.planets.b.pl_vec) ≈ collect(Float64, n.planets.b.pl_vec)
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
