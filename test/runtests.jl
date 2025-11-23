using Test
using Octofitter
using Distributions
using Random

# Set a fixed seed for reproducibility
Random.seed!(42)

@testset "Octofitter.jl" begin

    # ============================================================================
    # Test: Basic Relative Astrometry Fit (from rel-astrom.md)
    # ============================================================================
    @testset "Basic Relative Astrometry Fit" begin
        # Data from rel-astrom.md tutorial
        astrom_dat_1 = Table(;
            epoch= [50000,  50120, 50240, 50360, 50480, 50600, 50720, 50840,],
            ra   = [-505.764, -502.57, -498.209, -492.678, -485.977, -478.11, -469.08, -458.896,],
            dec  = [-66.9298, -37.4722, -7.92755, 21.6356, 51.1472,  80.5359,  109.729,  138.651,],
            σ_ra = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,],
            σ_dec = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,],
            cor =  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,]
        )
        astrom_like_1 = PlanetRelAstromLikelihood(astrom_dat_1, name="relastrom")

        # Sep/PA format data
        astrom_dat_2 = Table(
            epoch = [42000,],
            sep = [505.7637580573554,],
            pa = [deg2rad(24.1),],
            σ_sep = [70,],
            σ_pa = [deg2rad(10.2),],
        )
        astrom_like_2 = PlanetRelAstromLikelihood(astrom_dat_2, name="relastrom2")

        planet_1 = Planet(
            name="b",
            basis=Visual{KepOrbit},
            likelihoods=[astrom_like_1, astrom_like_2],
            variables=@variables begin
                plx = system.plx
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                a ~ Uniform(0, 100)
                e ~ Uniform(0.0, 0.5)
                i ~ Sine()
                ω ~ UniformCircular()
                Ω ~ UniformCircular()
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50420; M, e, a, i, ω, Ω)
            end
        )

        sys = System(
            name = "Tutoria",
            companions=[planet_1],
            likelihoods=[],
            variables=@variables begin
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )

        model = Octofitter.LogDensityModel(sys)
        @test model isa Octofitter.LogDensityModel

        # Sample from the posterior
        chain = octofit(model, iterations=100, adaptation=100)
        @test chain isa Octofitter.Chains
        @test size(chain, 1) > 0  # Has samples
    end

    # ============================================================================
    # Test: Proper Motion Anomaly / HGCA (from pma.md)
    # ============================================================================
    @testset "Proper Motion Anomaly (HGCA)" begin
        # Using HD 91312's Gaia ID
        hgca_like = HGCAInstantaneousLikelihood(gaia_id=756291174721509376, N_ave=1)

        planet_b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            variables=@variables begin
                a ~ LogUniform(0.1, 20)
                e ~ Uniform(0, 0.999)
                ω ~ UniformCircular()
                i ~ Sine()
                Ω ~ UniformCircular()

                mass = system.M_sec

                M = system.M
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 57423.0; M, e, a, i, ω, Ω)
            end
        )

        sys = System(
            name="HD91312_pma",
            companions=[planet_b],
            likelihoods=[hgca_like],
            variables=@variables begin
                M_pri ~ truncated(Normal(1.61, 0.1), lower=0.1)
                M_sec ~ LogUniform(0.5, 1000)
                M = M_pri + M_sec * Octofitter.mjup2msol

                plx ~ gaia_plx(gaia_id=756291174721509376)

                pmra ~ Normal(-137, 10)
                pmdec ~ Normal(2, 10)
            end
        )

        model_pma = Octofitter.LogDensityModel(sys)
        @test model_pma isa Octofitter.LogDensityModel

        # Sample from the posterior
        chain_pma = octofit(model_pma, iterations=100, adaptation=100)
        @test chain_pma isa Octofitter.Chains
        @test size(chain_pma, 1) > 0
    end

    # ============================================================================
    # Test: Hierarchical Co-Planar Model (from fit-coplanar.md)
    # ============================================================================
    @testset "Hierarchical Co-Planar Model" begin
        # HR8799 data from fit-coplanar.md
        astrom_dat_b = Table(;
            epoch = [53200.0, 54314.0, 54398.0, 54727.0, 55042.0, 55044.0, 55136.0, 55390.0, 55499.0, 55763.0, 56130.0, 56226.0, 56581.0, 56855.0, 58798.03906, 59453.245, 59454.231],
            ra    = [1471.0, 1504.0, 1500.0, 1516.0, 1526.0, 1531.0, 1524.0, 1532.0, 1535.0, 1541.0, 1545.0, 1549.0, 1545.0, 1560.0, 1611.002, 1622.924, 1622.872],
            dec   = [887.0, 837.0, 836.0, 818.0, 797.0, 794.0, 795.0, 783.0, 766.0, 762.0, 747.0, 743.0, 724.0, 725.0, 604.893, 570.534, 571.296],
            σ_ra  = [6.0, 3.0, 7.0, 4.0, 4.0, 7.0, 10.0, 5.0, 15.0, 5.0, 5.0, 4.0, 22.0, 13.0, 0.133, 0.32, 0.204],
            σ_dec = [6.0, 3.0, 7.0, 4.0, 4.0, 7.0, 10.0, 5.0, 15.0, 5.0, 5.0, 4.0, 22.0, 13.0, 0.199, 0.296, 0.446],
            cor   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.406, -0.905, -0.79]
        )

        astrom_dat_c = Table(;
            epoch = [53200.0, 54314.0, 54398.0, 54727.0, 55042.0, 55136.0, 55390.0, 55499.0, 55763.0, 56130.0, 56226.0, 56581.0, 56855.0],
            ra    = [-739.0, -683.0, -678.0, -663.0, -639.0, -636.0, -619.0, -607.0, -595.0, -578.0, -572.0, -542.0, -540.0],
            dec   = [612.0, 671.0, 678.0, 693.0, 712.0, 720.0, 728.0, 744.0, 747.0, 761.0, 768.0, 784.0, 799.0],
            σ_ra  = [6.0, 4.0, 7.0, 3.0, 4.0, 9.0, 4.0, 12.0, 4.0, 5.0, 3.0, 22.0, 12.0],
            σ_dec = [6.0, 4.0, 7.0, 3.0, 4.0, 9.0, 4.0, 12.0, 4.0, 5.0, 3.0, 22.0, 12.0],
            cor   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )

        astrom_b = PlanetRelAstromLikelihood(
            astrom_dat_b,
            name = "GPI_b",
            variables = @variables begin
                jitter = 0
                northangle = 0
                platescale = 1
            end
        )

        astrom_c = PlanetRelAstromLikelihood(
            astrom_dat_c,
            name = "GPI_c",
            variables = @variables begin
                jitter = 0
                northangle = 0
                platescale = 1
            end
        )

        # Exact co-planar model
        planet_b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            likelihoods=[astrom_b],
            variables=@variables begin
                e = 0.0
                ω = 0.0
                M_pri = system.M_pri
                M_b = system.M_b
                M_c = system.M_c
                M = M_pri + M_b * Octofitter.mjup2msol + M_c * Octofitter.mjup2msol
                mass = M_b

                i = system.i
                Ω = system.Ω

                P_mul ~ Normal(1, 0.1)
                P_nominal = system.P_nominal
                P = 2 * P_nominal * P_mul

                a = cbrt(M * P^2)
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 59454.231; M, e, a, i, ω, Ω)
            end
        )

        planet_c = Planet(
            name="c",
            basis=Visual{KepOrbit},
            likelihoods=[astrom_c],
            variables=@variables begin
                e = 0.0
                ω = 0.0
                M_pri = system.M_pri
                M_b = system.M_b
                M_c = system.M_c
                M = M_pri + M_b * Octofitter.mjup2msol
                mass = M_c

                i = system.i
                Ω = system.Ω

                P_mul ~ truncated(Normal(1, 0.1), lower=0.1)
                P_nominal = system.P_nominal
                P = P_nominal * P_mul

                a = cbrt(M * P^2)

                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 59454.231; M, e, a, i, ω, Ω)
            end
        )

        sys = System(
            name="HR8799_res_co",
            companions=[planet_b, planet_c],
            likelihoods=[],
            variables=@variables begin
                plx ~ gaia_plx(; gaia_id=2832463659640297472)
                M_pri ~ truncated(Normal(1.5, 0.02), lower=0.1)
                M_b ~ Uniform(0, 12)
                M_c ~ Uniform(0, 12)
                i ~ Sine()
                Ω ~ UniformCircular()
                P_nominal ~ Uniform(50, 300)
            end
        )

        model = Octofitter.LogDensityModel(sys)
        @test model isa Octofitter.LogDensityModel

        # Sample from the posterior
        results = octofit(model, iterations=100, adaptation=100)
        @test results isa Octofitter.Chains
        @test size(results, 1) > 0
    end

    # ============================================================================
    # Test: Hipparcos IAD Modeling (from hipparcos.md)
    # ============================================================================
    @testset "Hipparcos IAD Modeling" begin
        # Nielsen test - reproduce catalog values with zero mass planet
        hip_like = Octofitter.HipparcosIADLikelihood(
            hip_id=21547,
            renormalize=true,
            variables=@variables begin end
        )

        planet_b = Planet(
            name="b",
            basis=AbsoluteVisual{KepOrbit},
            variables=@variables begin
                mass = 0.0
                e = 0.0
                ω = 0.0
                a = 1.0
                i = 0.0
                Ω = 0.0
                tp = 0.0
            end
        )

        sys = System(
            name="c_Eri_straight_line",
            companions=[planet_b],
            likelihoods=[hip_like],
            variables=@variables begin
                M = 1.0
                rv = 0.0
                plx ~ Uniform(10, 100)
                pmra ~ Uniform(-100, 100)
                pmdec ~ Uniform(-100, 100)

                ra_hip_offset_mas ~ Normal(0, 10000)
                dec_hip_offset_mas ~ Normal(0, 10000)
                dec = $hip_like.hip_sol.dedeg + ra_hip_offset_mas / 60 / 60 / 1000
                ra = $hip_like.hip_sol.radeg + dec_hip_offset_mas / 60 / 60 / 1000 / cosd(dec)

                ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd
            end
        )

        model = Octofitter.LogDensityModel(sys)
        @test model isa Octofitter.LogDensityModel

        # Sample from the posterior
        chain = octofit(model, iterations=100, adaptation=100)
        @test chain isa Octofitter.Chains
        @test size(chain, 1) > 0
    end

    # ============================================================================
    # Test: Hipparcos with Planet Mass Constraint (from hipparcos.md)
    # ============================================================================
    @testset "Hipparcos Planet Mass Constraint" begin
        hip_like = Octofitter.HipparcosIADLikelihood(
            hip_id=21547,
            renormalize=true,
            variables=@variables begin end
        )

        astrom_dat = Table(;
            epoch = [57009.1, 57052.1, 57053.1, 57054.3, 57266.4, 57332.2, 57374.2, 57376.2, 57415.0, 57649.4, 57652.4, 57739.1, 58068.3, 58442.2],
            sep   = [454.24, 451.81, 456.8, 461.5, 455.1, 452.88, 455.91, 455.01, 454.46, 454.81, 451.43, 449.39, 447.54, 434.22],
            σ_sep = [1.88, 2.06, 2.57, 23.9, 2.23, 5.41, 6.23, 3.03, 6.03, 2.02, 2.67, 2.15, 3.02, 2.01],
            pa    = [2.98835, 2.96723, 2.97038, 2.97404, 2.91994, 2.89934, 2.89131, 2.89184, 2.8962, 2.82394, 2.82272, 2.79357, 2.70927, 2.61171],
            σ_pa  = [0.00401426, 0.00453786, 0.00523599, 0.0523599, 0.00453786, 0.00994838, 0.00994838, 0.00750492, 0.00890118, 0.00453786, 0.00541052, 0.00471239, 0.00680678, 0.00401426]
        )

        astrom_like1 = PlanetRelAstromLikelihood(
            astrom_dat,
            name="VLT_SPHERE",
            variables=@variables begin
                jitter = 0
                northangle = 0
                platescale = 1
            end
        )

        planet_b_mass = Planet(
            name="b",
            basis=AbsoluteVisual{KepOrbit},
            likelihoods=[astrom_like1],
            variables=@variables begin
                a ~ truncated(Normal(10, 1), lower=0.1)
                e ~ Uniform(0, 0.99)
                ω ~ Uniform(0, 2pi)
                i ~ Sine()
                Ω ~ Uniform(0, 2pi)
                θ ~ Uniform(0, 2pi)
                M = system.M
                tp = θ_at_epoch_to_tperi(θ, 58442.2; M, e, a, i, ω, Ω)
                mass = system.M_sec
            end
        )

        sys_mass = System(
            name="cEri",
            companions=[planet_b_mass],
            likelihoods=[hip_like],
            variables=@variables begin
                M_pri ~ truncated(Normal(1.75, 0.05), lower=0.03)
                M_sec ~ LogUniform(0.1, 100)
                M = M_pri + M_sec * Octofitter.mjup2msol

                rv = 12.60e3
                plx ~ Uniform(20, 40)
                pmra ~ Uniform(-100, 100)
                pmdec ~ Uniform(-100, 100)

                ra_hip_offset_mas ~ Normal(0, 1000)
                dec_hip_offset_mas ~ Normal(0, 1000)
                dec = $hip_like.hip_sol.dedeg + ra_hip_offset_mas / 60 / 60 / 1000
                ra = $hip_like.hip_sol.radeg + dec_hip_offset_mas / 60 / 60 / 1000 / cos(dec)

                ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd
            end
        )

        model = Octofitter.LogDensityModel(sys_mass)
        @test model isa Octofitter.LogDensityModel

        # Sample from the posterior
        chain = octofit(model, iterations=100, adaptation=100)
        @test chain isa Octofitter.Chains
        @test size(chain, 1) > 0
    end

    # ============================================================================
    # Test: Thiele-Innes Orbital Basis (from thiele-innes.md)
    # ============================================================================
    @testset "Thiele-Innes Orbital Basis" begin
        astrom_dat = Table(;
            epoch = [50000, 50120, 50240, 50360, 50480, 50600, 50720, 50840],
            ra    = [-505.7637580573554, -502.570356287689, -498.2089148883798, -492.67768482682357, -485.9770335870402, -478.1095526888573, -469.0801731788123, -458.89628893460525],
            dec   = [-66.92982418533026, -37.47217527025044, -7.927548139010479, 21.63557115669823, 51.147204404903704, 80.53589069730698, 109.72870493064629, 138.65128697876773],
            σ_ra  = [10, 10, 10, 10, 10, 10, 10, 10],
            σ_dec = [10, 10, 10, 10, 10, 10, 10, 10],
            cor   = [0, 0, 0, 0, 0, 0, 0, 0]
        )

        astrom_like = PlanetRelAstromLikelihood(
            astrom_dat,
            name = "GPI",
            variables = @variables begin
                jitter = 0
                northangle = 0
                platescale = 1
            end
        )

        planet_b = Planet(
            name="b",
            basis=ThieleInnesOrbit,
            likelihoods=[astrom_like],
            variables=@variables begin
                e ~ Uniform(0.0, 0.5)
                A ~ Normal(0, 1000)
                B ~ Normal(0, 1000)
                F ~ Normal(0, 1000)
                G ~ Normal(0, 1000)

                M = system.M
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000.0; system.plx, M, e, A, B, F, G)
            end
        )

        sys = System(
            name="TutoriaPrime",
            companions=[planet_b],
            likelihoods=[],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )

        model = Octofitter.LogDensityModel(sys)
        @test model isa Octofitter.LogDensityModel

        # Sample from the posterior (Thiele-Innes is typically very fast)
        results = octofit(model, iterations=200, adaptation=200)
        @test results isa Octofitter.Chains
        @test size(results, 1) > 0
    end

    # ============================================================================
    # Test: KDE Priors (from priors.md)
    # ============================================================================
    @testset "KDE Priors" begin
        # Create a smoothed KDE estimate from samples
        kde = Octofitter.KDEDist(randn(1000) .+ 10)

        planet_b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            likelihoods=[],
            variables=@variables begin
                a ~ kde
                e ~ Uniform(0.0, 0.99)
                i ~ Sine()
                M = system.M
                ω ~ UniformCircular()
                Ω ~ UniformCircular()
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000; M, e, a, i, ω, Ω)
            end
        )

        sys = System(
            name="Tutoria",
            companions=[planet_b],
            likelihoods=[],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )

        model = Octofitter.LogDensityModel(sys)
        @test model isa Octofitter.LogDensityModel

        chain = octofit(model, iterations=100, adaptation=100)
        @test chain isa Octofitter.Chains
        @test size(chain, 1) > 0
    end

    # ============================================================================
    # Test: Quick Start Example (from quick-start.md)
    # ============================================================================
    @testset "Quick Start Example" begin
        astrom_dat = Table(
            epoch = [50000, 50120, 50240],
            ra = [-505.7, -502.5, -498.2],
            dec = [-66.9, -37.4, -7.9],
            σ_ra = [10.0, 10.0, 10.0],
            σ_dec = [10.0, 10.0, 10.0],
            cor = [0.0, 0.0, 0.0]
        )
        astrom = PlanetRelAstromLikelihood(astrom_dat, name="GPI astrom")

        planet_b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            likelihoods=[astrom],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                a ~ Uniform(0, 100)
                e ~ Uniform(0.0, 0.5)
                i ~ Sine()
                ω ~ UniformCircular()
                Ω ~ UniformCircular()
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50000; M, e, a, i, ω, Ω)
            end
        )

        sys = System(
            name="HD1234",
            companions=[planet_b],
            likelihoods=[],
            variables=@variables begin
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )

        model = Octofitter.LogDensityModel(sys)
        @test model isa Octofitter.LogDensityModel

        chain = octofit(model, iterations=100, adaptation=100)
        @test chain isa Octofitter.Chains
        @test size(chain, 1) > 0

        # Test saving and loading chain
        Octofitter.savechain("test_chain.fits", chain)
        chain_loaded = Octofitter.loadchain("test_chain.fits")
        @test chain_loaded isa Octofitter.Chains
        rm("test_chain.fits")
    end

    # ============================================================================
    # Test: Astrometry with Instrument Calibration Variables (from rel-astrom.md)
    # ============================================================================
    @testset "Astrometry with Instrument Calibration" begin
        astrom_dat = Table(;
            epoch= [50000, 50120, 50240, 50360, 50480, 50600, 50720, 50840],
            ra   = [-505.764, -502.57, -498.209, -492.678, -485.977, -478.11, -469.08, -458.896],
            dec  = [-66.9298, -37.4722, -7.92755, 21.6356, 51.1472, 80.5359, 109.729, 138.651],
            σ_ra = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0],
            σ_dec = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0],
            cor =  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )

        # With calibration variables (jitter, northangle, platescale)
        astrom_like = PlanetRelAstromLikelihood(
            astrom_dat,
            name = "GPI astrom",
            variables = @variables begin
                jitter ~ Uniform(0, 10)
                northangle ~ Normal(0, deg2rad(1))
                platescale ~ truncated(Normal(1, 0.01), lower=0)
            end
        )

        planet_b = Planet(
            name="b",
            basis=Visual{KepOrbit},
            likelihoods=[astrom_like],
            variables=@variables begin
                plx = system.plx
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                a ~ Uniform(0, 100)
                e ~ Uniform(0.0, 0.5)
                i ~ Sine()
                ω ~ UniformCircular()
                Ω ~ UniformCircular()
                θ ~ UniformCircular()
                tp = θ_at_epoch_to_tperi(θ, 50420; M, e, a, i, ω, Ω)
            end
        )

        sys = System(
            name = "Tutoria_calib",
            companions=[planet_b],
            likelihoods=[],
            variables=@variables begin
                plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
            end
        )

        model = Octofitter.LogDensityModel(sys)
        @test model isa Octofitter.LogDensityModel

        chain = octofit(model, iterations=100, adaptation=100)
        @test chain isa Octofitter.Chains
        @test size(chain, 1) > 0
    end

end  # main testset
