using TypedTables
using Octofitter: Variables, Priors, Derived, FixedPosition
using CSV
using MCMCChains
using HDF5
using Random
using PairPlots, CairoMakie

@testset "Core Functionality" begin


    @testset "Variables and Priors" begin
        @testset "Priors struct" begin
            # Test basic prior creation
            priors = Priors(a = Uniform(0,1), e = Beta(1,2))
            @test priors isa Priors
            @test haskey(priors.priors, :a)
            @test priors.priors[:a] isa Uniform
            @test haskey(priors.priors, :e)
            @test priors.priors[:e] isa Beta
        end

        @testset "Derived struct" begin
            # Test basic derived variable creation
            derived = Derived(
                c = system -> system.a + system.b,
                d = (system, planet) -> planet.e * 2
            )
            @test derived isa Derived
            @test haskey(derived.variables, :c)
            @test haskey(derived.variables, :d)
            @test derived.variables[:c] isa Function
            @test derived.variables[:d] isa Function
        end

        @testset "UniformCircular" begin
            uc = UniformCircular()
            @test uc isa UniformCircular
            @test uc.domain ≈ 2π

            uc = UniformCircular(1.0)
            @test uc.domain ≈ 1.0
        end

        @testset "Variables constructor" begin
            vars = Variables(
                a = Uniform(0,1),
                b = Normal(0,1),
                c = system -> system.a + system.b
            )
            @test vars isa Tuple
            @test length(vars) >= 2  # Should have at least Priors and Derived
            @test first(vars) isa Priors
            @test vars[2] isa Derived
        end
    end

    @testset "Basic Planet Construction" begin
        # Create minimal test data
        astrom = PlanetRelAstromLikelihood(
            (epoch=5000.0, ra=100.0, dec=50.0, σ_ra=1.0, σ_dec=1.0)
        )

        # Test basic planet construction
        planet = Planet{Visual{KepOrbit}}(
            Priors(a = Uniform(0,1), e = Beta(1,2)),
            Derived(),
            (astrom,),
            name=:b
        )
        @test planet isa Planet
        @test planet.name == :b
        @test length(planet.observations) == 1
        @test first(planet.observations) === astrom
        @test planet.priors isa Priors

        # Test with derived variables
        planet = Planet{Visual{KepOrbit}}(
            Priors(a = Uniform(0,1), e = Beta(1,2)),
            Derived(tp = (sys, pl) -> pl.θ * pl.P),
            (astrom,),
            name=:b
        )
        @test !isnothing(planet.derived)
        @test planet.derived isa Derived
    end

    @testset "Fixed Position Orbit" begin
        # Test basic constructor
        fp = Octofitter.FixedPosition(1.0, 2.0, 3.0)
        @test fp.x == 1.0
        @test fp.y == 2.0
        @test fp.z == 3.0

        # Test keyword constructor
        fp = Octofitter.FixedPosition(x=1.0, y=2.0, z=3.0)
        @test fp.x == 1.0
        @test fp.y == 2.0
        @test fp.z == 3.0

        # Test properties
        @test period(fp) == Inf
        @test meanmotion(fp) == 0.0
        @test eccentricity(fp) == 0.0
        @test totalmass(fp) == 0.0
        @test semimajoraxis(fp) == 0.0
        @test periastron(fp) == 0.0
    end
    
    @testset "KepOrbit Models" begin

        # Create planet model with KepOrbit
        @planet b Visual{KepOrbit} begin
            a ~ truncated(Normal(10, 4), lower=0.1)
            e ~ Uniform(0.0, 0.5)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(system,b,50000)
        end 

        @system SimpleSystem begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
        end b

        # Test orbit type extraction
        @test Octofitter.orbittype(b) == Visual{KepOrbit}

        # Test element construction from parameters
        θ_system = (
            M = 1.2,
            plx = 50.0,
            planets = (
                b = (
                    a = 5.0,
                    e = 0.1,
                    i = 0.5,
                    ω = 1.0,
                    Ω = 2.0,
                    tp = 50000.0
                ),
            )
        )
        
        orbit = Octofitter.construct_elements(θ_system, :b)
        @test orbit isa Visual{KepOrbit{Float64},Float64}
        @test orbit.plx ≈ 50.0
        @test semimajoraxis(orbit) ≈ 5.0
    end

    @testset "ThieleInnes Models" begin

        @planet b ThieleInnesOrbit begin
            e ~ Uniform(0.0, 0.5)
            A ~ Normal(0, 10000)
            B ~ Normal(0, 10000)
            F ~ Normal(0, 10000)
            G ~ Normal(0, 10000)
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(system,b,50000.0)
        end 

        @system TISystem begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
        end b

        # Test orbit type extraction
        @test Octofitter.orbittype(b) == ThieleInnesOrbit

        # Test element construction from parameters
        θ_system = (
            M = 1.2,
            plx = 50.0,
            planets = (
                b = (
                    e = 0.1,
                    A = 100.0,
                    B = 200.0,
                    F = 300.0,
                    G = 400.0,
                    tp = 50000.0
                ),
            )
        )
        
        orbit = Octofitter.construct_elements(θ_system, :b)
        @test orbit isa ThieleInnesOrbit
        @test orbit.plx ≈ 50.0
        @test orbit.A ≈ 100.0
    end

    

    @testset "FixedPosition Models" begin

        @planet b Visual{Octofitter.FixedPosition} begin
            x ~ Uniform(-2000, 2000)
            y ~ Uniform(-2000, 2000)
        end 

        @system FixedSystem begin
            plx = 24.4620
        end b

        # Test orbit type extraction
        @test Octofitter.orbittype(b) == Visual{Octofitter.FixedPosition}

        # Test element construction
        θ_system = (
            plx = 24.4620,
            planets = (
                b = (
                    x = 2000.0,
                    y = 1.0,
                ),
            )
        )

        orbit = Octofitter.construct_elements(θ_system, :b)
        @test orbit isa Visual{Octofitter.FixedPosition{Float64}}
        @test posx(orbit,0) ≈ 2000.0
    end

end


@testset "Likelihood Objects" begin
    @testset "PlanetRelAstromLikelihood" begin
        # Test RA/Dec format
        data_radec = PlanetRelAstromLikelihood(
            (epoch=5000.0, ra=100.0, dec=50.0, σ_ra=1.0, σ_dec=1.0),
            (epoch=5100.0, ra=110.0, dec=55.0, σ_ra=1.0, σ_dec=1.0)
        )
        @test data_radec isa PlanetRelAstromLikelihood
        @test length(data_radec.table) == 2
        @test hasproperty(data_radec.table, :ra)

        # Test sep/PA format
        data_seppa = PlanetRelAstromLikelihood(
            (epoch=5000.0, sep=100.0, pa=1.0, σ_sep=1.0, σ_pa=0.1),
            (epoch=5100.0, sep=110.0, pa=1.1, σ_sep=1.0, σ_pa=0.1)
        )
        @test data_seppa isa PlanetRelAstromLikelihood
        @test length(data_seppa.table) == 2
        @test hasproperty(data_seppa.table, :sep)

        # Test that invalid column combinations throw errors
        @test_throws Exception PlanetRelAstromLikelihood(
            (epoch=5000.0, ra=100.0, pa=1.0, σ_ra=1.0, σ_pa=0.1)
        )

        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(data_seppa, 1:1)
        @test length(subset.table) == 1
    end

    @testset "PhotometryLikelihood" begin
        phot = PhotometryLikelihood(
            (band=:Z, phot=15.0, σ_phot=0.1),
            (band=:J, phot=14.0, σ_phot=0.2)
        )
        @test phot isa PhotometryLikelihood
        @test length(phot.table) == 2
        @test all([:band, :phot, :σ_phot] .∈ Ref(propertynames(phot.table)))

        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(phot, 1:1)
        @test length(subset.table) == 1
    end

    @testset "HGCALikelihood" begin
        # Using a known Gaia source ID (HD 91312)
        gaia_id = 756291174721509376
        
        # Test basic constructor
        hgca = HGCALikelihood(;gaia_id=gaia_id,)  
        @test hgca isa HGCALikelihood
        
        # Test epoch subsetting
        # TODO: not yet supported
        # subset = Octofitter.likeobj_from_epoch_subset(hgca, 1:2)
        # @test_broken length(subset.table) == 2
    end

end



# @testset "ImageLikelihood" begin
#     using AstroImages
    
#     # Create a simple test image
#     image = ones(10, 10)
#     image_centered = AstroImages.recenter(AstroImage(image))
    
#     # Create likelihood object
#     imglike = ImageLikelihood(
#         (band=:L, image=image_centered, platescale=10.0, epoch=50000.0),
#         (band=:L, image=image_centered, platescale=10.0, epoch=50100.0)
#     )
    
#     @test imglike isa ImageLikelihood
#     @test length(imglike.table) == 2
#     @test all([:band, :image, :platescale, :epoch] .∈ Ref(propertynames(imglike.table)))
    
#     # Test subsetting
#     subset = likeobj_from_epoch_subset(imglike, 1:1)
#     @test length(subset.table) == 1
# end

# @testset "InterferometryLikelihood" begin
#     # Create test data with closure phases
#     vislike = InterferometryLikelihood(
#         (
#             filename="test.oifits",  # Dummy filename
#             epoch=mjd("2023-06-01"),
#             spectrum_var=:contrast_F480M,
#             use_vis2=false
#         ),
#     )
    
#     @test vislike isa InterferometryLikelihood
#     @test hasproperty(vislike.table, :epoch)
    
#     # Test subsetting
#     subset = likeobj_from_epoch_subset(vislike, 1:1)
#     @test length(subset.table) == 1
# end

# @testset "LogLikelihoodMap" begin
#     using AstroImages
    
#     # Create a simple likelihood map
#     map_data = zeros(10, 10)
#     map_data[5,5] = 1.0  # Peak at center
#     map1 = AstroImage(map_data, (X(-5:4), Y(-5:4)))
    
#     # Create likelihood object
#     loglikemap = LogLikelihoodMap(
#         (
#             epoch=50000.0,
#             map=map1,
#             platescale=1.0
#         ),
#     )
    
#     @test loglikemap isa LogLikelihoodMap
#     @test length(loglikemap.table) == 1
#     @test all([:epoch, :map, :platescale] .∈ Ref(propertynames(loglikemap.table)))
    
#     # Test subsetting
#     subset = likeobj_from_epoch_subset(loglikemap, 1:1)
#     @test length(subset.table) == 1
# end

# @testset "GRAVITYWideKPLikelihood" begin
#     # Create test data
#     kplike = GRAVITYWideKPLikelihood(
#         (
#             filename="test.fits",
#             epoch=60676.0,
#             wavelength_min_meters=2.025e-6,
#             wavelength_max_meters=2.15e-6,
#             spectrum_var=:K_band_spectrum,
#             jitter=:kp_jit,
#             kp_Cy=:kp_Cy
#         ),
#     )
    
#     @test kplike isa GRAVITYWideKPLikelihood
#     @test hasproperty(kplike.table, :epoch)
    
#     # Test subsetting
#     subset = likeobj_from_epoch_subset(kplike, 1:1)
#     @test length(subset.table) == 1
# end



@testset "Prior Specifications" begin
    @testset "Standard Distribution Priors" begin
        # Test normal prior construction and bounds checking
        @test_nowarn @planet b Visual{KepOrbit} begin
            a ~ truncated(Normal(10, 2), lower=0.1)
            e ~ Beta(1, 2)  # Naturally bounded [0,1]
            i ~ Normal(1.0, 0.1)
            ω ~ Uniform(0, 2π)
            Ω ~ Uniform(0, 2π)
            θ ~ UniformCircular()
        end

        # Test LogUniform prior
        @test_nowarn @planet b Visual{KepOrbit} begin
            a ~ LogUniform(0.1, 100)
            e ~ Uniform(0, 0.99)
            i ~ Normal(1.0, 0.1)
            ω ~ Uniform(0, 2π)
            Ω ~ Uniform(0, 2π)
            θ ~ UniformCircular()
        end
    end

    @testset "Special Prior Distributions" begin
        @testset "Sine Distribution" begin
            sine = Sine()
            @test Distributions.minimum(sine) ≈ 0.0 + eps()
            @test Distributions.maximum(sine) ≈ π - eps()
            @test Distributions.insupport(sine, π/2)
            @test !Distributions.insupport(sine, -0.1)
            @test !Distributions.insupport(sine, π + 0.1)
            
            # Test PDF shape
            @test Distributions.pdf(sine, π/2) > Distributions.pdf(sine, 0.1)
            @test Distributions.pdf(sine, π/2) > Distributions.pdf(sine, π-0.1)
        end

        @testset "UniformCircular" begin
            uc = UniformCircular()
            @test uc.domain ≈ 2π
            
            uc_custom = UniformCircular(1.0)
            @test uc_custom.domain ≈ 1.0

            # Test in model context
            @test_nowarn @planet b Visual{KepOrbit} begin
                a ~ LogUniform(0.1, 100)
                e ~ Uniform(0, 0.99)
                i ~ Sine()
                ω ~ UniformCircular()
                Ω ~ UniformCircular(π)  # Custom domain
                θ ~ UniformCircular()
            end
        end
    end

    @testset "Observable-based Priors" begin
        # Create test data
        astrom_data = PlanetRelAstromLikelihood(
            (epoch = 50000, ra = 100.0, dec = 50.0, σ_ra = 1.0, σ_dec = 1.0),
            (epoch = 50100, ra = 110.0, dec = 55.0, σ_ra = 1.0, σ_dec = 1.0)
        )

        # Test O'Neil observable-based prior construction
        obs_prior = ObsPriorAstromONeil2019(astrom_data)
        @test obs_prior isa ObsPriorAstromONeil2019
        @test obs_prior.wrapped_like === astrom_data

        # Test in model context with period prior
        @test_nowarn @planet b Visual{KepOrbit} begin
            e ~ Uniform(0.0, 0.5)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            # Results sensitive to period prior
            P ~ LogUniform(0.1, 150)
            a = ∛(system.M * b.P^2)
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(system,b,50000)
        end astrom_data obs_prior

        # Test subsetting
        subset = Octofitter.likeobj_from_epoch_subset(obs_prior, 1:1)
        @test subset isa ObsPriorAstromONeil2019
    end

    @testset "KDE Priors" begin
        # Create sample data
        samples = randn(1000) .+ 10  # Normal(10,1) samples
        kde = Octofitter.KDEDist(samples)
        
        # Test KDE distribution properties
        @test kde isa Octofitter.KDEDist
        @test Distributions.minimum(kde) ≈ minimum(samples)
        @test Distributions.maximum(kde) ≈ maximum(samples)
        @test Distributions.insupport(kde, 10.0)
        @test !Distributions.insupport(kde, minimum(samples) - 1)

        # Test in model context
        @test_nowarn @planet b Visual{KepOrbit} begin
            a ~ kde  # Use KDE as prior
            e ~ Uniform(0, 0.99)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(system,b,50000)
        end
    end


    @testset "Prior Sampling" begin
        # Create a simple model
        astrom_data = PlanetRelAstromLikelihood(
            (epoch = 50000, ra = 100.0, dec = 50.0, σ_ra = 1.0, σ_dec = 1.0)
        )

        @planet b Visual{KepOrbit} begin
            a ~ LogUniform(0.1, 100)
            e ~ Uniform(0, 0.99)
            i ~ Sine()
            ω ~ UniformCircular()
            Ω ~ UniformCircular()
            θ ~ UniformCircular()
            tp = θ_at_epoch_to_tperi(system,b,50000)
        end astrom_data

        @system TestSystem begin
            M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
        end b

        # Test drawing from priors
        samples = Octofitter.sample_priors(TestSystem)
        @test length(samples) == 11  # Should match number of free parameters (after expanding UniformCircular)
        
        # Test multiple samples
        n_samples = 100
        samples_n = Octofitter.sample_priors(TestSystem, n_samples)
        @test length(samples_n) == n_samples
        @test all(length.(samples_n) .== 11)

        # Verify samples are within bounds
        model = Octofitter.LogDensityModel(TestSystem)
        arr2nt = Octofitter.make_arr2nt(TestSystem)
        
        for sample in samples_n
            params = arr2nt(sample)
            # Test key parameter bounds
            @test 0.1 ≤ params.planets.b.a ≤ 100
            @test 0 ≤ params.planets.b.e ≤ 0.99
            @test 0 ≤ params.planets.b.i ≤ π
            @test params.M ≥ 0.1
            @test params.plx ≥ 0.1
        end
    end
end


@testset "Data I/O" begin
    
    @testset "FITS Chain basic I/O" begin
        # Test handling of unicode characters in parameter names
        test_chain = Chains(
            randn(100, 5, 1),  # 100 samples, 5 parameters, 1 chain
            [:a, :e, :ω, :Ω, :θ]  # Unicode parameter names
        )
        
        # Save and reload with FITS
        Octofitter.savechain("test_unicode.fits", test_chain)
        loaded_chain = Octofitter.loadchain("test_unicode.fits")
        
        @test keys(test_chain) == keys(loaded_chain)
        for name in keys(test_chain)
            @test test_chain[name] ≈ loaded_chain[name]
        end
        
        # Clean up
        rm("test_unicode.fits")
    end

    # Create a test chain 
    @planet b Visual{KepOrbit} begin 
        a ~ LogUniform(0.1, 100)
        e ~ Uniform(0, 0.99)
        i ~ Sine()
        ω ~ UniformCircular() # Unicode omega
        Ω ~ UniformCircular() # Unicode capital omega  
        θ ~ UniformCircular() # Unicode theta
        tp = θ_at_epoch_to_tperi(system,b,50000)
    end

    @system TestSystem begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
    end b

    model = Octofitter.LogDensityModel(TestSystem)
    chain_original = octofit(model, iterations=100, adaptation=50)
    @testset "FITS Chain I/O" begin
        # Save and reload
        Octofitter.savechain("test_chain.fits", chain_original)
        chain_loaded = Octofitter.loadchain("test_chain.fits")
        
        # Test that all columns are preserved
        @test keys(chain_original) == keys(chain_loaded)
        
        # Test that values are preserved
        for name in keys(chain_original)
            @test chain_original[name] ≈ chain_loaded[name]
        end
    end


    # Clean up
    rm("test_chain.fits")

    @testset "HDF5 Chain I/O" begin
        # Save to HDF5 format
        Octofitter.savehdf5("test_chain.h5", model, chain_original)
        
        # Load back
        chain_loaded = Octofitter.loadhdf5("test_chain.h5")
        
        # Test that essential orbital parameters are preserved
        params = ["b_a", "b_e", "b_i", "b_ω", "b_Ω", "M", "plx"]
        for param in params
            @test median(chain_original[param]) ≈ median(chain_loaded[param]) rtol=1e-5
        end
        
        # Clean up
        rm("test_chain.h5")
    end

    @testset "CSV Astrometry Data I/O" begin
        # Create a test CSV file
        test_data = """
        epoch,ra,dec,σ_ra,σ_dec,cor
        50000.0,100.0,50.0,1.0,1.0,0.0
        50100.0,110.0,55.0,1.0,1.0,0.1
        50200.0,120.0,60.0,1.0,1.0,-0.1
        """
        
        open("test_astrometry.csv", "w") do io
            write(io, test_data)
        end
        
        # Load using CSV.read
        astrom_data = CSV.read("test_astrometry.csv", PlanetRelAstromLikelihood)
        
        # Test the loaded data
        @test length(astrom_data.table) == 3
        @test astrom_data.table.epoch[1] ≈ 50000.0
        @test astrom_data.table.ra[2] ≈ 110.0
        @test astrom_data.table.dec[3] ≈ 60.0
        @test astrom_data.table.cor[2] ≈ 0.1
        
        # Clean up
        rm("test_astrometry.csv")
    end

end


@testset "Plotting" begin
    # Create a test model and chain once to use across tests
    astrom_like = PlanetRelAstromLikelihood(
        (epoch = 50000, ra = -494.4, dec = -76.7, σ_ra =  12.6, σ_dec =  12.6, cor=  0.2),
        (epoch = 50120, ra = -495.0, dec = -44.9, σ_ra =  10.4, σ_dec =  10.4, cor=  0.5),
        (epoch = 50240, ra = -493.7, dec = -12.9, σ_ra =   9.9, σ_dec =   9.9, cor=  0.1),
        (epoch = 50360, ra = -490.4, dec =  19.1, σ_ra =   8.7, σ_dec =   8.7, cor= -0.8),
        (epoch = 50480, ra = -485.2, dec =  51.0, σ_ra =   8.0, σ_dec =   8.0, cor=  0.3),
        (epoch = 50600, ra = -478.1, dec =  82.8, σ_ra =   6.9, σ_dec =   6.9, cor= -0.0),
        (epoch = 50720, ra = -469.1, dec = 114.3, σ_ra =   5.8, σ_dec =   5.8, cor=  0.1),
        (epoch = 50840, ra = -458.3, dec = 145.3, σ_ra =   4.2, σ_dec =   4.2, cor= -0.2),
    )
    @planet b Visual{KepOrbit} begin 
        a ~ LogUniform(0.1, 100)
        e ~ Uniform(0, 0.99)
        i ~ Sine()
        ω ~ UniformCircular() # Unicode omega
        Ω ~ UniformCircular() # Unicode capital omega  
        θ ~ UniformCircular() # Unicode theta
        tp = θ_at_epoch_to_tperi(system,b,50000)
        mass = 10
    end astrom_like
    gaia_id = 756291174721509376
    hgca = HGCALikelihood(;gaia_id=gaia_id)  
    @system TestSystem begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
        pmra ~ Normal(-137, 10)
        pmdec ~ Normal(2,  10)    
    end hgca b

    model = Octofitter.LogDensityModel(TestSystem)
    Random.seed!(0)
    Octofitter.default_initializer!(model,nruns=1) # work around non-reproducibility bug
    chain = octofit(model, iterations=500, adaptation=500)
    
    # Test octoplot returns valid figure and handles options
    fig = octoplot(model, chain)
    @test fig isa Makie.Figure
    @test isfile("$(model.system.name)-plot-grid.png")
    
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
    @test length(fig_astrom.content) == 2

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
    @test length(fig_phys_orb.content) == 2


    fig_astrom_time = octoplot(
        model, chain,
        N = 1,
        show_astrom=false,
        show_physical_orbit=false,
        show_astrom_time=true,
        show_pma=false,
        show_mass=false,
        show_rv=false,
        show_relative_rv=false,
        show_hipparcos=false,
    )
    @test length(fig_astrom_time.content) == 3


    fig_hgca = octoplot(
        model, chain,
        N = 1,
        show_astrom=false,
        show_physical_orbit=false,
        show_astrom_time=false,
        show_pma=true,
        show_mass=false,
        show_rv=false,
        show_relative_rv=false,
        show_hipparcos=false,
    )
    @test length(fig_hgca.content) == 6


    fig_mass = octoplot(
        model, chain,
        N = 1,
        show_astrom=false,
        show_physical_orbit=false,
        show_astrom_time=false,
        show_pma=false,
        show_mass=true,
        show_rv=false,
        show_relative_rv=false,
        show_hipparcos=false,
    )
    @test length(fig_mass.content) == 4


    fig_rv = octoplot(
        model, chain,
        N = 1,
        show_astrom=false,
        show_physical_orbit=false,
        show_astrom_time=false,
        show_pma=false,
        show_mass=false,
        show_rv=true,
        show_relative_rv=false,
        show_hipparcos=false,
    )
    @test length(fig_rv.content) == 3


    fig_relative_rv = octoplot(
        model, chain,
        N = 1,
        show_astrom=false,
        show_physical_orbit=false,
        show_astrom_time=false,
        show_pma=false,
        show_mass=false,
        show_rv=false,
        show_relative_rv=true,
        show_hipparcos=false,
    )
    @test length(fig_relative_rv.content) == 3

    # Test corner plot
    fig_corner = octocorner(model, chain, small=true)
    @test fig_corner isa Makie.Figure
end




@testset "Multi-Planet Systems" begin
    # Create test data for two planets
    astrom_b = PlanetRelAstromLikelihood(
        (epoch=50000.0, ra=100.0, dec=50.0, σ_ra=1.0, σ_dec=1.0),
        (epoch=50100.0, ra=110.0, dec=55.0, σ_ra=1.0, σ_dec=1.0)
    )
    astrom_c = PlanetRelAstromLikelihood(
        (epoch=50000.0, ra=-200.0, dec=-100.0, σ_ra=1.0, σ_dec=1.0),
        (epoch=50100.0, ra=-210.0, dec=-110.0, σ_ra=1.0, σ_dec=1.0)
    )

    # Create two-planet model with resonance
    @planet b Visual{KepOrbit} begin
        a ~ LogUniform(0.1, 10)
        e ~ Uniform(0, 0.5)
        i = system.i  # Coplanar
        ω ~ UniformCircular()
        Ω = system.Ω  # Coplanar
        P = 2*system.P_nominal * b.P_mul
        P_mul ~ Normal(1, 0.1)
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(system,b,50000)
    end astrom_b

    @planet c Visual{KepOrbit} begin
        a ~ LogUniform(0.1, 10)
        e ~ Uniform(0, 0.5)
        i = system.i  # Coplanar
        ω ~ UniformCircular()
        Ω = system.Ω  # Coplanar
        P = system.P_nominal * c.P_mul
        P_mul ~ truncated(Normal(1, 0.1), lower=0.1)
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(system,c,50000)
    end astrom_c

    @system TwoPlanetSystem begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ truncated(Normal(50.0, 0.02), lower=0.1)
        i ~ Sine()  # Common inclination
        Ω ~ UniformCircular()  # Common node
        P_nominal ~ LogUniform(50, 300)  # Base period
    end b c

    model = Octofitter.LogDensityModel(TwoPlanetSystem)
    chain = octofit(model, iterations=100)

    # Test that period ratio is preserved
    period_ratio = chain[:b_P][:] ./ chain[:c_P][:]
    @test mean(period_ratio) ≈ 2.0 rtol=0.2

    # Test coplanarity is preserved
    @test all(chain[:b_i][:] .== chain[:c_i][:])
    @test all(chain[:b_Ω][:] .== chain[:c_Ω][:])
end

@testset "Joint Fitting" begin
    # Create test data of different types
    astrom = PlanetRelAstromLikelihood(
        (epoch=50000.0, ra=100.0, dec=50.0, σ_ra=1.0, σ_dec=1.0)
    )
    
    hgca = HGCALikelihood(;gaia_id=756291174721509376)

    @planet b Visual{KepOrbit} begin
        a ~ LogUniform(0.1, 10)
        e ~ Uniform(0, 0.5)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        mass ~ LogUniform(0.1, 100)
        θ ~ UniformCircular()
        tp = θ_at_epoch_to_tperi(system,b,50000)
    end astrom

    @system AstromSystem begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ gaia_plx(;gaia_id=756291174721509376)
        pmra ~ Normal(-975, 10)
        pmdec ~ Normal(20, 10)
        jitter ~ LogUniform(0.1, 100)
        rv0 ~ Normal(0, 100)
    end b

    @system JointSystem begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx ~ gaia_plx(;gaia_id=756291174721509376)
        pmra ~ Normal(-975, 10)
        pmdec ~ Normal(20, 10)
        jitter ~ LogUniform(0.1, 100)
        rv0 ~ Normal(0, 100)
    end hgca b

    # Test that the mass is better constrained with joint data
    # by comparing to chains from individual fits
    model_astrom = Octofitter.LogDensityModel(AstromSystem)
    model_astrom_pma = Octofitter.LogDensityModel(JointSystem)
    
    chain_astrom = octofit(model_astrom, iterations=100)
    chain_astrom_pma = octofit(model_astrom_pma, iterations=100)

    @test std(chain_astrom_pma[:b_mass]) < std(chain_astrom[:b_mass])
end

@testset "Cross Validation" begin
    # Create test model with multiple epochs
    astrom1 = PlanetRelAstromLikelihood(
        (epoch=t, ra=100.0+t/100, dec=50.0+t/200, σ_ra=1.0, σ_dec=1.0)
        for t in range(50000, 51000, length=10)
    )
    astrom2 = PlanetRelAstromLikelihood(
        (epoch=t, ra=100.0+t/100, dec=50.0+t/200, σ_ra=1.0, σ_dec=1.0)
        for t in range(50000, 51000, length=10)
    )

    @planet b Visual{KepOrbit} begin
        a  = 1.0
        e  = 0.0
        i  = 1.0
        ω  = 1.0
        Ω  = 1.0
        tp = 50000
    end astrom1 astrom2

    @system Sys begin
        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
        plx = 100.
    end b
    model = Octofitter.LogDensityModel(Sys)
    chain = octofit(model, iterations=100)

    # Test pointwise likelihood calculation
    like_mat = Octofitter.pointwise_like(model, chain)
    @test size(like_mat, 2) == 20  # One column per epoch
    @test all(isfinite, like_mat)

    # Test k-fold systems generation
    kfold_systems = Octofitter.generate_kfold_systems(model.system)
    @test length(kfold_systems) == 2  # One per dataset

    # Test per-epoch system generation
    per_epoch_systems = Octofitter.generate_system_per_epoch(model.system)
    @test length(per_epoch_systems) == 20
end