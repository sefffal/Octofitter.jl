using Distributions, Logging, ForwardDiff

# There is some runtime code generation, so the easiest way
# to make sure everything we need is precompiled is just to
# run a super (super) quick fit during precompilation.
using SnoopPrecompile

# define methods, types, etc

@precompile_setup begin
    # Silence logging during precompile (but don't precompile these logging calls)
    io = IOBuffer();
    # io = stdout
    logger = SimpleLogger(io)

    with_logger(logger) do
        @precompile_all_calls begin
        
            planet = DirectDetections.Planet{VisualOrbit}(
                Variables(
                    a = Uniform(1, 50),
                    e = Beta(1.2, 5),
                    τ = UniformCircular(1.0),
                    ω = UniformCircular(),
                    i = Sine(),
                    Ω = UniformCircular(),
                    mass = Uniform(0, 1),
                ),
                Astrometry(
                    (epoch=1234.0, ra=123., dec=123., σ_ra=12., σ_dec=34.),
                ),
                name=:test
            )

            system = System(
                Variables(
                    M = TruncatedNormal(1.0, 0.2, 0.1, Inf),
                    plx = TruncatedNormal(12.0, 0.2, 0.1, Inf),
                ),
                planet,
                name=:Test
            )
            show(io, "text/plain", planet)
            show(io, "text/plain", system)


            # This can't actually start since AdvancedHMC uses
            # Requires to load in ForwardDiff. It's not available
            # during precompile.

            output = DirectDetections.hmc(
                system, 0.85,
                adaptation = 5,
                iterations = 10,
                tree_depth = 5,
                initial_samples = 5,
                verbosity=4
            )

            # This would precompile the chain display,
            # but currently fails during precompile due to
            # an eval inside PrettyTables.jl
            # show(io, "text/plain", output)
            nothing

        end
    end
end
