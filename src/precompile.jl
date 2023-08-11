using Distributions, Logging, ForwardDiff

# There is some runtime code generation, so the easiest way
# to make sure everything we need is precompiled is just to
# run a super (super) quick fit during precompilation.
using PrecompileTools

# define methods, types, etc

@setup_workload  begin
    # Silence logging during precompile (but don't precompile these logging calls)
    io = IOBuffer();
    # io = stdout
    logger = SimpleLogger(io)

    with_logger(logger) do
        @compile_workload begin
        
            astrom = AstrometryLikelihood(
                (epoch=1234.0, ra=123., dec=123., σ_ra=12., σ_dec=34.),
            )
            @planet test Visual{KepOrbit} begin
                a ~ Uniform(1, 50)
                e ~ Beta(1.2, 5)
                τ ~ UniformCircular(1.0)
                ω ~ UniformCircular()
                i ~ Sine()
                Ω ~ UniformCircular()
                mass ~ Uniform(0, 1)
                computed1 = test.a
            end astrom

            @system Test begin
                M ~ TruncatedNormal(1.0, 0.2, 0.1, Inf)
                plx ~ TruncatedNormal(12.0, 0.2, 0.1, Inf)
            end test
            show(io, "text/plain", test)
            show(io, "text/plain", Test)

            model = Octofitter.LogDensityModel(Test)

            # This can't actually start since AdvancedHMC uses
            # Requires to load in ForwardDiff. It's not available
            # during precompile.

            output = Octofitter.advancedhmc(
                model, 0.85,
                adaptation = 5,
                iterations = 10,
                tree_depth = 5,
                initial_samples = 5,
                verbosity=0
            )
            # This would precompile the chain display (but crashes)
            # show(io, "text/plain", output)

            # We can't just precompile displaying the chain due to annoying
            # issues with MCMCChains.jl.
            # So lets just call some diagnostics manually and hope that 
            # speeds up most things.
            # MCMCChains.summarize(output)

            # Neither of these work, and both crash Julia.
            nothing
        end
    end
end
