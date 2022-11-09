using Distributions, Logging

# There is some runtime code generation, so the easiest way
# to make sure everything we need is precompiled is just to
# run a super (super) quick fit during precompilation.
let

    logger = SimpleLogger(stdout, Logging.Error)
    with_logger(logger) do
    
        planet = DirectDetections.Planet{VisualOrbit}(
            Variables(
                a = Uniform(1, 50),
                e = Beta(1.2, 5),
                τ = Uniform(0, 1),
                ω = Uniform(0, 2pi),
                i = Uniform(0, pi),
                Ω = Uniform(0, pi),
                mass = Uniform(0, 1),
            ),
            Astrometry(
                (epoch=1234.0, ra=123., dec=123., σ_ra=12., σ_dec=34.),
            ),
            name=:test
        )

        pma = ProperMotionAnom(
            (;
                ra_epoch=1234.,
                dec_epoch=1234.,
                pm_ra=12.,
                σ_pm_ra=.34,
                pm_dec=12.,
                σ_pm_dec=.34,
                dt = 365.
            ),
        )
        system = System(
            Variables(
                M = TruncatedNormal(1.0, 0.2, 0.1, Inf),
                plx = TruncatedNormal(12.0, 0.2, 0.1, Inf),
            ),
            pma,  
            planet,
            name=:Test
        )

        # This can't actually start since AdvancedHMC uses
        # Requires to load in ForwardDiff. It's not available
        # during precompile.

        # DirectDetections.hmc(
        #     system, 0.85,
        #     adaptation = 5,
        #     iterations = 10,
        #     tree_depth = 5,
        #     initial_samples = 5,
        # )

        # let's trigger some things manually though
        initial_θ_0 = sample_priors(system)
        ln_prior_transformed = make_ln_prior_transformed(system)
        arr2nt = DirectDetections.make_arr2nt(system) 

    end
    nothing
end