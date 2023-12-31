using Octofitter, Distributions, CSV
##
## Create a 1-planet model
@named b = Planet{Visual{KepOrbit}}(
    Variables(
        a    = Uniform(9, 100),
        e    = Uniform(0, 1),
        τ    = Uniform(0, 1),
        # i    = Sine(),
        i = (sys, pl) -> sys.i,
        mass = Uniform(0, 50),
        # mass = 5,
        Ω = (sys, pl) -> sys.Ω,
        ω    = Uniform(0, 2π),
    ),
    CSV.read("cEri-astrometry.csv", PlanetRelAstromLikelihood)
)

@named c = Planet{Visual{KepOrbit}}(
    Variables(
        a    = Uniform(1, 9),
        e    = 0,
        τ    = Uniform(0, 1),
        i = (sys, pl) -> sys.i,
        # mass = Uniform(0, 50),
        mass = TruncatedNormal(5, 5, 0, Inf),
        Ω = (sys, pl) -> sys.Ω,
        ω    = 0,
    ),
)

gaia_id = 3205095125321700480

@named cEri = System(
    Variables(
        plx    = gaia_plx(;gaia_id),
        M      = TruncatedNormal(1.75, 0.05, 0, Inf),
        i      = Sine(),
        Ω      = Uniform(0, 2π),
        pmra   = Normal(0, 200),
        pmdec  = Normal(0, 200),
        # rv     = Normal(0, 200),
        # jitter = TruncatedNormal(0, 5, 0, Inf)
    ),
    HGCALikelihood(;gaia_id),
    # RadialVelocityLikelihood(
    #     (epoch=57374.0, rv=400 - 112.918, σ_rv=10.657),
    #     (epoch=57376.0, rv=400 - 112.464, σ_rv=10.392),
    #     (epoch=57415.0, rv=400 - 110.406, σ_rv=10.188),
    #     (epoch=57649.0, rv=400 - 142.053, σ_rv=10.059),
    #     (epoch=57652.0, rv=400 - 141.521, σ_rv=10.465),
    #     (epoch=57739.0, rv=400 - 153.258, σ_rv=10.121),
    #     (epoch=58068.0, rv=400 - 187.509, σ_rv=10.041),
    #     (epoch=58442.0, rv=400 - 219.468, σ_rv=10.815),
    # ),
    b#, c
)


system = cEri

 # Choose parameter dimensionality and initial parameter value
initial_θ_0 = sample_priors(system)
D = length(initial_θ_0)

ln_prior_transformed = Octofitter.make_ln_prior_transformed(system)
# ln_prior = make_ln_prior(system)
arr2nt = Octofitter.make_arr2nt(system) 

priors_vec = Octofitter._list_priors(system)
Bijector_invlinkvec = Octofitter.make_Bijector_invlinkvec(priors_vec)

# Capture these variables in a let binding to improve performance
lnp = let system=system, ln_prior_transformed=ln_prior_transformed, arr2nt=arr2nt#, ln_prior=ln_prior
    function (θ_t)
        # Transform back from the unconstrained support to constrained support for the likelihood function
        θ = Bijector_invlinkvec(θ_t)
        # θ = θ_t
        θ_res = arr2nt(θ)
        ll = ln_prior_transformed(θ) + Octofitter.ln_like(system, θ_res)
        return ll
    end
end

##
θ = sample_priors(system)
θ_t = Octofitter.Bijectors.link.(priors_vec, θ)


##
lnp(θ_t)

##
using BenchmarkTools
@benchmark lnp(θ_t)

##
@code_warntype lnp(θ_t)

##
@code_warntype Bijector_invlinkvec(θ_t)
@benchmark Bijector_invlinkvec(θ_t)
@code_typed Bijector_invlinkvec(θ_t)

@code_warntype arr2nt(θ)

##
using JET
@report_call arr2nt(θ)
@report_opt arr2nt(θ)
@report_call Bijector_invlinkvec(θ_t)
@report_opt Bijector_invlinkvec(θ_t)
@report_call lnp(θ_t)
@report_opt lnp(θ_t)



##
@report_opt Bijector_invlinkvec(θ_t)
@report_opt ln_prior_transformed(θ)
@report_opt Octofitter.ln_like(system, θ_res)


@report_opt Bijector_invlinkvec(θ_t)
@report_opt arr2nt(θ)
@report_opt ln_prior_transformed(θ) + Octofitter.ln_like(system, θ_res)
