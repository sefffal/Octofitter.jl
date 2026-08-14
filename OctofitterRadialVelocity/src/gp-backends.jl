# ---------------------------------------------------
# Gaussian-process backends for `RadialVelocityObs`
#
# The likelihood moved to core Octofitter, which depends on neither
# AbstractGPs nor Celerite, so the backend-specific calls come in as methods
# on the three hooks core declares (`gp_condition`, `gp_ln_like`,
# `gp_predict`). v1 spelled the same fork as `gp isa Celerite.CeleriteGP`
# branches inside the likelihood body, which was only possible while the
# likelihood lived here.
#
# AbstractGPs needs nothing: core's defaults *are* the AbstractGPs shape
# (`gp(x, σ²)` then `logpdf`), which is plain duck typing and pulls in no
# package. Only Celerite, whose API is imperative — factorize into the GP
# object, then ask it for a likelihood — needs methods.
# ---------------------------------------------------

# `compute!` takes the *standard deviation* per point (it squares internally),
# while `RadialVelocityObs` carries variances, hence the `sqrt`. The factorized
# state lives in `gp`, so that is what the other two hooks receive.
function Octofitter.gp_condition(gp::Celerite.CeleriteGP, epochs, σ²)
    Celerite.compute!(gp, collect(epochs), sqrt.(σ²))
    return gp
end

Octofitter.gp_ln_like(gp::Celerite.CeleriteGP, residuals) =
    Celerite.log_likelihood(gp, collect(residuals))

# v1 called this as `Main.Celerite.predict`, which resolves only if the user
# happens to have a `Celerite` binding in `Main` — i.e. the cross-validation
# path threw `UndefVarError` for everyone else. Same call, correct module.
Octofitter.gp_predict(gp::Celerite.CeleriteGP, residuals, epochs) =
    Celerite.predict(gp, collect(residuals), collect(epochs); return_var=true)

# AbstractGPs prediction, which core declares but cannot implement (it does not
# depend on AbstractGPs) and v1 never implemented at all — its cross-validation
# path errored for every non-Celerite GP, and its plots reached for
# `Main.AbstractGPs` and silently drew nothing if the user had not imported it.
# It is two lines through the package's own API, and both `gp_predict` callers
# want it: cross-validating a GP-correlated `RadialVelocityObs`, and the
# correlated-noise band `Octofitter.noisemodel` hands the plot layer.
Octofitter.gp_predict(fx::AbstractGPs.FiniteGP, residuals, epochs) =
    AbstractGPs.mean_and_var(AbstractGPs.posterior(fx, collect(residuals)), collect(epochs))
