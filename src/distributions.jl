# This file contains additional custom distibutions using the interface of Distibutions.jl

"""
    Sine()

A custom univariate distribution.
The pdf is a sine function defined between 0 and π.
This is a common prior distribution used when fitting orbits to astrometry.

The full Distributions.jl interface is not yet defined for this distribution,
but the following methods work:
pdf, logpdf, minimum, maximum, insupport, mean, var, cdf, quantile
"""
struct Sine <: ContinuousUnivariateDistribution end
export Sine

function Distributions.pdf(d::Sine, x::Real)
    if 0 < x < π
        return sin(x)/2
    else
        return 0
    end
end
function Distributions.logpdf(d::Sine, x::Real)
    if 0 < x < π
        return log(sin(x)/2)
    else
        return -Inf
    end
end
Distributions.minimum(d::Sine) = 0.0+eps()
Distributions.maximum(d::Sine) = π-eps()
Distributions.insupport(d::Sine, x::Real) = 0 < x < π
Distributions.mean(d::Sine) = π/2

# See https://stats.libretexts.org/Bookshelves/Probability_Theory/Probability_Mathematical_Statistics_and_Stochastic_Processes_(Siegrist)/05%3A_Special_Distributions/5.27%3A_The_Sine_Distribution
Distributions.var(d::Sine) = 1/4 - 2/pi^2 
Distributions.cdf(d::Sine, x::Real)= 1/2*(1-cos(x))
Distributions.quantile(d::Sine, p::Real) = acos(1-2p)



"""
    UniformImproper()

A custom univariate distribution.
The pdf is exactly 1 between -Inf and Inf. This is an **improper**
distribution. This is useful for some change of variable transformations
where we want to apply no prior on a variable and instead apply it on
a transformed quantity. Use with caution.

The full Distributions.jl interface is not obeyed by this distribution,
but the following methods work:
pdf, logpdf, minimum, maximum, insupport, mean, var, cdf, quantile
"""
struct UniformImproper <: ContinuousUnivariateDistribution end
Distributions.pdf(d::UniformImproper, x::Real) = 1.0
Distributions.logpdf(d::UniformImproper, x::Real) = log(1.0)
# obviously these are not statistically correct since they cannot be
# defined for an improper distribution
Distributions.minimum(d::UniformImproper) = -Inf
Distributions.maximum(d::UniformImproper) = +Inf
Distributions.insupport(d::UniformImproper, x::Real) = isfinite(x)
Distributions.mean(d::UniformImproper) = 1.0
Distributions.var(d::UniformImproper) = Inf
# Distributions.cdf(d::UniformImproper, x::Real)= 1.0
Distributions.quantile(d::UniformImproper, p::Real) = p


"""
    kde = KDEDist(data)

A univariate distribution that obeys the Distributions.jl interface.
Uses KernelDensity.jl to create a 1D kernel density estimator using the provided
input data and optional bandwidth scale factor.

Appropriate to use as a prior in an Octofitter model.
"""
struct KDEDist{TKDE<:InterpKDE,TDat<:AbstractArray} <: ContinuousUnivariateDistribution
    ik::TKDE
    data::TDat
    bandwidth::Float64
    lower::Float64
    upper::Float64
end
function KDEDist(data; bandwidth=KernelDensity.default_bandwidth(data), lower=minimum(data), upper=maximum(data))
    T = eltype(data)
    k  = KernelDensity.kde(data; bandwidth, boundary=(lower,upper))
    ik = KernelDensity.InterpKDE(k)
    return KDEDist(ik, data, bandwidth, convert(T,lower), convert(T,upper))
end

Distributions.pdf(kded::KDEDist, x::Real) = pdf(kded.ik, x)
Distributions.logpdf(kded::KDEDist, x::Real) = log(pdf(kded.ik, x))
# obviously these are not statistically correct since they cannot be
# defined for an improper distribution
Distributions.minimum(kded::KDEDist) = kded.lower
Distributions.maximum(kded::KDEDist) = kded.upper
Distributions.insupport(kded::KDEDist, x::Real) = minimum(kded) <= x <= maximum(kded.data)
Distributions.mean(kded::KDEDist) = mean(kded.ik.kde.density)
Distributions.var(kded::KDEDist) = var(kded.ik.kde.density)
# Distributions.cdf(kded::KDEDist, x::Real)= 
Distributions.quantile(kded::KDEDist, p::Real) = quantile(kded.data, p) # TODO: this will be kind of choppy. better to define and use cdf.

# Define a random sampler from a kernel density estimate
# https://discourse.julialang.org/t/sample-from-kernel-density-estimator/50639/2
function Random.rand(rng::AbstractRNG, kded::KDEDist)
    while true
        val = rand(rng, Normal(rand(rng, kded.data), kded.bandwidth))
        if kded.lower < val < kded.upper
            return val
        end
    end
end

function Base.show(io::IO, mime::MIME"text/plain", @nospecialize p::KDEDist)
    println(io, "KDEDist kernel density estimate distribution")
end
function Base.show(io::IO, @nospecialize p::KDEDist)
    print(io, "KDEDist kernel density estimate distribution")
end
export KDEDist