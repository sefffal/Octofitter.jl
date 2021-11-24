

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
Distributions.minimum(d::Sine) = 0
Distributions.maximum(d::Sine) = π
Distributions.insupport(d::Sine, x::Real) = 0 < x < π
Distributions.mean(d::Sine) = π/2

# See https://stats.libretexts.org/Bookshelves/Probability_Theory/Probability_Mathematical_Statistics_and_Stochastic_Processes_(Siegrist)/05%3A_Special_Distributions/5.27%3A_The_Sine_Distribution
Distributions.var(d::Sine) = 1/4 - 2/pi^2 
Distributions.cdf(d::Sine, x::Real)= 1/2*(1−cos(x))
Distributions.quantile(d::Sine, p::Real) = acos(1−2p)
