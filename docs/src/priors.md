# [Priors](@id priors)

All parameters to your model must have a prior defined.
You may provide any continuous, univariate distribution from the Distributions.jl.
A few useful distributions include:

* `Normal`
* `Uniform`
* `LogNormal`
* `LogUniform`
* `TrucatedNormal`
* `VonMises`

This pacakge also defined the `Sine()` distribution for e.g. inclination priors.

The VonMise distribution is notable but not commonly used. It is the analog of a normal distribution defined on a circular domain (-π, +π). If you have a Gaussian prior on an angular parameter, a Von Mises distribution is probably more appropriate.

Behind the scenes, DirectDetections remaps your parameters to unconstrained domains using the Bijectors.jl (and corrects the priors accordingly). This is essential for good sampling efficiency with HMC based samplers.

This means that e.g. if you define the eccentricity prior as `e=Uniform(0,0.5)`, the sampler will actually generate values across the whole real line and transform them back into the `[0,0.5]` range before evaluating the orbit.
**It is therefore essential that your priors do not include invalid domains.**

For example, setting `a=Normal(3,2)` will result in poor sampling efficiency as sometimes negative values for semi-major axis will be drawn (especially if you're using the parallel tempered sampler).

Instead, for parameters like semi-major axis, eccentricity, parallax, and masses, you should truncate any distributions that have negative tails.
This can easily be accomplished with `TrauncatedNormal` or `Trunacted(dist, low, high)` for any arbitrary distribution.

