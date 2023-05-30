# Calling from Python

This page provides some guidance on how Octofitter can be used from Python.
Please note that this is not officially supported and performance will likely suffer. In general, we recommend you download Julia and copy-paste the examples as needed.

Still, it might be useful to embed Octofitter as part of a larger Python project or pipeline.
In those cases, you might consider using Octofitter via [pyjulia](https://pyjulia.readthedocs.io/en/stable/index.html).

The broad instructions are as follows:

1. Install `pyjulia`: `python3 -m pip install julia`. 
2. Configure the `pyjulia` installation by starting python and running `import julia; julia.install()`
3. Install Octofitter: `from julia import Pkg; Pkg.add(['Octofitter','Distributions'])`

Please note that these instructions have not been tested. Corrections would be welcomed.

From here, you should hopefully be able to use Octofitter from Python:
```python
from julia import Main

# Load packages
Main.eval("""
using Octofitter
using Distributions
""")

# Define model
Main.eval("""
astrom = AstrometryLikelihood(
    (epoch = 5000, ra = -505.7637580573554, dec = -66.92982418533026, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 5120, ra = -502.570356287689, dec = -37.47217527025044, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 5240, ra = -498.2089148883798, dec = -7.927548139010479, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 5360, ra = -492.67768482682357, dec = 21.63557115669823, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 5480, ra = -485.9770335870402, dec = 51.147204404903704, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 5600, ra = -478.1095526888573, dec = 80.53589069730698, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 5720, ra = -469.0801731788123, dec = 109.72870493064629, σ_ra = 10, σ_dec = 10, cor=0),
    (epoch = 5840, ra = -458.89628893460525, dec = 138.65128697876773, σ_ra = 10, σ_dec = 10, cor=0),
)

@planet B Visual{KepOrbit} begin
    a ~ truncated(Normal(10, 4), lower=0, upper=100)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    τ ~ UniformCircular(1.0)
end astrom

@system HD82134 begin
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end B

model = Octofitter.LogDensityModel(HD82134)

""")

# Sample
Main.eval("""
chain = Octofitter.advancedhmc(
    model, 0.85;
    adaptation =   500,
    iterations =  1000,
    verbosity = 4,
    tree_depth = 12
)
""")

# Save chain
Main.eval("""
Octofitter.savechain("mychain.fits", chain)
""")


```