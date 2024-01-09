# Calling from Python

This page provides some guidance on how Octofitter can be used from Python.
Please note that this is not officially supported and performance may suffer. In general, we recommend you download Julia and copy-paste the examples as needed.

Still, it might be useful to embed Octofitter as part of a larger Python project or pipeline.
In those cases, you might consider using Octofitter via [juliacall.py](https://pyjulia.readthedocs.io/en/stable/index.html).

The broad instructions are as follows:

### Step 1
Install the "JuliaCall" python package:
```bash
python3 -m pip install juliacall
```

### Step 3
From inside python, install the Octofitter and Distributions Julia packages:
```python
import juliacall
from juliacall import Main as jl
from juliacall import Pkg
Pkg.add(jl.map(jl.String, ['Octofitter','Distributions','Plots']))
```

!!! note You only need to run this step once to install everything. Don't repeat it each time you fit a model.

### Step 4 
Use Octofitter from inside python.

```python

# Import packages
from juliacall import Main as jl
jl.seval("using Octofitter")
jl.seval("using Distributions")

import numpy as np

# Now create a data table
jl.astrom = jl.PlanetRelAstromLikelihood(jl.Table(
    epoch = np.array([5000,5120,5240,5360,5480,5600,5720,5840]),
    ra = np.array([-505.7637580573554,-502.570356287689,-498.2089148883798,-492.67768482682357,-485.9770335870402,-478.1095526888573,-469.0801731788123,-458.89628893460525]),
    dec = np.array([-66.92982418533026,-37.47217527025044,-7.927548139010479,21.63557115669823,51.147204404903704,80.53589069730698,109.72870493064629,138.65128697876773]),
    σ_ra = np.array([10,10,10,10,10,10,10,10.0]),
	σ_dec = np.array([10,10,10,10,10,10,10,10]),
	cor= np.array([0,0,0,0,0,0,0,0.0])
))
# Print it out to the screen (optional)
jl.astrom._jl_display()

# Create a planet model
# If needed, you could pass these variables in from Python by setting e.g. `jl.e_max=1` and then writing `e ~ Uniform(0, e_max)` below
jl.seval("""
@planet B Visual{KepOrbit} begin
    a ~ truncated(Normal(10, 4), lower=0, upper=100)
    e ~ Uniform(0.0, 0.5)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    τ ~ UniformCircular(1.0)

    τ ~ UniformCircular(1.0)
    P = √(b.a^3/system.M)
    tp =  b.τ*b.P + 58849 # reference epoch for τ. Choose an MJD date near your data.
end astrom
""")
# Print it out to the screen (optional)
jl.B._jl_display()

# Create system model
system = jl.seval("""
@system HD82134 begin
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end B
""")
# Print it out to the screen (optional)
jl.B._jl_display()


model = jl.Octofitter.LogDensityModel(system)
# Print it out to the screen (optional)
model._jl_display()

# Sample
chain = jl.octofit(
    model, 0.85,
    adaptation =   500,
    iterations =  1000,
    verbosity = 4,
    max_depth = 12
)
# Display results (recommended)
chain._jl_display()

# Save chain to FITS file (optional)
jl.Octofitter.savechain("mychain.fits", chain)

# Plot chains (optional)
jl.seval("using Plots")
jl.Octofitter.plotchains(chain, jl.Symbol("B"), kind=jl.Symbol("astrometry"), 
color="B_a")
jl.Plots.savefig("orbits.png")
```
