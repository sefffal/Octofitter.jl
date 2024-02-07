# [Calling from Python](@id python)

This page provides some guidance on how Octofitter can be used from Python. 
Our general recomendation is to download Julia and copy-paste the examples as needed, but there may be cases where it useful to embed Octofitter within a larger Python project or pipeline

In those cases, you might consider using Octofitter via [juliacall.py](https://pyjulia.readthedocs.io/en/stable/index.html), as demonstrated in these examples.

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
Pkg.add(jl.map(jl.String, ['Octofitter','Distributions','Plots', 'CairoMakie', 'PairPlots']))
```

!!! note
    You only need to run this step once to install everything. Don't repeat it each time you fit a model.

### Step 4 
Use Octofitter from inside python.

```python
# Import packages
import juliacall
from juliacall import Main as jl
jl.seval("using Octofitter, Distributions, Plots, CairoMakie, PairPlots")
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

jl.seval("""
@planet B Visual{KepOrbit} begin
    a ~ LogUniform(0.1, 500)
    e ~ Uniform(0.0, 0.99)
    i ~ Sine()
    ω ~ UniformCircular()
    Ω ~ UniformCircular()
    θ ~ UniformCircular()
    tp = θ_at_epoch_to_tperi(system,b,50000) # use MJD epoch of your data here!!
end astrom
""")
jl.seval("""
@system HD82134 begin
    M ~ truncated(Normal(1.2, 0.1), lower=0)
    plx ~ truncated(Normal(50.0, 0.02), lower=0)
end B
""")

# Print it out to the screen (optional)
jl.HD82134._jl_display()

model = jl.Octofitter.LogDensityModel(jl.HD82134)
# Print it out to the screen (optional)
model._jl_display()

# Sample
chain = jl.octofit(model,)
# Display results summary (recommended)
chain._jl_display()

# Save chain to FITS file (optional)
jl.Octofitter.savechain("mychain.fits", chain)

# Make corner plot
fig = jl.octocorner(model, chain, small=True) # saved to "HD82134-pair-plot-small.png"

# Plot orbits
jl.octoplot(model, chain) # saved to "HD82134-plot-grid.png"
```
