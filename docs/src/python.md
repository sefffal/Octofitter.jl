# [Calling from Python](@id python)

This page provides some guidance on how Octofitter can be used from Python. 
Our general recomendation is to download Julia and copy-paste the examples as needed, but there may be cases where it useful to embed Octofitter within a larger Python project or pipeline.

In those cases, you might consider using [octofitterpy](https://github.com/sefffal/octofitterpy). This python package uses [juliacall.py](https://pyjulia.readthedocs.io/en/stable/index.html) to make some Octofitter functionality available in python.

Besides the model definition, most functions can be used the same in Python as in Julia. This [notebook](https://github.com/sefffal/octofitterpy/blob/master/examples/demo.ipynb)  provides some examples translated into Python.

See the [octofitterpy](https://github.com/sefffal/octofitterpy) site for installation instructions.

!!! warning "octofitterpy has not yet been updated for Octofitter v9"
    The v9 model surface is a breaking change: `Planet`/`companions=`/`basis=` were
    replaced by a flat list of `Body` nodes and observations that name their own
    `target` and `ref` (see [Migrating to Octofitter v9](@ref v9-migration)). Examples
    written against octofitterpy's current model-definition helpers will not work
    unchanged. Analysis functions that take a model and a chain — `octoplot`,
    `octocorner`, `construct_system`, the chain IO — are unaffected in shape.
