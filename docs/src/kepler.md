# Kepler Solver

The heart of this package is being able to take a set of Keplerian elements and output relative positions, velocities, etc.
For this, we use DirectOrbits.jl which adopts the same conventions as Orbitize!.

The Kepler solver used to go from mean anomaly to eccentric anomaly is a tweaked version copied from [AstroLib.jl](http://juliaastro.github.io/AstroLib.jl/stable/ref/#AstroLib.kepler_solver).

From AstroLib.jl:

> Many different numerical methods exist to solve Kepler's equation. This function implements the algorithm proposed in Markley (1995) Celestial Mechanics and Dynamical Astronomy, 63, 101 ([DOI:10.1007/BF00691917](http://dx.doi.org/10.1007/BF00691917)). This method is not iterative, requires only four transcendental function evaluations, and has been proved to be fast and efficient over the entire range of elliptic motion 0≤e≤10.

On my laptop, this solves for a single eccentric anomaly in just 47 ns.
Since it is implemented in pure Julia, there is no overhead from calling into a C or Cython compiled function and no need for vectorization.