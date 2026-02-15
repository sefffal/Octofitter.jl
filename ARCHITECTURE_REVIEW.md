# Octofitter Architecture Review: Simplification Opportunities

## Overview

This review identifies structural patterns in Octofitter.jl where the code architecture could be simplified. The findings are grouped by theme, ordered roughly by impact.

---

## 1. The "System Hierarchy Traversal" Pattern is Duplicated ~8 Times

**The core issue:** Octofitter's data model is a nested hierarchy: System -> (priors, observations, planets) -> planet -> (priors, observations). Nearly every code-generation function walks this same hierarchy independently, with copy-pasted loop structures.

**Files and functions affected:**

| Function | File | Lines | What it does over the hierarchy |
|----------|------|-------|---------------------------------|
| `make_arr2nt` | `variables.jl` | 723-1024 | Unrolls array -> named tuple |
| `make_ln_prior` | `variables.jl` | 1053-1164 | Unrolls prior log-density |
| `make_ln_prior_transformed` | `variables.jl` | 1170-1334 | Unrolls transformed prior log-density |
| `make_prior_sampler` | `variables.jl` | 1337-1409 | Unrolls prior sampling |
| `make_Bijector_invlinkvec` | `variables.jl` | 1414-1458 | Unrolls inverse-link transform |
| `_list_priors` | `variables.jl` | 656-695 | Collects flat list of priors |
| `result2mcmcchain` | `sampling.jl` | 414-498 | Flattens nested NT to matrix |
| `flatten_named_tuple` | `sampling.jl` | 760-833 | Flattens nested NT to flat NT |
| `mcmcchain2result` | `sampling.jl` | 506-752 | Reconstructs nested NT from chain |

Each function has the same four-level nested loop pattern:
```
for prior in system.priors ...
  for obs in system.observations ...
    for planet in system.planets ...
      for obs in planet.observations ...
```

**Concrete example of duplication** — extracting values from the hierarchy. Compare `result2mcmcchain` (sampling.jl:427-488) with `flatten_named_tuple` (sampling.jl:764-830): these are essentially the same code doing the same traversal with slightly different output formats (push to `Float64[]` vs push to `Pair{Symbol,Float64}[]`). The comment block at line 756-759 is even duplicated verbatim.

**Suggested simplification:** Extract a generic `foreach_prior(f, system)` iterator that calls `f(distribution, index)` for each prior in the canonical order. The six code-generation functions could each become a thin wrapper that provides different `f` callbacks. Similarly, a `foreach_variable(f, named_tuple)` could replace the duplicated flattening/reconstruction logic.

---

## 2. Prior Evaluation Code Generation is Heavily Duplicated

**`make_ln_prior` vs `make_ln_prior_transformed`** (variables.jl:1053-1334)

These two functions are ~280 lines each and share ~80% of their structure. The only meaningful differences are:
- `logpdf` vs `logpdf_with_trans`
- The `Beta` distribution special-case (in `make_ln_prior`) vs the "heal out-of-bounds" logic (in `make_ln_prior_transformed`)
- The extra `sampled` parameter in the transformed version
- The vector-valued distribution handling (SVector construction in transformed version)

Both walk the same system/observation/planet/planet-observation hierarchy. The prior evaluation expression construction at each level is copy-pasted with minor variations.

**Suggested simplification:** A single `_make_ln_prior_body(system; transformed=false)` that generates the expressions, parameterized by which `logpdf` variant to use and what error-handling strategy to apply.

---

## 3. `make_arr2nt` Has O(N^2) Variable Rebinding

**Location:** `variables.jl:761-784` (system derived), `893-914` (planet derived), `949-978` (planet observation derived)

For each derived variable, the code re-binds ALL previously available variables:
```julia
$([:($(k) = _prev.$k) for k in all_available_keys]...)
```

With N derived variables, this creates N*(N-1)/2 total binding expressions. For a model with 20 derived variables, that's 190 redundant bindings in the generated code. Each derived variable sees bindings for every variable before it, even if it only uses one.

**Suggested simplification:** Instead of re-binding everything, the derived expressions could just reference `_prev.varname` directly (or use a single `let` that destructures once). The `_prev` named tuple already contains everything needed.

---

## 4. Epoch-Gathering Logic is Duplicated Between `make_ln_like` and `generate_from_params`

**Location:** `likelihoods/system.jl:20-39` and `likelihoods/system.jl:240-262`

The exact same epoch-collection logic appears twice — once in `make_ln_like` and once in `generate_from_params`. Both build `all_epochs` and `epoch_start_index_mapping` with identical code:

```julia
all_epochs = Float64[]
epoch_start_index_mapping = Dict{Any,Int}()
j = 1
for obs in system.observations
    if hasproperty(obs, :table) && hasproperty(obs.table, :epoch)
        epoch_start_index_mapping[obs] = j
        j += length(obs.table.epoch)
        append!(all_epochs, obs.table.epoch)
    end
end
for i in 1:length(system.planets)
    for like in system.planets[i].observations
        ...identical...
    end
end
```

**Suggested simplification:** Extract a `gather_epochs(system)` helper that returns `(all_epochs, epoch_start_index_mapping)`.

---

## 5. `simulate()` Functions Bypass Context Structs, Taking 8-10 Loose Arguments

**Location:** `likelihoods/relative-astrometry.jl`, `likelihoods/g23h.jl`, `likelihoods/hgca.jl`, `likelihoods/hgca-linfit.jl`, `likelihoods/hipparcos.jl`

The codebase defines `SystemObservationContext` and `PlanetObservationContext` structs to bundle observation parameters. The `ln_like` methods correctly use these context types. However, the companion `simulate()` functions bypass these structs entirely and take the same parameters as 5-10 loose arguments:

```julia
# relative-astrometry.jl:104 — 10 arguments!
function simulate!(ra_model_buf, dec_model_buf, astrom::PlanetRelAstromObs,
    θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet,
    orbit_solutions_i_epoch_start)

# hgca-linfit.jl:289 — 6 arguments
function simulate(hgca_like::HGCAObs,
    θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)

# hipparcos.jl:543 — 6 arguments
function simulate(hiplike::HipparcosIADObs,
    θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
```

The argument group `(θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)` is repeated across at least 7 files. This is exactly the set of fields in `SystemObservationContext`. For planet observations, `(θ_system, θ_planet, θ_obs, orbits, orbit_solutions, i_planet, orbit_solutions_i_epoch_start)` matches `PlanetObservationContext`.

Some `simulate()` functions are also pure pass-throughs that just forward all arguments (e.g., `relative-astrometry.jl:145-150` allocates buffers then passes all 8 args to `simulate!`).

**Suggested simplification:** Refactor `simulate()` to accept the context struct, matching the `ln_like` signature. This eliminates the repeated argument lists and makes the API consistent.

---

## 6. `solutions_list` Tuple Passed Redundantly to Every Likelihood

**Location:** `likelihoods/system.jl:48, 81, 96, 158`

Every observation likelihood receives the *entire* `solutions_list` tuple (all orbit solutions for all planets), even though planet-specific likelihoods only need their own planet's solutions. The `PlanetObservationContext` already carries `i_planet` to index into this tuple, but the full tuple is still passed.

Similarly, `SystemObservationContext` (variables.jl:22-28) carries `orbit_solutions::NTuple{N,TSolutions}` — the full set for all planets — even when a system observation only cares about specific ones.

This isn't just a clarity issue — it means the generated code interpolates the full solutions list expression into every likelihood call expression, bloating the generated function body.

**Suggested simplification:** For `PlanetObservationContext`, pass only `orbit_solutions[i_planet]` instead of the full tuple, and remove the `i_planet` field. The individual likelihood functions that need cross-planet solutions (if any exist) could receive the full tuple as a special case.

---

## 6. `mcmcchain2result` Uses Fragile String-Based Key Matching

**Location:** `sampling.jl:506-752`

The `reform()` inner function (lines 602-744) reconstructs the nested named tuple from a flat MCMCChain by:
1. Using `startswith(string(kout), pk*"_")` to determine if a key belongs to a planet
2. Using `startswith(string(kout), pk*"_"*ok*"_")` to determine if it's a planet observation
3. Using regex replacement `r"^"*string(pk)*"_"` to strip prefixes

This is fragile — if a system variable name happens to start with a planet name followed by `_`, it would be misclassified. The function is also ~150 lines of manual Dict construction that mirrors the flattening logic in reverse.

**Suggested simplification:** Instead of string-matching to reverse the flattening, store the hierarchy structure (which keys belong to which level) alongside the chain, or use the model's `arr2nt` to reconstruct directly.

---

## 7. `ℓπcallback` Closure Uses an Unusual Pattern to Avoid Type Instability

**Location:** `logdensitymodel.jl:97-198`

The likelihood callback is constructed via an immediately-invoked function expression (IIFE):

```julia
ℓπcallback, ∇ℓπcallback = (function(arr2nt, system, ...)
    function ℓπcallback(θ_transformed, system=system, arr2nt=arr2nt, ...)
        ...
    end
    ...
    return ℓπcallback, ∇ℓπcallback
end)(arr2nt, system, Bijector_invlinkvec, ln_prior_transformed, ln_like_generated, D)
```

The inner `ℓπcallback` function also uses default arguments as a closure capture mechanism (`system=system, arr2nt=arr2nt, ...`), which is unconventional. A normal `let` block would be clearer:

```julia
ℓπcallback = let arr2nt=arr2nt, system=system, ...
    function(θ_transformed; sampled=true)
        ...
    end
end
```

This is a minor readability issue but the current pattern is confusing to readers.

---

## 8. Type Stability Diagnostics Are Repetitive

**Location:** `logdensitymodel.jl:200-226`

Five `Core.Compiler.return_type` calls check different components, each followed by a similar `@warn` block. The pattern:

```julia
out_type_X = Core.Compiler.return_type(X, typeof((args...)))
if !isconcretetype(out_type_X)
    @warn "... not type stable ..."
end
```

...is repeated for `ℓπcallback`, `∇ℓπcallback`, `arr2nt`, `ln_prior_transformed`, and `ln_like_generated`.

**Suggested simplification:** A helper `_check_type_stability(name, f, args...; warn_msg)` would reduce this to five one-liners.

---

## 9. The `@variables` Macro Output Format is a Raw Tuple

**Location:** `macros.jl:167-211`

The `@variables` macro returns a raw `Tuple` of `(Priors, Derived, likelihoods...)`. Whether likelihoods are present changes the tuple length, and downstream code has to check `length(vars) > 2` to determine this (see `_vcat_two_variables` at macros.jl:317-401).

**Suggested simplification:** A `VariablesBlock` struct wrapping `priors`, `derived`, and `likelihoods` fields would make the API self-documenting and remove the need for length-checking.

---

## 10. `ln_like` Methods Are Duplicated Across Context Types

**Location:** `variables.jl:313-416`

Several observation types define near-identical `ln_like` methods for both `SystemObservationContext` and `PlanetObservationContext`. The core logic is the same; only the parameter extraction differs:

**`UnitLengthPrior`** (lines 313-327) — two methods that differ only in how `x` and `y` are extracted:
```julia
# System version (line 314):
x = getproperty(ctx.θ_system, X)
# Planet version (line 322):
θ = merge(ctx.θ_system, ctx.θ_planet)
x = getproperty(θ, X)
```

**`UserLikelihood`** (lines 354-384) — two 14-line methods with identical logic after the first line of parameter extraction.

**`BlankLikelihood`** (lines 409-416) — two methods that are 100% identical in their bodies.

This pattern also appears in `generate_from_params` which has duplicate one-liner methods for both context types (line 350-351).

**Suggested simplification:** Add a `get_params(ctx::SystemObservationContext)` / `get_params(ctx::PlanetObservationContext)` helper that returns a merged parameter named tuple, then write each `ln_like` once against the unified parameter set. Alternatively, use `AbstractObservationContext` dispatch with a shared implementation.

---

## 11. Minor: Commented-Out Code and Duplicated Comments

Several files contain substantial blocks of commented-out code that add noise:

- `variables.jl:1040-1052` — commented-out alternative `make_ln_prior` implementation
- `variables.jl:835-837, 905-907, 969-971` — identical commented-out Distribution error checks in three places
- `likelihoods/system.jl:162-164` — commented-out validity check
- `logdensitymodel.jl:135-144` — commented-out debug logging
- `sampling.jl:756-759` — duplicated comment block (the same comment appears twice verbatim)

---

## Summary of Recommendations (by estimated impact)

1. **High impact:** Extract a shared hierarchy traversal pattern used by the ~8 functions that walk system/obs/planet/planet-obs. This would eliminate hundreds of lines of near-identical code in `variables.jl` and `sampling.jl`.

2. **High impact:** Unify `make_ln_prior` and `make_ln_prior_transformed` via a shared parameterized code generator.

3. **Medium impact:** Refactor `simulate()` functions to accept the existing context structs instead of 5-10 loose arguments. This affects ~8 functions across 5+ likelihood files.

4. **Medium impact:** Extract `gather_epochs(system)` to deduplicate epoch collection in `likelihoods/system.jl`.

5. **Medium impact:** Pass per-planet solutions to `PlanetObservationContext` instead of the full solutions tuple.

6. **Medium impact:** Replace string-based key reconstruction in `mcmcchain2result` with a structural approach.

7. **Low impact:** Fix O(N^2) variable rebinding in `make_arr2nt` derived variable expansion.

8. **Low impact:** Simplify the IIFE closure pattern in `logdensitymodel.jl` to a `let` block.

9. **Low impact:** Extract type-stability checking helper in `logdensitymodel.jl`.

10. **Low impact:** Add a `get_params(ctx)` helper to unify `ln_like` implementations across `SystemObservationContext` and `PlanetObservationContext`.

11. **Low impact:** Introduce a `VariablesBlock` struct for the `@variables` macro output.

12. **Low impact:** Remove commented-out code and duplicated comments.
