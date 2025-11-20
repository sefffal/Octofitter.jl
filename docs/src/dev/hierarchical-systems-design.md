# HierarchicalSystem Design Document

**Status:** Design Phase
**Author:** Architecture Planning
**Date:** 2025-11-20
**Version:** 1.0

---

## 1. Overview and Goals

### 1.1 Motivation

Octofitter currently models individual stellar systems with one or more planets. This design document proposes adding support for **hierarchical Bayesian models** where:

1. **Population-level hyperparameters** constrain multiple systems
2. **Nested stellar systems** (e.g., circumbinary planets, hierarchical triples) share parameters
3. **Multi-system fits** enable joint inference with shared priors

### 1.2 Design Goals

**Primary Goals:**
- Enable hierarchical models with hyperparameters at the top level
- Allow child systems to access parent-level variables
- Maintain full backward compatibility with existing `System` models
- Keep metaprogramming performance benefits (compile-time code generation)

**Non-Goals (MVP):**
- Multi-level nesting (only 2 levels: HierarchicalSystem → Systems)
- Cross-validation for multi-system models (single-system only)
- Specialized plotting for multi-system models

---

## 2. Use Cases

### 2.1 Population Study

Model multiple systems with shared hyperparameters:

```julia
population = HierarchicalSystem(
    name = :hot_jupiters,
    variables = @variables begin
        # Population hyperparameters
        mean_log_period ~ Normal(log(3), 0.5)
        sigma_log_period ~ LogNormal(0, 0.3)
        mean_ecc ~ Beta(1, 5)
    end,
    systems = (
        System(
            name = :HD189733,
            variables = @variables begin
                log_P ~ Normal(parent.mean_log_period, parent.sigma_log_period)
                e ~ Beta(1, 1)  # Can override population prior
            end,
            companions = (planet_b,),
            likelihoods = (rv_data_189733,)
        ),
        System(
            name = :HD209458,
            variables = @variables begin
                log_P ~ Normal(parent.mean_log_period, parent.sigma_log_period)
                e ~ Beta(1, 1)
            end,
            companions = (planet_b,),
            likelihoods = (rv_data_209458,)
        ),
    )
)
```

### 2.2 Circumbinary Planet

Model a planet orbiting a binary star:

```julia
# Binary system (inner)
binary = System(
    name = :binary,
    variables = @variables begin
        q ~ Uniform(0.3, 1.0)  # mass ratio
        M_secondary = parent.M_primary * q
        P_binary ~ LogUniform(1, 100)
    end,
    companions = (secondary_star,),
    likelihoods = (binary_rv,)
)

# Circumbinary planet (outer system)
primary = System(
    name = :primary,
    variables = @variables begin
        M_primary ~ Uniform(0.8, 1.5)
        plx ~ Uniform(10, 100)
    end,
    companions = (circumbinary_planet,),
    likelihoods = (photometry,)
)

model = HierarchicalSystem(
    name = :circumbinary_system,
    systems = (primary, binary)
)
```

### 2.3 Single System (Most Common)

Existing single-system models work unchanged:

```julia
system = System(
    name = :HD12345,
    variables = @variables begin
        M ~ Uniform(0.5, 2.0)
        plx ~ Uniform(10, 100)
    end,
    companions = (planet_b, planet_c),
    likelihoods = (rv, astrometry)
)

# Auto-wrapped internally
model = LogDensityModel(system)
```

---

## 3. Design Overview

### 3.1 Type Hierarchy

```
HierarchicalSystem (new)
  ├─ priors::Priors             # Hyperparameters
  ├─ derived::Derived           # Derived hyperparameters
  ├─ observations::Tuple        # Population-level data (optional)
  ├─ systems::NamedTuple        # One or more System objects
  └─ name::Symbol

      System (unchanged structure)
        ├─ priors::Priors         # System-level parameters
        ├─ derived::Derived       # Can access parent.* variables
        ├─ observations::Tuple
        ├─ planets::NamedTuple
        └─ name::Symbol

            Planet (unchanged)
              ├─ priors::Priors
              ├─ derived::Derived  # Can access system.* variables
              ├─ observations::Tuple
              └─ name::Symbol
```

### 3.2 Parameter Structure

Flat array ordering:
1. HierarchicalSystem priors
2. HierarchicalSystem observation priors
3. For each System in order:
   - System priors
   - System observation priors
   - For each Planet:
     - Planet priors
     - Planet observation priors

Nested NamedTuple structure (from `arr2nt`):
```julia
(
    # HierarchicalSystem level
    mean_mass = 1.05,
    sigma_mass = 0.15,
    observations = (;),

    # Systems
    systems = (
        HD1234 = (
            M = 1.12,
            plx = 45.2,
            observations = (rv = (...),),
            planets = (
                b = (a=..., e=..., observations=(...)),
            )
        ),
        HD5678 = (
            M = 0.98,
            plx = 87.3,
            observations = (rv = (...),),
            planets = (
                c = (a=..., e=..., observations=(...)),
            )
        ),
    )
)
```

---

## 4. Data Structures

### 4.1 HierarchicalSystem Struct

**Location:** `src/variables.jl`

```julia
"""
    HierarchicalSystem([derived,] priors, [observations,] systems..., name=:symbol)

Construct a hierarchical model containing one or more System objects.
The hierarchical level can define hyperparameters (priors) and derived variables
that are accessible to child systems via the `parent` binding.

# Arguments
- `priors::Priors`: Hyperparameters for the population/hierarchy
- `derived::Derived`: Derived hyperparameters (optional)
- `observations::Tuple{AbstractLikelihood}`: Population-level data (optional)
- `systems`: One or more System objects
- `name::Symbol`: Name of the hierarchical model

# Examples
```julia
# Population model
pop = HierarchicalSystem(
    name = :population,
    variables = @variables begin
        mean_period ~ Normal(10, 5)
        sigma_period ~ LogNormal(0, 0.5)
    end,
    systems = (system1, system2, system3)
)

# Single system (auto-wrapped)
sys = System(name=:HD1234, ...)
model = LogDensityModel(sys)  # Internally creates HierarchicalSystem
```
"""
struct HierarchicalSystem{TPriors<:Priors, TDerived<:Union{Derived,Nothing},
                          TObs<:Tuple, TSystems<:Union{NamedTuple,Tuple}}
    priors::TPriors
    derived::TDerived
    observations::TObs
    systems::TSystems
    name::Symbol
end
export HierarchicalSystem
```

### 4.2 Constructor

```julia
function HierarchicalSystem(;
    name::Union{Symbol,AbstractString},
    variables::Tuple=(Priors(), Derived()),
    systems::Union{Tuple,NamedTuple,AbstractVector},
    likelihoods::Tuple=()
)
    name = Symbol(name)

    # Parse variables
    if length(variables) >= 2
        priors, derived = variables[1:2]
        additional_likelihoods = length(variables) > 2 ? variables[3:end] : ()
    else
        error("variables must contain at least (Priors, Derived)")
    end

    # Validate systems
    if systems isa AbstractVector
        systems = Tuple(systems)
    end
    for s in systems
        if !(s isa System)
            error("All systems must be System objects, got $(typeof(s))")
        end
    end

    # Convert systems to NamedTuple if not already
    if systems isa Tuple && !isempty(systems) && hasfield(eltype(systems), :name)
        systems_nt = namedtuple(getproperty.(systems, :name), systems)
    else
        systems_nt = systems
    end

    # Validate no duplicate system names
    if systems_nt isa NamedTuple
        if length(keys(systems_nt)) != length(unique(keys(systems_nt)))
            error("Duplicate system names detected")
        end
    end

    # Combine likelihoods
    likes = (likelihoods..., additional_likelihoods...)

    # Check for duplicate observation names
    like_names = String[]
    for like in likes
        like_name = likelihoodname(like)
        if like_name in like_names
            error("HierarchicalSystem $name: Duplicate observation/likelihood name '$like_name'")
        end
        push!(like_names, like_name)
    end

    # Check for duplicate variables in priors and derived
    prior_keys = Set(keys(priors.priors))
    if !isnothing(derived)
        derived_keys = Set(keys(derived.variables))
        overlap = intersect(prior_keys, derived_keys)
        if !isempty(overlap)
            error("HierarchicalSystem $name: Variables $(collect(overlap)) defined as both prior and derived")
        end
    end

    return HierarchicalSystem{typeof(priors), typeof(derived), typeof(likes),
                             typeof(systems_nt)}(
        priors, derived, likes, systems_nt, name
    )
end
```

---

## 5. Variable Scoping and Access

### 5.1 Access Patterns

**Three scoping levels:**

1. **HierarchicalSystem variables** - Accessible to Systems via `parent.*`
2. **System variables** - Accessible to Planets via `system.*` (unchanged)
3. **Planet variables** - Local scope only (unchanged)

### 5.2 Derived Variable Evaluation

In `make_arr2nt`, derived variables are evaluated with appropriate context:

```julia
# HierarchicalSystem derived - no parent
hier_derived = @variables begin
    total_mass = M1 + M2  # Access hier-level priors only
end

# System derived - has parent
sys_derived = @variables begin
    M = parent.mean_mass * mass_ratio  # Access hier-level via parent
    a_au = a * parent.distance  # Can mix parent and local
end

# Planet derived - has system (unchanged)
planet_derived = @variables begin
    mass_msun = system.M * q  # Access system-level
end
```

### 5.3 Implementation in Generated Code

```julia
# In make_arr2nt generated function:
function (arr)
    # 1. HierarchicalSystem level
    hier0 = (; hyperprior1=arr[1], hyperprior2=arr[2])
    hier1 = let _prev=hier0
        (; _prev..., derived_var = expr_without_parent)
    end
    hier = hier1

    # 2. System level (accessing parent)
    sys1_0 = (; prior1=arr[3], prior2=arr[4])
    sys1_1 = let parent=hier, _prev=sys1_0
        (; _prev..., derived_var = expr_with_parent)
    end
    # ... planets ...
    sys1 = (; sys1_1..., planets=(...))

    # 3. Return
    return (; hier..., systems=(; sys1, sys2, ...))
end
```

---

## 6. Core Function Modifications

All functions in `src/variables.jl` and `src/likelihoods/system.jl` need to handle `HierarchicalSystem`.

### 6.1 `_list_priors(hiersys::HierarchicalSystem)`

**Purpose:** Return flat vector of all prior distributions

**Pattern:**
```julia
function _list_priors(hiersys::HierarchicalSystem)
    priors_vec = []

    # 1. HierarchicalSystem priors
    for prior_dist in values(hiersys.priors.priors)
        push!(priors_vec, prior_dist)
    end

    # 2. HierarchicalSystem observation priors
    for obs in hiersys.observations
        if hasproperty(obs, :priors)
            for prior_dist in values(obs.priors.priors)
                push!(priors_vec, prior_dist)
            end
        end
    end

    # 3. Loop over systems (call existing _list_priors(system::System))
    for sys in hiersys.systems
        append!(priors_vec, _list_priors(sys))
    end

    return map(identity, priors_vec)
end
```

**Note:** Existing `_list_priors(system::System)` remains unchanged

---

### 6.2 `make_prior_sampler(hiersys::HierarchicalSystem)`

**Purpose:** Generate function to sample from all priors

**Pattern:**
```julia
function make_prior_sampler(hiersys::HierarchicalSystem)
    prior_sample_expressions = Expr[]

    # 1. HierarchicalSystem priors
    for prior_dist in values(hiersys.priors.priors)
        push!(prior_sample_expressions, :(sample = rand(rng, $prior_dist)))
        for i in 1:length(prior_dist)
            push!(prior_sample_expressions, :(prior_samples = (prior_samples..., sample[$i])))
        end
    end

    # 2. HierarchicalSystem observation priors
    for obs in hiersys.observations
        if hasproperty(obs, :priors)
            for prior_dist in values(obs.priors.priors)
                # ... same pattern ...
            end
        end
    end

    # 3. Loop over systems - inline their prior sampling
    for sys in hiersys.systems
        # System priors
        for prior_dist in values(sys.priors.priors)
            # ... same pattern ...
        end

        # System observation priors
        for obs in sys.observations
            # ...
        end

        # Planets
        for planet in sys.planets
            for prior_dist in values(planet.priors.priors)
                # ...
            end
            # Planet observations
            for obs in planet.observations
                # ...
            end
        end
    end

    return @RuntimeGeneratedFunction(:(function (rng)
        prior_samples = ()
        @inbounds begin
            $(prior_sample_expressions...)
        end
        return prior_samples
    end))
end
```

**Implementation Note:** This is essentially the existing `make_prior_sampler` logic with an outer loop over systems.

---

### 6.3 `make_ln_prior_transformed(hiersys::HierarchicalSystem)`

**Purpose:** Generate function to compute log prior density

**Pattern:**
```julia
function make_ln_prior_transformed(hiersys::HierarchicalSystem)
    i = 0
    prior_evaluations = Expr[]

    # 1. HierarchicalSystem priors
    for prior_dist in values(hiersys.priors.priors)
        if length(prior_dist) > 1
            # Vector-valued
            samples = []
            for _ in 1:length(prior_dist)
                i += 1
                push!(samples, :(arr[$i]))
            end
            samples_expr = :(SVector($(samples...)))
        else
            i += 1
            samples_expr = :(arr[$i])
        end

        ex = :(
            p = logpdf_with_trans($prior_dist, $samples_expr, sampled);
            if !isfinite(p) && $(eltype(prior_dist) <: AbstractFloat)
                # Healing logic...
            end;
            lp += p
        )
        push!(prior_evaluations, ex)
    end

    # 2. HierarchicalSystem observation priors
    for obs in hiersys.observations
        if hasproperty(obs, :priors)
            for prior_dist in values(obs.priors.priors)
                # ... same pattern ...
            end
        end
    end

    # 3. Loop over systems - inline system prior evaluations
    for sys in hiersys.systems
        # System priors
        for prior_dist in values(sys.priors.priors)
            # ... same pattern ...
        end

        # System observation priors
        for obs in sys.observations
            # ...
        end

        # Planets
        for planet in sys.planets
            for prior_dist in values(planet.priors.priors)
                # ... same pattern ...
            end
            # Planet observations
            for obs in planet.observations
                # ...
            end
        end
    end

    return @RuntimeGeneratedFunction(:(function (arr, sampled)
        l = $i
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array (got $(length(arr)))")
        end
        lp = zero(first(arr))
        @inbounds begin
            $(prior_evaluations...)
        end
        return lp
    end))
end
```

---

### 6.4 `make_arr2nt(hiersys::HierarchicalSystem)`

**Purpose:** Convert flat parameter array to nested NamedTuple structure

**High-level structure:**
```julia
function make_arr2nt(hiersys::HierarchicalSystem)
    i = 0

    # ===== HIERARCHICAL LEVEL =====
    # 1. Hierarchical priors
    body_hier_priors = Expr[]
    for key in keys(hiersys.priors.priors)
        if length(hiersys.priors.priors[key]) > 1
            # Handle vector-valued distributions
            ex_is = []
            for _ in 1:length(hiersys.priors.priors[key])
                i += 1
                push!(ex_is, :(arr[$i]))
            end
            ex = :($key = ($(ex_is...),))
        else
            i += 1
            ex = :($key = arr[$i])
        end
        push!(body_hier_priors, ex)
    end

    # 2. Hierarchical derived variables
    body_hier_determ = Expr[]
    if !isnothing(hiersys.derived)
        for (j, (key, expr)) in enumerate(zip(keys(hiersys.derived.variables),
                                              values(hiersys.derived.variables)))
            captured_bindings = [:($(hiersys.derived.captured_names[k]) =
                                   $(hiersys.derived.captured_vals[k]))
                                for k in 1:length(hiersys.derived.captured_names)]

            prior_keys = collect(keys(hiersys.priors.priors))
            derived_keys_so_far = j == 1 ? Symbol[] :
                                  collect(keys(hiersys.derived.variables))[1:j-1]
            all_available = vcat(prior_keys, derived_keys_so_far)

            ex = :(
                $(Symbol("hier$j")) = let _prev=$(Symbol("hier$(j-1)"))
                    (; _prev..., $key = let $(captured_bindings...)
                        $([:($(k) = _prev.$k) for k in all_available]...)
                        result = $expr
                        result
                    end)
                end
            )
            push!(body_hier_determ, ex)
        end
        l = length(keys(hiersys.derived.variables))
        push!(body_hier_determ, :(hier = $(Symbol("hier$l"))))
    else
        push!(body_hier_determ, :(hier = hier0))
    end

    # 3. Hierarchical observations
    body_hier_obs = Expr[]
    for obs in hiersys.observations
        # Similar to system observation handling
        # Generate priors and derived for this observation
        # Store as named tuple entry
    end

    # ===== SYSTEM LEVEL (LOOP) =====
    body_systems = Expr[]
    for sys in hiersys.systems
        # System priors
        body_sys_priors = Expr[]
        for key in keys(sys.priors.priors)
            if length(sys.priors.priors[key]) > 1
                # Vector-valued
                ex_is = []
                for _ in 1:length(sys.priors.priors[key])
                    i += 1
                    push!(ex_is, :(arr[$i]))
                end
                ex = :($key = ($(ex_is...),))
            else
                i += 1
                ex = :($key = arr[$i])
            end
            push!(body_sys_priors, ex)
        end

        # System derived (with parent binding)
        j = 0
        body_sys_determ = Expr[]
        if !isnothing(sys.derived)
            for (key, expr) in zip(keys(sys.derived.variables),
                                   values(sys.derived.variables))
                captured_bindings = [:($(sys.derived.captured_names[k]) =
                                       $(sys.derived.captured_vals[k]))
                                    for k in 1:length(sys.derived.captured_names)]

                prior_keys = collect(keys(sys.priors.priors))
                derived_keys_so_far = j == 0 ? Symbol[] :
                                      collect(keys(sys.derived.variables))[1:j]
                all_available = vcat(prior_keys, derived_keys_so_far)

                ex = :(
                    $(Symbol("sys$(j+1)")) = let parent=hier, _prev=$(Symbol("sys$j"))
                        (; _prev..., $key = let $(captured_bindings...)
                            $([:($(k) = _prev.$k) for k in all_available]...)
                            result = $expr
                            result
                        end)
                    end
                )
                push!(body_sys_determ, ex)
                j += 1
            end
            push!(body_sys_determ, :(sys_resolved = $(Symbol("sys$j"))))
        else
            push!(body_sys_determ, :(sys_resolved = sys0))
        end

        # System observations (existing pattern, adapted)
        body_sys_obs = Expr[]
        # ... similar to current implementation ...

        # Planets (existing pattern, adapted)
        body_planets = Expr[]
        for planet in sys.planets
            # ... existing planet handling ...
            # Key: use 'system=sys_resolved' binding
        end

        # Assemble this system
        ex = :(
            $(sys.name) = begin
                sys0 = (;$(body_sys_priors...))
                $(body_sys_determ...)
                (;
                    sys_resolved...,
                    observations = (;$(body_sys_obs...)),
                    planets = (;$(body_planets...))
                )
            end
        )
        push!(body_systems, ex)
    end

    # ===== RETURN GENERATED FUNCTION =====
    return @RuntimeGeneratedFunction(:(function (arr)
        l = $i
        @boundscheck if length(arr) != l
            error("Expected exactly $l elements in array. Got ", length(arr))
        end

        # Hierarchical level
        hier0 = (;$(body_hier_priors...))
        $(body_hier_determ...)

        # Systems
        $(body_systems...)

        # Merge and return
        return (;
            hier...,
            observations = (;$(body_hier_obs...)),
            systems = (;$(body_systems...))
        )
    end))
end
```

**Key Points:**
- Hierarchical derived variables have no `parent` (or `parent=nothing`)
- System derived variables have `parent=hier`
- Planet derived variables have `system=sys_resolved` (unchanged pattern)
- Careful ordering of `i` increment to match `_list_priors`

---

### 6.5 `make_ln_like(hiersys::HierarchicalSystem, θ_hier)`

**Purpose:** Generate likelihood evaluation function

**High-level structure:**
```julia
function make_ln_like(hiersys::HierarchicalSystem, θ_hier)

    # ===== EPOCH COLLECTION =====
    # Collect all epochs from:
    # 1. HierarchicalSystem observations
    # 2. Each System's observations
    # 3. Each System's Planet observations

    all_epochs = Float64[]
    epoch_start_index_mapping = Dict{Any,Int}()
    j = 1

    # Hierarchical observations
    for obs in hiersys.observations
        if hasproperty(obs, :table) && hasproperty(obs.table, :epoch)
            epoch_start_index_mapping[obs] = j
            j += length(obs.table.epoch)
            append!(all_epochs, obs.table.epoch)
        end
    end

    # Loop over systems
    for sys in hiersys.systems
        # System observations
        for obs in sys.observations
            if hasproperty(obs, :table) && hasproperty(obs.table, :epoch)
                epoch_start_index_mapping[obs] = j
                j += length(obs.table.epoch)
                append!(all_epochs, obs.table.epoch)
            end
        end

        # Planet observations
        for planet in sys.planets
            for obs in planet.observations
                if hasproperty(obs, :table) && hasproperty(obs.table, :epoch)
                    epoch_start_index_mapping[obs] = j
                    j += length(obs.table.epoch)
                    append!(all_epochs, obs.table.epoch)
                end
            end
        end
    end

    # ===== CODE GENERATION =====
    planet_declarations = Expr[]
    planet_construction_exprs = Expr[]
    planet_orbit_solution_exprs = Expr[]
    planet_like_exprs = Expr[]

    planet_sol_keys = Symbol[]
    ll_counter = 0

    # Loop over systems to generate planet code
    for (sys_idx, sys) in enumerate(hiersys.systems)
        for (planet_idx, planet) in enumerate(sys.planets)
            OrbitType = _planet_orbit_type(planet)
            key = Symbol("planet_sys$(sys_idx)_$(planet_idx)")
            sols_key = Symbol("sols_$(key)")
            push!(planet_sol_keys, sols_key)

            # Declaration
            push!(planet_declarations, :($key = nothing))

            # Construction - merge hier, system, and planet variables
            planet_construction = quote
                $key = $(OrbitType)(;
                    merge(
                        θ_hier,  # Access hierarchical variables
                        θ_hier.systems[$(Meta.quot(sys.name))],  # System variables
                        θ_hier.systems[$(Meta.quot(sys.name))].planets[$planet_idx]  # Planet variables
                    )...
                )
            end
            push!(planet_construction_exprs, planet_construction)

            # Orbit solutions
            if isempty(all_epochs)
                orbit_sol_expr = quote
                    $sols_key = ()
                end
            else
                orbit_sol_expr = quote
                    epochs = @alloc(Float64, $(length(all_epochs)))
                    $((:(epochs[$k] = $(all_epochs[k])) for k in 1:length(all_epochs))...)

                    sol0 = orbitsolve($key, first(epochs))
                    $sols_key = @alloc(typeof(sol0), length(epochs))
                    $sols_key[begin] = sol0
                    $_kepsolve_all!(view($sols_key, 2:length(epochs)), $key,
                                   view(epochs, 2:length(epochs)))
                end
            end
            push!(planet_orbit_solution_exprs, orbit_sol_expr)

            # Planet observation likelihoods
            for (like_idx, like) in enumerate(planet.observations)
                i_epoch_start = get(epoch_start_index_mapping, like, 0)
                obs_name = normalizename(likelihoodname(like))

                expr = :(
                    $(Symbol("ll$(ll_counter+1)")) = $(Symbol("ll$ll_counter")) + ln_like(
                        hiersys.systems[$(Meta.quot(sys.name))].planets[$planet_idx].observations[$like_idx],
                        θ_hier,
                        θ_hier.systems[$(Meta.quot(sys.name))],
                        θ_hier.systems[$(Meta.quot(sys.name))].planets[$planet_idx],
                        hasproperty(θ_hier.systems[$(Meta.quot(sys.name))].planets[$planet_idx].observations,
                                   $(Meta.quot(obs_name))) ?
                            θ_hier.systems[$(Meta.quot(sys.name))].planets[$planet_idx].observations.$obs_name :
                            (;),
                        elems,
                        solutions_list,
                        $(length(planet_sol_keys)),  # This planet's index
                        $(i_epoch_start-1)
                    )
                )
                push!(planet_like_exprs, expr)
                ll_counter += 1
            end
        end
    end

    solutions_list = :(tuple($((:($sols_key) for sols_key in planet_sol_keys)...)))

    # System observation likelihoods
    system_like_exprs = Expr[]

    # Hierarchical observations
    for (obs_idx, obs) in enumerate(hiersys.observations)
        i_epoch_start = get(epoch_start_index_mapping, obs, 0)
        obs_name = normalizename(likelihoodname(obs))

        expr = :(
            $(Symbol("ll$(ll_counter+1)")) = $(Symbol("ll$ll_counter")) + ln_like(
                hiersys.observations[$obs_idx],
                θ_hier,
                hasproperty(θ_hier.observations, $(Meta.quot(obs_name))) ?
                    θ_hier.observations.$obs_name : (;),
                elems,
                $solutions_list,
                $(i_epoch_start-1)
            )
        )
        push!(system_like_exprs, expr)
        ll_counter += 1
    end

    # Each system's observations
    for sys in hiersys.systems
        for (obs_idx, obs) in enumerate(sys.observations)
            i_epoch_start = get(epoch_start_index_mapping, obs, 0)
            obs_name = normalizename(likelihoodname(obs))

            expr = :(
                $(Symbol("ll$(ll_counter+1)")) = $(Symbol("ll$ll_counter")) + ln_like(
                    hiersys.systems[$(Meta.quot(sys.name))].observations[$obs_idx],
                    θ_hier,
                    θ_hier.systems[$(Meta.quot(sys.name))],
                    hasproperty(θ_hier.systems[$(Meta.quot(sys.name))].observations,
                               $(Meta.quot(obs_name))) ?
                        θ_hier.systems[$(Meta.quot(sys.name))].observations.$obs_name : (;),
                    elems,
                    $solutions_list,
                    $(i_epoch_start-1)
                )
            )
            push!(system_like_exprs, expr)
            ll_counter += 1
        end
    end

    # Collect all planet keys for elem tuple
    planet_keys = [Symbol("planet_sys$(i)_$(j)")
                   for (i, sys) in enumerate(hiersys.systems)
                   for j in 1:length(sys.planets)]

    # ===== RETURN GENERATED FUNCTION =====
    return @RuntimeGeneratedFunction(:(function (hiersys::HierarchicalSystem, θ_hier)
        T = _system_number_type(θ_hier)
        ll0 = zero(T)

        # Declare planet variables
        $(planet_declarations...)

        # Try to construct orbits
        try
            $(planet_construction_exprs...)
        catch err
            return convert(T, -Inf)
        end

        ll_out = @no_escape begin
            # Construct elem tuple
            elems = tuple($(planet_keys...))

            # Solve orbits
            $(planet_orbit_solution_exprs...)

            # Evaluate likelihoods
            $(planet_like_exprs...)
            $(system_like_exprs...)

            $(Symbol("ll$ll_counter"))
        end

        return ll_out
    end))
end
```

**Key Design Points:**
1. Planets from different systems can orbit at same epochs (all collected together)
2. Planet construction merges hier + system + planet variables
3. Likelihood functions receive full `θ_hier` structure to access any level
4. Solution indexing must account for planets across multiple systems

---

## 7. Backward Compatibility

### 7.1 Auto-wrapping Strategy

**Single-system models automatically wrapped:**

```julia
# In src/logdensitymodel.jl

function LogDensityModel(system::System; kwargs...)
    # Auto-wrap in HierarchicalSystem
    hiersys = HierarchicalSystem(
        name = system.name,
        variables = (Priors(), Derived()),  # Empty hierarchical level
        systems = (system,),
        likelihoods = ()
    )
    return LogDensityModel(hiersys; kwargs...)
end

function LogDensityModel(hiersys::HierarchicalSystem; autodiff=nothing, verbosity=2, chunk_sizes=nothing)
    # Main constructor implementation
    # ... existing LogDensityModel logic, adapted for HierarchicalSystem ...
end
```

### 7.2 Storage in LogDensityModel

Update `LogDensityModel` struct:

```julia
mutable struct LogDensityModel{D,Tℓπ,T∇ℓπ,THier,TLink,TInvLink,TArr2nt,TPriSamp,ADType}
    const D::Int
    const ℓπcallback::Tℓπ
    const ∇ℓπcallback::T∇ℓπ
    const hierarchical_system::THier  # Changed from 'system::TSys'
    const link::TLink
    const invlink::TInvLink
    const arr2nt::TArr2nt
    const sample_priors::TPriSamp
    starting_points::Union{Nothing,Vector}
    # ...
end
```

### 7.3 Accessor Compatibility

Add convenience accessor for backward compatibility:

```julia
# For single-system models, allow model.system
function Base.getproperty(model::LogDensityModel, sym::Symbol)
    if sym === :system
        hiersys = getfield(model, :hierarchical_system)
        if length(hiersys.systems) == 1
            return first(hiersys.systems)
        else
            error("model.system is ambiguous for multi-system models. Use model.hierarchical_system.systems")
        end
    else
        return getfield(model, sym)
    end
end
```

### 7.4 Display Methods

Update display to show structure appropriately:

```julia
function Base.show(io::IO, mime::MIME"text/plain", model::LogDensityModel)
    hiersys = model.hierarchical_system
    n_systems = length(hiersys.systems)
    L = _count_epochs(hiersys)

    if n_systems == 1
        # Single system - use familiar display
        sys = first(hiersys.systems)
        println(io, "LogDensityModel for System $(sys.name) of dimension $(model.D) and $(L) epochs")
    else
        # Multi-system
        println(io, "LogDensityModel for HierarchicalSystem $(hiersys.name)")
        println(io, "  $(n_systems) systems, dimension $(model.D), $(L) epochs")
        for sys in hiersys.systems
            println(io, "    - $(sys.name)")
        end
    end
end
```

---

## 8. Impact on Existing Features

### 8.1 Features Requiring Updates

#### Cross-Validation (`src/cross-validation.jl`)

**Functions to update:**
- `prior_only_model`
- `_count_likeobj`
- `_count_epochs`

**Strategy:** Only support single-system models for cross-validation:

```julia
function prior_only_model(hiersys::HierarchicalSystem; exclude_all=false)
    if length(hiersys.systems) > 1
        error("prior_only_model: Cross-validation not supported for multi-system models")
    end

    # Reconstruct with observations stripped
    sys = first(hiersys.systems)
    newplanets = map(sys.planets) do planet
        # ... existing logic ...
    end

    newsys_obs = map(sys.observations) do obs
        if exclude_all || !_isprior(obs)
            return BlankLikelihood((obs.priors, obs.derived), likelihoodname(obs))
        end
    end
    newsys_obs = filter(!isnothing, newsys_obs)

    newsys = System(
        variables=(sys.priors, sys.derived),
        likelihoods=newsys_obs,
        companions=newplanets,
        name=sys.name
    )

    return HierarchicalSystem(
        variables=(hiersys.priors, hiersys.derived),
        likelihoods=(),
        systems=(newsys,),
        name=hiersys.name
    )
end

function generate_kfold_systems(hiersys::HierarchicalSystem)
    if length(hiersys.systems) > 1
        error("generate_kfold_systems: Cross-validation not supported for multi-system models")
    end
    # Delegate to existing logic for single system
    # ...
end

# Similar for other cross-validation functions
```

**Count functions:**
```julia
function _count_likeobj(hiersys::HierarchicalSystem)::Int
    likeobj_count = 0

    # Hierarchical observations
    for obs in hiersys.observations
        if !_isprior(obs)
            likeobj_count += 1
        end
    end

    # Loop over systems
    for sys in hiersys.systems
        likeobj_count += _count_likeobj(sys)
    end

    return likeobj_count
end

function _count_epochs(hiersys::HierarchicalSystem)::Int
    observation_count = 0

    # Hierarchical observations
    for obs in hiersys.observations
        if hasproperty(obs, :table)
            observation_count += Tables.rowcount(obs.table)
        else
            observation_count += 1
        end
    end

    # Loop over systems
    for sys in hiersys.systems
        observation_count += _count_epochs(sys)
    end

    return observation_count
end
```

#### Sampling Functions (`src/sampling.jl`)

**Functions to check:**
- `sample_priors` - Should work (calls `model.sample_priors`)
- `drawfrompriors` - Update to handle HierarchicalSystem

```julia
function drawfrompriors(hiersys::HierarchicalSystem)
    θ = sample_priors(hiersys)
    arr2nt = make_arr2nt(hiersys)
    θnt = arr2nt(θ)
    return θnt
end
```

#### Optimization (`src/optimization.jl`)

**No changes expected** - optimization works on `LogDensityModel` interface

#### Result Processing

**Functions accessing `.system` field:**
- Update to use `.hierarchical_system`
- Use accessor for backward compat
- Document that multi-system models require explicit access

### 8.2 Features NOT Requiring Changes

- **Samplers** (Pigeons, NUTS, etc.) - work via `LogDensityModel` interface
- **Orbit solvers** - unchanged
- **Individual likelihood types** - unchanged
- **Parameter transformations** (Bijectors) - unchanged

### 8.3 Plotting (Out of Scope)

Plotting functions will need updates but are **not included in MVP**. Document limitations:

```julia
# In plotmodel, etc.
function plotmodel(model::LogDensityModel)
    hiersys = model.hierarchical_system
    if length(hiersys.systems) > 1
        error("Plotting multi-system models not yet supported. Plot individual systems: plotmodel(model, :system_name)")
    end
    # Existing plotting logic...
end
```

---

## 9. Pigeons Extension

### 9.1 Updates Required

**File:** `ext/OctofitterPigeonsExt/OctofitterPigeonsExt.jl`

```julia
# Update non-sampleable prior check
function _has_non_sampleable_priors(model)
    hiersys = model.hierarchical_system
    ret = false

    # Check hierarchical observations
    ret |= any(Octofitter._isprior, hiersys.observations)

    # Check each system
    for sys in hiersys.systems
        ret |= any(Octofitter._isprior, sys.observations)
        for planet in sys.planets
            ret |= any(Octofitter._isprior, planet.observations)
        end
    end

    return ret
end

# Update Chains constructor
function MCMCChains.Chains(model::LogDensityModel, pt::Pigeons.PT, chain_num::Union{Nothing,Int}=nothing)
    hiersys = model.hierarchical_system
    ln_prior = Octofitter.make_ln_prior_transformed(hiersys)
    ln_like = Octofitter.make_ln_like(hiersys, model.arr2nt(model.sample_priors(Random.default_rng())))

    # ... rest unchanged, works on θ arrays ...
end
```

### 9.2 Reference Model

```julia
function Pigeons.default_reference(target::Octofitter.LogDensityModel)
    reference_hiersys = prior_only_model(target.hierarchical_system)
    reference = Octofitter.LogDensityModel(reference_hiersys; verbosity=0)
    reference.starting_points = target.starting_points
    return reference
end
```

---

## 10. Implementation Phases

### Phase 1: Core Structure (2-3 hours)

**Files:** `src/variables.jl`

- [ ] Define `HierarchicalSystem` struct
- [ ] Implement constructor with validation
- [ ] Update `_list_priors(::HierarchicalSystem)`
- [ ] Add display methods for `HierarchicalSystem`
- [ ] Write unit tests for construction

**Checkpoint:** Can construct HierarchicalSystem, list priors correctly

---

### Phase 2: Prior Functions (2-3 hours)

**Files:** `src/variables.jl`

- [ ] Implement `make_prior_sampler(::HierarchicalSystem)`
- [ ] Implement `make_ln_prior_transformed(::HierarchicalSystem)`
- [ ] Test prior sampling and evaluation
- [ ] Verify parameter ordering matches `_list_priors`

**Checkpoint:** Can sample from and evaluate priors for hierarchical models

---

### Phase 3: Parameter Structure (3-4 hours)

**Files:** `src/variables.jl`

- [ ] Implement `make_arr2nt(::HierarchicalSystem)`
- [ ] Handle hierarchical, system, and planet levels
- [ ] Implement `parent` binding for system derived variables
- [ ] Test nested NamedTuple structure
- [ ] Verify variable scoping works correctly

**Checkpoint:** Can convert flat array to full nested structure with proper scoping

---

### Phase 4: Likelihood (3-4 hours)

**Files:** `src/likelihoods/system.jl`

- [ ] Implement `make_ln_like(::HierarchicalSystem, θ_hier)`
- [ ] Handle epoch collection across hierarchy
- [ ] Generate orbit construction with merged variables
- [ ] Generate likelihood evaluation for all levels
- [ ] Test likelihood computation

**Checkpoint:** Full likelihood evaluation works

---

### Phase 5: LogDensityModel Integration (2-3 hours)

**Files:** `src/logdensitymodel.jl`

- [ ] Update `LogDensityModel` struct to store `hierarchical_system`
- [ ] Implement auto-wrapping for single `System`
- [ ] Add backward-compatible accessor for `.system`
- [ ] Update display methods
- [ ] Update type stability checks
- [ ] Test model construction

**Checkpoint:** Can create LogDensityModel from both System and HierarchicalSystem

---

### Phase 6: Utilities (2-3 hours)

**Files:** `src/cross-validation.jl`, `src/sampling.jl`

- [ ] Update `prior_only_model(::HierarchicalSystem)`
- [ ] Update `_count_likeobj(::HierarchicalSystem)`
- [ ] Update `_count_epochs(::HierarchicalSystem)`
- [ ] Add multi-system guards to cross-validation
- [ ] Update `drawfrompriors(::HierarchicalSystem)`
- [ ] Test utility functions

**Checkpoint:** Utilities work for single-system, error appropriately for multi-system

---

### Phase 7: Pigeons Extension (1-2 hours)

**Files:** `ext/OctofitterPigeonsExt/OctofitterPigeonsExt.jl`

- [ ] Update `_has_non_sampleable_priors`
- [ ] Update `default_reference`
- [ ] Update `Chains` constructor
- [ ] Test Pigeons sampling

**Checkpoint:** Pigeons sampler works with hierarchical models

---

### Phase 8: Documentation & Testing (2-3 hours)

**Files:** `docs/`, `test/`

- [ ] Write comprehensive integration tests
- [ ] Test single-system backward compatibility
- [ ] Test population model example
- [ ] Test circumbinary planet example
- [ ] Document usage in user guide
- [ ] Document limitations (no cross-validation for multi-system)
- [ ] Update migration guide

**Checkpoint:** Feature is tested and documented

---

**Total Estimated Time:** 17-25 hours

---

## 11. Testing Strategy

### 11.1 Unit Tests

**Test hierarchical construction:**
```julia
@testset "HierarchicalSystem construction" begin
    hier = HierarchicalSystem(
        name = :test,
        variables = @variables begin
            hyper ~ Normal(0, 1)
        end,
        systems = (
            System(name=:sys1, variables=@variables(M ~ Uniform(0.5, 1.5))),
            System(name=:sys2, variables=@variables(M ~ Uniform(0.5, 1.5)))
        )
    )

    @test length(hier.systems) == 2
    @test hier.name == :test
    @test length(hier.priors.priors) == 1
end
```

**Test prior ordering:**
```julia
@testset "Prior ordering" begin
    # Build hierarchical model
    # Sample priors
    # Check arr2nt produces correct structure
    # Verify parameters in correct positions
end
```

### 11.2 Integration Tests

**Single system (backward compatibility):**
```julia
@testset "Single system compatibility" begin
    # Build traditional System
    # Create LogDensityModel
    # Verify .system accessor works
    # Sample with Pigeons
    # Check output structure
end
```

**Population model:**
```julia
@testset "Population model" begin
    pop = HierarchicalSystem(
        name = :population,
        variables = @variables begin
            mean_period ~ Normal(10, 5)
            sigma_period ~ LogNormal(0, 0.5)
        end,
        systems = [
            System(
                name = Symbol("sys$i"),
                variables = @variables begin
                    log_P ~ Normal(log(parent.mean_period), parent.sigma_period)
                end,
                companions = (test_planet,)
            ) for i in 1:3
        ]
    )

    model = LogDensityModel(pop)

    # Test structure
    θ = sample_priors(Random.default_rng(), model)
    θ_nt = model.arr2nt(θ)
    @test haskey(θ_nt, :mean_period)
    @test haskey(θ_nt, :systems)
    @test length(θ_nt.systems) == 3

    # Test likelihood
    ll = model.ℓπcallback(model.link(θ))
    @test isfinite(ll)

    # Test sampling
    chain = octofit_pigeons(model; n_rounds=5, n_chains=4)
    @test size(chain.chain, 1) > 0
end
```

**Circumbinary planet:**
```julia
@testset "Circumbinary system" begin
    # Build nested system model
    # Verify planet can access parent.M_primary
    # Test orbit construction
    # Sample and verify
end
```

### 11.3 Performance Tests

```julia
@testset "Performance" begin
    # Single-system model
    # Check type stability of arr2nt
    # Check type stability of ℓπcallback
    # Verify no performance regression
end
```

---

## 12. Documentation Requirements

### 12.1 User Guide

**New section:** "Hierarchical Models"

Topics:
- What are hierarchical models?
- When to use them
- Population studies example
- Nested systems example
- Variable scoping (`parent.*` syntax)
- Limitations

### 12.2 API Documentation

Document:
- `HierarchicalSystem` constructor
- `parent` binding in derived variables
- Output structure from `arr2nt`
- Backward compatibility notes

### 12.3 Migration Guide

Add section on:
- `model.system` → `model.hierarchical_system.systems[1]` for explicit access
- Automatic wrapping behavior
- What breaks (if anything)

---

## 13. Future Extensions (Beyond MVP)

### 13.1 Multi-Level Nesting

Allow `HierarchicalSystem` to contain other `HierarchicalSystem` objects:
- Recursively handle arbitrary depth
- More complex variable scoping
- Requires recursive code generation

### 13.2 Cross-Validation for Multi-System

Support leave-one-system-out:
- Generate N models, each with N-1 systems
- Evaluate held-out system likelihood
- Useful for population model validation

### 13.3 Multi-System Plotting

- Joint corner plots across systems
- Population-level parameter visualization
- System comparison plots

### 13.4 System-Level Observations

Observations that constrain relationships between systems:
- Relative astrometry between wide companions
- Constraints on population statistics
- Hierarchical priors from external catalogs

---

## 14. Open Questions

### 14.1 Resolved

None yet - awaiting implementation start

### 14.2 To Resolve During Implementation

1. **Orbit construction:** Should planets in subsystems have access to both parent and local variables automatically, or require explicit merging?

2. **Observation likelihood signatures:** Do we need new `ln_like` methods, or can existing methods handle the full hierarchy by receiving `θ_hier`?

3. **Parameter names in chains:** Should multi-system chains namespace parameters like `systems.HD1234.M`, or keep flat like `HD1234_M`?

4. **Empty hierarchical level:** When auto-wrapping single systems, should we include the empty hierarchical level in the output structure, or flatten it away?

---

## 15. Success Criteria

### 15.1 Functional Requirements

- [ ] Can construct `HierarchicalSystem` with multiple systems
- [ ] Child systems can access parent variables via `parent.*`
- [ ] Single-system models work unchanged (backward compatible)
- [ ] Can sample from hierarchical priors
- [ ] Can evaluate likelihood for hierarchical models
- [ ] Pigeons sampler works with hierarchical models
- [ ] Parameter chains have correct structure

### 15.2 Performance Requirements

- [ ] No performance regression for single-system models
- [ ] Type-stable generated functions
- [ ] Compile times remain reasonable (<10s for typical models)

### 15.3 Code Quality Requirements

- [ ] All functions have docstrings
- [ ] Unit tests cover core functionality
- [ ] Integration tests demonstrate examples
- [ ] Code follows existing style conventions
- [ ] No breaking changes to public API

---

## 16. References

### 16.1 Related Code

- `src/variables.jl` - Core type definitions and code generation
- `src/likelihoods/system.jl` - Likelihood evaluation
- `src/logdensitymodel.jl` - Model interface
- `ext/OctofitterPigeonsExt/` - Sampler integration

### 16.2 Design Patterns

- Metaprogramming with `@RuntimeGeneratedFunction` for performance
- Nested named tuples for parameter structure
- Let-binding for variable scoping in derived variables
- Type parameters for compile-time specialization

---

**End of Design Document**
