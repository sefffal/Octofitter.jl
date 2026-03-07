#=
Gaia Non-Single Star (NSS) catalog integration for Octofitter.

Provides functions to query NSS orbital solutions and convert them into
starting points for Octofitter models via `initialize!`.

The NSS solutions are used as starting points only (not priors), because
NSS uncertainties may be overly optimistic.
=#

using HTTP, CSV

# ──────────────────────────────────────────────────────────────────────
# NSS catalog query
# ──────────────────────────────────────────────────────────────────────

"""
    query_nss(; gaia_id, catalog=:dr3)

Query the Gaia Non-Single Star (NSS) two-body orbit table for a given source ID.
Returns a named tuple of the NSS solution columns, or `nothing` if no solution is found.

Results are cached locally in `_gaia_nss_{catalog}/` directories.

# Arguments
- `gaia_id`: Gaia source ID (integer)
- `catalog`: `:dr3` or `:dr4` (default `:dr3`)
"""
function query_nss(; gaia_id, catalog=:dr3)
    if catalog == :dr3
        table_name = "gaiadr3.nss_two_body_orbit"
        cache_dir = "_gaia_nss_dr3"
    elseif catalog == :dr4
        table_name = "gaiadr4.nss_two_body_orbit"
        cache_dir = "_gaia_nss_dr4"
    else
        error("Unsupported catalog: $catalog. Use :dr3 or :dr4.")
    end

    fname = joinpath(cache_dir, "source-$gaia_id.csv")
    if !isfile(fname)
        @info "Querying NSS catalog at gea.esac.esa.int/tap-server" source_id=gaia_id catalog
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST" => "doQuery",
                "LANG" => "ADQL",
                "FORMAT" => "CSV",
                "QUERY" => "SELECT * FROM $table_name WHERE source_id=$gaia_id"
            ],
            cookies=false,
        )
        if resp.status != 200
            error("Error querying NSS catalog: HTTP $(resp.status)")
        end
        if !isdir(cache_dir)
            mkpath(cache_dir)
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end

    lines = split(strip(buf), "\n")
    if length(lines) < 2
        @warn "No NSS solution found for source_id=$gaia_id in $catalog"
        return nothing
    end

    # Parse CSV header and first data row
    header_line = lines[2]  # skip potential comment lines
    # Use CSV to parse properly (handles quoted fields etc.)
    tbl = CSV.read(IOBuffer(buf), Tables.columntable, normalizenames=true)

    if Tables.rowcount(tbl) == 0
        @warn "No NSS solution found for source_id=$gaia_id in $catalog"
        return nothing
    end

    # If multiple solutions exist, take the first one
    if Tables.rowcount(tbl) > 1
        @info "Multiple NSS solutions found for source_id=$gaia_id; using the first row"
    end

    # Build named tuple from first row
    cols = Tables.columnnames(tbl)
    vals = [Tables.getcolumn(tbl, col)[1] for col in cols]
    return NamedTuple{Tuple(cols)}(Tuple(vals))
end
export query_nss


# ──────────────────────────────────────────────────────────────────────
# NSS → Octofitter parameter mapping
# ──────────────────────────────────────────────────────────────────────

# NSS periastron reference epoch: JD = t_periastron + 2457389.0 (DR3 reference)
const _NSS_DR3_T_PERIASTRON_REF_JD = 2457389.0

"""
    nss_to_starting_point(nss_sol, model; planet_key=:b)

Convert an NSS orbital solution (as returned by [`query_nss`](@ref)) into a
named tuple suitable for passing to [`initialize!`](@ref) as fixed starting parameters.

This function inspects the model's planet parameters and maps NSS values to
whichever parameterization the user has chosen (Thiele-Innes or Campbell).

# Mapped parameters

For **Thiele-Innes** models (planet has `A`, `B`, `F`, `G`):
- `A`, `B`, `F`, `G` from NSS `a_thiele_innes`, `b_thiele_innes`, etc.

For **Campbell** models (planet has `a` or `P`, and `i`, `Ω`, `ω`):
- Converts NSS Thiele-Innes constants to Campbell elements
- Sets `a` (semi-major axis in AU) or `P` (period in days)
- Sets `i`, `Ω`, `ω`

Common parameters mapped in both cases:
- `e` (eccentricity)
- `tp` (periastron time in MJD, if the planet has `tp`)

# Arguments
- `nss_sol`: Named tuple from [`query_nss`](@ref)
- `model`: A `LogDensityModel`
- `planet_key`: Symbol identifying which planet to set (default `:b`)

# Returns
A named tuple like `(; planets=(; b=(; e=0.3, A=5.2, ...)))` ready for `initialize!`.
"""
function nss_to_starting_point(nss_sol, model; planet_key=:b)
    # Sample from the model to inspect the parameter structure
    sample = model.sample_priors(Random.default_rng())
    params_nt = model.arr2nt(sample)

    if !hasproperty(params_nt, :planets) || !hasproperty(params_nt.planets, planet_key)
        error("Model does not have a planet named :$planet_key")
    end

    planet_params = params_nt.planets[planet_key]
    # Only consider free (prior) parameters, not derived ones.
    # Derived parameters cannot be set via initialize! and will cause errors.
    planet_obj = nothing
    for p in model.system.planets
        if p.name == planet_key
            planet_obj = p
            break
        end
    end
    if isnothing(planet_obj)
        error("Could not find planet :$planet_key in system")
    end
    planet_param_names = Tuple(collect(keys(planet_obj.priors.priors)))

    mapped = Dict{Symbol, Any}()
    unmapped = String[]

    # Helper: set an angular parameter, handling UniformCircular parameterization.
    # If the model uses `Ω ~ UniformCircular()`, the free params are `Ωx` and `Ωy`,
    # not `Ω` itself. We detect this and set x=cos(angle), y=sin(angle).
    function _try_set_angle!(mapped, name::Symbol, value_rad)
        if name in planet_param_names
            # Direct free parameter (e.g. Ω ~ Uniform(0, 2π))
            mapped[name] = value_rad
            return true
        end
        namex = Symbol("$(name)x")
        namey = Symbol("$(name)y")
        if namex in planet_param_names && namey in planet_param_names
            # UniformCircular parameterization: set x,y components
            mapped[namex] = cos(value_rad)
            mapped[namey] = sin(value_rad)
            return true
        end
        return false
    end

    # --- Eccentricity ---
    if :e in planet_param_names && hasproperty(nss_sol, :eccentricity)
        e_val = nss_sol.eccentricity
        if !ismissing(e_val) && !isnothing(e_val) && isfinite(Float64(e_val))
            mapped[:e] = Float64(e_val)
        end
    end

    # --- Period ---
    nss_period = nothing
    if hasproperty(nss_sol, :period) && !ismissing(nss_sol.period) && !isnothing(nss_sol.period)
        nss_period = Float64(nss_sol.period)  # days
    end

    # --- Periastron time ---
    # tp is usually derived from θ ~ UniformCircular() in Octofitter models.
    # If tp is a direct free parameter, set it. Otherwise, we skip it and let
    # the optimizer find a good θ.
    nss_tp_mjd = nothing
    if hasproperty(nss_sol, :t_periastron) && !ismissing(nss_sol.t_periastron) && !isnothing(nss_sol.t_periastron)
        nss_tp_jd = Float64(nss_sol.t_periastron) + _NSS_DR3_T_PERIASTRON_REF_JD
        nss_tp_mjd = jd2mjd(nss_tp_jd)
    end

    if :tp in planet_param_names && !isnothing(nss_tp_mjd)
        mapped[:tp] = nss_tp_mjd
    end

    # --- Check for Thiele-Innes parameterization ---
    has_ti = all(s -> s in planet_param_names, (:A, :B, :F, :G))

    # --- Check for Campbell parameterization ---
    # Campbell angles may be free directly (i, Ω, ω) or via UniformCircular (Ωx/Ωy, ωx/ωy)
    function _has_angle(name::Symbol)
        name in planet_param_names ||
        (Symbol("$(name)x") in planet_param_names && Symbol("$(name)y") in planet_param_names)
    end
    has_campbell_angles = _has_angle(:i) && _has_angle(:Ω) && _has_angle(:ω)

    # Extract NSS Thiele-Innes constants
    nss_A = _nss_get_float(nss_sol, :a_thiele_innes)
    nss_B = _nss_get_float(nss_sol, :b_thiele_innes)
    nss_F = _nss_get_float(nss_sol, :f_thiele_innes)
    nss_G = _nss_get_float(nss_sol, :g_thiele_innes)
    has_nss_ti = !isnothing(nss_A) && !isnothing(nss_B) && !isnothing(nss_F) && !isnothing(nss_G)

    if has_ti && has_nss_ti
        # Direct Thiele-Innes mapping
        mapped[:A] = nss_A
        mapped[:B] = nss_B
        mapped[:F] = nss_F
        mapped[:G] = nss_G
        @info "Mapped NSS Thiele-Innes constants (A, B, F, G) to model"
    elseif has_campbell_angles && has_nss_ti
        # Convert Thiele-Innes → Campbell
        i, Ω, ω, α_mas = _ti_to_campbell(nss_A, nss_B, nss_F, nss_G)

        _try_set_angle!(mapped, :i, i)
        _try_set_angle!(mapped, :Ω, Ω)
        _try_set_angle!(mapped, :ω, ω)

        @info "Converted NSS Thiele-Innes to Campbell elements" i_deg=rad2deg(i) Ω_deg=rad2deg(Ω) ω_deg=rad2deg(ω)
    elseif !has_nss_ti
        push!(unmapped, "angular elements (NSS solution has no Thiele-Innes constants)")
    else
        push!(unmapped, "orbital geometry (model uses neither TI nor Campbell parameterization)")
    end

    # --- Semi-major axis or period ---
    if :a in planet_param_names && !isnothing(nss_period)
        M_rough = _estimate_system_mass(model)
        P_years = nss_period / PlanetOrbits.year2day_julian
        a_au = ∛(M_rough * P_years^2)
        mapped[:a] = a_au
        @info "Mapped NSS period to semi-major axis" P_days=nss_period a_AU=a_au M_est=M_rough
    elseif :P in planet_param_names && !isnothing(nss_period)
        mapped[:P] = nss_period
        @info "Mapped NSS period to model" P_days=nss_period
    elseif !isnothing(nss_period) && (:a in planet_param_names || :P in planet_param_names)
        # Already handled above
    elseif !isnothing(nss_period)
        push!(unmapped, "period (model planet has neither free `a` nor `P`)")
    end

    if isempty(mapped)
        @warn "Could not map any NSS parameters to model planet :$planet_key"
        return (;)
    end

    if !isempty(unmapped)
        @info "Some NSS parameters could not be mapped" unmapped
    end

    planet_nt = NamedTuple{Tuple(collect(keys(mapped)))}(Tuple(collect(values(mapped))))
    return (; planets=NamedTuple{(planet_key,)}((planet_nt,)))
end
export nss_to_starting_point


"""
    initialize_from_nss!(model; gaia_id, planet_key=:b, catalog=:dr3, kwargs...)

Convenience function that queries the NSS catalog for `gaia_id`, converts the
orbital solution to model parameters, and calls `initialize!` with those as
starting points.

The NSS parameters are used to anchor the global optimization search, but are
**not** used as priors. All model parameters remain free during sampling.

# Example
```julia
model = Octofitter.LogDensityModel(sys)
chain = initialize_from_nss!(model; gaia_id=4295745059252873600, planet_key=:b)
```

Any additional keyword arguments are forwarded to `initialize!`.
"""
function initialize_from_nss!(model; gaia_id, planet_key=:b, catalog=:dr3, kwargs...)
    nss_sol = query_nss(; gaia_id, catalog)
    if isnothing(nss_sol)
        @warn "No NSS solution found; falling back to default initialization"
        return initialize!(model; kwargs...)
    end

    fixed_params = nss_to_starting_point(nss_sol, model; planet_key)
    if isempty(propertynames(fixed_params))
        @warn "Could not map any NSS parameters; falling back to default initialization"
        return initialize!(model; kwargs...)
    end

    @info "Initializing from NSS solution" gaia_id catalog mapped_params=propertynames(fixed_params.planets[planet_key])
    return initialize!(model, fixed_params; kwargs...)
end
export initialize_from_nss!


# ──────────────────────────────────────────────────────────────────────
# NSS → Model + Chain for comparison plotting
# ──────────────────────────────────────────────────────────────────────

"""
    nss_model, nss_chain = nss_to_model_chain(nss_sol; plx=nothing, gaia_id=nothing, N=10_000)

Build a minimal Octofitter model and chain from an NSS orbital solution, suitable
for comparison plotting with your own posterior (e.g. via PairPlots or `octoplot`).

The returned `nss_chain` contains `N` draws from Normal distributions centred on the
NSS best-fit values with the NSS-reported uncertainties. The returned `nss_model` is
a one-planet `LogDensityModel` using a `ThieleInnesOrbit` parameterization.

The total system mass is derived automatically from the NSS period, Thiele-Innes
constants, and parallax via Kepler's third law, so you do not need to provide it.
Parallax is taken from the NSS table if available, otherwise looked up from Gaia DR3.

!!! warning "NSS uncertainties"
    NSS error bars may be overly optimistic. Use the returned chain for visual
    comparison only, not as a prior or ground truth.

# Arguments
- `nss_sol`: Named tuple from [`query_nss`](@ref)
- `plx`: Parallax in mas. If `nothing`, uses the NSS table value or queries Gaia DR3.
- `gaia_id`: Gaia source ID (used to look up parallax if not in `nss_sol`)
- `N`: Number of draws (default 10_000)

# Returns
`(nss_model, nss_chain)` — a `LogDensityModel` and `MCMCChains.Chains` object.

# Example
```julia
nss_sol = query_nss(gaia_id=4295745059252873600)
nss_model, nss_chain = nss_to_model_chain(nss_sol)

# Compare with your posterior in a pair plot
using PairPlots
pairplot(
    "Posterior" => chain,
    "NSS"       => nss_chain,
)

# Or overlay on octoplot
octoplot(nss_model, nss_chain)
```
"""
function nss_to_model_chain(nss_sol; plx=nothing, gaia_id=nothing, N=10_000)

    # Resolve parallax: try NSS table → Gaia DR3 query → error
    if isnothing(plx)
        nss_plx = _nss_get_float(nss_sol, :parallax)
        if !isnothing(nss_plx) && nss_plx > 0
            plx = nss_plx
        else
            src_id = !isnothing(gaia_id) ? gaia_id :
                     (hasproperty(nss_sol, :source_id) ? nss_sol.source_id : nothing)
            if isnothing(src_id)
                error("No parallax found in NSS solution. Provide `plx` or `gaia_id`.")
            end
            gaia_cat = _query_gaia_dr3(gaia_id=Int(src_id))
            plx = Float64(gaia_cat.parallax)
            @info "Looked up parallax from DR3" plx
        end
    end

    # Extract NSS values and errors
    function _val_err(sol, key_val, key_err)
        val = _nss_get_float(sol, key_val)
        err = _nss_get_float(sol, key_err)
        return val, err
    end

    e_val, e_err       = _val_err(nss_sol, :eccentricity, :eccentricity_error)
    P_val, P_err       = _val_err(nss_sol, :period, :period_error)
    tp_val, tp_err     = _val_err(nss_sol, :t_periastron, :t_periastron_error)
    A_val, A_err       = _val_err(nss_sol, :a_thiele_innes, :a_thiele_innes_error)
    B_val, B_err       = _val_err(nss_sol, :b_thiele_innes, :b_thiele_innes_error)
    F_val, F_err       = _val_err(nss_sol, :f_thiele_innes, :f_thiele_innes_error)
    G_val, G_err       = _val_err(nss_sol, :g_thiele_innes, :g_thiele_innes_error)

    has_ti = !isnothing(A_val) && !isnothing(B_val) && !isnothing(F_val) && !isnothing(G_val)
    if !has_ti
        error("NSS solution does not contain Thiele-Innes constants (a/b/f/g_thiele_innes). Cannot build model.")
    end

    # Convert tp to MJD
    tp_mjd = nothing
    if !isnothing(tp_val)
        tp_mjd = jd2mjd(tp_val + _NSS_DR3_T_PERIASTRON_REF_JD)
    end

    # Derive total system mass from NSS period + TI constants + parallax
    # via Kepler's third law: M = a³/P² (a in AU, P in years, M in Msol)
    if isnothing(P_val)
        error("NSS solution does not contain a period. Cannot derive system mass.")
    end
    u = (A_val^2 + B_val^2 + F_val^2 + G_val^2) / 2
    v = A_val * G_val - B_val * F_val
    α_mas = sqrt(u + sqrt((u + v) * (u - v)))
    a_AU = α_mas / plx
    P_years = P_val / PlanetOrbits.year2day_julian
    M = a_AU^3 / P_years^2

    # Build priors: Normal(value, error) for each parameter.
    # If error is missing, use a tiny spread (effectively fixed).
    _make_prior(val, err) = isnothing(err) || err ≤ 0 ? Normal(val, abs(val) * 1e-6 + 1e-10) : Normal(val, err)

    e_prior  = truncated(_make_prior(e_val, e_err), lower=0, upper=0.999)
    A_prior  = _make_prior(A_val, A_err)
    B_prior  = _make_prior(B_val, B_err)
    F_prior  = _make_prior(F_val, F_err)
    G_prior  = _make_prior(G_val, G_err)
    tp_prior = _make_prior(tp_mjd, isnothing(tp_err) ? nothing : tp_err)

    planet_b = Planet(
        name="b",
        basis=ThieleInnesOrbit,
        observations=[],
        variables=@variables begin
            e  ~ e_prior
            A  ~ A_prior
            B  ~ B_prior
            F  ~ F_prior
            G  ~ G_prior
            tp ~ tp_prior
        end
    )

    # Build system variables directly (can't use @variables with $ for derived
    # values inside a module-level function during precompilation)
    sys_priors = Priors(OrderedDict{Symbol,Distribution}())
    sys_derived = Derived(OrderedDict{Symbol,Any}(:M => M, :plx => plx), (), ())
    sys = System(
        name="NSS",
        companions=[planet_b],
        observations=[],
        variables=(sys_priors, sys_derived)
    )

    nss_model = LogDensityModel(sys, verbosity=0)

    # Draw N samples from priors and format as a chain
    rng = Random.Xoshiro(0)
    samples_nt = map(1:N) do _
        θ = nss_model.sample_priors(rng)
        nt = nss_model.arr2nt(θ)
        (; logpost=0.0, nt...)
    end

    nss_chain = result2mcmcchain(samples_nt, Dict(:internals => [:logpost]))

    return nss_model, nss_chain
end
export nss_to_model_chain


# ──────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────

"""Extract a Float64 value from an NSS named tuple, returning nothing if missing/NaN."""
function _nss_get_float(nss_sol, key::Symbol)
    if !hasproperty(nss_sol, key)
        return nothing
    end
    val = getproperty(nss_sol, key)
    if ismissing(val) || isnothing(val)
        return nothing
    end
    fval = Float64(val)
    if !isfinite(fval)
        return nothing
    end
    return fval
end

"""
Convert Thiele-Innes constants (A, B, F, G) in mas to Campbell elements (i, Ω, ω, α_mas).

Returns (i, Ω, ω, α_mas) where angles are in radians and α_mas is the angular
semi-major axis in milliarcseconds.

Uses the standard inversion formulas (e.g., Binnendijk 1960, Wright & Howard 2009).
"""
function _ti_to_campbell(A, B, F, G)
    # Intermediate quantities
    u = (A^2 + B^2 + F^2 + G^2) / 2
    v = A * G - B * F

    # Angular semi-major axis
    α = sqrt(u + sqrt((u + v) * (u - v)))

    # Inclination: cos(i) = v / (α²) but need to be careful with signs
    # v = α² * cos(i), so cos_i = v / α²
    cos_i = v / α^2
    cos_i = clamp(cos_i, -1.0, 1.0)
    i = acos(cos_i)

    # Ascending node and argument of periastron
    # From A, B, F, G definitions:
    # A = α * (cos(ω)cos(Ω) - sin(ω)sin(Ω)cos(i))
    # B = α * (cos(ω)sin(Ω) + sin(ω)cos(Ω)cos(i))
    # F = α * (-sin(ω)cos(Ω) - cos(ω)sin(Ω)cos(i))
    # G = α * (-sin(ω)sin(Ω) + cos(ω)cos(Ω)cos(i))
    #
    # A + G = α * (1 + cos(i)) * cos(ω + Ω)
    # A - G = α * (1 - cos(i)) * cos(ω - Ω)
    # B - F = α * (1 + cos(i)) * sin(ω + Ω)
    # B + F = α * (1 - cos(i)) * (-sin(ω - Ω))

    ω_plus_Ω = atan(B - F, A + G)
    ω_minus_Ω = atan(-(B + F), A - G)

    ω = (ω_plus_Ω + ω_minus_Ω) / 2
    Ω = (ω_plus_Ω - ω_minus_Ω) / 2

    # Normalize to [0, 2π)
    ω = mod(ω, 2π)
    Ω = mod(Ω, 2π)

    return i, Ω, ω, α
end

"""
Estimate system mass from model priors by sampling and taking the median.
"""
function _estimate_system_mass(model)
    rng = Random.Xoshiro(42)
    masses = Float64[]
    for _ in 1:1000
        sample = model.sample_priors(rng)
        nt = model.arr2nt(sample)
        if hasproperty(nt, :M)
            push!(masses, Float64(nt.M))
        end
    end
    if isempty(masses)
        @warn "No system mass parameter found; assuming 1.0 Msol"
        return 1.0
    end
    return median(masses)
end
