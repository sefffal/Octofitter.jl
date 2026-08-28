# ---------------------------------------------------
# Precomputed log-likelihood maps
#
# Unlike `ImageObs`, this one stays single-target: a logL surface is the
# output of somebody else's fit to one companion, so a map is
# one-companion-by-construction and combining several is the user's problem
# (add another `LogLikelihoodMapObs`).
#
# What changed from v1 is the front-end only: the epicyclic superposition
# loop and the private copy of the rotate-and-scale block are gone, replaced
# by `pixel_position` (see pixels.jl). The surface itself is read exactly as
# before.
# ---------------------------------------------------

const likemap_cols = (:map, :epoch, :platescale,)

"""
    LogLikelihoodMapObs(data; target, ref, name="likemap", variables=@variables begin end)

One or more maps of log-likelihood versus Δright-ascension and Δdeclination,
computed with some other tool, as a function of where one companion is.

    LogLikelihoodMapObs(likemap_dat; target=b, ref=A, name="GRAVITY")

`data` needs the columns $likemap_cols — one row per map — where `map` is an
`AstroImage` **recentred so that index `[0,0]` is `ref`**, `epoch` is MJD,
and `platescale` maps its pixels to mas.

`target` is the companion the map describes and `ref` the point it is
measured from; both take the usual grammar. The map is interpolated at the
modelled position of `target` relative to `ref` and the value added straight
to the log-density, so whatever normalization, jitter and correlated noise
the upstream tool applied is inherited verbatim.

An optional `fillvalue` column gives the value used outside the map and
wherever it is not finite; the minimum finite value of each map is used by
default.

# Variables
  - `platescale` — multiplier on the plate scale of every map [default 1].
  - `northangle` [rad] — rotation of the maps relative to true north
    [default 0].

# Example
```julia
likemap_dat = Table(;
    epoch      = [mjd("2016-01-01"), mjd("2017-01-01")],
    map        = [AstroImages.recenter(map1), AstroImages.recenter(map2)],
    platescale = [19.4, 19.4],
)

LogLikelihoodMapObs(likemap_dat; target=b, ref=A, name="GRAVITY")
```
"""
struct LogLikelihoodMapObs{TTable<:Table,TT,TR} <: Octofitter.AbstractObs
    table::TTable
    priors::Octofitter.Priors
    derived::Octofitter.Derived
    target::TT
    ref::TR
    name::String
end

function LogLikelihoodMapObs(
    observations;
    target,
    ref,
    name::String="likemap",
    variables::Tuple{Octofitter.Priors,Octofitter.Derived}=(Octofitter.@variables begin end)
)
    (priors, derived) = variables
    table = Table(observations)
    if !issubset(likemap_cols, Tables.columnnames(table))
        error("Expected columns $likemap_cols")
    end
    if any(>=(mjd("2050")), table.epoch) || any(<=(mjd("1950")), table.epoch)
        @warn "Epochs fell outside 1950–2050; the expected format is MJD. Double-check your input."
    end
    if !in(:fillvalue, Tables.columnnames(table))
        @info ":fillvalue not provided. Filling Inf/NaN with minimum log-like value"
        fillvalue = map(eachrow(table)) do row
            img = row[].map
            m = minimum(filter(isfinite, img))
            @info "row:" epoch = row[].epoch fillvalue = m
            return m
        end
        table = Table(table; fillvalue)
    end
    # Non-finite cells are replaced by the fill value *before* interpolation,
    # and the same value extrapolates outside the map, so a modelled position
    # anywhere on the plane gets a finite answer — a companion the map has
    # nothing to say about is neither favoured nor ruled out beyond whatever
    # floor the upstream tool's map already implies.
    mapinterp = map(eachrow(table)) do row
        img = row[].map
        data = collect(img)
        data[.!isfinite.(data)] .= row[].fillvalue
        linear_interpolation(parent.(dims(img)), data, extrapolation_bc=convert(eltype(img), row[].fillvalue))
    end
    table = Table(table; mapinterp)

    t, r = Octofitter.refspec(target), Octofitter.refspec(ref)
    return LogLikelihoodMapObs{typeof(table),typeof(t),typeof(r)}(
        table, priors, derived, t, r, String(name))
end

export LogLikelihoodMapObs

Octofitter.refspecs(obs::LogLikelihoodMapObs) = (obs.target, obs.ref)

Octofitter.likeobj_from_epoch_subset(obs::LogLikelihoodMapObs, inds) = LogLikelihoodMapObs(
    obs.table[inds]; target=obs.target, ref=obs.ref, name=obs.name,
    variables=(obs.priors, obs.derived))

"""
    simulate(obs::LogLikelihoodMapObs, ctx) -> (; x, y, epochs)

Where the target is predicted to land in each map, in that map's own pixel
coordinates. The plate scale and north angle are applied — the model is
carried onto the map's grid, never the other way around.
"""
function Octofitter.simulate(obs::LogLikelihoodMapObs, ctx::Octofitter.ObsContext)
    T = Octofitter._system_number_type(ctx.θ_system)
    platescale_mult, northangle = Octofitter.sky_calibration(ctx)
    target = Octofitter.ref(ctx, obs.target)
    reference = Octofitter.ref(ctx, obs.ref)
    L = length(obs.table.epoch)
    x = Vector{T}(undef, L)
    y = Vector{T}(undef, L)
    for i in 1:L
        sol = Octofitter.solutionat(ctx, i)
        ps = obs.table.platescale[i] * platescale_mult
        x[i], y[i] = pixel_position(sol, target, reference, ps, northangle)
    end
    return (; x, y, epochs=obs.table.epoch)
end

"""
    ln_like(obs::LogLikelihoodMapObs, ctx)

Sum of each map, interpolated at the modelled position of `target` relative
to `ref`.
"""
function Octofitter.ln_like(obs::LogLikelihoodMapObs, ctx::Octofitter.ObsContext)
    T = Octofitter._system_number_type(ctx.θ_system)
    tbl = obs.table
    platescale_mult, northangle = Octofitter.sky_calibration(ctx)
    target = Octofitter.ref(ctx, obs.target)
    reference = Octofitter.ref(ctx, obs.ref)

    ll = zero(T)
    for i in eachindex(tbl.epoch)
        sol = Octofitter.solutionat(ctx, i)
        ps = tbl.platescale[i] * platescale_mult
        x, y = pixel_position(sol, target, reference, ps, northangle)
        # The constructor already replaced every non-finite cell and set the
        # extrapolation boundary, so this is finite by construction; the
        # guard is for a map handed in with a non-finite `fillvalue`. (v1
        # wrote `if !isfinite(χ²) ll += fillvalue end; ll += χ²`, which added
        # the NaN as well as the fill value — dead code, since the same
        # constructor made the branch unreachable.)
        χ² = tbl.mapinterp[i](x, y)
        ll += isfinite(χ²) ? χ² : tbl.fillvalue[i]
    end
    return ll
end

"""
    generate_from_params(obs::LogLikelihoodMapObs, ctx; add_noise)

A replicate set of maps, each the measured surface translated so that its
peak sits at the modelled position of `target`.

A log-likelihood surface is a *reduction output*, not a measurement with a
noise model: there is no sampling distribution to draw a new one from. What
can be done — and what this does — is preserve the surface's measured shape
and put it where the model says the companion is, which is what a
posterior-predictive overlay of a localization map wants. `add_noise=true`
therefore has nothing to add and is an error rather than a silent no-op.
"""
function Octofitter.generate_from_params(obs::LogLikelihoodMapObs, ctx::Octofitter.ObsContext; add_noise)
    add_noise && error(
        "`generate_from_params(::LogLikelihoodMapObs; add_noise=true)` is not " *
        "available: a precomputed log-likelihood surface has no sampling " *
        "distribution to draw a replicate from. Pass add_noise=false for a " *
        "noiseless replicate (the measured surface, recentred on the modelled " *
        "position), or simulate the underlying data with the observation type " *
        "the map was reduced from.")
    sim = Octofitter.simulate(obs, ctx)
    maps = map(eachindex(obs.table.epoch)) do i
        img = obs.table.map[i]
        itp = obs.table.mapinterp[i]
        xs, ys = parent.(dims(img))
        # Peak of the surface as it stands, in map coordinates…
        data = collect(img)
        data[.!isfinite.(data)] .= obs.table.fillvalue[i]
        ipk = argmax(data)
        Δx = sim.x[i] - xs[ipk[1]]
        Δy = sim.y[i] - ys[ipk[2]]
        # …and the same surface, shifted so that peak lands on the model.
        out = copy(img)
        for (jx, x) in enumerate(xs), (jy, y) in enumerate(ys)
            out[jx, jy] = convert(eltype(img), itp(x - Δx, y - Δy))
        end
        return out
    end
    return LogLikelihoodMapObs(Table(obs.table; map=maps); target=obs.target, ref=obs.ref,
        name=obs.name, variables=(obs.priors, obs.derived))
end

# ---------------------------------------------------
# Plotting protocol
#
# A log-likelihood surface is a reduction output, not a measurement: there is
# no (data, model, σ) triple to be had, and no sampling distribution behind it.
# What *is* well defined per epoch is how far below that map's own maximum the
# modelled position falls — zero when the orbit passes exactly through the
# localization peak, growing as it misses. That is the honest 1-D reduction,
# and it is what the bespoke panel draws; there is no generic time-series
# panel here, because "residual over σ" would invent a σ the map never
# supplied.
# ---------------------------------------------------

Octofitter.plotchannels(::LogLikelihoodMapObs) = (
    Octofitter.PlotChannel(:lnlike, "log-likelihood below map peak", ""),
)

# The largest finite value of each map — the best any orbit could score there.
# Cached per observation: the maps do not change, and `residuals` is called
# once per posterior draw.
const _LIKEMAP_PEAKS = Base.IdDict{Any,Vector{Float64}}()
const _LIKEMAP_PEAKS_LOCK = ReentrantLock()

function _likemap_peaks(obs::LogLikelihoodMapObs)
    lock(_LIKEMAP_PEAKS_LOCK) do
        get!(_LIKEMAP_PEAKS, obs) do
            map(eachindex(obs.table.epoch)) do i
                data = collect(obs.table.map[i])
                fin = filter(isfinite, vec(data))
                isempty(fin) ? Float64(obs.table.fillvalue[i]) : Float64(maximum(fin))
            end
        end
    end
end

function Octofitter.residuals(obs::LogLikelihoodMapObs, ctx::Octofitter.ObsContext)
    sim = Octofitter.simulate(obs, ctx)
    tbl = obs.table
    peaks = _likemap_peaks(obs)
    L = length(tbl.epoch)
    model = Vector{Float64}(undef, L)
    for i in 1:L
        χ² = tbl.mapinterp[i](Float64(sim.x[i]), Float64(sim.y[i]))
        model[i] = isfinite(χ²) ? Float64(χ²) : Float64(tbl.fillvalue[i])
    end
    data = collect(Float64, peaks)
    # σ = 1 is a unit, not an uncertainty: it keeps the protocol's shape while
    # saying plainly that the map carries no error bar of its own.
    return (;
        lnlike=(; epoch=collect(Float64, tbl.epoch), data, model,
            resid=data .- model, σ=ones(L), σ_eff=ones(L), use=trues(L)),
    )
end

Octofitter.defaultpanels(obs::LogLikelihoodMapObs) =
    (:likemap => (gp, series) -> Octofitter.likemappanel!(gp, series, obs),)
