# ---------------------------------------------------
# orbitize! / HDF5 interoperability
#
# Reading and writing posteriors in the HDF5 layout used by orbitize! and by
# whereistheplanet.com, plus the whereistheplanet astrometry loader.
#
# The orbital-element conventions are the content of this file and are
# untouched: orbitize!'s standard basis is (sma, ecc, inc, aop, pan, tau, plx,
# mtot) with the same angle definitions Octofitter uses, which is why the
# conversion is a rename plus a τ↔tp change of phase variable.
#
# Two conventions on the Octofitter side that the conversion has to respect:
#
#   - **Masses are M⊙ throughout**, which is also what orbitize! stores, so
#     `m0`/`m1` land in `A_mass` / `b_mass` unscaled.
#
#   - **A body's dynamical mass comes from the hierarchy**, so there is no
#     per-companion total-mass column to synthesise. The total needed to turn
#     orbitize!'s `tau` into a `tp` is taken as `mtot` when the file supplies
#     it, and `m0 + m_i` otherwise — orbitize!'s own standard-basis definition
#     for the i-th companion.
# ---------------------------------------------------

using HDF5
using StringDistances
using MCMCChains: MCMCChains

"""
    Whereistheplanet_search("targetname")

Search for an orbit posterior and/or astrometry hosted on whereistheplanet.com by a given target name.
If not found, a list of similar target names will be reported.
"""
function Whereistheplanet_search(target, catalog=datadep"Whereistheplanet")

    dirpath = joinpath(catalog, "whereistheplanet-master", "data")
    fnames = readdir(dirpath, join=true)
    avail_targets = map(fnames) do fname
        m = match(r"post_(.+)\.hdf5", fname)
        return !isnothing(m) ? m.captures[1] : ""
    end
    fname_matched_i = findfirst(==(target), avail_targets)

    if isnothing(fname_matched_i)
        avail_filt = avail_targets[avail_targets.!=""]
        similarity = evaluate.(Ref(Levenshtein()), target, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3, end)]]
        @error "No results were found for the target $target."
        @info "Here are a list of similar and available target names" closest_3
        error()
    end

    return fnames[fname_matched_i]

end

"""
    Whereistheplanet_astrom("targetname"; target=:b, ref=:A)

Load relative astrometry hosted on whereistheplanet.com by a given target
name, as a vector of [`RelAstromObs`](@ref). If the name is not found, a list
of similar target names is reported.

orbitize! stores separation/position-angle and RA/Dec rows in one table;
Octofitter needs one likelihood object per format, so up to two are returned:

    seppa, radec = Octofitter.Whereistheplanet_astrom("51erib"; target=b, ref=A)

`target` and `ref` are the model references the astrometry measures — the
companion and the host in the usual case. They take the full v9 grammar
(a `Body`, a `Symbol`, `Barycentre(…)`, `Photocentre(…)`), because
whereistheplanet's "object 1" is only a companion by convention.

The returned objects are named `"<name>_seppa"` and `"<name>_radec"`; both can
go into the same `System`, and their names must differ for that to be legal.
"""
function Whereistheplanet_astrom(targetname, catalog=datadep"Whereistheplanet";
                                 object=1, target=:b, ref=:A,
                                 name="whereistheplanet")

    fname = Whereistheplanet_search(targetname, catalog)
    return h5open(fname, "r") do f
        records = read(f["data"])

        records = filter(row -> row.object == object, records)

        # Group observations by type
        seppa = filter(row -> row.quant_type == "seppa", records)
        radec = filter(row -> row.quant_type == "radec", records)

        out = RelAstromObs[]
        if length(seppa) > 0
            dat = map(seppa) do row
                cor = row.quant12_corr
                if !isfinite(cor)
                    cor = 0.0
                end
                (; row.epoch, sep=row.quant1, σ_sep=row.quant1_err,
                   pa=deg2rad(row.quant2), σ_pa=deg2rad(row.quant2_err), cor)
            end
            push!(out, RelAstromObs(dat; target, ref, name="$(name)_seppa"))
        end
        if length(radec) > 0
            dat = map(radec) do row
                cor = row.quant12_corr
                if !isfinite(cor)
                    cor = 0.0
                end
                (; row.epoch, ra=row.quant1, σ_ra=row.quant1_err,
                   dec=row.quant2, σ_dec=row.quant2_err, cor)
            end
            push!(out, RelAstromObs(dat; target, ref, name="$(name)_radec"))
        end
        # v1 built this vector and then never returned it — the `return out`
        # was commented out, so the documented
        # `astrom1, astrom2 = Whereistheplanet_astrom(…)` destructured
        # `nothing`. It returns now.
        return out
    end
end

# orbitize!'s standard basis, per companion index, → the v2 chain spelling.
# `tau` is kept alongside the `tp` derived from it: it is the quantity actually
# stored, and dropping it would make the import lossy.
const _ORBITIZE_ELEMENTS = ("sma" => "a", "ecc" => "e", "inc" => "i",
                            "aop" => "ω", "pan" => "Ω", "tau" => "τ")

"""
    loadhdf5("fname.h5")
    loadhdf5("51erib")
    loadhdf5(fname, numchains; host=:A, bodynames=("b","c","d","e"))

Load an orbitize! posterior from an HDF5 file and convert it into an
Octofitter chain. Both tools use the same orbital-element conventions, so this
is a rename plus one change of phase variable (orbitize!'s `tau`, a fraction
of a period past a reference epoch, becomes Octofitter's `tp`).

Passing a name rather than a filename looks the target up on
whereistheplanet.com.

`numchains` interprets the stored array as that many chains concatenated
together.

`host` and `bodynames` give the v9 body names the columns are written under —
`m0` becomes `<host>_mass` and `sma1, ecc1, …` become `<bodynames[1]>_a,
<bodynames[1]>_e, …`. Masses are solar masses, as orbitize! stores them and as
v9 uses throughout.
"""
function loadhdf5(fname_or_targetname, numchains=1; colnames=nothing,
                  host::Union{Symbol,AbstractString}=:A,
                  bodynames=("b", "c", "d", "e"))
    if !(occursin(".hdf5", fname_or_targetname) || occursin(".h5", fname_or_targetname))
        fname = Whereistheplanet_search(fname_or_targetname)
    else
        fname = fname_or_targetname
    end
    host = Symbol(host)
    planet_keys = String.(collect(bodynames))
    return h5open(fname, "r") do f
        # Standard orbitize basis assumed: semi-major axis (sma), eccentricity
        # (ecc), inclination (inc), argument of periastron (aop), position
        # angle of the nodes (pan), epoch of periastron expressed as a fraction
        # of the period past a reference epoch (tau), parallax (plx) and total
        # system mass (mtot).
        arr = transpose(read(f["post"]))
        if isnothing(colnames)
            try
                colnames = read_attribute(f, "parameter_labels")
            catch
                @warn "`parameter_labels` not present, will have to fall back on a guess. You can also provide a vector of colname strings with the argument `colnames`"
                colnames = ["sma1", "ecc1", "inc1", "aop1", "pan1", "tau1", "plx", "mtot"]
            end
        end
        colnames = String.(collect(colnames))
        size(arr, 2) == length(colnames) || error(
            "The posterior in $fname has $(size(arr,2)) columns but $(length(colnames)) " *
            "column names. Pass `colnames=` explicitly.")

        # Reshape to (iterations, parameters, chains) up front and stay there.
        # v1 stacked into 3D only when `numchains > 1` and then `hcat`ed the
        # derived `tp` columns on, which cannot work for a 3D array — so the
        # multi-chain path errored.
        nrow = size(arr, 1)
        numchains >= 1 || error("`numchains` must be at least 1")
        nrow % numchains == 0 || error(
            "The posterior has $nrow samples, which is not divisible into $numchains chains.")
        lenchain = nrow ÷ numchains
        cube = Array{Float64}(undef, lenchain, size(arr, 2), numchains)
        for c in 1:numchains
            cube[:, :, c] .= @view arr[(c-1)*lenchain+1:c*lenchain, :]
        end

        num_planets = count(i -> "sma$i" in colnames, 1:length(planet_keys))
        num_planets > 0 || error(
            "No `sma<i>` column found in $fname; this does not look like an orbitize! " *
            "posterior in the standard basis.")

        renames = Dict{String,String}()
        for i in 1:num_planets, (orb, oct) in _ORBITIZE_ELEMENTS
            renames["$orb$i"] = "$(planet_keys[i])_$oct"
        end
        for i in 1:num_planets
            renames["m$i"] = "$(planet_keys[i])_mass"
        end
        renames["m0"] = "$(host)_mass"
        # `mtot` keeps its name: v2 has no system-level `M`, the masses live on
        # the bodies, and silently renaming a total to a body's mass is exactly
        # the kind of quiet mis-association this refactor is removing.
        outnames = [get(renames, c, c) for c in colnames]

        colidx = Dict(nm => i for (i, nm) in enumerate(outnames))
        getcol(nm) = view(cube, :, colidx[nm], :)

        # Total mass per companion, for the τ → tp conversion. orbitize!'s
        # standard basis defines the i-th companion's two-body total as
        # `m0 + m_i`; `mtot` is what a single-companion fit stores instead.
        function totalmass(i)
            mi = "$(planet_keys[i])_mass"
            if haskey(colidx, "$(host)_mass") && haskey(colidx, mi)
                return getcol("$(host)_mass") .+ getcol(mi)
            elseif haskey(colidx, "mtot")
                return getcol("mtot")
            else
                error("$fname has neither `mtot` nor `m0`+`m$i`, so the total mass needed " *
                      "to convert orbitize!'s `tau` into an epoch of periastron is unknown.")
            end
        end

        tau_ref_epoch = 58849
        if haskey(attrs(f), "tau_ref_epoch")
            tau_ref_epoch = attrs(f)["tau_ref_epoch"]
        end

        extra = Array{Float64}(undef, lenchain, num_planets, numchains)
        extranames = String[]
        for i in 1:num_planets
            a = getcol("$(planet_keys[i])_a")
            τ = getcol("$(planet_keys[i])_τ")
            M = totalmass(i)
            period_days = @. √(a^3 / M) * PlanetOrbits.kepler_year_to_julian_day_conversion_factor
            extra[:, i, :] .= @. τ * period_days + tau_ref_epoch
            push!(extranames, "$(planet_keys[i])_tp")
        end

        chn = MCMCChains.Chains(cat(cube, extra; dims=2),
                                Symbol.([outnames; extranames]))

        # Read additional attributes in and convert to named tuple
        metadata = [Symbol(k) => v for (k, v) in attrs(f)]
        return MCMCChains.setinfo(chn, NamedTuple(metadata))
    end
end

"""
    savehdf5("filename.hdf5", model, chain)
    savehdf5("filename.hdf5", model, chain, :b)

Save an Octofitter chain in the HDF5 layout used by orbitize! and by
whereistheplanet.com.

Only the eight columns of orbitize!'s standard basis are written — `sma`,
`ecc`, `inc`, `aop`, `pan`, `tau`, `plx`, `mtot` — so this exports **one
companion's** visual orbit and nothing else. No data are exported.

The companion defaults to the model's first non-root body. `mtot` is the total
mass of the bodies its orbit binds — `<body>_mass` plus the mass of everything
the model places it about. For the usual `about=A` companion that is
`A_mass + b_mass`; for a Jacobi chain (`about=(A, b)`) it is the sum over the
whole interior, which is the same convention orbitize! uses for its own
multi-planet fits.
"""
function savehdf5(fname::AbstractString, model, chain::Chains,
                  body::Symbol=_first_companion(_system_of(model)))
    sys = _system_of(model)
    interior = _interior_of(sys, body)

    tau_ref_epoch = 58849

    get3(nm) = haskey(chain, nm) ? vec(collect(chain[:, nm, :])) : nothing
    require(nm, what) = begin
        v = get3(nm)
        isnothing(v) && error(
            "The chain has no `$nm` column, which is needed for $what. orbitize! export " *
            "only supports models that carry $what explicitly.")
        v
    end

    # The gravitating mass of the row: the companion plus everything its orbit is
    # about. With a single host that is `A_mass + b_mass`; with a Jacobi chain
    # (`about=(A, b)`) it is the sum over the interior, which is what orbitize!'s
    # own multi-planet basis means by `mtot`.
    m_body = require(Symbol(body, "_mass"), "the companion mass")
    mtot = m_body .+ sum(require(Symbol(h, "_mass"), "the mass of interior body :$h")
                         for h in interior)

    sma = let a = get3(Symbol(body, "_a"))
        if !isnothing(a)
            a
        else
            P = require(Symbol(body, "_P"),
                        "the semi-major axis (neither `$(body)_a` nor `$(body)_P` is present)")
            # Periods are days, matching `PlanetOrbits.period`.
            @. cbrt((P / PlanetOrbits.kepler_year_to_julian_day_conversion_factor)^2 * mtot)
        end
    end
    period_days = @. √(sma^3 / mtot) * PlanetOrbits.kepler_year_to_julian_day_conversion_factor
    # `tp` is only a chain column when the model sampled it. The phase spelling
    # the docs recommend is `θ` + `epoch`, which produces no `tp` — so rather
    # than refusing, rebuild each draw's orbit and read `tp` off the row, which
    # is where the conversion already lives. Costs one system build per draw and
    # only happens on the fallback path.
    tp = let t = get3(Symbol(body, "_tp"))
        isnothing(t) ? _tp_from_chain(model, chain, body) : t
    end
    tau = @. mod((tp - tau_ref_epoch) / period_days, 1)

    ecc = require(Symbol(body, "_e"), "eccentricity")
    inc = require(Symbol(body, "_i"), "inclination")
    aop = require(Symbol(body, "_ω"), "the argument of periastron")
    pan = require(Symbol(body, "_Ω"), "the position angle of the nodes")
    plx = require(:plx, "parallax")

    return h5open(fname, "w") do f
        dat = transpose(hcat(sma, ecc, inc, aop, pan, tau, plx, mtot))

        labels = ["sma1", "ecc1", "inc1", "aop1", "pan1", "tau1", "plx", "mtot"]
        f["col_names"] = labels
        # `loadhdf5` reads the `parameter_labels` *attribute*, and so does
        # orbitize!; writing only `col_names` means the export cannot be read
        # back without falling through to a guess, with a warning.
        attrs(f)["parameter_labels"] = labels
        attrs(f)["tau_ref_epoch"] = tau_ref_epoch
        attrs(f)["sampler_name"] = "Octofitter"

        dset = create_dataset(f, "post", Float32, size(dat))
        write(dset, convert(Matrix{Float32}, dat))
        return fname
    end
end

_system_of(model) = hasproperty(model, :system) ? model.system : model

"""First body placed by a single-body-exterior row: the "planet" of a simple model."""
function _first_companion(sys::System)
    for (_, ext, _) in sys.rows
        length(ext) == 1 && return ext[1]
    end
    error("System $(sys.name) has no body placed by a single-body orbit, so there is no " *
          "companion to export. Name one explicitly.")
end

"""
    _tp_from_chain(model, chain, body) -> Vector{Float64}

Epoch of periastron [MJD] for `body`, one entry per draw, rebuilt from the
chain when it carries no `<body>_tp` column.

A model that parametrizes phase as `θ` + `epoch` (the v9 spelling that
replaced `θ_at_epoch_to_tperi`) never produces a `tp` column, but `tp` is
determined by the elements — `PlanetOrbits.Row` computes it in its
constructor — so it is recoverable rather than missing.
"""
function _tp_from_chain(model, chain::Chains, body::Symbol)
    sys = _system_of(model)
    k = findfirst(r -> r[2] == (body,), sys.rows)
    isnothing(k) && error(
        "Cannot recover the epoch of periastron for :$body: no orbit places it on its own.")
    systems = construct_system(model, chain, :)
    return Float64[PlanetOrbits.periastron(s, k) for s in systems]
end

"""
    _interior_of(sys, body) -> Tuple{Vararg{Symbol}}

The bodies whose barycentre `body`'s own orbit is about — one name for an
astrocentric companion (`about=A`), several for a Jacobi chain
(`about=(A, b)`). Both export: orbitize! parametrizes its multi-planet fits in
Jacobi coordinates too, so the interior total is exactly what its `mtot` means.
"""
function _interior_of(sys::System, body::Symbol)
    for (_, ext, int) in sys.rows
        ext == (body,) && return int
    end
    error("System $(sys.name) has no orbit placing body :$body on its own. Its bodies are " *
          "$(join(sys.bodynames, ", ")).")
end
