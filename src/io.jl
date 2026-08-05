# ---------------------------------------------------
# Chain and model serialization   (agent H)
#
# `savechain` / `loadchain`.
#
# Chains are stored as a FITS binary table so that the samples and the run's
# metadata travel in one file. Unicode parameter names (`b_Ω`, `c_ω`) are not
# legal FITS column names, so they are latexified on the way out and restored
# on the way back in — that part is unchanged from v1 and deliberately so:
# people have saved chains with it.
#
# What *is* new is that the file says which model surface produced it. v1
# stored no structural metadata at all — only `chain.info` verbatim — so a
# chain's shape lived entirely in its column names, and those changed with the
# flattening: v1 wrote a planet-owned observation's variable as
# `<planet>_<obs>_<var>` (`b_GPI_jitter`), v2 writes `<obs>_<var>`
# (`GPI_jitter`), because observations are no longer owned by a body. Loading
# is not where that bites — a chain is just named columns, and they round-trip
# either way. It bites when a v1 chain meets a v2 *model*:
# `mcmcchain2result` looks each parameter up by name and substitutes `missing`
# for anything it cannot find, so the mismatch is silent. Hence two additions:
#
#   - a format stamp (`OCTOFMT`), so `loadchain` can say plainly that a file
#     predates the flattening rather than leaving the user to discover it;
#   - `checkchain(model, chain)`, which is where the actual model↔chain
#     agreement is asserted, and which `loadchain(fname; model)` calls for you.
#
# The section map (`:parameters` vs `:internals`) is now saved too. v1 dropped
# it, so every reloaded chain came back with `tree_depth` and `numerical_error`
# masquerading as model parameters — visible in `describe`, in corner plots,
# and in anything that iterates `keys(chain)`.
# ---------------------------------------------------

using Latexify, MathTeXEngine

using FITSIO

"""
Version of the Octofitter FITS chain layout written by [`savechain`](@ref).

  - **1** — v1 (implicit; those files carry no stamp at all). Two HDUs: an
    empty primary carrying `chain-info` cards, and the sample table. Column
    names follow v1's three-level naming for planet-owned observations.
  - **2** — v2. Adds `OCTOFMT`/`OCTOVER` cards to the primary header and a
    third HDU holding the chain's section map.

A v2 reader loads a v1 file (with a warning); a v1 reader loads a v2 file,
since the extra HDU and cards are additive.
"""
const CHAIN_FORMAT_VERSION = 2

# --- header value marshalling -------------------------------------------------

# FITS header values are numbers, booleans or short ASCII strings. `chain.info`
# is none of those in general: a NUTS run stores `samples_transformed` (one
# entry per draw), an adaptor, and a metric.
#
# v1 handled that with `filter(isascii, string(x)[1:min(255,end)])`, which has
# two problems. It builds the *entire* string first — for a 10k-draw chain
# that is tens of megabytes of formatting to keep 255 bytes of it — and
# `String` slicing is by byte, so a multi-byte character straddling byte 255
# throws `StringIndexError` from inside `savechain`. Summarising containers by
# type and size avoids both. Nothing reads these cards back as values; they
# exist so a human can see what produced the file.
_hdrval(x::Number) = x
_hdrval(x::Bool) = x
_hdrval(::Nothing) = nothing
_hdrval(x::AbstractString) = _ascii_trunc(String(x))
_hdrval(x::Symbol) = _ascii_trunc(String(x))
_hdrval(@nospecialize x) = _ascii_trunc(_hdrsummary(x))

_hdrsummary(x::AbstractArray) = string(typeof(x), size(x))
function _hdrsummary(@nospecialize x)
    io = IOBuffer()
    try
        show(IOContext(io, :limit => true, :compact => true), x)
    catch
        print(io, typeof(x))
    end
    return String(take!(io))
end

function _ascii_trunc(s::AbstractString, n::Int=255)
    t = filter(isascii, s)
    return length(t) <= n ? t : first(t, n)
end

# Numeric vectors are expanded across `<key>_1`, `<key>_2`, … cards, which is
# how v1 round-tripped them and how `loadchain` reassembles them. The cap is
# new: an unlucky `info` entry would otherwise write one card per element and
# produce a header larger than the data.
const _HDR_ARRAY_MAX = 1024

"""
    savechain("saved-chain.fits", chain)

Save an `MCMCChains.Chains` to a FITS binary table, together with the run
metadata in `chain.info` and the chain's section map.

Unicode parameter names are stored as their LaTeX spellings and restored by
[`loadchain`](@ref).
"""
function savechain(fname, chain::MCMCChains.Chains)
    tbl = Table(chain)
    FITS(fname, "w") do fits
        # Format stamp first, so it is near the top of the header when a human
        # runs `fitsheader` on the file.
        ks = String["OCTOFMT", "OCTOVER"]
        vals = Any[CHAIN_FORMAT_VERSION, OCTO_VERSION_STR]
        comments = String["octo-chainfmt", "octo-chainfmt"]

        info = chain.info
        for k in keys(info)
            v = info[k]
            if v isa AbstractArray{<:Number}
                n = length(v)
                if n > _HDR_ARRAY_MAX
                    @warn "chain.info[:$k] has $n elements; storing a summary instead of expanding it into header cards." maxlog = 1
                    push!(ks, string(k)); push!(vals, _hdrval(v)); push!(comments, "chain-info")
                    continue
                end
                # The marker tells `loadchain` to go looking for the elements.
                push!(ks, string(k)); push!(vals, "ARRAY"); push!(comments, "chain-info")
                for (i, el) in enumerate(v)
                    push!(ks, "$(k)_$i"); push!(vals, el); push!(comments, "chain-info")
                end
            else
                push!(ks, string(k)); push!(vals, _hdrval(v)); push!(comments, "chain-info")
            end
        end

        write(fits, zeros(1, 1), header=FITSHeader(ks, vals, comments))

        # Take the titles from the table itself rather than assembling them
        # from `keys(chain)` separately: the two orders have to agree, and
        # deriving both from one source is the only way to be sure they do.
        coltitles = collect(propertynames(tbl))
        write(
            fits,
            _coltitle_serialize.(string.(coltitles)),
            [collect(getproperty(tbl, k)) for k in coltitles];
            hdutype=FITSIO.TableHDU,
        )

        # Section map. `Chains` keeps it as a NamedTuple of name vectors;
        # anything not listed here loads back as `:parameters`, which is the
        # default, so only the other sections need storing.
        secnames = String[]
        secmembers = String[]
        for section in keys(chain.name_map)
            section === :parameters && continue
            for nm in chain.name_map[section]
                push!(secnames, _coltitle_serialize(string(nm)))
                push!(secmembers, String(section))
            end
        end
        if !isempty(secnames)
            write(fits, ["name", "section"], Any[secnames, secmembers];
                hdutype=FITSIO.TableHDU, name="OCTOSECTIONS")
        end
    end
    return fname
end

_coltitle_serialize(title::AbstractString) =
    all(isascii, title) ? String(title) : String(latexify(title))

function coltitle_restorer(title::AbstractString)
    replacements = Pair{String,Char}[]
    for m in eachmatch(r"(\\[^\W_]+)", title)
        if !haskey(MathTeXEngine.command_definitions, m.captures[1])
            @error "Could not restore column" title m.captures[1]
            error()
        end
        push!(
            replacements,
            m.captures[1] =>
                Char(MathTeXEngine.command_definitions[m.captures[1]][1].args[1])
        )
    end
    replace(title, replacements..., '$' => "", '{' => "", '}' => "", "\\" => "")
end

"""
    loadchain("saved-chain.fits")
    loadchain("saved-chain.fits"; model)

Load an `MCMCChains.Chains` from a FITS binary table written by
[`savechain`](@ref).

Pass `model` to have the loaded chain checked against that model's parameter
list ([`checkchain`](@ref)) before it is returned. That is the check worth
doing: a chain and a model that disagree about a parameter's name do not
error, they silently produce `missing`.

Files written by Octofitter v1 load, with a warning — their column names use
v1's `<planet>_<observation>_<variable>` spelling for observations that hung
off a companion, which no v2 model will match.
"""
function loadchain(fname; model=nothing)
    chn = FITS(fname) do fits
        tbl = Table(fits[2])
        coltitles_serialized = propertynames(tbl)
        coltitles_deserialized = map(coltitle_restorer, string.(coltitles_serialized))
        out = map(coltitles_serialized, coltitles_deserialized) do colin, colout
            Symbol(colout) => getproperty(tbl, colin)
        end
        input = Table(; out...)
        # Now we have the original table back with proper column names.

        # Convert back into MCMCChains. The saved table is one row per
        # (iteration, chain) pair; rebuild the (iters × params × chains) cube
        # by position in the *sorted unique* values rather than by using the
        # iteration number as an index directly. v1 did the latter, which
        # silently required chains numbered 1:n and iterations numbered 1:m —
        # a chain that had been thinned or sliced (`chain[501:end]`) wrote
        # fine and then threw `BoundsError` on load.
        iters = sort(unique(input.iteration))
        chainnos = sort(unique(input.chain))
        iterpos = Dict(v => i for (i, v) in enumerate(iters))
        chainpos = Dict(v => i for (i, v) in enumerate(chainnos))
        cols = setdiff(propertynames(input), (:iteration, :chain))
        T = typeof(sum(collect(first(input))))
        data = zeros(T, length(iters), length(cols), length(chainnos))
        seen = falses(length(iters), length(chainnos))
        for (j, key) in enumerate(cols)
            col = getproperty(input, key)
            for i in eachindex(input)
                data[iterpos[input.iteration[i]], j, chainpos[input.chain[i]]] = col[i]
            end
        end
        for i in eachindex(input)
            seen[iterpos[input.iteration[i]], chainpos[input.chain[i]]] = true
        end
        all(seen) || error(
            "$fname does not hold a rectangular chain: $(count(!, seen)) of " *
            "$(length(seen)) (iteration, chain) slots are missing. It may have been " *
            "written by something other than `savechain`, or truncated.")

        sectionmap = _read_sectionmap(fits, cols)
        chn = MCMCChains.Chains(data, cols, sectionmap)
        if iters != collect(1:length(iters))
            chn = MCMCChains.setrange(chn, iters)
        end

        # And restore metadata.
        h = read_header(fits[1])
        meta = Dict{Symbol,Any}()
        for k in keys(h)
            if get_comment(h, k) == "chain-info"
                if h[k] != "ARRAY"
                    meta[Symbol(lowercase(k))] = h[k]
                else
                    arrdat = []
                    i = 0
                    while true
                        i += 1
                        ki = "$(k)_$i"
                        if haskey(h, ki)
                            push!(arrdat, h[ki])
                        else
                            break
                        end
                    end
                    meta[Symbol(lowercase(k))] = identity.(arrdat)
                end
            end
        end
        fmt = haskey(h, "OCTOFMT") ? Int(h["OCTOFMT"]) : 1
        meta[:chainformat] = fmt
        if fmt < CHAIN_FORMAT_VERSION
            @warn """
            $fname was written by Octofitter v1 (no chain-format stamp).

            Its samples load unchanged, but its column names follow the v1 model
            surface: an observation attached to a companion was named
            `<planet>_<observation>_<variable>` (`b_GPI_jitter`), where v2 names it
            `<observation>_<variable>` (`GPI_jitter`) because observations are no
            longer owned by a body. Rebuilding parameters from this chain against a
            v2 model will silently yield `missing` for every such column — pass
            `model=` to `loadchain`, or call `Octofitter.checkchain(model, chain)`,
            to have that checked.
            """
        end
        return MCMCChains.setinfo(chn, NamedTuple(meta))
    end
    if model !== nothing
        checkchain(model, chn)
    end
    return chn
end

# v1 files have two HDUs; the section table is HDU 3 when present.
function _read_sectionmap(fits, cols)
    sectionmap = Dict{Symbol,Vector{Symbol}}()
    length(fits) >= 3 || return sectionmap
    hdu = fits[3]
    sectbl = Table(hdu)
    (hasproperty(sectbl, :name) && hasproperty(sectbl, :section)) || return sectionmap
    known = Set(cols)
    for (nm, sec) in zip(sectbl.name, sectbl.section)
        key = Symbol(coltitle_restorer(strip(String(nm))))
        key in known || continue
        push!(get!(sectionmap, Symbol(strip(String(sec))), Symbol[]), key)
    end
    return sectionmap
end

# --- model ↔ chain agreement --------------------------------------------------

"""
    checkchain(model, chain; strict=true)

Assert that `chain` carries every free parameter `model` expects, and explain
the mismatch if it does not.

Worth doing explicitly because the failure it catches is quiet:
`mcmcchain2result` looks each parameter up by name and yields `missing` for
anything absent, so a chain from a *different* model — or from Octofitter v1,
whose planet-owned observations were named `<planet>_<obs>_<var>` — flows
through plotting and post-prediction producing nonsense rather than an error.

With `strict=false` the mismatch is a warning instead. Returns `chain`.
"""
function checkchain(model, chain::MCMCChains.Chains; strict::Bool=true)
    sys = hasproperty(model, :system) ? model.system : model
    present = Set(names(chain))
    missingnames = Symbol[]
    for nm in list_parameter_names(sys)
        # A vector-valued prior is flattened to `<name>_1`, `<name>_2`, … on
        # the way into a chain, so its base name is legitimately absent.
        (nm in present || Symbol(nm, "_1") in present) || push!(missingnames, nm)
    end
    isempty(missingnames) && return chain

    # If the chain has a column whose name *ends* with the one we wanted, it is
    # almost certainly a v1 file: `b_GPI_jitter` vs `GPI_jitter`.
    hints = String[]
    for nm in missingnames
        suffix = "_" * string(nm)
        for have in names(chain)
            if endswith(string(have), suffix)
                push!(hints, "$have  →  $nm")
                break
            end
        end
    end
    msg = "The chain does not match this model. Missing parameter(s): " *
          join(missingnames, ", ") * "."
    if !isempty(hints)
        msg *= "\n\nThe chain does have similarly-named columns:\n  " *
               join(hints, "\n  ") *
               "\n\nThat is v1's three-level `<planet>_<observation>_<variable>` naming. " *
               "v2 observations are not owned by a body, so they are named " *
               "`<observation>_<variable>`. Rename the columns, or re-run the fit."
    end
    strict ? error(msg) : @warn msg
    return chain
end

# Deliberately not exported, matching v1 and every documented call site: these
# are spelled `Octofitter.savechain(…)` / `Octofitter.loadchain(…)`, and
# `save`/`load` are names users' own scripts often own.
