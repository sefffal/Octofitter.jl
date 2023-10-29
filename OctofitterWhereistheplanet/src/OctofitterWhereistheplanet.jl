module OctofitterWhereistheplanet

# using TypedTables
using Octofitter
using HDF5
using DataDeps
using StringDistances
using MCMCChains: MCMCChains

function search(target, catalog=datadep"Whereistheplanet")

    dirpath = joinpath(catalog, "whereistheplanet-master", "data")
    fnames = readdir(dirpath,join=true)
    avail_targets = map(fnames) do fname
        m = match(r"post_(.+)\.hdf5", fname)
        return !isnothing(m) ? m.captures[1] : ""
    end
    fname_matched_i = findfirst(==(target), avail_targets)

    if isnothing(fname_matched_i)
        avail_filt = avail_targets[avail_targets.!=""] 
        similarity = evaluate.(Ref(Levenshtein()), target, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3,end)]]
        @error "No results were found for the target $target."
        @info  "Here are a list of similar and available target names" closest_3
        error()
    end

    return fnames[fname_matched_i]

end

function astrom(target, catalog=datadep"Whereistheplanet"; object=1)

    fname = search(target)
    return h5open(fname, "r") do f
        records = read(f["data"])

        records = filter(row->row.object==object, records)

        # Group observations by type
        seppa = filter(row->row.quant_type=="seppa", records)
        radec = filter(row->row.quant_type=="radec", records)

        out = AstrometryLikelihood[]
        if length(seppa) > 0
            astrom_seppa = AstrometryLikelihood(map(seppa) do row
                cor=row.quant12_corr
                if !isfinite(cor)
                    cor = 0.0
                end
                (;row.epoch, sep=row.quant1, σ_sep=row.quant1_err, pa=row.quant2, σ_pa=row.quant2_err, cor)
            end)
            push!(out, astrom_seppa)
        end
        if length(radec) > 0
            astrom_radec = AstrometryLikelihood(map(radec) do row
                cor=row.quant12_corr
                if !isfinite(cor)
                    cor = 0.0
                end
                (;row.epoch, ra=row.quant1, σ_ra=row.quant1_err, dec=row.quant2, σ_dec=row.quant2_err, cor)
            end)
            push!(out, astrom_radec)
        end

        # return(out...)
    end
end


"""
Load an Orbitize! posterior from an HDF5 file and convert it into
an Octofitter-compatible chains format.
Both tools use the same orbit conventions so this is fairly straightforward.

If you pass `numchains` as a second argument, the array will be interpretted
as coming from multiple chains concatenated together, 
"""
function posterior(fname_or_targetname, numchains=1)
    if !occursin(".hdf5", fname_or_targetname)
        fname = search(fname_or_targetname)
    else
        fname = fname_or_targetname
    end
    return h5open(fname, "r") do f
        # Standard orbitize basis assumed:
        # semi-major axis (sma), eccentricity (ecc), inclination (inc), argument of periastron (aop), position angle of the nodes (pan), epoch of periastron expressed as a fraction of the period past a reference epoch (tau), parallax (plx) and total system mass (mtot)
        arr = transpose(read(f["post"]))
        if numchains > 1
            lenchain = size(arr,1)÷numchains
            arr = stack(
                (
                    arr[(chain_no-1)*lenchain+1:chain_no*lenchain,:]
                    for chain_no in 1:numchains
                ),
            )
        end
        chn =  MCMCChains.Chains(
            arr,
            [
                :b_a,
                :b_e,
                :b_i,
                :b_ω,
                :b_Ω,
                :b_τ,
                :plx,
                :M
            ],
        )
        # Read additional attributes in and convert to named tuple
        metadata = [Symbol(k)=>v for (k,v) in attrs(f)]
        return MCMCChains.setinfo(chn, NamedTuple(metadata))
    end
end

function __init__()

    register(DataDep("Whereistheplanet",
        """
        Dataset:     Planet astrometry and orbit fits from whereistheplanet.com
        Author:      Wang et al.
        License:     BSD-3 Clause
        Website:     https://github.com/semaphoreP/whereistheplanet

        File size: 10MiB
        """,
        "https://github.com/semaphoreP/whereistheplanet/archive/refs/heads/master.zip",
        # "c02e7c601dc94d7acd0c58398b518038b036d1507f790f3419b574b39d515197",
        post_fetch_method=unpack
    ))




    return
end
end
