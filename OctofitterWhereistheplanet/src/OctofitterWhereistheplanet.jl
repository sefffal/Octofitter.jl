module OctofitterWhereistheplanet

# using TypedTables
using Octofitter
using HDF5
using DataDeps

function Whereistheplanet_search(target, catalog=datadep"Whereistheplanet")

    dirpath = joinpath(catalog, "whereistheplanet-master", "data")
    fnames = readdir(dirpath,join=true)
    fname_matched_i = findfirst(fnames) do fname
        m = match(r"post_(.+)\.hdf5", fname)
        !isnothing(m) && m.captures[1] == target
    end

    return fnames[fname_matched_i]

end

function Whereistheplanet_astrom(target, catalog=datadep"Whereistheplanet")

    fname = Whereistheplanet_search(target)
    return h5open(fname, "r") do f
        records = read(f["data"])

        # Group observations by type
        seppa = filter(row->row.quant_type=="seppa", records)
        radec = filter(row->row.quant_type=="radec", records)

        out = Astrometry[]
        if length(seppa) > 0
            astrom_seppa = Astrometry(map(seppa) do row
                cor=row.quant12_corr
                if !isfinite(cor)
                    cor = 0.0
                end
                (;row.epoch, sep=row.quant1, σ_sep=row.quant1_err, pa=row.quant2, σ_pa=row.quant2_err, cor)
            end)
            push!(out, astrom_seppa)
        end
        if length(radec) > 0
            astrom_radec = Astrometry(map(radec) do row
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




    
    # rvbank = CSV.read(joinpath(catalog, "HARPS_RVBank_v1.csv"), Table)

    # target_rows = findall(==(target), rvbank.target)
    # return rvbank[target_rows]
   
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
        "c02e7c601dc94d7acd0c58398b518038b036d1507f790f3419b574b39d515197",
        post_fetch_method=unpack
    ))




    return
end
end
