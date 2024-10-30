using CSV
using StringDistances


# TODO: make a search function. there is a file with a list of all stars.
function CES_lc_rvs(target, catalog=datadep"CES_rvs"; inst_idx::Int=1)
    fname = joinpath(catalog, "lc", target*".dat")
    tbl = open(fname, read=true) do f
        map(readlines(f)) do line
            epoch_bjd = parse(Float64, line[1:16])
            epoch = jd2mjd(epoch_bjd)
            rv = parse(Float64, line[18:27])
            σ_rv = parse(Float64, line[29:34])
            return (;epoch, inst_idx=1, rv, σ_rv)
        end
    end

    return tbl
end


function CES_vlc_rvs(target, catalog=datadep"CES_rvs"; inst_idx::Int=1)
    fname = joinpath(catalog, "vlc", target*".dat")
    tbl = open(fname, read=true) do f
        map(readlines(f)) do line
            epoch_bjd = parse(Float64, line[1:16])
            epoch = jd2mjd(epoch_bjd)
            rv = parse(Float64, line[18:27])
            σ_rv = parse(Float64, line[29:34])
            return (;epoch, inst_idx=1, rv, σ_rv)
        end
    end

    return tbl
end

# --------------------------------------------------------------------------------
# Bytes Format Units   Label     Explanations
# --------------------------------------------------------------------------------
# 1- 16  F16.8 d       BJD       Barycentric Julian date
# 18- 27  F10.2 m/s     RV        Barycentric radial velocity (1)
# 29- 34  F6.2  m/s   e_RV        Uncertainty on RV