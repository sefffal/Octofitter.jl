using CSV
using StringDistances


function HARPS_RVBank_observations(target, catalog=datadep"HARPS_RVBank")

    
    rvbank = CSV.read(joinpath(catalog, "HARPS_RVBank_ver02.csv"), Table)

    target_matched_i = findfirst(==(target), rvbank.target)

    if isnothing(target_matched_i)
        avail_filt = unique(rvbank.target[rvbank.target.!=""])
        similarity = evaluate.(Ref(Levenshtein()), target, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3,end)]]
        @error "No results were found for the target $target."
        @info  "Here are a list of similar and available target names" closest_3
        error()
    end
    return rvbank[target_matched_i]
   
end

function HARPS_RVBank_rvs(target, catalog=datadep"HARPS_RVBank"; inst_idx::Int=1)

    rvbank = CSV.read(joinpath(catalog, "HARPS_RVBank_ver02.csv"), Table)

    target_rows = findall(==(target), rvbank.target)
    table = rvbank[target_rows]

    return StarAbsoluteRVLikelihood(Table(;
        epoch=OctofitterRadialVelocity.jd2mjd.(table.BJD),
        inst_idx=fill(inst_idx,size(table,1)),
        rv=-table.RV_mlc_nzp,
        σ_rv=table.e_RV_mlc_nzp,
    ))
end