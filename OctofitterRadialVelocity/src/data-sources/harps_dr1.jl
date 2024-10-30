
function HARPS_DR1_rvs(target_id; catalog=datadep"ESOHARPS_DR1_rvs", inst_idx::Int=1)
    table = HARPS_DR1_rvs_observations(target_id; catalog)
    Table(;
        epoch = jd2mjd.(table.drs_bjd),
        inst_idx = fill(inst_idx,size(table,1)),
        rv = table.drs_ccf_rvc,
        σ_rv = table.drs_dvrms,
    )
end

function HARPS_DR1_rvs_observations(target_id; catalog=datadep"ESOHARPS_DR1_rvs", inst_idx::Int=1)
    rvcat = FITS(joinpath(catalog, "ADP.2023-12-04T15:16:53.464.fits"),"r") do fits
        Table(fits[2])
    end
    target_matched_i = findall(==(target_id), rvcat.main_id_simbad)

    if isempty(target_matched_i)
        avail_filt = unique(rvcat.main_id_simbad[rvcat.main_id_simbad.!=""])
        similarity = evaluate.(Ref(Levenshtein()), target_id, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3,end)]]
        @error "No results were found for the target $target_id."
        @info  "Here are a list of similar and available target names" closest_3
        error()
    end
    return rvcat[target_matched_i]
   
end