
function HARPS_DR1_rvs(target_main_id_simbad; catalog=datadep"ESOHARPS_DR1_rvs", inst_idx::Int=1)
#     rvcat = FITS(joinpath(catalog, "ADP.2023-12-04T15:16:53.464.fits"),"r") do fits
#         Table(fits[2])
#     end
# return rvcat
#     target_rows = findall(==(target_main_id_simbad), rvcat.main_id_simbad)
#     table = rvcat[target_rows]
    table = HARPS_DR1_rvs_observations(target_main_id_simbad; catalog)
    Table(;
        epoch = table.mjd_obs,
        inst_idx = fill(inst_idx,size(table,1)),
        rv = table.drs_ccf_rvc,
        σ_rv = table.drs_dvrms,
    )
    # target_rows = findall(==(target), rvbank.target)
    # table = rvbank[target_rows]

    # return StarAbsoluteRVLikelihood(Table(;
    #     epoch=OctofitterRadialVelocity.jd2mjd.(table.BJD),
    #     inst_idx=fill(inst_idx,size(table,1)),
    #     rv=-table.RV_mlc_nzp,
    #     σ_rv=table.e_RV_mlc_nzp,
    # ))
end

function HARPS_DR1_rvs_observations(target_main_id_simbad; catalog=datadep"ESOHARPS_DR1_rvs", inst_idx::Int=1)
    rvcat = FITS(joinpath(catalog, "ADP.2023-12-04T15:16:53.464.fits"),"r") do fits
        Table(fits[2])
    end
    target_matched_i = findall(==(target_main_id_simbad), rvcat.main_id_simbad)

    if isnothing(target_matched_i)
        avail_filt = unique(rvcat.main_id_simbad[rvcat.main_id_simbad.!=""])
        similarity = evaluate.(Ref(Levenshtein()), target_main_id_simbad, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3,end)]]
        @error "No results were found for the target $target."
        @info  "Here are a list of similar and available target names" closest_3
        error()
    end
    return rvcat[target_matched_i]
   
end