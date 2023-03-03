



function HIRES_search(target, catalog=datadep"HIRES_rvs")

    fnames = readdir(catalog,join=true)
    fnames_search = readdir(catalog,join=false)
    target_names = map(fnames_search) do fname
        obj = first(split(first(split(fname, '_')), '.'))
    end

    @show fnames_search
    target_matched_i = findfirst(==(target), target_names)

    if isnothing(target_matched_i)
        avail_filt = target_names[target_names.!=""] 
        similarity = evaluate.(Ref(Levenshtein()), target, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3,end)]]
        @error "No results were found for the target $target."
        @info  "Here are a list of similar and available target names" closest_3
        error()
    end
    return fnames[target_matched_i]


    
    # rvbank = CSV.read(joinpath(catalog, "HARPS_RVBank_v1.csv"), Table)

    # target_rows = findall(==(target), rvbank.target)
    # return rvbank[target_rows]
   
end
