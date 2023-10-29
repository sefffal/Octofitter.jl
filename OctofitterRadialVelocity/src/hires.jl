



function HIRES_search(target, catalog=datadep"HIRES_rvs")

    fnames = readdir(catalog,join=true)
    fnames_search = readdir(catalog,join=false)
    target_names = map(fnames_search) do fname
        obj = first(split(first(split(fname, '_')), '.'))
    end

    target_matched_i = findfirst(==(target), target_names)

    if isnothing(target_matched_i)
        avail_filt = unique(target_names[target_names.!=""])
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



function HIRES_rvs(target, catalog=datadep"HIRES_rvs"; inst_idx::Int=1)

    fname = HIRES_search(target)

    # Julian Date	Velocity (m/s)	Uncertainty (m/s)	S_value	Halpha	median photons per pixel	exposure time (seconds)
    
    # We have to load in all of the papers and then filter down.
    df_all = Table(
        Tables.rowtable(
        (open(fname,"r") do f
            for _ in 1:17
                readline(f)
            end
            rows = NamedTuple[]
            for line in eachline(f)
                jd = tryparse(Float64, line[1:13])
                if isnothing(jd)
                    continue
                end
                mjd = jd2mjd(jd)
                rv = parse(Float64, line[14:24])#,
                σ_rv=parse(Float64, line[25:31])#,
                push!(rows, (;jd,epoch=mjd,rv,σ_rv,inst_idx))
            end
            return rows
    end)))
    return RadialVelocityLikelihood(df_all)
end