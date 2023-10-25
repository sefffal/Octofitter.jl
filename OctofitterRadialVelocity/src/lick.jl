


function Lick_rvs(target, catalog=datadep"Lick_rvs")

    # We have to load in all of the papers and then filter down.
    df_lick_all = Table(
        Tables.rowtable(
        (open(joinpath(catalog,"apjs488421t2_mrt.txt"),"r") do f
            for _ in 1:17
                readline(f)
            end
            rows = NamedTuple[]
            for line in eachline(f)
                # --------------------------------------------------------------------------------
                # Bytes Format Units   Label  Explanations
                # --------------------------------------------------------------------------------
                #     1-  7 A7     ---     Name   Star identfier
                #     9- 19 F11.5  d       JD     ? Julian Date of the observation; - 2440000 (1)
                #    21- 29 F9.2   m/s     RVel   ? Radial velocity (1)
                #    31- 36 F6.2   m/s   e_RVel   ? Error in RVel (1)
                #    38- 41 I4     ---     SNR    Signal-to-Noise Ratio
                #    43- 44 I2     ---     Dewar  CCD dewar used
                star=strip(line[1:7])#,
                jd=tryparse(Float64, line[9:19])#,
                if isnothing(jd)
                    continue
                end
                jd += 2440000
                mjd=jd2mjd(jd)#,
                rv=parse(Float64, line[21:29])#,
                σ_rv=parse(Float64, line[31:36])#,
                dewar=parse(Int, line[43:44])
                push!(rows, (;target=star,jd,epoch=mjd,rv,σ_rv,dewar))
            end
            return rows
    end)))
    target_names = df_lick_all.target
    target_matched_ii = findall(==(target), target_names)

    if isempty(target_matched_ii)
        avail_filt = unique(target_names[target_names.!=""])
        similarity = evaluate.(Ref(Levenshtein()), target, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3,end)]]
        @error "No results were found for the target $target."
        @info  "Here are a list of similar and available target names" closest_3
        error()
    end
    return RadialVelocityLikelihood(df_lick_all[target_matched_ii,:])


    
    # rvbank = CSV.read(joinpath(catalog, "HARPS_RVBank_v1.csv"), Table)

    # target_rows = findall(==(target), rvbank.target)
    # return rvbank[target_rows]
   
end