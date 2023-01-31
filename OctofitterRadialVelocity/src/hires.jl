



function HIRES_search(target, catalog=datadep"HIRES_rvs")

    fname = findfirst(readdir(catalog,join=true)) do fname
        obj = first(split(fname, '_'))
        target == obj
    end



    
    # rvbank = CSV.read(joinpath(catalog, "HARPS_RVBank_v1.csv"), Table)

    # target_rows = findall(==(target), rvbank.target)
    # return rvbank[target_rows]
   
end