
# Logic to read each record array by age, and load using CSV
function _load_bhac15_models(fname="BHAC15-models.dat")
    lines = readlines(fname)
    colnames = [
        :M_Ms
        :Teff
        :L_Ls
        :g
        :R_Rs
        :Li_Li0
        :Mj
        :Mh
        :Mk  
    ]

    record_start_i = 0
    record_stop_i = 0
    age_Gyr = 0.0
    records = FlexTable[]
    for i in eachindex(lines)
        if contains(lines[i],"t (Gyr) = ",)
            age_Gyr = parse(Float64, split(lines[i],'=')[end])
            record_start_i = i + 4
        end
        if i < record_start_i
            continue
        end
        if contains(lines[i],r"!---",)
            record_stop_i = i-1
            continue
        end
        if record_stop_i < i && record_stop_i > 0
            str = join(lines[record_start_i:record_stop_i], "\n")
            io = IOBuffer(str,)
            df0 = CSV.read(io, DataFrame, header=0, normalizenames=true, delim=' ', ignorerepeated=true, skipto=4, comment="!")
            df = FlexTable(df0, colnames)
            df.age_Gyr = age_Gyr
            push!(records, df)
            record_stop_i = 0
            record_start_i = 0
        end
    end
    return records
end

"""
    itp = bhac15_mass_age_interpolator()


Create a function mapping (age_Myr, mass_Mjup) -> absolute K band magnitude using the
BHAC15 model grids.
    
"""
function bhac15_mass_age_interpolator(fname="BHAC15-models.dat")
    records = _load_bhac15_models(fname)

    # Need list grid of (age_myr X mass_jup) -> (K Mag) 
    dfall = reduce(vcat, records)

    agemyr = dfall.age_Gyr .* 1000
    mmjup = dfall.M_Ms ./ mjup2msol


    points = [ 
        log.(agemyr) log.(mmjup)
    ]
    samples = dfall.Mk

    sitplog = RBFInterpolator(points, samples, 0.1)

    sitp = (agemyr, mmjup) -> sitplog(log(agemyr), log(mmjup))
    # return (;sitp, agemyr, mmjup, teffk)
    
    # Now build a fast linear interpolator over this grid.
    minmmjup, maxmmjup = extrema(mmjup)
    minagemyr, maxagemyr = extrema(agemyr)
    agemyrrange = range(minagemyr, maxagemyr,length=2000)
    mmjuprange = range(minmmjup, maxmmjup,length=500)
    gridded = sitp.(agemyrrange, mmjuprange')

    # Now build linear interpolator
    litp = LinearInterpolation((agemyrrange, mmjuprange), gridded, extrapolation_bc=NaN)
    function model_interpolator(agemyr, mmjup)
        if minagemyr <= agemyr <= maxagemyr && minmmjup < mmjup < maxmmjup
            return litp(agemyr, mmjup)
        else
            return NaN
        end
    end

    return model_interpolator
end
