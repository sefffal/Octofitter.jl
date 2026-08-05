# ---------------------------------------------------
# BHAC15 isochrone tables   (agent H)
#
# The second flux-table source behind `flux_<band>` body variables. Same shape
# as `sonora.jl`, same two v2 accommodations: masses in are M⊙ by default, and
# magnitudes out have to be converted before they can be a linear
# `flux_<band>`. See the header of `sonora.jl` for the reasoning.
# ---------------------------------------------------

# Logic to read each record array by age, and load using CSV
function _load_bhac15_models(fname)
    lines = readlines(fname)
    local colnames

    record_start_i = 0
    record_stop_i = 0
    age_Gyr = 0.0
    records = FlexTable[]
    for i in eachindex(lines)
        if contains(lines[i], "t (Gyr) = ",)
            age_Gyr = parse(Float64, split(lines[i], '=')[end])
            record_start_i = i + 4
        end
        if i == record_start_i - 2
            colnames = map(match -> Symbol(match.captures[1]), eachmatch(r"([\w\/]+)", lines[i]))
        end
        if i < record_start_i
            continue
        end
        if contains(lines[i], r"!---",)
            record_stop_i = i - 1
            continue
        end
        if record_stop_i < i && record_stop_i > 0
            str = join(lines[record_start_i:record_stop_i], "\n")
            io = IOBuffer(str,)
            df = CSV.read(
                io,
                FlexTable,
                header=colnames,
                normalizenames=true,
                delim=' ',
                ignorerepeated=true,
                skipto=4,
                comment="!"
            )
            df.age_Gyr = fill(age_Gyr, size(df, 1))
            push!(records, df)
            record_stop_i = 0
            record_start_i = 0
        end
    end
    return records
end

# The mass column is spelled `M/Ms` in the file. Whether that survives to a
# property name depends on CSV.jl's `normalizenames`, which differs across
# versions, so look it up by any of its plausible spellings rather than
# hard-coding one and failing at the last line of a two-minute load.
function _bhac_masscol(tbl)
    for nm in (Symbol("M/Ms"), :M_Ms, :MMs, :M_Ms_)
        hasproperty(tbl, nm) && return getproperty(tbl, nm)
    end
    error("Could not find the BHAC15 mass column in the loaded table; its columns are " *
          "$(join(propertynames(tbl), ", ")).")
end

"""
    itp = bhac15_mass_age_interpolator(; key=:G)
    itp = bhac15_mass_age_interpolator("BHAC15_iso.GAIA"; key=:G)

Create a function mapping `(age_Myr, mass)` -> absolute magnitude in the
column named by `key`, using the BHAC15 model grids. With no filename the
`BHAC15_GAIA` DataDep is used (and downloaded on first use).

    itp = Octofitter.bhac15_mass_age_interpolator(key=:G)
    itp(15.0, 0.08)               # mass in M⊙, like every mass in v2

`key` names a column of the isochrone file (`:Teff`, `:G`, `:G_BP`, … for the
GAIA tables). Out-of-grid inputs give `NaN` rather than an extrapolation.
`mass_unit` selects the input unit: `:Msol` (default), `:Mjup` (the v1
convention), or `:Mearth`.

Like the Sonora interpolators, this returns an **absolute magnitude**. A body's
`flux_<band>` variable is a linear flux, so convert before assigning it:

    flux_G = 10^(-0.4 * \$itp(system.age, mass))
"""
function bhac15_mass_age_interpolator(fname=datadep"BHAC15_GAIA"; key,
                                      mass_unit::Symbol=:Msol)
    mscale = _mass_to_mjup(mass_unit)
    records = _load_bhac15_models(fname)

    # Need list grid of (age_myr X mass_jup) -> (magnitude)
    dfall = reduce(vcat, records)

    agemyr = dfall.age_Gyr .* 1000
    mmjup = _bhac_masscol(dfall) ./ PlanetOrbits.mjup

    points = [
        log.(agemyr) log.(mmjup)
    ]
    samples = getproperty(dfall, key)

    sitplog = RBFInterpolator(points, samples, 0.1)

    sitp = (agemyr, mmjup) -> sitplog(log(agemyr), log(mmjup))

    # Now build a fast linear interpolator over this grid.
    minmmjup, maxmmjup = extrema(mmjup)
    minagemyr, maxagemyr = extrema(agemyr)
    agemyrrange = range(minagemyr, maxagemyr, length=2000)
    mmjuprange = range(minmmjup, maxmmjup, length=500)
    gridded = sitp.(agemyrrange, mmjuprange')

    # Now build linear interpolator
    litp = LinearInterpolation((agemyrrange, mmjuprange), gridded, extrapolation_bc=NaN)
    function model_interpolator(agemyr, mass)
        m = mass * mscale
        if minagemyr <= agemyr <= maxagemyr && minmmjup < m < maxmmjup
            return litp(agemyr, m)
        else
            return convert(promote_type(typeof(agemyr), typeof(m)), NaN)
        end
    end

    return model_interpolator
end
