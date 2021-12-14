using DelimitedFiles
using NamedTupleTools
using BasicInterpolators
using Interpolations

"""
    sonora_photometry_interpolator(:Keck_L′, [metalicity="+0.0"])

Given a supported photometric band and [M/H] metalicity (default=solar),
return a function of temperature (K) and mass (M_jup) that gives the 
absolute magnitude of the planet in that bandpass.

Supported bands:
:MKO_Y, :MKO_Z, :MKO_J, :MKO_H, :MKO_K, :MKO_L′, :MKO_M′, :TwoMASS_J, :TwoMASS_H, :TwoMASS_Ks, :Keck_Ks, :Keck_L′, :Keck_Ms, :SDSS_g′, :SDSS_r′, :SDSS_i′, :SDSS_z′, :IRAC_36, :IRAC_45, :IRAC_57, :IRAC_79, :WISE_W1, :WISE_W2, :WISE_W3, :WISE_W4

Supported metalicities:
"+0.0", "-0.5", "+0.5"
"""
function sonora_photometry_interpolator(band, metalicity="+0.0";catalog=datadep"SonoraBobcatEvoPhot")

    # ## Load Sonora cooling track 
    # valid_lines = filter(readlines(joinpath(datadep"SonoraBobcatEvoPhot", "evolution_tables", "evo_tables+0.0", "nc+0.0_co1.0_age"))) do line
    #     length(line) > 10
    # end
    # headers = split(valid_lines[1], r"  +")[2:7]
    # headers = lowercase.(replace.(headers, r"\W"=>""))
    # grid= mapreduce(vcat, valid_lines[2:end]) do line
    #     permutedims(parse.(Float64, split(line, r"  +")[2:7]))
    # end
    ##
    # cooling_grid = DataFrame(grid, headers)


    #  Load Sonora magnitude table
    mag_table = load_table(joinpath(catalog, "photometry_tables", "mag_table"*metalicity))

    # We first use BasicInterpolators to grid the sparsely sampled models
    # Then we use Interpolations to make a fast, autodiff compatible 
    # linear interpolator of the resulting data.
    points = [ 
        mag_table.Teff./10 mag_table.mass
    ]
    samples = collect(mag_table[band])

    sitp = RBFInterpolator(points, samples, 2)
    minmass, maxmass = extrema(mag_table.mass)
    minT, maxT = extrema(mag_table.Teff)
    temp_mass_to_abs_mag = (teffk, mass_mjup) -> sitp(teffk/10, mass_mjup)
    teff = range(minT, maxT,length=200)
    mass = range(minmass, maxmass,length=200)
    gridded = temp_mass_to_abs_mag.(teff, mass')

    # Now build linear interpolator
    litp = LinearInterpolation((teff./10, mass), gridded, extrapolation_bc=NaN)
    function model_interpolator(teffk, mass)
        if minT <= teffk <= maxT && minmass < mass < maxmass
            return litp(teffk/10, mass)
        else
            return NaN
        end
    end
        
    return model_interpolator
end
export sonora_photometry_interpolator


function sonora_cooling_interpolator(metalicity="+0.0";catalog=datadep"SonoraBobcatEvoPhot")

    # Load Sonora cooling track 
    # valid_lines = filter(readlines("./sonora/evolution_tables/evo_tables+0.0/nc+0.0_co1.0_age")) do line
    valid_lines = filter(readlines(joinpath(catalog, "evolution_tables", "evo_tables"*metalicity, "nc$(metalicity)_co1.0_age"))) do line
        length(line) > 10
    end
    headers = split(valid_lines[1], r"  +")[2:7]
    headers = lowercase.(replace.(headers, r"\W"=>""))
    grid= mapreduce(vcat, valid_lines[2:end]) do line
        permutedims(parse.(Float64, split(line, r"  +")[2:7]))
    end
    # cooling_grid = DataFrame(grid, headers)
    return grid, headers

    iteffk = findfirst(==("teffk"), headers)
    immsun = findfirst(==("teffk"), headers)
    teffk = grid[iteffk]
    mmsun = grid[immsun]

    points = [ 
        teffk/10 mmsun
    ]
    samples = collect(mag_table[band])

    sitp = RBFInterpolator(points, samples, 2)
    litp = LinearInterpolation(teffk, grid[], extrapolation_bc=NaN)


    # Want (age, mass) -> teffk. This can be combined with the others to get specific fluxes.

end
export sonora_cooling_interpolator


function load_table(fname)

    headers = open(fname, lock = false, read = true) do f
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        h1 = readline(f)
        h2 = readline(f)

        headers_1 = strip.(split(h1, '|'))
        headers_2 = strip.(split(h2, r"  +"))

        # Starting after 6
        headers_2_1_indices =
            [2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7]

        combined = vcat(
            headers_2[1:6],
            headers_1[headers_2_1_indices] .* '_' .* headers_2[7:end],
        )
        # Normalize headers
        combined = replace.(combined, '\'' => '′')
        combined = replace.(combined, ' ' => "")
        combined = replace.(combined, "/" => '_')
        combined = replace.(combined, "2MASS" => "TwoMASS")
        combined = replace.(combined, r"[^\w′_]" => "")
        return combined
    end

    data = readdlm(fname, String, skipstart = 10, header = false)
    data = [
        try
            typeof(d) <: AbstractString ? parse(Float64, replace(d, '*' => "")) : d
        catch err
            NaN
        end
        for d in data
    ]
    # Return simple table
    return namedtuple(headers, eachcol(data))
end