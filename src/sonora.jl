# ---------------------------------------------------
# Sonora Bobcat evolutionary/photometry tables
#
# One of the two flux-table sources that feed `PhotometryObs` and the
# `flux_<band>` body variables. Under the flux/band unification a body's
# `flux_<band>` is what an interpolator like this produces, so these are the
# mass ↔ photometry bridge:
#
#     b = Body(name="b", about=A, variables=@variables begin
#         mass   ~ LogUniform(0.5mjup, 50mjup)          # M⊙, like every mass
#         tempK   = $cooling($(system).age, mass)
#         flux_H  = 10^(-0.4 * $absmag_H(tempK, mass))  # linear, see below
#         …
#     end)
#
# Two things that composition needs, and that v1 did not provide:
#
#  1. **Masses are M⊙.** v1's interpolators took Jupiter masses, because a v1
#     `Planet`'s `mass` variable was in Jupiter masses. v2 has one mass unit,
#     so the default input unit is M⊙ and `mass_unit=:Mjup` recovers the v1
#     calling convention exactly.
#
#  2. **These return absolute magnitudes, and `flux_<band>` must be linear.**
#     Photocentres are flux-weighted means (`PlanetOrbits.photocentre` builds
#     `w_j ∝ f_j`), so feeding a magnitude straight into `flux_H` would weight
#     bodies by a logarithm. The conversion is the one line above; the
#     docstrings say so at every entry point.
# ---------------------------------------------------

using DelimitedFiles
using NamedTupleTools
using BasicInterpolators
using Interpolations

"""Scale factor taking a mass in `unit` to the Jupiter masses these grids are tabulated in."""
function _mass_to_mjup(unit::Symbol)
    unit === :Msol && return 1 / PlanetOrbits.mjup
    unit === :Mjup && return 1.0
    unit === :Mearth && return PlanetOrbits.mearth / PlanetOrbits.mjup
    error("`mass_unit` must be :Msol (the v2 default), :Mjup (the v1 convention), or :Mearth; got :$unit")
end

"""
    sonora_photometry_interpolator(:Keck_L′, [metalicity="+0.0"])

Given a supported photometric band and [M/H] metalicity (default=solar),
return a function of temperature (K) and mass that gives the **absolute
magnitude** of the object in that bandpass.

    absmag_L = sonora_photometry_interpolator(:Keck_L′)
    absmag_L(1200.0, 12mjup)      # mass in M⊙, like every mass in v9

Out-of-grid inputs give `NaN` rather than an extrapolation.

# Feeding a `flux_<band>` variable
The result is a magnitude. A body's `flux_<band>` variable is a *linear* flux —
photocentres weight bodies by it — so convert:

    flux_L = 10^(-0.4 * \$absmag_L(tempK, mass))

or, if the host's flux is pinned to 1.0 so that companion fluxes are contrast
ratios, take the magnitude difference to the host first.

# Mass units
`mass_unit` selects how the second argument is interpreted: `:Msol` (default,
matching v9 body `mass` variables), `:Mjup` (what v8 passed), or `:Mearth`.

Supported bands:
:MKO_Y, :MKO_Z, :MKO_J, :MKO_H, :MKO_K, :MKO_L′, :MKO_M′, :TwoMASS_J, :TwoMASS_H, :TwoMASS_Ks, :Keck_Ks, :Keck_L′, :Keck_Ms, :SDSS_g′, :SDSS_r′, :SDSS_i′, :SDSS_z′, :IRAC_36, :IRAC_45, :IRAC_57, :IRAC_79, :WISE_W1, :WISE_W2, :WISE_W3, :WISE_W4

Supported metalicities:
"+0.0", "-0.5", "+0.5"
"""
function sonora_photometry_interpolator(band, metalicity="+0.0"; jwst=false,
                                        catalog=datadep"SonoraBobcatEvoPhot",
                                        mass_unit::Symbol=:Msol)

    mscale = _mass_to_mjup(mass_unit)

    #  Load Sonora magnitude table
    if jwst
        mag_table = load_table(joinpath(catalog, "photometry_tables", "mag_table_JWST" * metalicity); jwst)
    else
        mag_table = load_table(joinpath(catalog, "photometry_tables", "mag_table" * metalicity); jwst)
    end
    # We first use BasicInterpolators to grid the sparsely sampled models
    # Then we use Interpolations to make a fast, autodiff compatible
    # linear interpolator of the resulting data.
    points = [
        mag_table.Teff ./ 10 mag_table.mass
    ]
    if !(hasproperty(mag_table, band))
        error("not a valid band: $(keys(mag_table))")
    end
    samples = collect(mag_table[band])

    sitp = RBFInterpolator(points, samples, 2)
    minmass, maxmass = extrema(mag_table.mass)
    minT, maxT = extrema(mag_table.Teff)
    temp_mass_to_abs_mag = (teffk, mass_mjup) -> sitp(teffk / 10, mass_mjup)
    teff = range(minT, maxT, length=200)
    mass = range(minmass, maxmass, length=200)
    gridded = temp_mass_to_abs_mag.(teff, mass')

    # Now build linear interpolator
    litp = LinearInterpolation((teff ./ 10, mass), gridded, extrapolation_bc=NaN)
    function model_interpolator(teffk, mass)
        m = mass * mscale
        if minT <= teffk <= maxT && minmass < m < maxmass
            return litp(teffk / 10, m)
        else
            # Typed rather than a bare `Float64` NaN: this is called from inside
            # `@variables` blocks under ForwardDiff, and a `Union{Float64,Dual}`
            # return infects everything downstream of it in the hot loop.
            return convert(promote_type(typeof(teffk), typeof(m)), NaN)
        end
    end

    return model_interpolator
end
export sonora_photometry_interpolator


"""
    itp = sonora_cooling_interpolator()

Create a function mapping (age_Myr, mass) -> temp_K using Sonora Bobcat
cooling model grids.

    cooling = sonora_cooling_interpolator()
    cooling(15.0, 12mjup)         # mass in M⊙, like every mass in v9

Out-of-grid inputs give `NaN` rather than an extrapolation. `mass_unit`
selects the input unit: `:Msol` (default), `:Mjup` (the v8 convention), or
`:Mearth`.
"""
function sonora_cooling_interpolator(metalicity="+0.0";
                                     catalog=datadep"SonoraBobcatEvoPhot",
                                     mass_unit::Symbol=:Msol)

    mscale = _mass_to_mjup(mass_unit)

    # Load Sonora cooling track
    valid_lines = filter(readlines(joinpath(catalog, "evolution_tables", "evo_tables" * metalicity, "nc$(metalicity)_co1.0_age"))) do line
        length(line) > 10
    end
    headers = split(valid_lines[1], r"  +")[2:7]
    headers = lowercase.(replace.(headers, r"\W" => ""))
    grid = mapreduce(vcat, valid_lines[2:end]) do line
        permutedims(parse.(Float64, split(line, r"  +")[2:7]))
    end

    iteffk = findfirst(==("teffk"), headers)
    immsun = findfirst(==("mmsun"), headers)
    iagegyr = findfirst(==("agegyr"), headers)
    agegyr = grid[:, iagegyr]
    agemyr = agegyr .* 1e3
    mmsun = grid[:, immsun]
    mmjup = mmsun ./ PlanetOrbits.mjup
    teffk = grid[:, iteffk]

    points = [
        log.(agemyr) mmjup
    ]
    samples = collect(teffk)

    sitplog = RBFInterpolator(points, samples, 0.5)

    sitp = (agemyr, mmjup) -> sitplog(log(agemyr), mmjup)

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
export sonora_cooling_interpolator


"""
    Octofitter.load_table(fname; jwst=false)

Read one of the **Sonora Bobcat** evolution/photometry grid files. It parses
that grid's exact fixed layout — eight lines of preamble, then a two-line
split header whose columns are recombined — so it is not a general table
loader and will fail on an ordinary CSV with a `BoundsError`. For a CSV use
`CSV.read(fname, Table)`.

`jwst=true` selects the wider JWST filter set's column grouping.
"""
function load_table(fname; jwst=false)

    headers = open(fname, lock=false, read=true) do f
        for _ in 1:8
            readline(f)
        end
        h1 = readline(f)
        h2 = readline(f)

        headers_1 = strip.(split(h1, '|'))
        headers_2 = strip.(split(h2, r"  +"))

        # Starting after 6
        if jwst
            headers_2_1_indices =
                [
                    fill(2, 29)
                    fill(3, 14)
                ]
        else
            headers_2_1_indices =
                [2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7]
        end

        combined = vcat(
            headers_2[1:6],
            headers_1[headers_2_1_indices] .* '_' .* headers_2[7:end],
        )
        # Normalize headers
        combined = replace.(combined, '\'' => '′')
        combined = replace.(combined, ' ' => "")
        combined = replace.(combined, "/" => '_')
        combined = replace.(combined, "2MASS" => "TwoMASS")
        combined = replace.(combined, "JWST" => "")
        combined = replace.(combined, "NIRCamNIRCamNIRCamNIRCam" => "NIRCam")
        combined = replace.(combined, "MIRIMIRI" => "MIRI")
        combined = replace.(combined, r"[^\w′_]" => "")
        return combined
    end

    data = readdlm(fname, String, skipstart=10, header=false)
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
