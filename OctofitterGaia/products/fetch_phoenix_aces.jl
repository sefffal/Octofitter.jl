# Fetch PHOENIX-ACES HiRes spectra (Husser et al. 2013, A&A 553, A6) and reduce
# each one to the two pieces the node products are built from:
#
#   sed   mean flux density in 0.1 nm bins over 300–1100 nm (8000 bins)
#   rvs   the native-sampling cut-out 842–875 nm
#
# Flux is the surface flux density per unit wavelength, erg s⁻¹ cm⁻² cm⁻¹, on
# vacuum wavelengths, as distributed. The 6 MB originals are not kept.
#
# One small file per spectrum goes to `<out>/staging/aces/<Z>/`, so the run can
# be interrupted and resumed. `build_sed_nodes.jl` consolidates them.
#
# Usage:
#   julia --project=. fetch_phoenix_aces.jl <out dir> [ntasks=8] [Z-1.0 Z-0.5 Z-0.0 Z+0.5]

using Downloads, FITSIO, HDF5

const BASE = "https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS"
const LIB = "PHOENIX-ACES-AGSS-COND-2011"

const SED_LO, SED_HI, SED_STEP = 3000.0, 11000.0, 1.0     # Å
const RVS_LO, RVS_HI = 8420.0, 8750.0                       # Å

function fetch(url, path; tries=5)
    for k in 1:tries
        try
            Downloads.download(url, path; timeout=600)
            return path
        catch err
            k == tries && rethrow()
            sleep(5k)
        end
    end
end

function listing(zdir)
    html = read(fetch("$BASE/$LIB/$(replace(zdir, "+" => "%2B"))/", tempname()), String)
    return sort(unique(m.captures[1] for m in eachmatch(r"href=\"(lte[^\"]*HiRes\.fits)\"", html)))
end

# Bin-averaging weights: each native sample contributes flux × Δλ to its bin.
struct Reducer
    bin::Vector{Int32}        # 0 outside the SED range
    dλ::Vector{Float64}
    norm::Vector{Float64}     # Σ Δλ per bin
    rvs::UnitRange{Int}
end

function Reducer(wave)
    n = length(wave)
    dλ = similar(wave)
    dλ[1] = wave[2] - wave[1]
    dλ[n] = wave[n] - wave[n-1]
    for i in 2:n-1
        dλ[i] = (wave[i+1] - wave[i-1]) / 2
    end
    nbin = round(Int, (SED_HI - SED_LO) / SED_STEP)
    bin = zeros(Int32, n)
    norm = zeros(nbin)
    for i in 1:n
        (SED_LO <= wave[i] < SED_HI) || continue
        b = floor(Int, (wave[i] - SED_LO) / SED_STEP) + 1
        bin[i] = b
        norm[b] += dλ[i]
    end
    all(>(0), norm) || error("an SED bin holds no native sample")
    return Reducer(bin, dλ, norm, searchsortedfirst(wave, RVS_LO):searchsortedlast(wave, RVS_HI))
end

function reduce_spectrum(r::Reducer, flux)
    sed = zeros(length(r.norm))
    @inbounds for i in eachindex(flux)
        b = r.bin[i]
        b == 0 && continue
        sed[b] += flux[i] * r.dλ[i]
    end
    sed ./= r.norm
    return Float32.(sed), Float32.(flux[r.rvs])
end

function main(out, ntasks, zdirs)
    stage = joinpath(out, "staging", "aces")
    mkpath(stage)
    wavefile = joinpath(stage, "WAVE_$LIB.fits")
    isfile(wavefile) || fetch("$BASE/WAVE_$LIB.fits", wavefile)
    wave = FITS(f -> Float64.(read(f[1])), wavefile)
    r = Reducer(wave)
    h5open(joinpath(stage, "grid.h5"), "w") do f
        f["sed_edges_angstrom"] = collect(SED_LO:SED_STEP:SED_HI)
        f["rvs_wave_angstrom"] = wave[r.rvs]
    end
    for zdir in zdirs
        files = listing(zdir)
        dir = joinpath(stage, zdir)
        mkpath(dir)
        todo = filter(f -> !isfile(joinpath(dir, replace(f, ".fits" => ".h5"))), files)
        println("$zdir: $(length(files)) spectra, $(length(todo)) to fetch"); flush(stdout)
        done = Threads.Atomic{Int}(0)
        asyncmap(todo; ntasks) do name
            tmp = tempname() * ".fits"
            fetch("$BASE/$LIB/$(replace(zdir, "+" => "%2B"))/$name", tmp)
            flux, hdr = FITS(tmp) do f
                read(f[1]), read_header(f[1])
            end
            rm(tmp)
            length(flux) == length(wave) || error("$name: $(length(flux)) samples, wavelength grid has $(length(wave))")
            sed, rvs = reduce_spectrum(r, flux)
            part = joinpath(dir, replace(name, ".fits" => ".h5.part"))
            h5open(part, "w") do f
                f["sed"] = sed
                f["rvs"] = rvs
                for k in ("PHXTEFF", "PHXLOGG", "PHXM_H", "PHXALPHA", "PHXREFF", "PHXLUM", "PHXMASS", "PHXXI_L")
                    haskey(hdr, k) && (attributes(f)[k] = Float64(hdr[k]))
                end
                haskey(hdr, "BUNIT") && (attributes(f)["BUNIT"] = String(hdr["BUNIT"]))
            end
            mv(part, replace(part, ".part" => ""); force=true)
            n = Threads.atomic_add!(done, 1) + 1
            n % 50 == 0 && (println("  $zdir $n / $(length(todo))"); flush(stdout))
        end
    end
    println("done")
end

out = ARGS[1]
ntasks = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
zdirs = length(ARGS) >= 3 ? ARGS[3:end] : ["Z-1.0", "Z-0.5", "Z-0.0", "Z+0.5"]
main(out, ntasks, zdirs)
