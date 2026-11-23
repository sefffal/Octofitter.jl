# Fetch BOSZ 2024 spectra (Mészáros et al. 2024, A&A 688, A197; the files
# recomputed in September 2025) from MAST and reduce each one to the staging
# layout described in `README.md`.
#
# What is taken: solar-scaled composition ([α/M] = [C/M] = 0), the R = 50 000
# resampled files. BOSZ is on air wavelengths and gives the Eddington flux H
# in erg s⁻¹ cm⁻² Å⁻¹; the staging files are on vacuum wavelengths with the
# surface flux density F = 4πH in erg s⁻¹ cm⁻² cm⁻¹.
#
# One (Teff, log g, [M/H]) node can exist several times over:
#
#   atmosphere        MARCS (2800–8000 K) and ATLAS9 (7500–16 000 K) overlap
#                     at 7500–8000 K. MARCS is the product there.
#   microturbulence   0, 1, 2 and 4 km/s. The product is 1 km/s at
#                     log g ≥ 3.5 and 2 km/s below.
#
# The product set goes to `<out>/staging/bosz/<m dir>/`. The other choice of
# each (ATLAS9 in the overlap; 2 km/s for dwarfs and 1 km/s for giants) goes
# to `<out>/staging/bosz-alt/<m dir>/`, which is a staging directory of the
# same layout, kept for measuring what those choices are worth.
#
# Usage:
#   julia --project=. fetch_bosz2024.jl <out dir> [ntasks=8] [m+0.00 m-0.50 ...]

using Downloads, CodecZlib, HDF5
include(joinpath(@__DIR__, "common.jl"))

const BASE = "https://archive.stsci.edu/hlsps/bosz/bosz2024"
const RES = "r50000"

const SED_LO, SED_HI, SED_STEP = 3000.0, 11000.0, 1.0     # Å, vacuum
const RVS_LO, RVS_HI = 8420.0, 8750.0                       # Å, vacuum

# Air to vacuum, Edlén (1966), the IAU standard that spectrum synthesis codes
# apply above 200 nm: n − 1 = 10⁻⁶ [64.328 + 29498.1/(146 − σ²) + 255.4/(41 − σ²)]
# with σ = 1/λ_vac in µm⁻¹. The formula is in the vacuum wavelength, so it is
# iterated. The modern (Ciddor 1996) index differs by about 1e-8 here: 3 m/s.
function air_to_vacuum(λair)
    λ = λair
    for _ in 1:4
        σ2 = (1e4 / λ)^2
        n = 1 + 1e-6 * (64.328 + 29498.1 / (146 - σ2) + 255.4 / (41 - σ2))
        λ = λair * n
    end
    return λ
end

product_vmicro(logg) = logg >= 3.5 ? 1 : 2

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

function listing(mdir)
    html = read(fetch("$BASE/$RES/$(replace(mdir, "+" => "%2B"))/", tempname()), String)
    pat = Regex("href=\"(bosz2024_([a-z]+)_t(\\d+)_g([+-][\\d.]+)_m([+-][\\d.]+)_a\\+0\\.00_c\\+0\\.00_v([12])_$(RES)_resam\\.txt\\.gz)\"")
    rows = [(; name=String(m.captures[1]), atmos=String(m.captures[2]), teff=parse(Float64, m.captures[3]),
        logg=parse(Float64, m.captures[4]), mh=parse(Float64, m.captures[5]), vmicro=parse(Int, m.captures[6]))
            for m in eachmatch(pat, html)]
    isempty(rows) && error("no spectra listed for $mdir")
    return unique(r -> r.name, rows)
end

# The product is MARCS where both families exist, at the product's microturbulence.
function is_product(r, rows)
    r.vmicro == product_vmicro(r.logg) || return false
    r.atmos == "ap" || return true
    return !any(o -> o.atmos != "ap" && o.teff == r.teff && o.logg == r.logg && o.vmicro == r.vmicro, rows)
end

function read_flux(path, n)
    flux = Vector{Float64}(undef, n)
    i = 0
    open(GzipDecompressorStream, path) do io
        for line in eachline(io)
            i += 1
            i <= n || error("$path has more than $n rows")
            flux[i] = parse(Float64, SubString(line, 1, findfirst(' ', line) - 1))
        end
    end
    i == n || error("$path has $i rows, the wavelength grid has $n")
    return flux
end

function main(out, ntasks, mdirs)
    stage, alt = joinpath(out, "staging", "bosz"), joinpath(out, "staging", "bosz-alt")
    mkpath(stage)
    mkpath(alt)
    wavefile = joinpath(stage, "bosz2024_wave_$RES.txt")
    isfile(wavefile) || fetch("$BASE/wavelength_grids/bosz2024_wave_$RES.txt", wavefile)
    wave = air_to_vacuum.(parse.(Float64, eachline(wavefile)))
    r = Reducer(wave, SED_LO, SED_HI, SED_STEP, RVS_LO, RVS_HI)
    for dir in (stage, alt)
        h5open(joinpath(dir, "grid.h5"), "w") do f
            f["sed_edges_angstrom"] = collect(SED_LO:SED_STEP:SED_HI)
            f["rvs_wave_angstrom"] = wave[r.rvs]
            a = attributes(f)
            a["library"] = "BOSZ 2024 (R = 50 000, solar-scaled, 2025 recomputation)"
            a["reference"] = "Mészáros et al. 2024, A&A 688, A197"
            a["source"] = BASE
            a["resolving_power"] = 50000.0
        end
    end
    for mdir in mdirs
        rows = listing(mdir)
        dest(row) = joinpath(is_product(row, rows) ? stage : alt, mdir, replace(row.name, ".txt.gz" => ".h5"))
        foreach(d -> mkpath(joinpath(d, mdir)), (stage, alt))
        todo = filter(row -> !isfile(dest(row)), rows)
        println("$mdir: $(length(rows)) spectra ($(count(row -> is_product(row, rows), rows)) in the product), $(length(todo)) to fetch"); flush(stdout)
        done = Threads.Atomic{Int}(0)
        asyncmap(todo; ntasks) do row
            tmp = tempname() * ".txt.gz"
            fetch("$BASE/$RES/$(replace(mdir, "+" => "%2B"))/$(replace(row.name, "+" => "%2B"))", tmp)
            flux = read_flux(tmp, length(wave)) .* (4π * 1e8)
            rm(tmp)
            sed, rvs = reduce_spectrum(r, flux)
            part = dest(row) * ".part"
            h5open(part, "w") do f
                f["sed"] = sed
                f["rvs"] = rvs
                a = attributes(f)
                a["PHXTEFF"] = row.teff
                a["PHXLOGG"] = row.logg
                a["PHXM_H"] = row.mh
                a["PHXALPHA"] = 0.0
                a["PHXXI_L"] = Float64(row.vmicro)
                a["ATMOS"] = row.atmos
            end
            mv(part, dest(row); force=true)
            n = Threads.atomic_add!(done, 1) + 1
            n % 100 == 0 && (println("  $mdir $n / $(length(todo))"); flush(stdout))
        end
    end
    println("done")
end

out = ARGS[1]
ntasks = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 8
mdirs = length(ARGS) >= 3 ? ARGS[3:end] : ["m+0.00", "m-0.50", "m+0.50", "m-1.00", "m-0.25", "m+0.25", "m-0.75"]
main(out, ntasks, mdirs)
