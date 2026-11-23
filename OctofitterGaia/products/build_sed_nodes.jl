# Build the low-resolution SED node product from the reduced library spectra
# written by a fetch script (`fetch_phoenix_aces.jl`).
#
# Output: one HDF5 file that `OctofitterGaia.SpectralLibrary` loads,
#
#   teff [K], logg [cgs], mh [dex]          node axes
#   edges_nm                                 bin edges, vacuum
#   log_flux[bin, teff, logg, mh]            log10 of the bin-averaged surface flux
#                                            density, erg s⁻¹ cm⁻² cm⁻¹
#   present[teff, logg, mh]                  1 where the library has that model
#
# The node grid is the full tensor product of the axes. Libraries are not
# rectangular (the lowest gravities do not exist for hot stars), so the
# missing nodes are filled by extrapolation along log g, then Teff. A filled
# node is never interpolated *to*: `present` marks the real ones and a lookup
# is valid only inside a cell whose eight corners are real. The filled values
# exist so that the cubic stencil has something smooth to take its end
# slopes from.
#
# Usage:
#   julia --project=. build_sed_nodes.jl <staging dir> <out.h5> [lo_nm=320] [hi_nm=1100] [step_nm=5]

using HDF5

function load_staging(stage)
    grid = h5open(joinpath(stage, "grid.h5")) do f
        (; edges=read(f["sed_edges_angstrom"]) ./ 10)
    end
    rows = NamedTuple[]
    for zdir in filter(d -> isdir(joinpath(stage, d)), readdir(stage))
        for name in filter(endswith(".h5"), readdir(joinpath(stage, zdir)))
            h5open(joinpath(stage, zdir, name)) do f
                a = attributes(f)
                alpha = haskey(a, "PHXALPHA") ? read(a["PHXALPHA"]) : 0.0
                push!(rows, (; teff=read(a["PHXTEFF"]), logg=read(a["PHXLOGG"]),
                    mh=read(a["PHXM_H"]), alpha, sed=Float64.(read(f["sed"]))))
            end
        end
    end
    return grid, rows
end

# Average the fine staging bins into the product bins. Both are uniform in
# wavelength and the product edges are a subset of the staging edges.
function rebin(sed, fine_edges, edges)
    out = zeros(length(edges) - 1)
    for b in eachindex(out)
        i0 = findfirst(e -> isapprox(e, edges[b]; atol=1e-6), fine_edges)
        i1 = findfirst(e -> isapprox(e, edges[b+1]; atol=1e-6), fine_edges)
        (i0 === nothing || i1 === nothing) && error("product bin edge $(edges[b]) nm is not a staging bin edge")
        out[b] = sum(@view sed[i0:i1-1]) / (i1 - i0)
    end
    return out
end

# Fill the missing nodes of `a[:, i, j, k]`, in place: two-point linear
# extrapolation along axis 3 (log g) where possible, then along axis 2
# (Teff), repeated until nothing is missing.
function fill_missing!(a, have)
    have = copy(have)
    nT, ng, nm = size(have)
    while !all(have)
        progress = false
        for axis in (3, 2)
            snapshot = copy(have)
            for k in 1:nm, j in 1:ng, i in 1:nT
                snapshot[i, j, k] && continue
                step(d) = axis == 3 ? (i, j + d, k) : (i + d, j, k)
                inb(I) = 1 <= I[1] <= nT && 1 <= I[2] <= ng
                for dir in (1, -1)
                    p1, p2 = step(dir), step(2dir)
                    (inb(p1) && inb(p2) && snapshot[p1...] && snapshot[p2...]) || continue
                    @views a[:, i, j, k] .= 2 .* a[:, p1...] .- a[:, p2...]
                    have[i, j, k] = true
                    progress = true
                    break
                end
            end
            progress && break
        end
        progress || error("could not fill every missing node")
    end
    return a
end

function main(stage, out, lo, hi, step)
    grid, rows = load_staging(stage)
    rows = filter(r -> r.alpha == 0, rows)
    teff = sort(unique(r.teff for r in rows))
    logg = sort(unique(r.logg for r in rows))
    mh = sort(unique(r.mh for r in rows))
    edges = collect(lo:step:hi)
    nb = length(edges) - 1
    log_flux = fill(NaN, nb, length(teff), length(logg), length(mh))
    present = fill(false, length(teff), length(logg), length(mh))
    for r in rows
        i, j, k = searchsortedfirst(teff, r.teff), searchsortedfirst(logg, r.logg), searchsortedfirst(mh, r.mh)
        present[i, j, k] && error("duplicate node $(r.teff) $(r.logg) $(r.mh)")
        f = rebin(r.sed, grid.edges, edges)
        all(>(0), f) || error("non-positive flux at node $(r.teff) $(r.logg) $(r.mh)")
        log_flux[:, i, j, k] .= log10.(f)
        present[i, j, k] = true
    end
    nfilled = count(!, present)
    fill_missing!(log_flux, present)
    h5open(out, "w") do f
        attributes(f)["library"] = "PHOENIX-ACES-AGSS-COND-2011 (HiRes)"
        attributes(f)["reference"] = "Husser et al. 2013, A&A 553, A6"
        attributes(f)["source"] = "https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS"
        attributes(f)["flux_unit"] = "log10 of surface flux density per unit wavelength, erg s-1 cm-2 cm-1, bin average"
        attributes(f)["wavelength"] = "vacuum, nm"
        attributes(f)["alpha"] = 0.0
        f["teff"] = teff
        f["logg"] = logg
        f["mh"] = mh
        f["edges_nm"] = edges
        f["log_flux", chunk=(nb, length(teff), 1, 1), deflate=4] = Float32.(log_flux)
        f["present"] = UInt8.(present)
    end
    println("wrote $out: $(count(present)) models, $nfilled filled nodes, $nb bins $(lo)–$(hi) nm; ",
        "Teff $(extrema(teff)), log g $(extrema(logg)), [M/H] $mh")
end

stage, out = ARGS[1], ARGS[2]
lo = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 320.0
hi = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1100.0
step = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 5.0
main(stage, out, lo, hi, step)
