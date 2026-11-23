# Build the RVS-window node product from the reduced library spectra written
# by a fetch script (`fetch_phoenix_aces.jl`).
#
# Each model's cut-out is brought to a Gaussian line-spread function of
# resolving power R (libraries distributed already broadened are convolved
# with the difference in quadrature) and put on a grid uniform in ln λ. The log fluxes of all models are then compressed by a principal
# component analysis: log10 F(λ) = mean(λ) + Σ_i w_i e_i(λ), and what is
# interpolated between models is the handful of weights w_i.
#
# Output: one HDF5 file that `OctofitterGaia.RVSLibrary` loads,
#
#   teff, logg, mh                     node axes
#   ln_wave0, dln_wave, n_wave         the wavelength grid (vacuum, λ in nm)
#   mean[wave]                         mean log10 surface flux density, erg s⁻¹ cm⁻² cm⁻¹
#   basis[wave, component]             orthonormal eigenspectra
#   weights[component, teff, logg, mh]
#   present[teff, logg, mh]            1 where the library has that model
#   rms_wave[wave], rms_node[teff, logg, mh]
#                                      rms of the truncation residual in log10 flux,
#                                      over real models and over wavelength
#
# Missing nodes of the tensor grid are filled as in `build_sed_nodes.jl`.
#
# Usage:
#   julia --project=. build_rvs_nodes.jl <staging dir> <out.h5> [ncomp=40] [R=11500] [step_kms=2] [lo_nm=845] [hi_nm=871.5]

using HDF5, LinearAlgebra, Statistics
include(joinpath(@__DIR__, "common.jl"))

function main(stage, out, ncomp, R, step, lo, hi)
    wave = h5read(joinpath(stage, "grid.h5"), "rvs_wave_angstrom") ./ 10     # nm
    meta = h5open(f -> Dict(k => read(attributes(f)[k]) for k in keys(attributes(f))), joinpath(stage, "grid.h5"))
    R_native = meta["resolving_power"]
    dlnw = step / C_KMS
    lnw0 = log(lo)
    n = floor(Int, (log(hi) - lnw0) / dlnw) + 1
    rows = NamedTuple[]
    for zdir in filter(d -> isdir(joinpath(stage, d)), readdir(stage))
        for name in filter(endswith(".h5"), readdir(joinpath(stage, zdir)))
            h5open(joinpath(stage, zdir, name)) do f
                a = attributes(f)
                (haskey(a, "PHXALPHA") ? read(a["PHXALPHA"]) : 0.0) == 0 || return
                flux = lsf_resample(wave, Float64.(read(f["rvs"])), lnw0, dlnw, n; R, R_native)
                all(>(0), flux) || error("non-positive flux in $name")
                push!(rows, (; teff=read(a["PHXTEFF"]), logg=read(a["PHXLOGG"]), mh=read(a["PHXM_H"]), lf=log10.(flux)))
            end
        end
    end
    teff = sort(unique(r.teff for r in rows))
    logg = sort(unique(r.logg for r in rows))
    mh = sort(unique(r.mh for r in rows))

    X = reduce(hcat, (r.lf for r in rows))                 # [wave, model]
    μ = vec(mean(X; dims=2))
    F = svd(X .- μ)
    basis = F.U[:, 1:ncomp]
    W = basis' * (X .- μ)                                  # [component, model]
    resid = X .- μ .- basis * W
    rms_wave = vec(sqrt.(mean(abs2, resid; dims=2)))
    rms_model = vec(sqrt.(mean(abs2, resid; dims=1)))
    println("singular values (relative): ", round.(F.S[[1, 5, 10, 20, 30, 40, 60, 80]] ./ F.S[1]; sigdigits=3))
    for k in (10, 20, 30, 40, 60)
        k <= ncomp || continue
        r = X .- μ .- F.U[:, 1:k] * (F.U[:, 1:k]' * (X .- μ))
        println("  $k components: rms residual $(round(sqrt(mean(abs2, r)) * log(10) * 100; sigdigits=3)) % of flux, ",
            "worst model $(round(maximum(sqrt.(mean(abs2, r; dims=1))) * log(10) * 100; sigdigits=3)) %, ",
            "worst sample $(round(maximum(abs, r) * log(10) * 100; sigdigits=3)) %")
    end

    weights = fill(NaN, ncomp, length(teff), length(logg), length(mh))
    present = fill(false, length(teff), length(logg), length(mh))
    rms_node = fill(NaN, length(teff), length(logg), length(mh))
    for (c, r) in enumerate(rows)
        i, j, k = searchsortedfirst(teff, r.teff), searchsortedfirst(logg, r.logg), searchsortedfirst(mh, r.mh)
        present[i, j, k] && error("duplicate node $(r.teff) $(r.logg) $(r.mh)")
        weights[:, i, j, k] .= W[:, c]
        rms_node[i, j, k] = rms_model[c]
        present[i, j, k] = true
    end
    nfilled = count(!, present)
    fill_missing!(weights, present)

    h5open(out, "w") do f
        foreach(k -> attributes(f)[k] = meta[k], ("library", "reference", "source"))
        attributes(f)["flux_unit"] = "log10 of surface flux density per unit wavelength, erg s-1 cm-2 cm-1"
        attributes(f)["wavelength"] = "vacuum, nm; ln(wave) = ln_wave0 + i * dln_wave, i = 0 … n_wave - 1"
        attributes(f)["lsf"] = "Gaussian, FWHM = wavelength / resolving_power"
        attributes(f)["resolving_power"] = R
        attributes(f)["alpha"] = 0.0
        f["teff"] = teff
        f["logg"] = logg
        f["mh"] = mh
        f["ln_wave0"] = lnw0
        f["dln_wave"] = dlnw
        f["n_wave"] = n
        f["mean"] = μ
        f["basis"] = Float32.(basis)
        f["weights"] = Float32.(weights)
        f["present"] = UInt8.(present)
        f["rms_wave"] = Float32.(rms_wave)
        f["rms_node"] = Float32.(rms_node)
    end
    println("wrote $out: $(length(rows)) models, $nfilled filled nodes, $n samples of $(step) km/s over $(lo)–$(hi) nm, $ncomp components")
end

stage, out = ARGS[1], ARGS[2]
arg(i, T, d) = length(ARGS) >= i ? parse(T, ARGS[i]) : d
main(stage, out, arg(3, Int, 40), arg(4, Float64, 11500.0), arg(5, Float64, 2.0), arg(6, Float64, 845.0), arg(7, Float64, 871.5))
