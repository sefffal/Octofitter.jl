# Shared by the node-product builders.

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

# Reduction of one library spectrum, sampled at `wave` [Å, vacuum], to the two
# staging pieces: the mean flux density in bins of `sed_step` over
# `sed_lo`–`sed_hi`, and the native samples in `rvs_lo`–`rvs_hi`. Each native sample contributes flux × Δλ to its bin.
struct Reducer
    bin::Vector{Int32}        # 0 outside the SED range
    dλ::Vector{Float64}
    norm::Vector{Float64}     # Σ Δλ per bin
    rvs::UnitRange{Int}
end

function Reducer(wave, sed_lo, sed_hi, sed_step, rvs_lo, rvs_hi)
    n = length(wave)
    dλ = similar(wave)
    dλ[1] = wave[2] - wave[1]
    dλ[n] = wave[n] - wave[n-1]
    for i in 2:n-1
        dλ[i] = (wave[i+1] - wave[i-1]) / 2
    end
    nbin = round(Int, (sed_hi - sed_lo) / sed_step)
    bin = zeros(Int32, n)
    norm = zeros(nbin)
    for i in 1:n
        (sed_lo <= wave[i] < sed_hi) || continue
        b = floor(Int, (wave[i] - sed_lo) / sed_step) + 1
        bin[i] = b
        norm[b] += dλ[i]
    end
    all(>(0), norm) || error("an SED bin holds no native sample")
    return Reducer(bin, dλ, norm, searchsortedfirst(wave, rvs_lo):searchsortedlast(wave, rvs_hi))
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

const C_KMS = 299792.458

# Convolve a spectrum sampled at `wave` with a Gaussian of resolving power R
# (FWHM = λ/R) and return it on the grid ln λ = lnw0 + i·dlnw. The work is
# done on a grid `over` times finer, so the output samples are averages over
# their own width.
function lsf_resample(wave, flux, lnw0, dlnw, n; R, R_native=Inf, over=4)
    h = dlnw / over
    R_native > R || error("the library's resolving power $R_native is not above the requested $R")
    σ = sqrt(1 / R^2 - 1 / R_native^2) / (2sqrt(2log(2)))      # in ln λ
    half = ceil(Int, 5σ / h)
    kernel = [exp(-0.5 * (k * h / σ)^2) for k in -half:half]
    kernel ./= sum(kernel)
    nfine = n * over + 2half
    x0 = lnw0 - dlnw / 2 + h / 2 - half * h             # first fine sample
    fine = Vector{Float64}(undef, nfine)
    lnwave = log.(wave)
    j = 1
    for i in 1:nfine
        x = x0 + (i - 1) * h
        while j < length(lnwave) - 1 && lnwave[j+1] < x
            j += 1
        end
        (lnwave[j] <= x <= lnwave[j+1]) || error("the cut-out does not cover the requested range plus the LSF wings")
        t = (x - lnwave[j]) / (lnwave[j+1] - lnwave[j])
        fine[i] = (1 - t) * flux[j] + t * flux[j+1]
    end
    out = zeros(n)
    for i in 1:n, o in 1:over
        c = (i - 1) * over + o + half                   # centre in `fine`
        acc = 0.0
        @inbounds for k in -half:half
            acc += kernel[k+half+1] * fine[c+k]
        end
        out[i] += acc / over
    end
    return out
end
