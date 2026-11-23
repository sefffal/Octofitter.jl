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
