# ---------------------------------------------------
# Image block, part 2: the likelihood of epoch images
#
# One transit of a source is up to nine AF windows. The samples y of window
# q are modelled as
#
#     y = μ(θ) + B_q(θ) a_q + T_q(θ) c + ε,    ε ~ N(0, diag(μ + σ_read² + (ϵ × sources)²))
#
# with μ the scene (background plus every point source through the device's
# chromatic image basis) and the columns of B_q and T_q the directions in
# which the scene can be wrong without the system being different: the
# background level, the flux scale, a shift along and across scan, and the
# device's shape modes. The coefficients a_q are the window's own, c is
# shared by the windows of the transit (the attitude, the across-scan rate
# and the focus are the same for all nine strips). Both are Gaussian with
# fixed widths and are marginalised analytically, so the covariance of a
# transit is diagonal plus nested low rank and never formed: each window
# reduces to the Gram matrix of a handful of vectors, the window's
# coefficients are eliminated from it, and what is left is summed over the
# windows and solved for the transit.
#
# The prior width of the shift is what the reduced astrometry calls excess
# noise; here it is a term of the pixel model.
#
# A transit is either like that or bad: the same model with a broad extra
# variance on every sample. The two states are summed.
# ---------------------------------------------------

"""
    ImageInstrument(bases; kwargs...)
    ImageInstrument(path; devices = nothing, kwargs...)

What the image block assumes about the instrument: the chromatic image
basis of every device (a `Dict` from device name to [`ImageBasis`](@ref), or
an instrument product to read them from) and the widths of everything that
is marginalised. Nothing here is sampled; the widths are placeholders to be
calibrated on the population.

  - `read_noise = 4.48` [e⁻]: per sample.
  - `model_error = 0.0`: what the image model is not trusted beyond, as a
    fraction of the source electrons in a sample; its square times those
    electrons squared is added to the variance of every sample, without
    correlation between samples. Photon noise is `1 / sqrt(electrons)`, so
    this is what limits a bright star.
  - `window_flux_sigma = 0.005`, `transit_flux_sigma = 0.02`: fractional
    flux scale of one window, and of all the windows of a transit.
  - `window_shift_sigma = (0.05, 5.0)`, `transit_shift_sigma = (0.1, 5.0)`
    [mas, along and across scan]: displacement of the whole scene in one
    window (the geometric calibration of a device) and in a transit (the
    attitude).
  - `window_shape_sigma = 0.0`, `transit_shape_sigma = 1.0`: coefficients of
    the devices' shape modes, in units of the calibration's scatter.
  - `nshapes`: how many shape modes to use; by default all that every device
    has.
  - `bad_sigma = 50.0` [e⁻], `bad_fraction = 0.05`: the bad state adds
    `bad_sigma² + (bad_fraction × source electrons)²` to the variance of
    every sample.

A width of zero removes that mode. The background's width is per window,
from the data.
"""
struct ImageInstrument
    bases::Dict{String,ImageBasis}
    read_noise::Float64
    model_error::Float64
    window_flux_sigma::Float64
    window_shift_sigma::NTuple{2,Float64}
    window_shape_sigma::Float64
    transit_flux_sigma::Float64
    transit_shift_sigma::NTuple{2,Float64}
    transit_shape_sigma::Float64
    nshapes::Int
    bad_sigma::Float64
    bad_fraction::Float64
end

function ImageInstrument(bases::AbstractDict; read_noise=4.48, model_error=0.0, window_flux_sigma=0.005, window_shift_sigma=(0.05, 5.0),
    window_shape_sigma=0.0, transit_flux_sigma=0.02, transit_shift_sigma=(0.1, 5.0), transit_shape_sigma=1.0,
    nshapes=minimum(OctofitterGaia.nshapes, values(bases)), bad_sigma=50.0, bad_fraction=0.05)
    isempty(bases) && throw(ArgumentError("an instrument needs at least one device"))
    0 <= nshapes <= minimum(OctofitterGaia.nshapes, values(bases)) ||
        throw(ArgumentError("nshapes = $nshapes, but a device has only $(minimum(OctofitterGaia.nshapes, values(bases))) shape modes"))
    allequal(nterms(b) for b in values(bases)) || throw(ArgumentError("the devices' bases differ in chromatic degree"))
    allequal(b.wavenumber_range for b in values(bases)) || throw(ArgumentError("the devices' bases differ in wavenumber range"))
    (nshapes == 0 || allequal(size(b.modes, 2) for b in values(bases))) ||
        throw(ArgumentError("the devices' shape modes differ in chromatic degree"))
    all(>=(0), (read_noise, model_error, window_flux_sigma, window_shift_sigma..., window_shape_sigma, transit_flux_sigma,
        transit_shift_sigma..., transit_shape_sigma, bad_fraction)) || throw(ArgumentError("widths must not be negative"))
    bad_sigma > 0 || throw(ArgumentError("bad_sigma must be positive"))
    return ImageInstrument(Dict{String,ImageBasis}(String(k) => v for (k, v) in bases), read_noise, model_error, window_flux_sigma,
        Float64.(Tuple(window_shift_sigma)), window_shape_sigma, transit_flux_sigma, Float64.(Tuple(transit_shift_sigma)),
        transit_shape_sigma, nshapes, bad_sigma, bad_fraction)
end

function ImageInstrument(path::AbstractString; devices=nothing, kwargs...)
    names = devices === nothing ? h5open(f -> [k for k in keys(f) if f[k] isa HDF5.Group && haskey(f[k], "terms")], path) : devices
    return ImageInstrument(Dict(String(n) => ImageBasis(path, String(n)) for n in names); kwargs...)
end

"A copy of the instrument with some widths changed."
ImageInstrument(inst::ImageInstrument; kwargs...) = ImageInstrument(inst.bases;
    (k => getfield(inst, k) for k in fieldnames(ImageInstrument) if k != :bases)..., kwargs...)

function Base.show(io::IO, inst::ImageInstrument)
    print(io, "ImageInstrument(", length(inst.bases), " devices, ", inst.nshapes, " shape modes)")
end

"Degree of the chromatic expansion of the instrument's bases."
chromatic_degree(inst::ImageInstrument) = nterms(first(values(inst.bases))) - 1

# One window: where its samples are, what it was exposed with.
struct _ImageWindow
    basis::Int                  # index into `EpochImages.bases`
    n_al::Int
    n_ac::Int
    w0::Float64                 # centre of sample (0, 0), mas
    z0::Float64
    dw::Float64                 # sample size
    dz::Float64
    cw::Float64                 # shear: along-scan shift per across-scan index, and the converse
    vz::Float64
    exposure::Float64           # s
    background::Float64         # e⁻ per sample
    background_sigma::Float64
    offset::Int                 # of its samples in `EpochImages.samples`
end

const IMAGE_COLUMNS = (:epoch, :transit, :scan_pos_angle, :parallax_factor_al, :parallax_factor_ac, :device, :exposure,
    :samples, :sample_pos_al, :sample_pos_ac, :sample_size_al, :sample_size_ac, :background, :background_error)

"""
    EpochImages(table, instrument; p_bad = 1e-3)

The epoch images of one source, ready for the likelihood. `table` has one
row per **window** (one CCD of one transit):

  - `epoch` [MJD] and `transit`, any identifier that is the same for the
    windows of one field-of-view transit; rows of a transit must be adjacent;
  - `scan_pos_angle` [**degrees**], of this CCD's crossing; `parallax_factor_al`,
    `parallax_factor_ac`, of the transit;
  - `device`: the name of the window's image basis in `instrument` (field of
    view, CCD strip, gate, window class);
  - `exposure` [s]: of the gate that was used;
  - `samples`: a matrix `[along scan, across scan]` of electrons, de-biased,
    background **not** subtracted, `NaN` where a sample is not to be used;
  - `sample_pos_al`, `sample_pos_ac` [mas]: local-plane coordinates of the
    centre of sample (0, 0), counted from the source's tangent point;
    `sample_size_al`, `sample_size_ac` [mas]; and optionally `shear_term_al`,
    `shear_term_ac` [mas]. Sample (i, j), from zero, is at
    `w = w₀ + i Δw + j c_w`, `z = z₀ + j Δz + i v_z`;
  - `background`, `background_error` [**e⁻ per sample**]: the published
    values are per pixel and second, so multiply by the exposure and the
    pixels per sample.

`p_bad` is the prior probability that a transit is bad, one number or one
per transit.
"""
struct EpochImages{M,NM,NT}             # shape modes, their chromatic terms, the chromatic terms: sizes the inner loops are compiled for
    instrument::ImageInstrument
    bases::Vector{ImageBasis}
    windows::Vector{_ImageWindow}
    samples::Vector{Float64}
    transit_windows::Vector{UnitRange{Int}}
    f_al::Vector{Float64}               # per transit
    f_ac::Vector{Float64}
    p_bad::Vector{Float64}
    window_modes::Vector{Int}           # which of the 4 + nshapes mode vectors each level marginalises
    transit_modes::Vector{Int}
end

ntransits(d::EpochImages) = length(d.transit_windows)
nwindows(d::EpochImages) = length(d.windows)
nwindows(d::EpochImages, t::Integer) = length(d.transit_windows[t])
nshapes(d::EpochImages) = d.instrument.nshapes

function Base.show(io::IO, d::EpochImages)
    print(io, "EpochImages(", ntransits(d), " transits, ", nwindows(d), " windows, ", count(!isnan, d.samples), " valid samples)")
end

# Prior width of mode `a` (1 background, 2 flux, 3 and 4 shift along and across scan, then the shape modes) at each level.
@inline _window_sigma(inst::ImageInstrument, win::_ImageWindow, a) =
    a == 1 ? win.background_sigma : a == 2 ? inst.window_flux_sigma : a <= 4 ? inst.window_shift_sigma[a-2] : inst.window_shape_sigma
@inline _transit_sigma(inst::ImageInstrument, a) =
    a == 1 ? 0.0 : a == 2 ? inst.transit_flux_sigma : a <= 4 ? inst.transit_shift_sigma[a-2] : inst.transit_shape_sigma

function EpochImages(table, instrument::ImageInstrument; p_bad=1e-3)
    cols = Octofitter.Tables.columnnames(table)
    issubset(IMAGE_COLUMNS, cols) || error("EpochImages: the table needs the columns $IMAGE_COLUMNS")
    nrow = length(table.epoch)
    nrow > 0 || throw(ArgumentError("no windows"))
    shear_al = :shear_term_al in cols ? table.shear_term_al : zeros(nrow)
    shear_ac = :shear_term_ac in cols ? table.shear_term_ac : zeros(nrow)
    names = unique(String.(table.device))
    for n in names
        haskey(instrument.bases, n) || throw(ArgumentError("the instrument has no device \"$n\"; it has $(join(sort!(collect(keys(instrument.bases))), ", "))"))
    end
    bases = [instrument.bases[n] for n in names]
    windows = _ImageWindow[]
    samples = Float64[]
    transit_windows = UnitRange{Int}[]
    f_al, f_ac = Float64[], Float64[]
    first_row = 1
    for r in 1:nrow
        s = table.samples[r]
        s isa AbstractMatrix || throw(ArgumentError("row $r: `samples` must be a matrix [along scan, across scan]"))
        table.exposure[r] > 0 || throw(ArgumentError("row $r: the exposure must be positive"))
        table.background_error[r] > 0 || throw(ArgumentError("row $r: the background error must be positive"))
        count(isfinite, s) > 4 + instrument.nshapes || throw(ArgumentError("row $r has $(count(isfinite, s)) valid samples; drop it"))
        push!(windows, _ImageWindow(findfirst(==(String(table.device[r])), names), size(s, 1), size(s, 2),
            table.sample_pos_al[r], table.sample_pos_ac[r], table.sample_size_al[r], table.sample_size_ac[r],
            shear_al[r], shear_ac[r], table.exposure[r], table.background[r], table.background_error[r], length(samples)))
        append!(samples, (isfinite(v) ? Float64(v) : NaN for v in s))
        if r == nrow || table.transit[r+1] != table.transit[r]
            table.transit[r] in view(table.transit, 1:first_row-1) &&
                throw(ArgumentError("the rows of transit $(table.transit[r]) are not adjacent"))
            rows = first_row:r
            (allequal(view(table.parallax_factor_al, rows)) && allequal(view(table.parallax_factor_ac, rows))) ||
                throw(ArgumentError("the parallax factors differ within transit $(table.transit[r])"))
            push!(transit_windows, rows)
            push!(f_al, table.parallax_factor_al[r])
            push!(f_ac, table.parallax_factor_ac[r])
            first_row = r + 1
        end
    end
    nt = length(transit_windows)
    p_bad = p_bad isa Number ? fill(Float64(p_bad), nt) : collect(Float64, p_bad)
    length(p_bad) == nt || throw(DimensionMismatch("p_bad needs one entry per transit ($nt)"))
    all(p -> 0 <= p < 1, p_bad) || throw(ArgumentError("p_bad must be in [0, 1)"))
    nmode = 4 + instrument.nshapes
    window_modes = [a for a in 1:nmode if a == 1 || _window_sigma(instrument, first(windows), a) > 0]
    transit_modes = [a for a in 1:nmode if _transit_sigma(instrument, a) > 0]
    return EpochImages{instrument.nshapes,size(first(bases).modes, 2),nterms(first(bases))}(instrument, bases, windows, samples,
        transit_windows, f_al, f_ac, p_bad, window_modes, transit_modes)
end

# ---------------------------------------------------
# The scene in one sample
# ---------------------------------------------------

# Fill `u[1:4+M]` for sample (i, j) of a window: 1, then the electrons of the
# sources, their derivatives with respect to a shift of the sample along and
# across scan, and their change per unit of each shape mode. `ws`, `zs`,
# `electrons` and `cvecs` are the sources' positions, electrons and
# (1, moments...) in this window.
@inline function _sample_vector!(u::AbstractVector, basis::ImageBasis, win::_ImageWindow, i, j,
    ws::AbstractVector, zs::AbstractVector, electrons::AbstractVector, cvecs::AbstractVector, ::Val{M}, ::Val{NM}) where {M,NM}
    w = win.w0 + i * win.dw + j * win.cw
    z = win.z0 + j * win.dz + i * win.vz
    T = eltype(u)
    f = fw = fz = zero(T)
    shapes = ntuple(_ -> zero(T), Val(M))
    @inbounds for k in eachindex(ws)
        (pf, pw, pz), ps = _point(basis, cvecs[k], w - ws[k], z - zs[k], Val(M), Val(NM))
        e = electrons[k]
        f += e * pf
        fw += e * pw
        fz += e * pz
        shapes = shapes .+ e .* ps
    end
    @inbounds u[1], u[2], u[3], u[4] = one(T), f, fw, fz
    for m in 1:M
        @inbounds u[4+m] = shapes[m]
    end
    return u
end

@inline _static_sizes(::EpochImages{M,NM,NT}) where {M,NM,NT} = (Val(M), Val(NM), Val(NT))

# Expected electrons of a sample and its variance in the good state, from
# the sample's vector: background plus sources, the star-level shape
# coefficients applied.
@inline function _sample_mean(u::AbstractVector, background, shape)
    m = @inbounds u[2]
    for a in eachindex(shape)
        m += shape[a] * @inbounds u[4+a]
    end
    return background + m, m
end

# ---------------------------------------------------
# Likelihood
# ---------------------------------------------------

# Eliminate the window's coefficients from its Gram matrix `G` (of the 4 + M
# mode vectors and the residual) and add what is left, for the transit's
# modes and the residual, to `H`. Returns log det of the window's posterior
# precision minus that of its prior precision, or -Inf.
function _eliminate_window!(H::AbstractMatrix, G::AbstractMatrix, A::AbstractMatrix, V::AbstractMatrix,
    d::EpochImages, win::_ImageWindow, transit_modes)
    iw = d.window_modes
    Kw, Kt, P = length(iw), length(transit_modes), size(G, 1)
    inst = d.instrument
    logprec = 0.0
    @inbounds for b in 1:Kw
        for a in b:Kw
            A[a, b] = G[iw[a], iw[b]]
        end
        p = 1 / _window_sigma(inst, win, iw[b])^2
        A[b, b] += p
        logprec += log(p)
    end
    logdet = _cholesky_small!(view(A, 1:Kw, 1:Kw))
    logdet == -Inf && return logdet
    @inbounds for c in 1:Kt+1
        jc = c <= Kt ? transit_modes[c] : P
        col = view(V, 1:Kw, c)
        for a in 1:Kw
            col[a] = G[iw[a], jc]
        end
        _forward_small!(col, view(A, 1:Kw, 1:Kw))
        for c2 in 1:c
            j2 = c2 <= Kt ? transit_modes[c2] : P
            acc = G[jc, j2]
            for a in 1:Kw
                acc -= V[a, c] * V[a, c2]
            end
            H[c, c2] += acc
        end
    end
    return logdet - logprec
end

# −2 log density of a transit, up to the constant, from what the windows left in `H`.
function _transit_solve!(H::AbstractMatrix, sigma_t)
    Kt = size(H, 1) - 1
    logprec = zero(eltype(H))
    @inbounds for b in 1:Kt
        p = 1 / sigma_t[b]^2
        H[b, b] += p
        logprec += log(p)
    end
    A = view(H, 1:Kt, 1:Kt)
    logdet = _cholesky_small!(A)
    logdet == -Inf && return oftype(logdet, Inf)
    h = view(H, Kt + 1, 1:Kt)
    _forward_small!(h, A)
    return H[Kt+1, Kt+1] - sum(abs2, h) + logdet - logprec
end

"""
    image_transit_loglike(data, t, w, z, rate, moments; shape = (), flux_scale = 1,
                          jitter = nothing, buf) -> (ll, ll_good)

Log-likelihood of the windows of transit `t` of `data`, given the scene:
`w[k, q]`, `z[k, q]` [mas] the local-plane position of point source `k` in
the transit's `q`-th window, along and across scan; `rate[k]` [e⁻ s⁻¹] its G
photon rate; `moments[n, k]` its chromatic moments (see
[`chromatic_moments!`](@ref)).

`shape` holds coefficients of the devices' shape modes that apply to every
window of the source (the star-level part, sampled); `flux_scale` multiplies
every rate; `jitter = (along, across)` [mas] replaces the instrument's
`transit_shift_sigma`.

`ll` is the sum over the good and bad states, `ll_good` the good state
alone. Works with dual numbers in any argument.
"""
function image_transit_loglike(d::EpochImages, t::Integer, w::AbstractMatrix, z::AbstractMatrix, rate::AbstractVector,
    moments::AbstractMatrix; shape=(), flux_scale=1.0, jitter=nothing, buf=Bumper.default_buffer())
    inst = d.instrument
    M = inst.nshapes
    K = length(rate)
    wins = d.transit_windows[t]
    size(w) == size(z) == (K, length(wins)) || throw(DimensionMismatch("w and z must be [source, window] = ($K, $(length(wins)))"))
    size(moments) == (chromatic_degree(inst), K) || throw(DimensionMismatch("moments must be [degree, source] = ($(chromatic_degree(inst)), $K)"))
    length(shape) <= M || throw(DimensionMismatch("$(length(shape)) shape coefficients for $M shape modes"))
    T = promote_type(eltype(w), eltype(z), eltype(rate), eltype(moments), map(typeof, Tuple(shape))..., typeof(flux_scale),
        jitter === nothing ? Float64 : promote_type(map(typeof, jitter)...), Float64)
    return _image_transit_loglike(T, d, t, w, z, rate, moments, shape, flux_scale, jitter, buf, _static_sizes(d)...)
end

function _image_transit_loglike(::Type{T}, d::EpochImages, t, w, z, rate, moments, shape, flux_scale, jitter, buf,
    ::Val{M}, ::Val{NM}, ::Val{NT}) where {T,M,NM,NT}
    inst = d.instrument
    K = length(rate)
    wins = d.transit_windows[t]
    it = d.transit_modes
    Kw, Kt, P = length(d.window_modes), length(it), 5 + M
    ron2, me2 = inst.read_noise^2, inst.model_error^2
    ll_good, ll_bad = @no_escape buf begin
        G, Gb = @alloc(T, P, P), @alloc(T, P, P)
        A, V = @alloc(T, Kw, Kw), @alloc(T, Kw, Kt + 1)
        H, Hb = @alloc(T, Kt + 1, Kt + 1), @alloc(T, Kt + 1, Kt + 1)
        u, electrons = @alloc(T, P), @alloc(T, K)
        cvecs = @alloc(NTuple{NT,eltype(moments)}, K)
        for k in 1:K
            cvecs[k] = _cvec(view(moments, :, k), Val(NT))
        end
        sigma_t = @alloc(T, Kt)
        for (c, a) in enumerate(it)
            sigma_t[c] = (jitter !== nothing && (a == 3 || a == 4)) ? jitter[a-2] : _transit_sigma(inst, a)
        end
        fill!(H, 0)
        fill!(Hb, 0)
        acc_good = acc_bad = zero(T)
        nvalid = 0
        ok = true
        for (q, iwin) in enumerate(wins)
            win = d.windows[iwin]
            basis = d.bases[win.basis]
            for k in 1:K
                electrons[k] = flux_scale * rate[k] * win.exposure
            end
            ws, zs = view(w, :, q), view(z, :, q)
            fill!(G, 0)
            fill!(Gb, 0)
            uP = view(u, 1:P-1)
            @inbounds for j in 0:win.n_ac-1, i in 0:win.n_al-1
                y = d.samples[win.offset+1+i+j*win.n_al]
                isnan(y) && continue
                nvalid += 1
                _sample_vector!(uP, basis, win, i, j, ws, zs, electrons, cvecs, Val(M), Val(NM))
                μ, m = _sample_mean(u, win.background, shape)
                var = max(μ, zero(μ)) + ron2 + me2 * m^2
                varb = var + inst.bad_sigma^2 + (inst.bad_fraction * m)^2
                u[P] = y - μ
                acc_good += log(var)
                acc_bad += log(varb)
                iv, ivb = 1 / var, 1 / varb
                for b in 1:P
                    ub, ubb = u[b] * iv, u[b] * ivb
                    for a in b:P
                        G[a, b] += u[a] * ub
                        Gb[a, b] += u[a] * ubb
                    end
                end
            end
            @inbounds for b in 1:P, a in b+1:P
                G[b, a] = G[a, b]
                Gb[b, a] = Gb[a, b]
            end
            lg = _eliminate_window!(H, G, A, V, d, win, it)
            lb = _eliminate_window!(Hb, Gb, A, V, d, win, it)
            ok &= lg != -Inf && lb != -Inf
            acc_good += lg
            acc_bad += lb
        end
        if ok
            acc_good += _transit_solve!(H, sigma_t)
            acc_bad += _transit_solve!(Hb, sigma_t)
        end
        c = nvalid * log(2π)
        ok ? (-(acc_good + c) / 2, -(acc_bad + c) / 2) : (T(-Inf), T(-Inf))
    end
    p = d.p_bad[t]
    ll = p == 0 ? ll_good : _logaddexp(log1p(-p) + ll_good, log(p) + ll_bad)
    return ll, ll_good
end

# ---------------------------------------------------
# The mean and draws from the model
# ---------------------------------------------------

# Columns of one window, sample by sample: U[sample, 1:4+M] the mode vectors,
# `mean` and `variance` of the good state. Samples in the order of the data.
function _window_design(d::EpochImages, iwin::Integer, ws, zs, rate, moments, shape, flux_scale)
    inst = d.instrument
    win = d.windows[iwin]
    basis = d.bases[win.basis]
    n, P1 = win.n_al * win.n_ac, 4 + inst.nshapes
    U = zeros(n, P1)
    mean, variance = zeros(n), zeros(n)
    electrons = [flux_scale * r * win.exposure for r in rate]
    vM, vNM, vNT = _static_sizes(d)
    cvecs = [_cvec(view(moments, :, k), vNT) for k in axes(moments, 2)]
    for j in 0:win.n_ac-1, i in 0:win.n_al-1
        s = 1 + i + j * win.n_al
        u = view(U, s, :)
        _sample_vector!(u, basis, win, i, j, ws, zs, electrons, cvecs, vM, vNM)
        mean[s], m = _sample_mean(u, win.background, shape)
        variance[s] = max(mean[s], 0.0) + inst.read_noise^2 + (inst.model_error * m)^2
    end
    return U, mean, variance
end

"""
    image_transit_model(data, t, w, z, rate, moments; shape = (), flux_scale = 1) -> Vector{Matrix}

The expected electrons in each window of transit `t`, at the prior mean of
the marginalised modes; `NaN` where the window has no data. Arguments as
for [`image_transit_loglike`](@ref).
"""
function image_transit_model(d::EpochImages, t::Integer, w::AbstractMatrix, z::AbstractMatrix, rate::AbstractVector,
    moments::AbstractMatrix; shape=(), flux_scale=1.0)
    return _image_transit_draw(nothing, d, t, w, z, rate, moments, shape, flux_scale, nothing)
end

"""
    image_transit_simulate(rng, data, t, w, z, rate, moments; shape = (), flux_scale = 1,
                           jitter = nothing, add_noise = true) -> Vector{Matrix}

Draw the windows of transit `t` from the good state of the model: the
transit's and the windows' mode coefficients from their priors, noise with
the model's variance. With `add_noise = false`, the same as
[`image_transit_model`](@ref).
"""
function image_transit_simulate(rng, d::EpochImages, t::Integer, w::AbstractMatrix, z::AbstractMatrix, rate::AbstractVector,
    moments::AbstractMatrix; shape=(), flux_scale=1.0, jitter=nothing, add_noise::Bool=true)
    return _image_transit_draw(add_noise ? rng : nothing, d, t, w, z, rate, moments, shape, flux_scale, jitter)
end

function _image_transit_draw(rng, d::EpochImages, t, w, z, rate, moments, shape, flux_scale, jitter)
    inst = d.instrument
    c = zeros(4 + inst.nshapes)
    if rng !== nothing
        for a in d.transit_modes
            σ = (jitter !== nothing && (a == 3 || a == 4)) ? jitter[a-2] : _transit_sigma(inst, a)
            c[a] = σ * randn(rng)
        end
    end
    return map(enumerate(d.transit_windows[t])) do (q, iwin)
        win = d.windows[iwin]
        U, mean, variance = _window_design(d, iwin, view(w, :, q), view(z, :, q), rate, moments, shape, flux_scale)
        y = copy(mean)
        if rng !== nothing
            a = copy(c)
            for m in d.window_modes
                a[m] += _window_sigma(inst, win, m) * randn(rng)
            end
            y .+= U * a .+ sqrt.(variance) .* randn(rng, length(y))
        end
        data = view(d.samples, win.offset+1:win.offset+length(y))
        reshape(ifelse.(isnan.(data), NaN, y), win.n_al, win.n_ac)
    end
end
