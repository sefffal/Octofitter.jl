# ---------------------------------------------------
# MIST stellar evolution: (initial mass, age, [Fe/H]) → (Teff, log g, R)
#
# The tables are MIST v1.2 isochrones (Choi et al. 2016, Dotter 2016),
# v/vcrit = 0.4, solar-scaled. An isochrone is a curve through
# (mass, Teff, log g, L) indexed by the equivalent evolutionary point (EEP).
# Each isochrone already carries MIST's own interpolation in mass; what is
# added here is the interpolation between ages and metallicities, and the
# question is which point of one isochrone corresponds to which point of its
# neighbour. There are two answers, each right in one regime:
#
#   - Evolved stars: the same EEP. Stars at the same EEP are at the same
#     stage of their evolution, whereas at fixed mass a star can be on the
#     main sequence at one age node and a giant at the next (the whole
#     post-main-sequence of an old isochrone spans about 0.02 dex of mass).
#   - Unevolved stars: the same mass. Below about 0.35 M☉ a star never
#     reaches the intermediate-age main sequence, its EEPs are spread
#     between the zero-age main sequence and wherever its track happened to
#     stop, and the numbering is not comparable between metallicities: at
#     10 Gyr a 0.2 M☉ star is at EEP 255 for [Fe/H] = 0, 220 for +0.25 and
#     244 for +0.5, with all of 0.1–0.27 M☉ inside one EEP at +0.25.
#
# So the lookup evaluates both and blends them on the EEP itself: fixed mass
# up to the intermediate-age main sequence (EEP 353), fixed EEP from the
# terminal-age main sequence (EEP 454) on, a smooth step in between. On a
# node of the age and metallicity axes both reduce to that node's isochrone,
# so the result is continuous everywhere.
#
# Only theoretical quantities are used by the model. MIST's own Gaia
# magnitudes are kept in the grid file for one purpose, the check that
# magnitudes integrated from the spectral library reproduce them.
# ---------------------------------------------------

# IAU 2015 nominal solar values (resolution B3), which are also MESA's.
const TEFF_SUN = 5772.0                 # K
const LOGG_SUN = log10(1.3271244e26 / 6.957e10^2)   # cgs; GM☉ / R☉²

# EEPs kept: everything before the post-AGB (primary EEP 1409). A star past
# the AGB tip is a remnant, which is a different emitter. The cut is on the
# EEP and not on MIST's phase flag, because the flag is not monotonic along
# an isochrone (Wolf–Rayet points sit among the core-helium-burning ones).
const MIST_LAST_EEP = 1408
# The tables number EEPs from 0; array index = EEP + 1.
const MIST_EEP_IAMS = 353.0     # intermediate-age main sequence
const MIST_EEP_TAMS = 454.0     # terminal-age main sequence

"""
    MISTGrid

MIST v1.2 isochrones on a regular `EEP × log age × [Fe/H]` grid. Call it to
evolve a star:

    evol = MISTGrid()
    evol(mass, logage, feh)   # → (; teff, logg, radius, mass_now, mh, logL, eep)

with `mass` the initial mass [M☉], `logage` log₁₀ of the age in years and
`feh` the initial [Fe/H]. `teff` is in K, `logg` in cgs, `radius` in R☉,
`mass_now` the present mass [M☉], and `mh` the present surface [Fe/H], which
is what a model atmosphere should be looked up at. Outside the grid, or for
a star that has evolved past the AGB tip at that age, every field is `NaN`.

Interpolation between isochrones is cubic (Catmull–Rom) in `log age` and
`[Fe/H]`, at fixed mass on the main sequence and at fixed EEP beyond it;
along an isochrone it is linear in log mass between adjacent EEPs. The
result is continuous, and works with `ForwardDiff` duals.

With no argument the grid is read from the `MIST_v1.2_vvcrit0.4` DataDep,
which downloads the isochrones (160 MB) and builds the grid on first use.
"""
struct MISTGrid
    feh::Vector{Float64}
    log_age::Vector{Float64}
    eep_lo::Matrix{Int}            # [age, feh] first EEP index present
    eep_hi::Matrix{Int}            # [age, feh] last EEP index present
    log_mass::Array{Float64,3}     # [EEP + 1, age, feh] log10 initial mass
    log_teff::Array{Float32,3}
    logg::Array{Float32,3}
    log_L::Array{Float32,3}
    feh_surf::Array{Float32,3}
    path::String
end

function Base.show(io::IO, g::MISTGrid)
    print(io, "MISTGrid(", size(g.log_mass, 1), " EEPs × ", length(g.log_age), " ages × ",
        length(g.feh), " [Fe/H]; log age ", extrema(g.log_age), ", [Fe/H] ", extrema(g.feh), ")")
end

MISTGrid() = MISTGrid(joinpath(datadep"MIST_v1.2_vvcrit0.4", MIST_GRID_FILE))

function MISTGrid(path::AbstractString)
    h5open(path, "r") do f
        MISTGrid(read(f["feh"]), read(f["log_age"]),
            Int.(read(f["eep_lo"])), Int.(read(f["eep_hi"])),
            read(f["log_mass"]), read(f["log_teff"]), read(f["logg"]), read(f["log_L"]),
            read(f["feh_surf"]), String(path))
    end
end

"""
    mist_magnitudes(grid) -> (; MG, MBP, MRP, phase)

MIST's own Gaia EDR3 absolute magnitudes and phase flags on the grid's
`[EEP + 1, age, feh]` nodes (`NaN` and `typemin(Int8)` where an isochrone has
no such EEP). The model does not use them; they exist so that magnitudes
integrated from the spectral library can be checked against MIST's.
"""
function mist_magnitudes(g::MISTGrid)
    h5open(g.path, "r") do f
        (; MG=read(f["MG"]), MBP=read(f["MBP"]), MRP=read(f["MRP"]), phase=read(f["phase"]))
    end
end

# Catmull–Rom stencil on a sorted, possibly non-uniform axis: four node
# indices and their weights, interpolating and C¹. Tangents are the
# three-point parabola derivatives, one-sided at the ends of the axis, where
# the unused index is clamped and carries zero weight. `nothing` outside.
@inline function _stencil(nodes::Vector{Float64}, x)
    n = length(nodes)
    (x >= nodes[1] && x <= nodes[n]) || return nothing
    i = min(searchsortedlast(nodes, x), n - 1)
    @inbounds x1, x2 = nodes[i], nodes[i+1]
    d1 = x2 - x1
    t = (x - x1) / d1
    h00 = (2t - 3) * t^2 + 1
    h10 = ((t - 2) * t + 1) * t
    h01 = (3 - 2t) * t^2
    h11 = (t - 1) * t^2
    w0, w1, w2, w3 = zero(t), h00, h01, zero(t)
    if i > 1
        d0 = x1 - @inbounds nodes[i-1]
        w0 -= h10 * d1 * d1 / (d0 * (d0 + d1))
        w1 += h10 * (d1 - d0) / d0
        w2 += h10 * d0 / (d0 + d1)
    else
        w1 -= h10
        w2 += h10
    end
    if i < n - 1
        d2 = @inbounds(nodes[i+2]) - x2
        w1 -= h11 * d2 / (d1 + d2)
        w2 += h11 * (d2 - d1) / d2
        w3 += h11 * d1 * d1 / (d2 * (d1 + d2))
    else
        w1 -= h11
        w2 += h11
    end
    return (max(i - 1, 1), i, i + 1, min(i + 2, n)), (w0, w1, w2, w3)
end

# Position of log mass `lm` along one stored isochrone: the index of the
# lower EEP and the fraction towards the next. An index of 0 means the mass
# is below the isochrone's first point, -1 that it is above its last (the
# star is dead at that age).
@inline function _along(g::MISTGrid, k, j, lm)
    @inbounds lo, hi = g.eep_lo[k, j], g.eep_hi[k, j]
    @inbounds lm < g.log_mass[lo, k, j] && return 0, zero(lm)
    @inbounds lm > g.log_mass[hi, k, j] && return -1, zero(lm)
    a, b = lo, hi
    @inbounds while b - a > 1
        m = (a + b) >>> 1
        if g.log_mass[m, k, j] <= lm
            a = m
        else
            b = m
        end
    end
    @inbounds la, lb = g.log_mass[a, k, j], g.log_mass[b, k, j]
    return a, (lb > la ? (lm - la) / (lb - la) : zero(lm))
end

function (g::MISTGrid)(mass, logage, feh)
    T = float(promote_type(typeof(mass), typeof(logage), typeof(feh)))
    dead = (; teff=T(NaN), logg=T(NaN), radius=T(NaN), mass_now=T(NaN), mh=T(NaN),
        logL=T(NaN), eep=T(NaN))
    sa = _stencil(g.log_age, logage)
    sf = _stencil(g.feh, feh)
    (sa === nothing || sf === nothing || !(mass > 0)) && return dead
    ks, was = sa
    js, wfs = sf
    lm = log10(mass)

    # --- fixed mass: each neighbouring isochrone evaluated at this mass
    alive = true
    estar = T(-Inf)                       # the most evolved of the neighbours
    m_teff = m_logg = m_logL = m_mh = m_eep = zero(T)
    for ia in 1:4, jf in 1:4
        (iszero(was[ia]) || iszero(wfs[jf])) && continue
        k, j = ks[ia], js[jf]
        e, s = _along(g, k, j, lm)
        e == 0 && return dead
        if e == -1
            alive = false
            continue
        end
        w = was[ia] * wfs[jf]
        @inbounds begin
            m_teff += w * ((1 - s) * g.log_teff[e, k, j] + s * g.log_teff[e+1, k, j])
            m_logg += w * ((1 - s) * g.logg[e, k, j] + s * g.logg[e+1, k, j])
            m_logL += w * ((1 - s) * g.log_L[e, k, j] + s * g.log_L[e+1, k, j])
            m_mh += w * ((1 - s) * g.feh_surf[e, k, j] + s * g.feh_surf[e+1, k, j])
        end
        m_eep += w * (e - 1 + s)
        estar = max(estar, e - 1 + s)
    end
    wE = if alive
        t = clamp((estar - MIST_EEP_IAMS) / (MIST_EEP_TAMS - MIST_EEP_IAMS), 0, 1)
        t^2 * (3 - 2t)
    else
        one(T)
    end

    # --- fixed EEP: the neighbouring isochrones blended EEP by EEP, and the
    # EEP then found from the mass
    e_teff = e_logg = e_logL = e_mh = e_eep = zero(T)
    if !iszero(wE)
        lo, hi = 1, size(g.log_mass, 1)
        for ia in 1:4, jf in 1:4
            (iszero(was[ia]) || iszero(wfs[jf])) && continue
            @inbounds lo = max(lo, g.eep_lo[ks[ia], js[jf]])
            @inbounds hi = min(hi, g.eep_hi[ks[ia], js[jf]])
        end
        lo < hi || return dead
        blend(arr, e) = begin
            acc = zero(T)
            for ia in 1:4, jf in 1:4
                # a zero weight marks a clamped index, whose isochrone may
                # not reach this EEP
                (iszero(was[ia]) || iszero(wfs[jf])) && continue
                @inbounds acc += was[ia] * wfs[jf] * arr[e, ks[ia], js[jf]]
            end
            acc
        end
        (lm >= blend(g.log_mass, lo) && lm <= blend(g.log_mass, hi)) || return dead
        a, b = lo, hi
        while b - a > 1
            m = (a + b) >>> 1
            if blend(g.log_mass, m) <= lm
                a = m
            else
                b = m
            end
        end
        la, lb = blend(g.log_mass, a), blend(g.log_mass, b)
        s = lb > la ? (lm - la) / (lb - la) : zero(T)
        e_teff = (1 - s) * blend(g.log_teff, a) + s * blend(g.log_teff, b)
        e_logg = (1 - s) * blend(g.logg, a) + s * blend(g.logg, b)
        e_logL = (1 - s) * blend(g.log_L, a) + s * blend(g.log_L, b)
        e_mh = (1 - s) * blend(g.feh_surf, a) + s * blend(g.feh_surf, b)
        e_eep = a - 1 + s
    end

    mix(m, e) = iszero(wE) ? m : isone(wE) ? e : (1 - wE) * m + wE * e
    log_teff, logg, logL = mix(m_teff, e_teff), mix(m_logg, e_logg), mix(m_logL, e_logL)
    teff = exp10(log_teff)
    radius = exp10(logL / 2) * (TEFF_SUN / teff)^2
    mass_now = exp10(logg - LOGG_SUN) * radius^2
    return (; teff=T(teff), logg=T(logg), radius=T(radius), mass_now=T(mass_now),
        mh=T(mix(m_mh, e_mh)), logL=T(logL), eep=T(mix(m_eep, e_eep)))
end

# ---------------------------------------------------
# Building the grid file from the distributed isochrones
# ---------------------------------------------------

const MIST_GRID_FILE = "mist_v1.2_vvcrit0.4_eep_grid.h5"

# Columns of the UBV(RI)c + 2MASS + Kepler + Hipparcos + Gaia `.iso.cmd` files
const _MIST_COLS = (eep=1, log_age=2, initial_mass=3, log_teff=5, logg=6, log_L=7,
    feh_init=8, feh_surf=9, MG=31, MBP=32, MRP=33, phase=34)

function _read_iso_cmd(path)
    rows = NTuple{length(_MIST_COLS),Float64}[]
    names = nothing
    for line in eachline(path)
        if startswith(line, '#')
            startswith(line, "# EEP") && (names = split(line)[2:end])
            continue
        end
        isempty(strip(line)) && continue
        f = split(line)
        push!(rows, map(c -> parse(Float64, f[c]), values(_MIST_COLS)))
    end
    # Guard the column numbers against a different photometric-system file.
    expect = (eep="EEP", log_age="log10_isochrone_age_yr", initial_mass="initial_mass",
        log_teff="log_Teff", logg="log_g", log_L="log_L", feh_init="[Fe/H]_init",
        feh_surf="[Fe/H]", MG="Gaia_G_EDR3", MBP="Gaia_BP_EDR3", MRP="Gaia_RP_EDR3",
        phase="phase")
    for (k, c) in pairs(_MIST_COLS)
        names[c] == expect[k] || error("$path: column $c is $(names[c]), expected $(expect[k])")
    end
    return rows
end

"""
    build_mist_grid(dir, out) -> report

Read every `*.iso.cmd` file in `dir` (the unpacked
`MIST_v1.2_vvcrit0.4_UBVRIplus` distribution) and write the regular
`EEP × log age × [Fe/H]` grid file `out` that [`MISTGrid`](@ref) loads.

Returns a named tuple describing what was built, including
`mass_inversions`, the number of adjacent EEP pairs along an isochrone whose
initial mass decreases (the lookup assumes there are none), and `n_filled`,
the number of EEPs an isochrone skipped that were filled by interpolation.
"""
function build_mist_grid(dir::AbstractString, out::AbstractString)
    files = sort(filter(endswith(".iso.cmd"), readdir(dir; join=true)))
    isempty(files) && error("no .iso.cmd files in $dir")
    col = map(i -> i, NamedTuple{keys(_MIST_COLS)}(1:length(_MIST_COLS)))
    tables = map(_read_iso_cmd, files)
    fehs = [t[1][col.feh_init] for t in tables]
    order = sortperm(fehs)
    tables, fehs = tables[order], fehs[order]
    ages = sort(unique(round(r[col.log_age]; digits=4) for t in tables for r in t))
    ne = MIST_LAST_EEP + 1
    na, nf = length(ages), length(fehs)

    cube(T) = fill(T(NaN), ne, na, nf)
    log_mass = cube(Float64)
    log_teff, logg, log_L, feh_surf = cube(Float32), cube(Float32), cube(Float32), cube(Float32)
    MG, MBP, MRP = cube(Float32), cube(Float32), cube(Float32)
    phase = fill(typemin(Int8), ne, na, nf)
    eep_lo = fill(typemax(Int32), na, nf)
    eep_hi = fill(Int32(0), na, nf)

    for (j, t) in enumerate(tables)
        all(r -> r[col.feh_init] == fehs[j], t) || error("$(files[order[j]]): mixed [Fe/H]")
    end
    for (j, t) in enumerate(tables), r in t
        r[col.eep] <= MIST_LAST_EEP || continue
        e = Int(r[col.eep]) + 1
        k = searchsortedfirst(ages, round(r[col.log_age]; digits=4))
        isnan(log_mass[e, k, j]) || error("duplicate EEP $(e - 1) at log age $(ages[k]), [Fe/H] $(fehs[j])")
        log_mass[e, k, j] = log10(r[col.initial_mass])
        log_teff[e, k, j] = r[col.log_teff]
        logg[e, k, j] = r[col.logg]
        log_L[e, k, j] = r[col.log_L]
        feh_surf[e, k, j] = r[col.feh_surf]
        MG[e, k, j], MBP[e, k, j], MRP[e, k, j] = r[col.MG], r[col.MBP], r[col.MRP]
        phase[e, k, j] = Int8(r[col.phase])
        eep_lo[k, j] = min(eep_lo[k, j], e)
        eep_hi[k, j] = max(eep_hi[k, j], e)
    end

    # A few isochrones skip EEPs (at [Fe/H] = +0.25 below 0.31 M☉ and at -3.5
    # near 1 M☉). The lookup needs every EEP between the first and the last, so
    # the skipped ones are filled by linear interpolation in EEP and counted.
    filled = fill(false, ne, na, nf)
    inversions = 0
    for j in 1:nf, k in 1:na
        eep_hi[k, j] == 0 && error("no isochrone at log age $(ages[k]), [Fe/H] $(fehs[j])")
        run = eep_lo[k, j]:eep_hi[k, j]
        present = [e for e in run if !isnan(log_mass[e, k, j])]
        for (a, b) in zip(present[1:end-1], present[2:end]), e in a+1:b-1
            t = (e - a) / (b - a)
            for arr in (log_mass, log_teff, logg, log_L, feh_surf, MG, MBP, MRP)
                arr[e, k, j] = (1 - t) * arr[a, k, j] + t * arr[b, k, j]
            end
            phase[e, k, j] = phase[a, k, j]
            filled[e, k, j] = true
        end
        inversions += count(e -> log_mass[e+1, k, j] < log_mass[e, k, j], run[1:end-1])
    end

    h5open(out, "w") do f
        attributes(f)["source"] = "MIST v1.2, v/vcrit = 0.4, [a/Fe] = 0; UBV(RI)c + 2MASS + Kepler + Hipparcos + Gaia isochrones"
        attributes(f)["reference"] = "Choi et al. 2016, ApJ 823, 102; Dotter 2016, ApJS 222, 8"
        attributes(f)["last_eep"] = MIST_LAST_EEP
        f["feh"] = fehs
        f["log_age"] = ages
        f["eep_lo"] = eep_lo
        f["eep_hi"] = eep_hi
        for (name, a) in (("log_mass", log_mass), ("log_teff", log_teff), ("logg", logg),
            ("log_L", log_L), ("feh_surf", feh_surf), ("MG", MG), ("MBP", MBP), ("MRP", MRP),
            ("phase", phase), ("filled", UInt8.(filled)))
            f[name, chunk=(ne, 1, 1), deflate=4] = a
        end
    end
    return (; feh=fehs, log_age=extrema(ages), n_eep=ne, mass_inversions=inversions,
        n_points=count(!isnan, log_mass), n_filled=count(filled),
        filled_isochrones=[(fehs[j], ages[k]) for j in 1:nf for k in 1:na if any(@view filled[:, k, j])])
end

# DataDeps post-fetch: unpack the distributed archive, build the grid file,
# and drop the 1.1 GB of text it was built from.
function _mist_post_fetch(archive)
    dir = dirname(abspath(archive))
    p7zip_jll.p7zip() do exe
        run(pipeline(`$exe x $archive -y -so`, `$exe x -si -y -ttar -o$dir`))
    end
    src = joinpath(dir, "MIST_v1.2_vvcrit0.4_UBVRIplus")
    build_mist_grid(src, joinpath(dir, MIST_GRID_FILE))
    rm(src; recursive=true)
    rm(archive)
end
