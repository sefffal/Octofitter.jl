# ---------------------------------------------------
# BHAC15 evolution of very-low-mass stars and brown dwarfs
#
# Baraffe, Homeier, Allard & Chabrier (2015): isochrones for 0.01–1.4 M☉ at
# solar metallicity. Used below the 0.1 M☉ floor of MIST (see
# `SeamedEvolution`). Only the effective temperature and the radius are
# taken; the surface gravity follows from the mass and the radius, which is
# more precise than the two decimals the table prints.
# ---------------------------------------------------

"""
    BHAC15Grid()
    BHAC15Grid(path)

The BHAC15 isochrones on their `mass × age` nodes. Call it to evolve an
object of solar metallicity:

    low = BHAC15Grid()
    low(mass, logage)   # → (; teff, logg, radius)

with `mass` in M☉ and `logage` log₁₀ of the age in years; `teff` in K, `logg`
in cgs, `radius` in R☉. Interpolation is cubic (Catmull–Rom) in log mass and
log age. Every field is `NaN` outside the isochrones: the tables stop where
the atmosphere models do, at lower masses for younger ages (0.07 M☉ at
5 Gyr, 0.05 M☉ at 1 Gyr).

With no argument the file is Octofitter's `BHAC15_GAIA` DataDep.
"""
struct BHAC15Grid
    log_mass::Vector{Float64}
    log_age::Vector{Float64}
    log_teff::Matrix{Float64}        # [mass, age]
    log_radius::Matrix{Float64}
    present::Matrix{Bool}
end

BHAC15Grid() = BHAC15Grid(joinpath(datadep"BHAC15_GAIA", "BHAC15_iso.GAIA"))

Base.show(io::IO, g::BHAC15Grid) = print(io, "BHAC15Grid(", count(g.present), " models, ",
    round(exp10(g.log_mass[1]); digits=3), "–", round(exp10(g.log_mass[end]); digits=2), " M☉, log age ",
    round(g.log_age[1]; digits=2), "–", round(g.log_age[end]; digits=2), ")")

function BHAC15Grid(path::AbstractString)
    rows = NamedTuple{(:age, :mass, :teff, :radius),NTuple{4,Float64}}[]
    age = NaN
    for line in eachline(path)
        m = match(r"t \(Gyr\)\s*=\s*([\d.]+)", line)
        if m !== nothing
            age = parse(Float64, m.captures[1])
            continue
        end
        tok = split(line)
        (length(tok) >= 5 && !startswith(tok[1], '!')) || continue
        x = tryparse.(Float64, tok[1:5])
        any(isnothing, x) && continue
        push!(rows, (; age, mass=x[1], teff=x[2], radius=x[5]))
    end
    isempty(rows) && error("no BHAC15 isochrone found in $path")
    masses = sort(unique(r.mass for r in rows))
    ages = sort(unique(r.age for r in rows))
    log_teff = fill(NaN, length(masses), length(ages))
    log_radius = fill(NaN, length(masses), length(ages))
    for r in rows
        i, k = searchsortedfirst(masses, r.mass), searchsortedfirst(ages, r.age)
        log_teff[i, k] = log10(r.teff)
        log_radius[i, k] = log10(r.radius)
    end
    present = isfinite.(log_teff)
    log_mass = log10.(masses)
    # Each isochrone covers a contiguous range of masses. The nodes beyond
    # its ends are filled by linear extrapolation in log mass: never
    # interpolated to (a lookup needs four real corners), they give the
    # cubic stencil its end slopes.
    for k in eachindex(ages), arr in (log_teff, log_radius)
        have = findall(view(present, :, k))
        length(have) >= 2 || continue
        for (edge, next, range) in ((have[1], have[2], 1:have[1]-1), (have[end], have[end-1], have[end]+1:length(masses)))
            slope = (arr[edge, k] - arr[next, k]) / (log_mass[edge] - log_mass[next])
            for i in range
                arr[i, k] = arr[edge, k] + slope * (log_mass[i] - log_mass[edge])
            end
        end
    end
    return BHAC15Grid(log_mass, log10.(ages) .+ 9, log_teff, log_radius, present)
end

function (g::BHAC15Grid)(mass, logage)
    T = float(promote_type(typeof(mass), typeof(logage)))
    dead = (; teff=T(NaN), logg=T(NaN), radius=T(NaN))
    mass > 0 || return dead
    lm = log10(mass)
    sm = _stencil(g.log_mass, lm)
    sa = _stencil(g.log_age, logage)
    (sm === nothing || sa === nothing) && return dead
    (im, wm, _), (ia, wa, _) = sm, sa
    @inbounds for b in (2, 3), a in (2, 3)
        g.present[im[a], ia[b]] || return dead
    end
    log_teff = log_radius = zero(T)
    @inbounds for b in 1:4, a in 1:4
        w = wm[a] * wa[b]
        iszero(w) && continue
        log_teff += w * g.log_teff[im[a], ia[b]]
        log_radius += w * g.log_radius[im[a], ia[b]]
    end
    return (; teff=T(exp10(log_teff)), logg=T(LOGG_SUN + lm - 2log_radius), radius=T(exp10(log_radius)))
end
