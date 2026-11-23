# ---------------------------------------------------
# Joining evolution models in mass
#
# MIST stops at 0.1 M☉, which is where the companions of interest start.
# BHAC15 continues to the bottom of the main sequence and below, but only at
# solar metallicity, while at 0.1 M☉ a MIST star changes by 300 K between
# [Fe/H] = -0.5 and 0. Switching tables at a mass would put a step in
# temperature and radius there at every metallicity but one.
#
# So below the join the low-mass model supplies the *dependence on mass* and
# the stellar model the level: the difference between the two at the top of
# the overlap, at the star's age and metallicity, is carried down unchanged.
# Inside the overlap the two are blended with a weight that rises smoothly.
# ---------------------------------------------------

"""
    SeamedEvolution(low, high; mass = (0.10, 0.15))

Two evolution models used as one: `high` (a [`MISTGrid`](@ref)) above
`mass[2]`, and below it `low` (a [`BHAC15Grid`](@ref), solar metallicity
only) shifted in log Teff and log radius by `high − low` at `mass[2]`, for
the same age and [Fe/H]. Between `mass[1]` and `mass[2]` the shifted `low`
and `high` are blended with the weight `t² (3 - 2t)`, `t` the fraction of the
way across in log mass, so the result is continuous with a continuous
first derivative in mass.

Called like the stellar model, `evol(mass, logage, feh)`, with the same
fields. Below `mass[1]` the present mass is the initial mass, `mh` is that
of the star at the top of the overlap, and `eep` is `NaN`.

Metallicity below the join is therefore not modelled: it is the offset of a
`mass[2]` star, applied to every lower mass.
"""
struct SeamedEvolution{L,H}
    low::L
    high::H
    mass_lo::Float64
    mass_hi::Float64
    function SeamedEvolution(low::L, high::H; mass=(0.10, 0.15)) where {L,H}
        0 < mass[1] < mass[2] || throw(ArgumentError("the seam needs 0 < mass[1] < mass[2]"))
        new{L,H}(low, high, mass[1], mass[2])
    end
end

Base.show(io::IO, e::SeamedEvolution) = print(io, "SeamedEvolution(", e.low, " below ", e.mass_lo, " M☉, ", e.high, " above ", e.mass_hi, " M☉)")

function (e::SeamedEvolution)(mass, logage, feh)
    mass >= e.mass_hi && return e.high(mass, logage, feh)
    top = e.high(e.mass_hi, logage, feh)
    T = float(promote_type(typeof(mass), typeof(top.teff)))
    low_top = e.low(e.mass_hi, logage)
    low = e.low(mass, logage)
    # the low-mass model, moved to pass through the stellar model at the top of the overlap
    log_teff = log10(low.teff) + log10(top.teff) - log10(low_top.teff)
    log_radius = log10(low.radius) + log10(top.radius) - log10(low_top.radius)
    logg = low.logg + top.logg - low_top.logg
    eep = T(NaN)
    if mass > e.mass_lo
        star = e.high(mass, logage, feh)
        t = (log(mass) - log(e.mass_lo)) / (log(e.mass_hi) - log(e.mass_lo))
        w = t^2 * (3 - 2t)
        log_teff = (1 - w) * log_teff + w * log10(star.teff)
        log_radius = (1 - w) * log_radius + w * log10(star.radius)
        logg = (1 - w) * logg + w * star.logg
        eep = T(star.eep)
    end
    teff, radius = exp10(log_teff), exp10(log_radius)
    return (; teff=T(teff), logg=T(logg), radius=T(radius), mass_now=T(exp10(logg - LOGG_SUN) * radius^2),
        mh=T(top.mh), logL=T(2log_radius + 4log10(teff / TEFF_SUN)), eep)
end
