# ---------------------------------------------------
# From the model's absolute frame to a transit's local plane
#
# Every DR4 epoch product places things in local-plane coordinates counted
# from one tangent point per source, (ra0, dec0). The model has an absolute
# frame, the barycentre's DR4 position, proper motion, parallax and radial
# velocity, which PlanetOrbits propagates in three dimensions, and bodies
# that move about the barycentre. This file is the only place where one
# becomes the other, and it does it the way `GaiaDR4AstromObs` does: the
# body's excursion from the barycentre is composed with the barycentre's
# propagated direction on the sphere and projected about the tangent point,
# and the parallax that multiplies the transit's parallax factors is the
# body's own at that epoch.
#
# Conventions that the DR4 documentation has not yet fixed live here too.
# ---------------------------------------------------

"Unit vector and local east and north axes of a tangent point [deg], as `transit_position` uses them."
tangent_basis(ra0, dec0) = Octofitter._dr4_tangent_basis(ra0, dec0)

"""
    scan_coordinates(ξ, η, sinψ, cosψ) -> (w, z)

Along- and across-scan components [mas] of a local-plane offset (`ξ` toward
increasing right ascension, `η` toward increasing declination) at scan
position angle ψ. Along scan is the convention of `centroid_pos_al`. The
across-scan sign is GaiaTwin's and is **not yet confirmed** against the DR4
documentation.
"""
@inline scan_coordinates(ξ, η, sinψ, cosψ) = (ξ * sinψ + η * cosψ, ξ * cosψ - η * sinψ)

"""
    transit_position(sol, body, barycentre, basis, sinψ, cosψ, f_al, f_ac) -> (w, z)

Local-plane position [mas] of `body` in one transit, along and across scan,
counted from the tangent point whose [`tangent_basis`](@ref) is `basis`.
`sol` is the solved trajectory at the transit's epoch, of a system with an
absolute frame; `f_al` and `f_ac` are the transit's parallax factors.
"""
@inline function transit_position(sol, body, barycentre, basis, sinψ, cosψ, f_al, f_ac)
    ξ, η = Octofitter._dr4_project(PlanetOrbits.frame_ra(sol), PlanetOrbits.frame_dec(sol),
        PlanetOrbits.raoff(sol, body, barycentre), PlanetOrbits.decoff(sol, body, barycentre), basis)
    plx = Octofitter._dr4_source_plx(sol, body, barycentre)
    w, z = scan_coordinates(ξ, η, sinψ, cosψ)
    return w + plx * f_al, z + plx * f_ac
end

"""
    transit_velocity(sol, body, barycentre, einstein) -> km/s

Barycentric line-of-sight velocity of `body` at one transit: the frame's
propagated radial velocity plus the body's motion about the barycentre,
with the relativistic terms of the orbit when `einstein` is true (the
system's `einstein_rv` setting).
"""
@inline function transit_velocity(sol, body, barycentre, einstein::Bool)
    orbital = einstein ? PlanetOrbits.radvel(sol, body, barycentre) :
              PlanetOrbits.velz(sol, body, barycentre) * Octofitter._AU_PER_YEAR_TO_M_PER_S
    return (PlanetOrbits.frame_rv(sol) + orbital) / 1000
end

"""
    catalogue_track(epoch, scan_pos_angle, parallax_factor_al, parallax_factor_ac;
                    parallax, pmra, pmdec, ref_epoch, Δra = 0.0, Δdec = 0.0) -> (; al, ac)

Where a five-parameter astrometric solution puts a source in each transit,
in local-plane coordinates [mas]: the position the RVS pipeline assumed
when it set the wavelength scale of a source without an orbital solution.
`Δra`, `Δdec` [mas] are the catalogue position's offsets from the tangent
point (zero when the tangent point is the catalogue position), `pmra`,
`pmdec` in mas/yr, `parallax` in mas, epochs in MJD and the scan angle in
degrees.
"""
function catalogue_track(epoch, scan_pos_angle, parallax_factor_al, parallax_factor_ac;
    parallax, pmra, pmdec, ref_epoch, Δra=0.0, Δdec=0.0)
    al, ac = similar(epoch, Float64), similar(epoch, Float64)
    for i in eachindex(epoch)
        Δt = (epoch[i] - ref_epoch) / PlanetOrbits.year2day_julian
        w, z = scan_coordinates(Δra + pmra * Δt, Δdec + pmdec * Δt, sind(scan_pos_angle[i]), cosd(scan_pos_angle[i]))
        al[i] = w + parallax * parallax_factor_al[i]
        ac[i] = z + parallax * parallax_factor_ac[i]
    end
    return (; al, ac)
end
