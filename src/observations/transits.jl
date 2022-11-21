using Transits

const light_curve_cols = (:epoch, :phot, :σ_phot)

struct LightCurve3{TLimbDark, TTable<:Table} <: AbstractObs
    limbdark::TLimbDark
    table::TTable
    function LightCurve3(ld, observations...)
        table = Table(observations...)
        if !issubset(light_curve_cols, Tables.columnnames(table))
            error("Ecpected columns $light_curve_cols")
        end
        return new{typeof(ld), typeof(table)}(ld, table)
    end
end
LightCurve3(observations::NamedTuple...) = LightCurve3(observations)
export LightCurve3


"""
Transit likelihood. Uses Transits.jl QuadLimbDark.
"""
function ln_like(lc::LightCurve3, θ_system, elements)
    T = Float64
    ll = zero(T)

    params = tuple()
    if hasproperty(θ_system, :u1)
        params = tuple(θ_system.u1)
        if hasproperty(θ_system, :u2)
            params = tuple(θ_system.u1, θ_system.u2)
            if hasproperty(θ_system, :u3)
                params = tuple(θ_system.u1, θ_system.u2, θ_system.u3)
                if hasproperty(θ_system, :u4)
                    params = tuple(θ_system.u1, θ_system.u2, θ_system.u3, θ_system.u4)
                end
            end
        end
    end

    ld = lc.limbdark(SVector(params))

    Rₛₜₐᵣ = θ_system.R

    for i in eachindex(lc.table.epoch)

        # TODO: Handle multiple transiting bodies at once
        phot_star = zero(T)
        # for j in eachindex(elements)
        j = firstindex(elements)
        θ_planet = θ_system.planets[j]
        orbit = elements[j]
        t = lc.table.epoch[i]
        r = θ_planet.r

        phot_star = transit_depth(orbit, t, r, Rₛₜₐᵣ, ld)
        resid = phot_star - lc.table.phot[i]
        σ² = lc.table.σ_phot[i]^2
        χ² = -0.5resid^2 / σ² - log(sqrt(2π * σ²))
        ll += χ²
    end

    return ll
end

# R_star in meters
function transit_depth(orbit, t,  r, Rₛₜₐᵣ, ld=QuadLimbDark(Float64[]))
    soln = orbitsolve(orbit, t)

    x = PlanetOrbits.posx(soln)*DirectDetections.PlanetOrbits.au2m
    y = PlanetOrbits.posy(soln)*DirectDetections.PlanetOrbits.au2m
    z = PlanetOrbits.posz(soln)*DirectDetections.PlanetOrbits.au2m

    
    # TODO: at the moment this only supports VisualOrbit
    cosi = orbit.cosi


    T = promote_type(typeof(t), typeof(r), typeof(cosi))

    # From Transits.jl.
    # Ensure we are in front of the star.
    # if y < 0
        # one(T)
    # else
        # bₜᵣₐₙ = sqrt(x^2 + z^2)/Rₛₜₐᵣ
        bₜᵣₐₙ = sqrt(z^2 + y^2)/Rₛₜₐᵣ
        convert(T, Transits.compute(ld, bₜᵣₐₙ, r))
    # end

end
