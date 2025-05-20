
using Octofitter.StaticArrays

# We want this to work with both Campbell and Thiele-Innes parameterizations
export θ_at_epoch_to_tperi
function θ_at_epoch_to_tperi(system,planet,theta_epoch)
    (;M,e,θ) = merge(system,planet)
    if  hasproperty(planet, :B) &&
        hasproperty(planet, :G) &&
        hasproperty(planet, :A) &&
        hasproperty(planet, :F)
        (;B,G,A,F) = planet
        (;plx) = system

        # Math in order to get semi-major axis -> mean motion and period (already in TI constructor)
        u = (A^2 + B^2 + F^2 + G^2)/2
        v = A*G - B * F
        α = sqrt(u + sqrt((u+v)*(u-v)))
        a = α/plx
    elseif hasproperty(planet, :i) &&
           hasproperty(planet, :Ω) &&
           hasproperty(planet, :ω)

        if hasproperty(planet, :a)
            (;a) = planet
        elseif hasproperty(planet, :P)
            (; P) = planet
            a = ∛(system.M * b.P^2)*PlanetOrbits.kepler_year_to_julian_day_conversion_factor
        else
            error("Your planet must specify either i, Ω, ω, and a/P, or B, G, A, F.")
        end

        (;i, Ω, ω) = planet
        A = ( cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i))
        B = ( sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i))
        F = (-cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i))
        G = (-sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i))
    else
        error("Your planet must specify either i, Ω, ω, and a/P, or B, G, A, F.")
    end

    # Thiele-Innes transformation matrix
    T = @SArray [
        A F
        B G
    ]
    x_over_r, y_over_r = T \ @SArray[
        cos(θ)
        sin(θ)
    ]
    # Note: this is missing the radius factor but we don't need it for true anomaly.
    # r = sqrt((sol.x*sol.elem.B+sol.y*sol.elem.G)^2 + (sol.x*sol.elem.A + sol.y*sol.elem.F)^2 )

    # True anomaly is now straightforward to calculate in the deprojected coordinate space
    ν = atan(y_over_r,x_over_r)

    # Mean anomaly (see Wikipedia page)
    MA = atan(-sqrt(1-e^2)*sin(ν), -e-cos(ν))+π-e*sqrt(1-e^2)*sin(ν)/(1+e*cos(ν))

    period_days = √(a^3/M)*PlanetOrbits.kepler_year_to_julian_day_conversion_factor
    period_yrs = period_days/PlanetOrbits.year2day_julian
    n = 2π/period_yrs # mean motion

    # Get epoch of periastron passage
    tₚ = theta_epoch - MA/n*PlanetOrbits.year2day_julian
    return tₚ
end

function θ_ωθ_at_epoch_to_ω_tperi(system,planet,theta_epoch)
    (;M,e,θ) = merge(system,planet)
   
    if hasproperty(planet, :a)
        (;a) = planet
    elseif hasproperty(planet, :P)
        (; P) = planet
        a = ∛(system.M * P^2)*PlanetOrbits.kepler_year_to_julian_day_conversion_factor
    else
        error("Your planet must specify either i, Ω, ω, and a/P")
    end

    # # Old version
    # (;i, Ω, ω) = planet
    # A = ( cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i))
    # B = ( sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i))
    # F = (-cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i))
    # G = (-sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i))

    # New version
    (;i, Ω, ωθ) = planet
    # ω = atan(tan(ωθ) / cos(i)) - Ω
    ω = atan(-tan(Ω - ωθ) / cos(i))

    # note that 
    # ωθ = atan2(B, A)
    
    # # start with calculation for position angle at mean anomaly = 0 and work backwards
    # MA = atan(-sqrt(1-e^2)*sin(ν), -e-cos(ν))+π-e*sqrt(1-e^2)*sin(ν)/(1+e*cos(ν))
    # 0 ==  atan(-sqrt(1-e^2)*sin(ν), -e-cos(ν))+π-e*sqrt(1-e^2)*sin(ν)/(1+e*cos(ν))
    # e*sqrt(1-e^2)*sin(ν)/(1+e*cos(ν)) ==  atan(-sqrt(1-e^2)*sin(ν), -e-cos(ν))+π
    # # so  if mean anomaly = 0, then true anomaly = 0

    # 0 = ν = atan(y_over_r,x_over_r)
    # tan(0) = y_over_r/ x_over_r
    # -> y_over_r = 0.
    # x_over_r, undefined?

    # [x_over_r, 0] = T \ @SArray[
    #     cos(ωθ)
    #     sin(ωθ)
    # ]
    # T = [x_over_r, 0] * @SArray[
    #     cos(ωθ)
    #     sin(ωθ)
    # ]
    A = ( cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i))
    B = ( sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i))
    F = (-cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i))
    G = (-sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i))


    T = @SArray [
        A F
        B G
    ]

    # Note: this is missing the radius factor but we don't need it for true anomaly.
    # r = sqrt((sol.x*sol.elem.B+sol.y*sol.elem.G)^2 + (sol.x*sol.elem.A + sol.y*sol.elem.F)^2 )

    # True anomaly is now straightforward to calculate in the deprojected coordinate space



    # Thiele-Innes transformation matrix
    T = @SArray [
        A F
        B G
    ]
    x_over_r, y_over_r = T \ @SArray[
        cos(θ)
        sin(θ)
    ]
    # Note: this is missing the radius factor but we don't need it for true anomaly.
    # r = sqrt((sol.x*sol.elem.B+sol.y*sol.elem.G)^2 + (sol.x*sol.elem.A + sol.y*sol.elem.F)^2 )

    # True anomaly is now straightforward to calculate in the deprojected coordinate space
    ν = atan(y_over_r,x_over_r)

    # Mean anomaly (see Wikipedia page)
    MA = atan(-sqrt(1-e^2)*sin(ν), -e-cos(ν))+π-e*sqrt(1-e^2)*sin(ν)/(1+e*cos(ν))

    period_days = √(a^3/M)*PlanetOrbits.kepler_year_to_julian_day_conversion_factor
    period_yrs = period_days/PlanetOrbits.year2day_julian
    n = 2π/period_yrs # mean motion

    # Get epoch of periastron passage
    tₚ = theta_epoch - MA/n*PlanetOrbits.year2day_julian
    return tₚ, ω
end

function θ_sep_at_epoch_to_tperi_sma(system,planet,theta_epoch)
    (;M,plx,e,θ,e,i, Ω, ω, θ, sep) = merge(system,planet)
    A = ( cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i))
    B = ( sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i))
    F = (-cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i))
    G = (-sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i))

    # Thiele-Innes transformation matrix
    T = @SArray [
        A F
        B G
    ]
    x_over_r, y_over_r = T \ @SArray[
        cos(θ)
        sin(θ)
    ]
    # Note: this is missing the radius factor but we don't need it for true anomaly.
    # r = sqrt((sol.x*sol.elem.B+sol.y*sol.elem.G)^2 + (sol.x*sol.elem.A + sol.y*sol.elem.F)^2 )
    #  sqrt(r^2* (x*B+y*G)^2 + r^2 * (x*A + y*F)^2)
    #  sqrt(r^2* (x*B+y*G)^2 + r^2 * (x*A + y*F)^2)
    # r * sqrt(r^2* (x*B+y*G)^2 + r^2 * (x*A + y*F)^2)
    # r_unit = sqrt(
    #     (x_over_r*B+y_over_r*G)^2 + (x_over_r*A + y_over_r*F)^2
    # )


    # True anomaly is now straightforward to calculate in the deprojected coordinate space
    ν = atan(y_over_r,x_over_r)

    # Mean anomaly (see Wikipedia page)
    if !(0 < e  < 1)
        # @show planet.sqrtecosω
        # @show planet.sqrtesinω
        # @show planet.ω
        # @show e
        return NaN, NaN
    end
    MA = atan(-sqrt(1-e^2)*sin(ν), -e-cos(ν))+π-e*sqrt(1-e^2)*sin(ν)/(1+e*cos(ν))

    # r = a*(1 - e^2)/(1 + ecosν)
    # # r_unit = (1 - e^2)/(1 + e*cos(ν))
    # r_unit = (1 - e^2)/(1 + e*cos(ν))

    # # unitless projection vectors for case of a=1
    # xcart_unit = r_unit*(cos(ν+ω)*sin(Ω) + sin(ν+ω)*cos(i)*cos(Ω))
    # ycart_unit = r_unit*(cos(ν+ω)*cos(Ω) - sin(ν+ω)*cos(i)*sin(Ω))
    # zcart_unit = r_unit*(sin(ν+ω)*sin(i))

    
    # # # scale a until projected sep matches variable
    dist = 1000/plx * PlanetOrbits.pc2au
    cart2angle = PlanetOrbits.rad2as*1e3/dist
    # sep_unit = sqrt(xcart_unit^2 + ycart_unit^2 + zcart_unit^2)*cart2angle
    # a = sep/sep_unit
    # # @show a

    # So we had the direction, now we have the mean anomaly.
    # We know the separation of a 1AU orbit at that mean anomaly 
    # Then we scan the orbit SMA in AU to match
    
    
    # Compute eccentric anomaly 
    # Ideally we will adjust our actual orbit parameterization to prevent
    # calculating this twice
    # EA = kepler_solver(MA, eccentricity(elem))
    
    # Calculate true anomaly
    # ν = 2*atan(elem.parent.ν_fact*tan(EA/2))
 
    #  return orbitsolve_ν(elem, ν, EA)


    sinν_ω, cosν_ω = sincos(ω + ν)
    sinΩ, cosΩ = sincos(Ω)
    ecosν = e*cos(ν)
    cosi = cos(i)

    # p = a*oneminusesq # semi-latus rectum [AU]
    # r = p/(1 + ecosν)
    # xcart = r*(cosν_ω*sinΩ + sinν_ω*cosi*cosΩ) # [AU]
    # ycart = r*(cosν_ω*cosΩ - sinν_ω*cosi*sinΩ) # [AU]
    # # Now we have sep in 
    # sep_au = sqrt(xcart^2 + ycart^2)
    # sep_mas = sep_au*cart2angle

    ####

    sep_au = sep/cart2angle
    # sep_au = sqrt(
    #     (r*(cosν_ω*sinΩ + sinν_ω*cosi*cosΩ))^2 + 
    #     (r*(cosν_ω*cosΩ - sinν_ω*cosi*sinΩ))^2
    # )
    # sep_au = r*sqrt(
    #     (cosν_ω*sinΩ + sinν_ω*cosi*cosΩ)^2 +
    #     (cosν_ω*cosΩ - sinν_ω*cosi*sinΩ)^2
    # )
    r = sep_au / sqrt(
        (cosν_ω*sinΩ + sinν_ω*cosi*cosΩ)^2 +
        (cosν_ω*cosΩ - sinν_ω*cosi*sinΩ)^2
    )
    # r = p/(1 + ecosν)
    p = r*(1+ecosν)
    oneminusesq = (1 - e^2)
    # p = a*oneminusesq # semi-latus rectum [AU]
    # p = a*oneminusesq # semi-latus rectum [AU]
    a = p/oneminusesq

    # @show a




    period_days = √(a^3/M)*PlanetOrbits.kepler_year_to_julian_day_conversion_factor
    period_yrs = period_days/PlanetOrbits.year2day_julian
    n = 2π/period_yrs # mean motion

    # Get epoch of periastron passage
    tₚ = theta_epoch - MA/n*PlanetOrbits.year2day_julian
    return tₚ, a
end
