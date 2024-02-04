
using Octofitter.StaticArrays

# We want this to work with both Campbell and Thiele-Innes parameterizations
export θ_at_epoch_to_tperi
function θ_at_epoch_to_tperi(system,planet,theta_epoch)
    (;M) = system
    (;e,θ) = planet
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
            a = ∛(system.M * b.P^2)
        else
            error("Your planet must specify either i, Ω, ω, and a/P, or B, G, A, F.")
        end

        (;a, i, Ω, ω) = planet
        A = a*( cos(Ω)*cos(ω)-sin(Ω)*sin(ω)*cos(i))
        B = a*( sin(Ω)*cos(ω)+cos(Ω)*sin(ω)*cos(i))
        F = a*(-cos(Ω)*sin(ω)-sin(Ω)*cos(ω)*cos(i))
        G = a*(-sin(Ω)*sin(ω)+cos(Ω)*cos(ω)*cos(i))
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

    n = 2π/√(a^3/M) 

    # Get epoch of periastron passage
    tₚ = theta_epoch - MA/n*PlanetOrbits.year2day
    return tₚ
end
