using Dates
using LinearAlgebra
using HORIZONS: HORIZONS

const hipparcos_catalog_epoch_mjd = years2mjd(1991.15)

###########################
# Stores data:
const hip_iad_cols = (:iorb, :epoch, :parf, :cpsi, :spsi, :res, :sres)
struct HipparcosIADLikelihood{THipSol,TIADTable<:Table} <: AbstractLikelihood
    hip_sol::THipSol
    table::TIADTable
    function HipparcosIADLikelihood(hip_sol, observations...)
        iad_table = Table(observations...)
        if !issubset(hip_iad_cols, Tables.columnnames(iad_table))
            error("Expected columns $hip_iad_cols")
        end
        if hip_sol.isol_n != 5
            @warn "Only 5 parameter solutions are currently implemented"
        end
        return new{typeof(hip_sol),typeof(iad_table)}(hip_sol, iad_table)
    end
end

##########################
# Loads data:
"""
    HipparcosIADLikelihood(;hipparcos_id)
"""
function HipparcosIADLikelihood(; hip_id, catalog=(datadep"Hipparcos_IAD"))

    # TODO: copy in file search functionality from OctofitterRadialVelocity
    file = @sprintf("H%06d.d", hip_id)
    fname = joinpath(catalog, "ResRec_JavaTool_2014", file[1:4], file)

    lines = readlines(fname)

    # HIP    MCE    NRES NC isol_n SCE  F2     F1
    # 27321  27251  111  1  5      0    -1.63  0 
    hip, mce, nres, nc, isol_n, sce, f2, f1 = parse.(Float64, split(lines[7])[2:end])

    # Hp      B-V    VarAnn NOB NR
    # 3.9077  0.171  0      111 0  
    hp, b_m_v, varann, nob, nr = parse.(Float64, split(lines[9])[2:end])

    (
        radeg,
        dedeg,
        plx,
        pm_ra,
        pm_de,
        e_ra,
        e_de,
        e_plx,
        e_pmra,
        e_pmde,
        dpmra,
        dpmde,
        e_dpmra,
        e_dpmde,
        ddpmra,
        ddpmde,
        e_ddpmra,
        e_ddpmde,
        upsra,
        upsde,
        e_upsra,
        e_upsde,
        var
    ) = tryparse.(Float64, split(lines[11])[2:end])
    hip_sol = (;
        hip, mce, nres, nc, isol_n, sce, f2, f1,
        hp, b_m_v, varann, nob, nr,
        radeg, dedeg, plx, pm_ra, pm_de, e_ra, e_de, e_plx,
        e_pmra, e_pmde, dpmra, dpmde, e_dpmra, e_dpmde, ddpmra, ddpmde,
        e_ddpmra, e_ddpmde, upsra, upsde, e_upsra, e_upsde, var
    )

    iad_table_rows = NamedTuple[]
    for line in lines[13:end]
        if startswith(line, '#')
            continue
        end
        iorb, epoch, parf, cpsi, spsi, res, sres = split(line)
        push!(iad_table_rows, (;
            iorb=parse(Int, iorb),
            epoch_yrs=parse(Float64, epoch),
            parf=parse(Float64, parf),
            cpsi=parse(Float64, cpsi),
            spsi=parse(Float64, spsi),
            res=parse(Float64, res),
            sres=parse(Float64, sres),
        ))
    end

    iad_table = FlexTable(iad_table_rows)
    # transform epoch to MJD
    iad_table.epoch = years2mjd.(1991.25 .+ iad_table.epoch_yrs)

    # Remove rejected scans, if any
    tokeep = iad_table.sres .> 0

    if !all(tokeep)
        @warn "rejected scans present" mean(tokeep)
    end
    iad_table = iad_table[tokeep, :]

    # Query the HORIZONS system for the Earth-Moon barycentre position at each specific epoch
    # of the dataset. 
    # This incurs a separate query for each epoch.
    @info "Querying Earth Geocentre position from HORIZONS"
    earth_pos_vel = map(iad_table.epoch) do epoch
        print(".")
        fname = tempname()
        t = DateTime(mjd2date(epoch))
        HORIZONS.vec_tbl("Geocenter", t, t + Hour(1), Day(1); FILENAME=fname, CENTER="@ssb", REF_PLANE="FRAME", OUT_UNITS="AU-D", CSV_FORMAT=true, VEC_TABLE=2)
        lines = readlines(fname)
        rm(fname)
        i = 0
        for line in lines
            i += 1
            if startswith(line, raw"$$SOE")
                break
            end
        end
        record = split(rstrip(lines[i+1], ','), ", ")
        x, y, z, vx, vy, vz = parse.(Float64, record[3:8])
        return (; x, y, z, vx, vy, vz)
    end
    println("done")
    table = FlexTable(eachcol(iad_table)..., eachcol(earth_pos_vel)...)
    
    ## Nielsen eq 10
    # Un-re-normalize uncertainties
    D = length(table.epoch) - 6
    G = f2
    f = (G√(2 / 9D) + 1 - (2 / 9D))^(3 / 2)
    table.sres_renorm = table.sres .* f
    ##

    # Nielsen et al. Beta Pic modelling Eq. 1 & 2
    # We start by calculating the sky-path (proper motion, parallax) given
    # the best Hipparcos solution. This is a simplified model that works
    # entirely in the tangent plane about the catalog position (α₀, δ₀).
    # The Hipparcos IAD catalog then gives us scan directions and residuals
    # versus this model---basically we are undoing the Hipparcos model
    # to get at the "raw"(-er) data.
    μα✱ = hip_sol.pm_ra # mas/yr -- /cos(δ₀)
    μδ = hip_sol.pm_de  # mas/yr
    α₀ = hip_sol.radeg  # deg
    δ₀ = hip_sol.dedeg  # deg

    # NB: the α✱ variables are multiplied by cos(δ₀)

    # Calculate sky coordinate Delta at each scan epoch from the catalog position
    # using the eath's motion ("ephemeris"), x,y,z in AU.
    table.Δα✱ = @. hip_sol.plx * (
        table.x * sind(α₀) -
        table.y * cosd(α₀)
    ) + (table.epoch_yrs) * μα✱
    table.Δδ = @. hip_sol.plx * (
        table.x * cosd(α₀) * sind(δ₀) +
        table.y * sind(α₀) * cosd(δ₀) -
        table.z * cosd(δ₀)
    ) + (table.epoch_yrs) * μδ
    
    # Nielsen Eq. 3: abcissa point
    table.α✱ₐ = @. table.res * table.cpsi + table.Δα✱
    table.δₐ = @. table.res * table.spsi + table.Δδ

    # These two points (picked ± 1 mas) specify a line--- the Hipparcos `sres`
    # is then then the 
    table.α✱ₘ = collect(eachrow(@. [-1, 1]' * table.spsi + table.α✱ₐ))
    table.δₘ = collect(eachrow(@. [1, -1]' * table.cpsi + table.δₐ))

    return HipparcosIADLikelihood(hip_sol, table)
end
export HipparcosIADLikelihood


##########################
# Computes likelihood
function ln_like(
    hiplike::HipparcosIADLikelihood,
    θ_system,
    orbits,
    num_epochs::Val{L}=Val(length(hiplike.table.epochs))
) where L
    ll = 0.0

    # All planets inthe system have orbits defined with the same ra, dec, and proper motion,
    # since these are properties of the system.
    first_planet_orbit = first(orbits)
    if length(orbits) > 1
        for i in eachindex(orbits)[2:end]
            if orbits[i].ra != first_planet_orbit.ra ||
               orbits[i].dec != first_planet_orbit.dec ||
               orbits[i].pmra != first_planet_orbit.rpma ||
               orbits[i].pmdec != first_planet_orbit.ra pmdec
                error("Planet orbits do not have matching ra, dec, pmpra, and pmdec.")
            end
        end
    end
    
    # Solution at the reference epoch (ignoring perturbation)
    sol₀ = orbitsolve(first_planet_orbit, θ_system.ref_epoch)

    for i in eachindex(hiplike.table.epoch)

        first_planet_sol = orbitsolve(first_planet_orbit, hiplike.table.epoch[i])

        # Rigorous propagation of coordinates
        α✱_model = first_planet_sol.compensated.parallax2 * (
            hiplike.table.x[i] * sind(first_planet_sol.compensated.ra2) -
            hiplike.table.y[i] * cosd(first_planet_sol.compensated.ra2)
        ) + (first_planet_sol.compensated.ra2 - sol₀.compensated.ra2)*cosd(sol₀.compensated.dec2)*60*60*1000 

        δ_model = first_planet_sol.compensated.parallax2 * (
            hiplike.table.x[i] * cosd(first_planet_sol.compensated.ra2) * sind(first_planet_sol.compensated.dec2) +
            hiplike.table.y[i] * sind(first_planet_sol.compensated.ra2) * cosd(first_planet_sol.compensated.dec2) -
            hiplike.table.z[i] * cosd(first_planet_sol.compensated.dec2)
        ) + (first_planet_sol.compensated.dec2 - sol₀.compensated.dec2)*60*60*1000



        # # Simplified tangent-plane model (matches Hipparcos)
        # α✱_model = first_planet_orbit.plx * (
        #     hiplike.table.x[i] * sind(first_planet_orbit.ra) -
        #     hiplike.table.y[i] * cosd(first_planet_orbit.ra)
        # ) + (first_planet_sol.compensated.epoch2 - sol₀.compensated.epoch2)*first_planet_orbit.pmra

        # δ_model = first_planet_orbit.plx * (
        #     hiplike.table.x[i] * cosd(first_planet_orbit.ra) * sind(first_planet_orbit.dec) +
        #     hiplike.table.y[i] * sind(first_planet_orbit.ra) * cosd(first_planet_orbit.dec) -
        #     hiplike.table.z[i] * cosd(first_planet_orbit.dec)
        # ) + (first_planet_sol.compensated.epoch2 - sol₀.compensated.epoch2)*first_planet_orbit.pmdec

        α✱_model_with_perturbation = α✱_model
        δ_model_with_perturbation = δ_model

        # Add perturbations from all planets
        for planet_i in eachindex(orbits)
            orbit = orbits[planet_i]
            if orbit == first_planet_orbit
                sol = first_planet_sol
            else
                sol = orbitsolve(orbit, hiplike.table.epoch[i])
            end
            # Add perturbation from planet
            α✱_model_with_perturbation += raoff(sol, θ_system.planets[planet_i].mass*Octofitter.mjup2msol)
            δ_model_with_perturbation += decoff(sol, θ_system.planets[planet_i].mass*Octofitter.mjup2msol)
        end

        # (X,Y) point of star, relative to re-epoch (ra,dec) position
        point = @SVector [
            α✱_model_with_perturbation,
            δ_model_with_perturbation
        ]
        # Two points defining a line along which the star's position was measured
        line_point_1 = @SVector [hiplike.table.α✱ₘ[i][1], hiplike.table.δₘ[i][1]]
        line_point_2 = @SVector [hiplike.table.α✱ₘ[i][2], hiplike.table.δₘ[i][2]]
        # Distance from model star Delta position to this line where it was measured 
        # by the satellite
        resid = distance_point_to_line(point, line_point_1, line_point_2)
        # Re-normalized uncertainty on the location of this line (not angle of line, no 
        # uncertainties considered there).
        σ = hiplike.table.sres_renorm[i]

        # Add likelihood at this epoch.
        ll += logpdf(Normal(0, σ), resid)

    end

    return ll
end

function simulate(hiplike::HipparcosIADLikelihood, θ_system, orbits)
    α✱_model_with_perturbation_out = zeros(length(hiplike.table.epoch))
    δ_model_with_perturbation_out = zeros(length(hiplike.table.epoch))
    resid_out = zeros(length(hiplike.table.epoch))


    # All planets inthe system have orbits defined with the same ra, dec, and proper motion,
    # since these are properties of the system.
    first_planet_orbit = first(orbits)
    if length(orbits) > 1
        for i in eachindex(orbits)[2:end]
            if orbits[i].ra != first_planet_orbit.ra ||
               orbits[i].dec != first_planet_orbit.dec ||
               orbits[i].pmra != first_planet_orbit.rpma ||
               orbits[i].pmdec != first_planet_orbit.ra pmdec
                error("Planet orbits do not have matching ra, dec, pmpra, and pmdec.")
            end
        end
    end
    
    # Solution at the catalog epoch (ignoring perturbation)
    sol₀ = orbitsolve(first_planet_orbit, θ_system.ref_epoch)

    for i in eachindex(hiplike.table.epoch)

        first_planet_sol = orbitsolve(first_planet_orbit, hiplike.table.epoch[i])

        # Rigorous propagation of coordinates
        α✱_model = first_planet_sol.compensated.parallax2 * (
            hiplike.table.x[i] * sind(first_planet_sol.compensated.ra2) -
            hiplike.table.y[i] * cosd(first_planet_sol.compensated.ra2)
        ) + (first_planet_sol.compensated.ra2 - sol₀.compensated.ra2)*cosd(sol₀.compensated.dec2)*60*60*1000 

        δ_model = first_planet_sol.compensated.parallax2 * (
            hiplike.table.x[i] * cosd(first_planet_sol.compensated.ra2) * sind(first_planet_sol.compensated.dec2) +
            hiplike.table.y[i] * sind(first_planet_sol.compensated.ra2) * cosd(first_planet_sol.compensated.dec2) -
            hiplike.table.z[i] * cosd(first_planet_sol.compensated.dec2)
        ) + (first_planet_sol.compensated.dec2 - sol₀.compensated.dec2)*60*60*1000



        # # Simplified tangent-plane model (matches Hipparcos)
        # α✱_model = first_planet_orbit.plx * (
        #     hiplike.table.x[i] * sind(first_planet_orbit.ra) -
        #     hiplike.table.y[i] * cosd(first_planet_orbit.ra)
        # ) + (first_planet_sol.compensated.epoch2 - sol₀.compensated.epoch2)*first_planet_orbit.pmra

        # δ_model = first_planet_orbit.plx * (
        #     hiplike.table.x[i] * cosd(first_planet_orbit.ra) * sind(first_planet_orbit.dec) +
        #     hiplike.table.y[i] * sind(first_planet_orbit.ra) * cosd(first_planet_orbit.dec) -
        #     hiplike.table.z[i] * cosd(first_planet_orbit.dec)
        # ) + (first_planet_sol.compensated.epoch2 - sol₀.compensated.epoch2)*first_planet_orbit.pmdec



        α✱_model_with_perturbation = α✱_model
        δ_model_with_perturbation = δ_model

        # Add perturbations from all planets
        for planet_i in eachindex(orbits)
            orbit = orbits[planet_i]
            if orbit == first_planet_orbit
                sol = first_planet_sol
            else
                sol = orbitsolve(orbit, hiplike.table.epoch[i])
            end
            # Add perturbation from planet
            α✱_model_with_perturbation += raoff(sol, θ_system.planets[planet_i].mass*Octofitter.mjup2msol)
            δ_model_with_perturbation += decoff(sol, θ_system.planets[planet_i].mass*Octofitter.mjup2msol)
        end

        # (X,Y) point of star, relative to re-epoch (ra,dec) position
        point = @SVector [
            α✱_model_with_perturbation,
            δ_model_with_perturbation
        ]
        # Two points defining a line along which the star's position was measured
        line_point_1 = @SVector [hiplike.table.α✱ₘ[i][1], hiplike.table.δₘ[i][1]]
        line_point_2 = @SVector [hiplike.table.α✱ₘ[i][2], hiplike.table.δₘ[i][2]]
        # Distance from model star Delta position to this line where it was measured 
        # by the satellite
        resid = distance_point_to_line(point, line_point_1, line_point_2)


        α✱_model_with_perturbation_out[i] = α✱_model_with_perturbation
        δ_model_with_perturbation_out[i]  = δ_model_with_perturbation
        resid_out[i] = resid
    end


    # for i in eachindex(hiplike.table.epoch), planet_i in eachindex(orbits)
    #     orbit = orbits[planet_i]
    #     sol = orbitsolve(orbit, hiplike.table.epoch[i])
    #     α✱_model = sol.compensated.parallax2 * (
    #         hiplike.table.x[i] * sind(sol.compensated.ra2) -
    #         hiplike.table.y[i] * cosd(sol.compensated.ra2)
    #     ) + (sol.compensated.ra2 - orbit.ra)*cosd(orbit.dec)*60*60*1000 
    #     δ_model = sol.compensated.parallax2 * (
    #         hiplike.table.x[i] * cosd(sol.compensated.ra2) * sind(sol.compensated.dec2) +
    #         hiplike.table.y[i] * sind(sol.compensated.ra2) * cosd(sol.compensated.dec2) -
    #         hiplike.table.z[i] * cosd(sol.compensated.dec2)
    #     ) + (sol.compensated.dec2 - orbit.dec)*60*60*1000
    #     α✱_model_with_perturbation[i,planet_i] = α✱_model + raoff(sol, θ_system.planets[planet_i].mass)
    #     δ_model_with_perturbation[i,planet_i] = δ_model + decoff(sol, θ_system.planets[planet_i].mass)
    #     point = @SVector [
    #         α✱_model_with_perturbation[i,planet_i],
    #         δ_model_with_perturbation[i,planet_i]
    #     ]
    #     line_point_1 = @SVector [hiplike.table.α✱ₘ[i][1], hiplike.table.δₘ[i][1]]
    #     line_point_2 = @SVector [hiplike.table.α✱ₘ[i][2], hiplike.table.δₘ[i][2]]
    #     resid[i,planet_i] = distance_point_to_line(point, line_point_1, line_point_2)
    # end
    return Table(;
        α✱_model_with_perturbation=α✱_model_with_perturbation_out,
        δ_model_with_perturbation=δ_model_with_perturbation_out,
        resid=resid_out
    )
end





# For fitting, sigma error  is sres
# I guess we think of the distance from the line as the residual
"""
    distance_point_to_line(point, line_point_1, line_point_2)

Given three points (each a vector with two values), calculate the distance
from `point` to the line defined by `line_point_1` and `line_point_2`.
"""
function distance_point_to_line(point, line_point_1, line_point_2)
    r₀ = point
    r₁ = line_point_1
    r₂ = line_point_2
    d = abs(
        (r₂[1] - r₁[1]) * (r₁[2] - r₀[2]) - (r₁[1] - r₀[1]) * (r₂[2] - r₁[2])
    ) / norm(r₂ - r₁)
end