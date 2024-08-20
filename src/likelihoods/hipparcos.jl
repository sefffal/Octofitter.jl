using Dates
using LinearAlgebra
using HORIZONS: HORIZONS

# Table 1.2.3 from the 1997 book, vol 1.
# RV increasing away from the observer
const hipparcos_stellar_rvs_used_for_cat_solutions = Table([
    (hip_id=439   , rv_kms= +22.9   ),
    (hip_id=3829  , rv_kms= -38.0   ),
    (hip_id=5336  , rv_kms= -98.1   ),
    (hip_id=15510 , rv_kms= +86.7   ),
    (hip_id=19849 , rv_kms= -42.7   ),
    (hip_id=24186 , rv_kms= +245.5  ),
    (hip_id=26857 , rv_kms= +105.6  ),
    (hip_id=54035 , rv_kms= -84.3   ),
    (hip_id=54211 , rv_kms= +68.8   ),
    (hip_id=57939 , rv_kms= -99.1   ),
    (hip_id=70890 , rv_kms= -16.0   ),
    (hip_id=71681 , rv_kms= -18.1   ),
    (hip_id=71683 , rv_kms= -26.2   ),
    (hip_id=74234 , rv_kms= +308.0  ),
    (hip_id=74235 , rv_kms= +294.3  ),
    (hip_id=86990 , rv_kms= -115.0  ),
    (hip_id=87937 , rv_kms= -111.0  ),
    (hip_id=99461 , rv_kms= -129.8  ),
    (hip_id=104214, rv_kms= -64.8   ),
    (hip_id=104217, rv_kms= -64.3   ),
    (hip_id=108870, rv_kms= -40.4   ),
])

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
        return new{typeof(hip_sol),typeof(iad_table)}(hip_sol, iad_table)
    end
end

##########################
# Loads data:
"""
    HipparcosIADLikelihood(;hipparcos_id)
"""
function HipparcosIADLikelihood(; hip_id, catalog=(datadep"Hipparcos_IAD"), renormalize=true)

    file = @sprintf("H%06d.d", hip_id)
    fname = joinpath(catalog, "ResRec_JavaTool_2014", file[1:4], file)

    lines = readlines(fname)

    # See table 1 of https://arxiv.org/pdf/1108.4971 for units

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

    if isol_n ∉ (5,7,9)
        error("Only stars with solution types 5, 7, or 9 are currently supported. This sol is $isol_n. If you need solution type 1, please open an issue on GitHub and we would be happy to add it.")
    end

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
    # iad_table.epoch = years2mjd.(hipparcos_catalog_epoch_decimal_year .+ iad_table.epoch_yrs)
    iad_table.epoch = hipparcos_catalog_epoch_mjd .+ iad_table.epoch_yrs*julian_year

    # Remove rejected scans, if any
    iad_table.reject = iad_table.sres .<= 0
    if any(iad_table.reject)
        @warn "rejected scans present" count(iad_table.reject)
    end

    earth_pos_vel = geocentre_position_query.(iad_table.epoch)
    table = FlexTable(eachcol(iad_table)..., eachcol(earth_pos_vel)...)
    
    ## Nielsen eq 10
    # Un-re-normalize uncertainties
    if renormalize
        D = length(table.epoch) - 6
        G = f2
        f = (G√(2 / 9D) + 1 - (2 / 9D))^(3 / 2)
        table.sres_renorm = table.sres .* f
        @info "renormalizing hipparcos uncertainties according to Nielsen et al (2020). Pass `renormalize=false` to disable."
    else
        table.sres_renorm = table.sres
        @info "Not renormalizing hipparcos uncertainties according to Nielsen et al (2020). Pass `renormalize=true` to enable."
    end
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

    # For 21 stars that are nearby and have high radial velocity, the Hipparcos team
    # did account for changing parallax. 
    # We have to reproduce that when undoing their model.
    # Note that we can use a different RV in our actual fits if we want to, but we have to
    # use their exact same RVs to recover the scan data correctly.
    # As a result, parallax is a function of time

    I = findfirst(==(hip_id), hipparcos_stellar_rvs_used_for_cat_solutions.hip_id)
    if isnothing(I)
        rv_kms = 0.
    else
        rv_kms = hipparcos_stellar_rvs_used_for_cat_solutions.rv_kms[I]
        @info "This star was flagged by the Hipparcos team as having perspective acceleration effects greater than 0.1 mas/yr. They (and we) correct for it using the following epoch RV" rv_kms
    end

    plx0 = hip_sol.plx
    dist0 = 1000/plx0
    δdist_pc_δt_sec = rv_kms / IAU_pc2km / PlanetOrbits.sec2day
    Δdist_pc = δdist_pc_δt_sec .* (table.epoch .- hipparcos_catalog_epoch_mjd)
    dist1 = dist0 .+ Δdist_pc
    table.plx_vs_time = 1000 ./ dist1
    # for use in likelihood function
    table.rv_kms = fill(rv_kms, length(table.plx_vs_time))

    # Calculate sky coordinate Delta at each scan epoch from the catalog position
    # using the eath's motion ("ephemeris"), x,y,z in AU.
    table.Δα✱ = @. table.plx_vs_time * (
        table.x * sind(α₀) -
        table.y * cosd(α₀)
    # Both of these calculations are correct, that's how we know the dates are computed correctly:
    # ) + (table.epoch_yrs) * μα✱
    ) + (table.epoch - hipparcos_catalog_epoch_mjd)/julian_year * μα✱
    table.Δδ = @. table.plx_vs_time * (
        table.x * cosd(α₀) * sind(δ₀) +
        table.y * sind(α₀) * sind(δ₀) -
        table.z * cosd(δ₀)
    # Both of these calculations are correct, that's how we know the dates are computed correctly:
    # ) + (table.epoch_yrs) * μδ
    ) + (table.epoch - hipparcos_catalog_epoch_mjd)/julian_year * μδ

         
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
    T = _system_number_type(θ_system)
    ll = zero(T)

    hip_model = simulate(hiplike, θ_system, orbits)
    for i in eachindex(hip_model.resid)
        if hiplike.table.reject[i]
            continue
        end
        ll+= logpdf(Normal(0, hiplike.table.sres_renorm[i]),hip_model.resid[i])
    end

    return ll
end

function simulate(hiplike::HipparcosIADLikelihood, θ_system, orbits)

    T = _system_number_type(θ_system)
    α✱_model_with_perturbation_out = zeros(T,length(hiplike.table.epoch))
    δ_model_with_perturbation_out = zeros(T,length(hiplike.table.epoch))
    resid_out = zeros(T,length(hiplike.table.epoch))

    # All planets inthe system have orbits defined with the same ra, dec, and proper motion,
    # since these are properties of the system.
    orbit = first(orbits)
    if length(orbits) > 1
        for i in eachindex(orbits)[2:end]
            if orbits[i].ra != orbit.ra ||
               orbits[i].dec != orbit.dec ||
               orbits[i].pmra != orbit.rpma ||
               orbits[i].pmdec != orbit.ra pmdec
                error("Planet orbits do not have matching ra, dec, pmpra, and pmdec.")
            end
        end
    end

    for i in eachindex(hiplike.table.epoch)

        orbitsol_hip_epoch = orbitsolve(orbit, hiplike.table.epoch[i])

        # Simplified tangent-plane model (matches Hipparcos)
        # The parallax versus time is considered for 21 high RV, nearby stars.
        # This matches the Hipparcos reduction.
        # For other stars, rv_kms will be zero and this isn't considered.
        delta_t_days = orbitsol_hip_epoch.t - hipparcos_catalog_epoch_mjd
        # Proper motion: each individual year is either 365 or 366 days long. Need to be careful!
        # the formal definition is to use Julian years, 365.25000000 exactly.
        delta_time_julian_year = delta_t_days/julian_year
        plx0 = orbit.plx
        if hiplike.table.rv_kms != 0
            dist0 = 1000/plx0
            δdist_pc_δt_sec = hiplike.table.rv_kms[i] / IAU_pc2km / PlanetOrbits.sec2day
            Δdist_pc = δdist_pc_δt_sec * (hiplike.table.epoch[i] - hipparcos_catalog_epoch_mjd)
            dist1 = dist0 + Δdist_pc
            plx_at_epoch = 1000 / dist1
        else
            plx_at_epoch = plx0
        end
        # TODO: could hoist the trig into the constructor for speed. It doesn't change in this simplified
        # model.
        α✱_model = (θ_system.ra-hiplike.hip_sol.radeg)*60*60*1000*cosd(hiplike.hip_sol.dedeg) + plx_at_epoch * (
            hiplike.table.x[i] * sind(hiplike.hip_sol.radeg) -
            hiplike.table.y[i] * cosd(hiplike.hip_sol.radeg)
        ) + delta_time_julian_year * orbit.pmra
        δ_model = (θ_system.dec-hiplike.hip_sol.dedeg)*60*60*1000 + plx_at_epoch * (
            hiplike.table.x[i] * cosd(hiplike.hip_sol.radeg) * sind(hiplike.hip_sol.dedeg) +
            hiplike.table.y[i] * sind(hiplike.hip_sol.radeg) * sind(hiplike.hip_sol.dedeg) -
            hiplike.table.z[i] * cosd(hiplike.hip_sol.dedeg)
        ) + delta_time_julian_year * orbit.pmdec

        α✱_perturbation = zero(T)
        δ_perturbation = zero(T)

        # Add perturbations from all planets
        for planet_i in eachindex(orbits)
            orbit = orbits[planet_i]
            if orbit == orbit
                sol = orbitsol_hip_epoch
            else
                sol = orbitsolve(orbit, hiplike.table.epoch[i])
            end
            # Add perturbation from planet
            α✱_perturbation += raoff(sol, θ_system.planets[planet_i].mass*Octofitter.mjup2msol)
            δ_perturbation += decoff(sol, θ_system.planets[planet_i].mass*Octofitter.mjup2msol)
        end
        # coordinates in degrees of RA and Dec
        α✱_model_with_perturbation = α✱_model + α✱_perturbation/cosd(hiplike.hip_sol.dedeg)
        δ_model_with_perturbation = δ_model + δ_perturbation

        # Calculate differences in milliarcseconds by mapping into a local tangent plane,
        # with the local coordinate defined by the Hipparcos solution at this *epoch* (not 
        # at the best solution)
        point = @SVector [
            α✱_model_with_perturbation,
            δ_model_with_perturbation,
        ]

        # Two points defining a line along which the star's position was measured
        line_point_1 = @SVector [
            hiplike.table.α✱ₘ[i][1],
            hiplike.table.δₘ[i][1],
        ]
        line_point_2 = @SVector [
            hiplike.table.α✱ₘ[i][2],
            hiplike.table.δₘ[i][2],
        ]
        # Distance from model star Delta position to this line where it was measured 
        # by the satellite
        resid = distance_point_to_line(point, line_point_1, line_point_2) # degrees

        α✱_model_with_perturbation_out[i] = α✱_model_with_perturbation
        δ_model_with_perturbation_out[i]  = δ_model_with_perturbation
        resid_out[i] = resid
    end

    # So we can add the model ra and dec to each of our model outputs
    # And then compare to the data outputs + the 1991.25 position

    return (;
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


function idenfity_rejectable_scans(hiplike, θ_system, orbits)
    hip_model = simulate(hiplike, θ_system, orbits)
    for i in eachindex(hip_model.resid)
        print("$i\t", (hip_model.resid[i]/hiplike.table.sres_renorm[i])^2)
        if hiplike.table.reject[i]
            println("\t-already rejected")
            continue
        end
        println()
    end
end