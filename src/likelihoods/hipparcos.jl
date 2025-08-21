using Dates
using LinearAlgebra
using Combinatorics

# Table 1.2.3 from the 1997 book, vol 1.
# RV increasing away from the observer
# These are only applicable to data from the 1997 reduction, 
# not the van Leeuwan reduction.
const hipparcos_stellar_rvs_used_for_cat_solutions = Table([
    (hip_id=439, rv_kms=+22.9),
    (hip_id=3829, rv_kms=-38.0),
    (hip_id=5336, rv_kms=-98.1),
    (hip_id=15510, rv_kms=+86.7),
    (hip_id=19849, rv_kms=-42.7),
    (hip_id=24186, rv_kms=+245.5),
    (hip_id=26857, rv_kms=+105.6),
    (hip_id=54035, rv_kms=-84.3),
    (hip_id=54211, rv_kms=+68.8),
    (hip_id=57939, rv_kms=-99.1),
    (hip_id=70890, rv_kms=-16.0),
    (hip_id=71681, rv_kms=-18.1),
    (hip_id=71683, rv_kms=-26.2),
    (hip_id=74234, rv_kms=+308.0),
    (hip_id=74235, rv_kms=+294.3),
    (hip_id=86990, rv_kms=-115.0),
    (hip_id=87937, rv_kms=-111.0),
    (hip_id=99461, rv_kms=-129.8),
    (hip_id=104214, rv_kms=-64.8),
    (hip_id=104217, rv_kms=-64.3),
    (hip_id=108870, rv_kms=-40.4),
])

# Table F.1 of van Leeuwan, column 1.
hipparcos_catalog_orbit_parameters_used = [
    2, 171, 443, 518, 677, 999, 1242, 1349, 1674, 2170, 2237, 2237, 2552, 2762,
    2912, 2941, 3504, 4809, 4849, 5249, 5300, 5336, 5531, 5685, 5778, 5842,
    6486, 6564, 6867, 7078, 7213, 7372, 7580, 7918, 7981, 8514, 8833, 8882,
    8903, 8922, 9236, 9480, 9500, 38538, 40167, 40239, 40326, 41261, 41426,
    41820, 42075, 42430, 42455, 42805, 43109, 43671, 44248, 44471, 44676,
    45075, 45170, 45527, 45571, 45617, 46396, 46404, 46454, 46651, 46706,
    46893, 47479, 48348, 49166, 49841, 50805, 51233, 51384, 51885, 51986, 52085,
    52271, 52419, 53240, 53423, 53763, 54061, 54155, 54204, 54632, 54677, 54977,
    55016, 55203, 55266, 55425, 55642, 56290, 56528, 56675, 56731, 57363, 57565,
    57791, 57994, 58113, 58590, 58799, 59459, 59468, 59750, 59780, 59816, 59856,
    60129, 60299, 60994, 61724, 61880, 61932, 61941, 62124, 62145, 62371, 62409,
    63406, 63503, 63613, 63742, 64241, 65026, 65135, 65203, 65417, 65420, 65783,
    66438, 66458, 66640, 67234, 67422, 67483, 67927, 68682, 68756, 69112, 69176,
    69283, 69879, 70327, 70857, 70973, 71094, 71141, 71469, 71510, 71683, 71729,
    71914, 72217, 72479, 72659, 72848, 73182, 73199, 73507, 73695, 73787, 74392,
    74893, 75312, 75379, 75389, 75411, 75415, 75508, 75695, 75695, 75949, 76031,
    76267, 76382, 76466, 76734, 76752, 76852, 76952, 77541, 77760, 78459, 78662,
    78727, 78918, 79101, 80166, 80346, 80725, 80816, 81023, 81126, 81693, 81726,
    81754, 82020, 82817, 82860, 83575, 83838, 83895, 84012, 84123, 84140, 84179,
    84425, 84709, 84924, 84949, 85019, 85106, 85141, 85333, 85582, 85667, 85683,
    85727, 85749, 85846, 86036, 86221, 86254, 86373, 86400, 86722, 87204, 87428,
    87655, 87861, 87895, 87899, 88404, 88436, 88601, 88637, 88715, 88728, 88745,
    88932, 88964, 89507, 89808, 89937, 90659, 91009, 91196, 91394, 91395, 91751,
    92122, 92175, 92512, 92677, 92757, 92818, 93017, 93244, 93383, 93506, 93574,
    93864, 93995, 94076, 94252, 94349, 94521, 94739, 95066, 95477, 95501, 95769,
    95995, 96302, 96683, 96807, 97016, 97222, 97237, 97837, 98001, 98416, 99312,
    99376, 99473, 99668, 99675, 99848, 99965, 100259, 100345, 101093, 101235,
    101382, 101750, 101769, 101955, 101958, 102782, 103519, 103546, 103655,
    104019, 104486, 104788, 104858, 104887, 105200, 105431, 105712, 105881,
    105947, 105969, 106711, 106811, 106972, 106985, 107354, 107522, 107522,
    107788, 107818, 108084, 108431, 108431, 108478, 109554, 110102, 110108,
    110130, 110893, 111062, 111170, 111200, 111314, 111528, 111685, 111805,
    111965, 111965, 111974, 112158, 112915, 113048, 113445, 113718, 113860,
    113996, 114222, 114273, 114421, 114922, 115031, 115126, 116233, 116436,
    116849, 117570, 117607, 117666, 118169, 118266
]

###########################
# Stores data:
const hip_iad_cols = (:iorb, :epoch, :parf, :cosϕ, :sinϕ, :res, :sres)
struct HipparcosIADLikelihood{THipSol,TIADTable<:Table,TDist,TFact} <: AbstractLikelihood
    hip_sol::THipSol
    table::TIADTable
    priors::Priors
    derived::Derived
    name::String
    # precomputed MvNormal distribution
    dist::TDist
    A_prepared_4::TFact
    A_prepared_5::TFact
end

# Add likelihoodname method
likelihoodname(::HipparcosIADLikelihood) = "Hipparcos IAD"

function likeobj_from_epoch_subset(obs::HipparcosIADLikelihood, obs_inds)
    return HipparcosIADLikelihood(obs.hip_sol, obs.table[obs_inds, :, 1], obs.priors, obs.derived, obs.name, obs.dist, obs.A_prepared_4, obs.A_prepared_5)
end

"""
    HipparcosIADLikelihood(;
        hip_id,
        variables=@variables begin
            fluxratio ~ Product([Uniform(0, 1), Uniform(0, 1)])  # one entry for each companion
        end
    )

Load the Hipparcos IAD likelihood.
By default, this fetches and catches the extracted Java Tool edition of the
van Leeuwan reduction. 

The `fluxratio` variable should be a Product distribution containing the flux ratio of each companion
in the same order as the planets in the system.

Additional arguments:
* `catalog`: path to the data directory. By default, will fetch from online and cache.
* `renormalize=true`: renormalize the uncertainties according to Nielsen et al. (2020)
* `attempt_correction=true`: perform ad-hoc correction of any corrupted scans according to G. Brandt et al. (2021)
* `is_van_leeuwen=true`: set to false if using e.g. the 1997 reduction.
    This impacts how the residuals are reconstructed. The 1997 version applied a time varying 
    parallax to high RV stars, while the van Leeuwan reduction did not. The van Leeuwan reduction
    residuals are relative to the 7 or 9 parameter model (when applicable) while the 1997 version
    is always relative to the 5 parameter solution, even when higher order terms are reported.  

"""
function HipparcosIADLikelihood(;
        hip_id,
        catalog=(datadep"Hipparcos_IAD"),
        renormalize=true,
        attempt_correction=true,
        is_van_leeuwen=true,
        ref_epoch_ra=nothing,
        ref_epoch_dec=nothing,
        variables::Tuple{Priors,Derived}=(@variables begin;end)
    )
    (priors,derived)=variables

    if hip_id ∈ hipparcos_catalog_orbit_parameters_used
        @warn "This object was originally fit using an external orbit solution. The reconstructed IAD may not be correct."
    end

    if isnothing(ref_epoch_ra)
        ref_epoch_ra = meta_gaia_DR3.ref_epoch_mjd
    end
    if isnothing(ref_epoch_dec)
        ref_epoch_dec = meta_gaia_DR3.ref_epoch_mjd
    end

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
        @warn("Only stars with solution types 5, 7, and 9 are currently supported. This sol is $isol_n. If you need another solution type, please open an issue on GitHub and we would be happy to add it.")
    end

    iad_table_rows = NamedTuple[]
    for line in lines[13:end]
        if startswith(line, '#')
            continue
        end
        iorb, epoch, parf, cosϕ, sinϕ, res, sres = split(line)
        push!(iad_table_rows, (;
            iorb=parse(Int, iorb),
            epoch_yrs=parse(Float64, epoch),
            parf=parse(Float64, parf),
            cosϕ=parse(Float64, cosϕ),
            sinϕ=parse(Float64, sinϕ),
            res=parse(Float64, res),
            sres=parse(Float64, sres),
        ))
    end

    iad_table = FlexTable(iad_table_rows)

    # Remove scans that were flagged as rejected in the original reduction, if any
    iad_table.reject = iad_table.sres .<= 0
    if any(iad_table.reject)
        @warn "rejected scans present" count(iad_table.reject)
    end

    ## Nielsen eq 10
    # Un-re-normalize uncertainties
    if renormalize
        D = length(iad_table.sres) - isol_n
        G = f2
        f = (G√(2 / 9D) + 1 - (2 / 9D))^(3 / 2)
        iad_table.sres_renorm = iad_table.sres .* f
        @info "renormalizing hipparcos uncertainties according to Nielsen et al (2020). Pass `renormalize=false` to disable."
    else
        iad_table.sres_renorm = iad_table.sres
        @info "Not renormalizing hipparcos uncertainties according to Nielsen et al (2020). Pass `renormalize=true` to enable."
    end

    if attempt_correction
        iad_table, was_corrected = correct_iad_corruption(iad_table)
        if was_corrected
            @info "Attempted correction of corrupted Hipparcos IAD data using algorithm from appendix A of G. Brandt et al., 2021."
        end
    end

    # transform epoch to MJD
    iad_table.epoch = hipparcos_catalog_epoch_mjd .+ iad_table.epoch_yrs * julian_year


    # Query earth position vs time using SPICE ephemeris (offline)
    earth_pos_vel = geocentre_position_query.(iad_table.epoch)
    table = FlexTable(eachcol(iad_table)..., eachcol(earth_pos_vel)...)



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

    # Not currently used
    dist = MvNormal([α₀,δ₀,hip_sol.plx,μα✱,μδ], [hip_sol.e_ra, hip_sol.e_de, hip_sol.e_plx, hip_sol.e_pmra,  hip_sol.e_pmde,] )

    # For 21 stars that are nearby and have high radial velocity, the Hipparcos team
    # did account for changing parallax. 
    # We have to reproduce that when undoing their model.
    # Note that we can use a different RV in our actual fits if we want to, but we have to
    # use their exact same RVs to recover the scan data correctly.
    # As a result, parallax is a function of time.

    # This is only true of the 1997 reduction. Tests show that van Leeuwan did not
    # use this RV table, despite it's use in previous reductions.
    # For best results, we ignore the RV when reconstructing but 
    # use the RV when fitting--this reduced the chi^2 beyond what was
    # achieved in the official Hipparcos reduction for high proper motion,
    # high RV stars.

    if is_van_leeuwen
        rv_kms = 0.
    else
        I = findfirst(==(hip_id), hipparcos_stellar_rvs_used_for_cat_solutions.hip_id)
        if isnothing(I)
            rv_kms = 0.
        else
            rv_kms = hipparcos_stellar_rvs_used_for_cat_solutions.rv_kms[I]
            @info "This star was flagged by the Hipparcos team as having perspective acceleration effects greater than 0.1 mas/yr. They modelled it using the following barycentric RV. We use the same RV to back out the IAD residuals" rv_kms
        end
    end
    table.rv_kms = fill(rv_kms, length(table.epoch_yrs))

    # # NB: the α✱ indicates the variable are multiplied by cos(δ)

    # ####################
    # # This version of the code uses a rigorous non-linear model of proper motion, under the assumption
    # # that Hipparcos used an iterative solution for high proper motion stars
    o = orbit(
        a=0,M=1,i=0,Ω=0,ω=0,e=0.,
        plx=hip_sol.plx,
        ra=hip_sol.radeg,
        dec=hip_sol.dedeg,
        rv=rv_kms*1e3,
        pmra=hip_sol.pm_ra,
        pmdec=hip_sol.pm_de,
        ref_epoch=hipparcos_catalog_epoch_mjd
    )
    sols = orbitsolve.(o, table.epoch)
    α = map(s->s.compensated.ra2, sols)
    δ = map(s->s.compensated.dec2, sols)
    plx_vs_time = map(s->s.compensated.parallax2, sols)
    table.Δα✱ = @. (α - hip_sol.radeg)*60*60*1000 *cosd.(δ) + plx_vs_time * (
        table.x * sind(α) -
        table.y * cosd(α)
    )
    table.Δδ = @. (δ - hip_sol.dedeg)*60*60*1000  + plx_vs_time * (
        table.x * cosd(α) * sind(δ) +
        table.y * sind(α) * sind(δ) -
        table.z * cosd(δ)
    )

    ##################
    # This version of the code uses a simplifed tangent-plane approximation, under the assumption
    # that Hipparcos used a purely linear model when they fit (or least, to report the IAD residuals against)
    # Calculate sky coordinate Delta at each scan epoch from the catalog position
    # using the eath's motion ("ephemeris"), x,y,z in AU.
    plx0 = hip_sol.plx
    dist0 = 1000/plx0
    δdist_pc_δt_sec = rv_kms / IAU_pc2km / PlanetOrbits.sec2day
    Δdist_pc = δdist_pc_δt_sec .* (table.epoch .- hipparcos_catalog_epoch_mjd)
    dist1 = dist0 .+ Δdist_pc
    table.plx_vs_time = 1000 ./ dist1
    table.Δα✱ = @. table.plx_vs_time * (
        table.x * sind(α₀) -
        table.y * cosd(α₀)
        # table.x * sind(α) -
        # table.y * cosd(α)
    # Both of these calculations are correct, that's how we know the dates are computed correctly:
    ) + (table.epoch - hipparcos_catalog_epoch_mjd)/julian_year * μα✱
    table.Δδ = @. table.plx_vs_time * (
        table.x * cosd(α₀) * sind(δ₀) +
        table.y * sind(α₀) * sind(δ₀) -
        table.z * cosd(δ₀)
        # table.x * cosd(α) * sind(δ) +
        # table.y * sind(α) * sind(δ) -
        # table.z * cosd(δ)
    # Both of these calculations are correct, that's how we know the dates are computed correctly:
    ) + (table.epoch - hipparcos_catalog_epoch_mjd)/julian_year * μδ


    # Nielsen Eq. 3: abcissa point
    table.α✱ₐ = @. table.res * table.cosϕ + table.Δα✱
    table.δₐ = @. table.res * table.sinϕ + table.Δδ

    # These two points (picked ± 1 mas) specify a line--- the Hipparcos `sres`
    # is then then the 
    table.α✱ₘ = collect(eachrow(@. [-1, 1]' * table.sinϕ + table.α✱ₐ))
    table.δₘ = collect(eachrow(@. [1, -1]' * table.cosϕ + table.δₐ))

    table.scanAngle_rad = @. atan(table.sinϕ, table.cosϕ)
    
    table.parallaxFactorAlongScan = @. (
        (table.x * sind(α₀) - table.y * cosd(α₀)) * table.cosϕ +
        (table.x * cosd(α₀) * sind(δ₀) + table.y * sind(α₀) * sind(δ₀) - table.z * cosd(δ₀)) * table.sinϕ
    )

    table = Table(table)

    # Prepare some matrices for linear system solves
    A_prepared_4 = prepare_A_4param(table, ref_epoch_ra, ref_epoch_dec)
    A_prepared_5 = prepare_A_5param(table, ref_epoch_ra, ref_epoch_dec)


    return HipparcosIADLikelihood(hip_sol, table, priors, derived, "Hipparcos IAD", dist, A_prepared_4, A_prepared_5)
end
export HipparcosIADLikelihood


"""
    detect_corruption(iad_table)

Check for corruption in Hipparcos IAD data by looking for duplicate AL errors
at the end of the table. Returns the number of corrupted entries detected.
"""
function detect_corruption(iad_table)
    # Need at least 4 measurements
    length(iad_table.sres_renorm) < 4 && return 0
    
    # Look at last 4 measurements
    last_idx = length(iad_table.sres_renorm)
    last_four_idx = last_idx-3:last_idx
    
    # Check if all from same orbit
    last_four_orbits = iad_table.iorb[last_four_idx]
    if !all(x -> x == last_four_orbits[1], last_four_orbits)
        return 0
    end
    
    # Check for any suspicious repeats in SRES
    last_four_sres = iad_table.sres_renorm[last_four_idx]
    
    # First ≈ last value
    if abs(last_four_sres[1] - last_four_sres[4]) < 0.0001
        # Middle values approximately equal
        if abs(last_four_sres[2] - last_four_sres[3]) < 0.1
            @debug "Detected corruption" first_val=last_four_sres[1] middle_vals=last_four_sres[2:3] last_val=last_four_sres[4]
            return 3
        end
    end
    
    return 0
end

"""
    find_best_correction(iad_table, n_corrupt)

Find the best correction pattern for corrupted IAD data by trying different
combinations of orbit rejections. Returns the indices to remove and the achieved
chi-squared value.
"""
function find_best_correction(iad_table, n_corrupt)
    # This is a translation of htof.parse.find_epochs_to_reject_java from G. Mirek Brandt

    possible_rejects = 1:length(iad_table.epoch_yrs)
    resid_reject_idx = [length(iad_table.epoch_yrs) - i + 1 for i in 1:n_corrupt]  # always reject the repeated observations.
    # need to iterate over popping orbit combinations
    orbits_to_keep = trues(length(iad_table.epoch_yrs))
    residuals_to_keep = trues(length(iad_table.epoch_yrs))
    residuals_to_keep[resid_reject_idx] .= false

    # @show iad_table.res
    # @show (iad_table.sres_renorm .^ 2)
    residual_factors = (iad_table.res ./ (iad_table.sres_renorm .^ 2))[residuals_to_keep]
    # @show residual_factors
    dt = iad_table.epoch_yrs


    # np.array([data.parallax_factors.values, sin_scan, cos_scan, dt * sin_scan, dt * cos_scan]).T
    _orbit_factors = [ iad_table.parf iad_table.cosϕ iad_table.sinϕ dt .* iad_table.cosϕ dt .* iad_table.sinϕ]

    # we should be able to do the orbit reject calculation fairly easily in memory.
    # for 100 choose 3 we have like 250,000 combinations of orbits -- we sghould be able to
    # do those in 10,000 orbit chunks in memory and gain a factor of 10,000 speed up.
    candidate_orbit_rejects = Vector{Int}[]
    candidate_orbit_chisquared_partials = Float64[]
    for orbit_to_reject in Combinatorics.combinations(possible_rejects, n_corrupt)
        orbits_to_keep[orbit_to_reject] .= false
        # now we want to try a variety of deleting orbits and sliding the other orbits
        # upward to fill the vacancy.
        # this pops the orbits out and shifts all the orbits after:
        orbit_factors = @view _orbit_factors[orbits_to_keep,:]
        # this simultaneously deletes one of the residuals, assigns the remaining residuals to the
        # shifted orbits, and calculates the chi2 partials vector per orbit:
        chi2_vector = (2 .* residual_factors .* orbit_factors)'
        # sum the square of the chi2 partials to decide for whether or not it is a stationary point.
        sum_chisquared_partials = sqrt(sum(sum(chi2_vector, dims=2).^2))
        push!(candidate_orbit_rejects, orbit_to_reject)
        push!(candidate_orbit_chisquared_partials, sum_chisquared_partials)
        # reset for the next loop:
        orbits_to_keep .= true
    end
    orbit_reject_idx = candidate_orbit_rejects[argmin(candidate_orbit_chisquared_partials)]

    @debug "Final solution" indices=orbit_reject_idx chi2=minimum(candidate_orbit_chisquared_partials)
    return orbit_reject_idx, minimum(candidate_orbit_chisquared_partials)
end

"""
    correct_iad_corruption(iad_table)

Apply the ad-hoc correction algorithm to potentially corrupted IAD data.
Modifies the input table in-place if corruption is detected and correctable.
Returns the (possibly) corrected table, and a flag indicating if the correction was applied.
"""
function correct_iad_corruption(iad_table)
    n_corrupt = detect_corruption(iad_table)
    
    if n_corrupt == 0
        return iad_table, false
    end
    
    @info "Detected $n_corrupt corrupted entries. Attempting correction following Appendix A of G. M. Brandt et al. (2021)"

    indices_to_remove, chi2 = find_best_correction(iad_table, n_corrupt)
    
    if isempty(indices_to_remove)
        @warn "Could not find valid correction pattern"
        return false
    end
    
    if chi2 > 0.5
        @warn "Correction achieved but chi² ($chi2) is high, treat results with caution"
    end
    
    # Remove corrupted entries
    mask = trues(length(iad_table.epoch_yrs))
    mask[indices_to_remove] .= false
    
    # Create new table with corrected data
    # We remove the rows from "L" that we determined above, while removing `n_corrupt` rows from the end of "R".
    corrected_table = FlexTable(;
        # L:
        iad_table[mask].iorb,
        iad_table[mask].epoch_yrs,
        iad_table[mask].parf,
        iad_table[mask].cosϕ,
        iad_table[mask].sinϕ ,
        # R:
        iad_table[1:end-n_corrupt].res ,
        iad_table[1:end-n_corrupt].sres,
        iad_table[1:end-n_corrupt].reject ,
        iad_table[1:end-n_corrupt].sres_renorm,
    )
    
    @info "Successfully applied correction" removed_indices=indices_to_remove chi2=chi2
    return corrected_table, true
end


##########################
# Computes likelihood
function ln_like(
    hiplike::HipparcosIADLikelihood,
    θ_system,
    θ_obs,
    orbits,
    orbit_solutions,
    sol_start_i
)
    T = _system_number_type(θ_system)
    ll = zero(T)

    hip_model = simulate(hiplike, θ_system, θ_obs, orbits, orbit_solutions, sol_start_i)
    for i in eachindex(hip_model.resid)
        if hiplike.table.reject[i]
            continue
        end
        ll += logpdf(Normal(0, hiplike.table.sres_renorm[i]), hip_model.resid[i])
    end

    return ll
end

function simulate(hiplike::HipparcosIADLikelihood, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)

    T = _system_number_type(θ_system)
    α✱_model_with_perturbation_out = zeros(T, length(hiplike.table.epoch))
    δ_model_with_perturbation_out = zeros(T, length(hiplike.table.epoch))
    resid_out = zeros(T, length(hiplike.table.epoch))

    # All planets inthe system have orbits defined with the same ra, dec, and proper motion,
    # since these are properties of the system.
    orbit = first(orbits)
    if length(orbits) > 1
        for i in eachindex(orbits)[2:end]
            if orbits[i].ra != orbit.ra ||
               orbits[i].dec != orbit.dec ||
               orbits[i].pmra != orbit.rpma ||
               orbits[i].pmdec != orbit.ra
                pmdec
                error("Planet orbits do not have matching ra, dec, pmpra, and pmdec.")
            end
        end
    end

    # Pre-compute perturbations from all planets for all epochs using standardized function
    α✱_perturbations_total = zeros(T, length(hiplike.table.epoch))
    δ_perturbations_total = zeros(T, length(hiplike.table.epoch))
    
    for planet_i in eachindex(orbits)
        planet_mass_msol = θ_system.planets[planet_i].mass * Octofitter.mjup2msol
        fluxratio = hasproperty(θ_obs, :fluxratio) ? θ_obs.fluxratio[planet_i] : zero(T)
        
        # Create temporary arrays for this planet's contribution
        Δα_mas = zeros(T, length(hiplike.table.epoch))
        Δδ_mas = zeros(T, length(hiplike.table.epoch))
        
        _simulate_skypath_perturbations!(
            Δα_mas, Δδ_mas,
            hiplike.table, orbits[planet_i],
            planet_mass_msol, fluxratio,
            orbit_solutions[planet_i], orbit_solutions_i_epoch_start[planet_i], T
        )
        
        # Add this planet's contribution to total perturbations
        α✱_perturbations_total .+= Δα_mas
        δ_perturbations_total .+= Δδ_mas
    end

    for i in eachindex(hiplike.table.epoch)

        orbitsol_hip_epoch = first(orbit_solutions)[i+orbit_solutions_i_epoch_start]

        ##############
        # Non-linear version
        α = orbitsol_hip_epoch.compensated.ra2
        δ = orbitsol_hip_epoch.compensated.dec2
        plx_at_epoch = orbitsol_hip_epoch.compensated.parallax2
        α✱_model = (α - hiplike.hip_sol.radeg).*60*60*1000 *cosd(δ) + plx_at_epoch * (
            hiplike.table.x[i] * sind(α) -
            hiplike.table.y[i] * cosd(α)
        )
        δ_model = (δ - hiplike.hip_sol.dedeg).*60*60*1000  + plx_at_epoch * (
            hiplike.table.x[i] * cosd(α) * sind(δ) +
            hiplike.table.y[i] * sind(α) * sind(δ) -
            hiplike.table.z[i] * cosd(δ)
        )


        # #############
        # # Linear version
        # # TODO: could hoist the trig into the constructor for speed. It doesn't change in this simplified
        # # model.

        # # Simplified tangent-plane model (matches Hipparcos)
        # # The parallax versus time is considered for 21 high RV, nearby stars.
        # # This matches the Hipparcos reduction.
        # # For other stars, rv_kms will be zero and this isn't considered.
        # delta_t_days = orbitsol_hip_epoch.t - hipparcos_catalog_epoch_mjd
        # # Proper motion: each individual year is either 365 or 366 days long. Need to be careful!
        # # the formal definition is to use Julian years, 365.25000000 exactly.
        # delta_time_julian_year = delta_t_days / julian_year
        # plx0 = orbit.plx
        # if orbit.rv != 0
        #     dist0 = 1000 / plx0
        #     δdist_pc_δt_sec = orbit.rv*1000 / IAU_pc2km / PlanetOrbits.sec2day
        #     Δdist_pc = δdist_pc_δt_sec * (hiplike.table.epoch[i] - hipparcos_catalog_epoch_mjd)
        #     dist1 = dist0 + Δdist_pc
        #     plx_at_epoch = 1000 / dist1
        # else
        #     plx_at_epoch = plx0
        # end
        # # α = orbitsol_hip_epoch.compensated.ra2
        # # δ = orbitsol_hip_epoch.compensated.dec2
        # α✱_model = (θ_system.ra - hiplike.hip_sol.radeg) * 60 * 60 * 1000 * cosd(hiplike.hip_sol.dedeg) + plx_at_epoch * (
        #                hiplike.table.x[i] * sind(hiplike.hip_sol.radeg) -
        #                hiplike.table.y[i] * cosd(hiplike.hip_sol.radeg)
        #     # hiplike.table.x[i] * sind(α) -
        #     # hiplike.table.y[i] * cosd(α)
        # ) + delta_time_julian_year * orbit.pmra
        # δ_model = (θ_system.dec - hiplike.hip_sol.dedeg) * 60 * 60 * 1000 + plx_at_epoch * (
        #               hiplike.table.x[i] * cosd(hiplike.hip_sol.radeg) * sind(hiplike.hip_sol.dedeg) +
        #               hiplike.table.y[i] * sind(hiplike.hip_sol.radeg) * sind(hiplike.hip_sol.dedeg) -
        #               hiplike.table.z[i] * cosd(hiplike.hip_sol.dedeg)
        #             # hiplike.table.x[i] * cosd(α) * sind(δ) +
        #             # hiplike.table.y[i] * sind(α) * sind(δ) -
        #             # hiplike.table.z[i] * cosd(δ)
        # ) + delta_time_julian_year * orbit.pmdec


        # Use pre-computed perturbations for this epoch
        α✱_perturbation = α✱_perturbations_total[i]
        δ_perturbation = δ_perturbations_total[i]
        # coordinates in degrees of RA*cos(dec) and Dec
        α✱_model_with_perturbation = α✱_model + α✱_perturbation #/ cosd(hiplike.hip_sol.dedeg)
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
        δ_model_with_perturbation_out[i] = δ_model_with_perturbation
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
