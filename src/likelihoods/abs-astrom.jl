

"""
    G23HObs(; gaia_id=nothing, hip_id=nothing, include_rv=true, ...)
    GaiaHipparcosUEVAJointObs(; gaia_id=nothing, hip_id=nothing, include_rv=true, ...)

A likelihood for joint Gaia-Hipparcos astrometry using the G23H catalog, including
proper motion accelerations and unit weight error variance analysis (UEVA).

# Arguments
- `gaia_id`: Gaia DR3 source ID (provide either this or `hip_id`)
- `hip_id`: Hipparcos ID (provide either this or `gaia_id`)
- `catalog`: Path to G23H catalog file, or a loaded DataFrame/Table.
  Defaults to the automatically downloaded G23H catalog via DataDeps.
- `include_rv`: Whether to include Gaia RV variability constraints (default: true)
- `ueva_mode`: Either `:RUWE` (default) or `:EAN` for astrometric excess noise modeling
- `freeze_epochs`: If true, fix Gaia observation epochs for faster sampling (default: false)
- `variables`: Optional custom priors (defaults are set from catalog values)

# Examples
```julia
# Using Gaia DR3 source ID
absastrom = G23HObs(gaia_id=756291174721509376)

# Using Hipparcos ID (automatically resolved to Gaia ID)
absastrom = G23HObs(hip_id=21547)
```

The G23H catalog (~14 GB) is automatically downloaded on first use.

# Variable Priors
Default priors are set automatically from the catalog. Custom priors can be specified:
```julia
variables=@variables begin
    fluxratio ~ Product([LogUniform(1e-6, 1e-1) for _ in 1:N_planets])  # For luminous companions
    σ_att ~ LogUniform(0.01, 1.0)     # Attitude error (mas)
    σ_AL ~ LogUniform(0.01, 1.0)      # Along-scan error (mas)
    σ_calib ~ LogUniform(0.01, 1.0)   # Calibration error (mas)
end
```
"""
struct GaiaHipparcosUEVAJointObs{TTable,TTableH,TTableG,TCat,THip} <: AbstractObs
    table::TTable
    priors::Priors
    derived::Derived
    hip_table::TTableH
    gaia_table::TTableG
    catalog::TCat
    hip_sol::THip
    A_prepared_5_hip::Matrix{Float64}
    A_prepared_5_dr2::Matrix{Float64}
    A_prepared_5_dr3::Matrix{Float64}
    include_iad::Bool
    ueva_mode::Symbol
end


function likelihoodname(like::GaiaHipparcosUEVAJointObs)
    return "GaiaHipparcosUEVA"
end

function GaiaHipparcosUEVAJointObs(;
        gaia_id=nothing,
        hip_id=nothing,
        scanlaw_table=nothing,
        catalog=joinpath(datadep"G23H_Catalog", "G23H-v1.0.feather"),
        variables::Union{Nothing,Tuple{Priors,Derived}}=nothing,
        include_rv=true,
        ueva_mode::Symbol=:RUWE,
        freeze_epochs=false
    )
    include_iad=false

    # Validate that exactly one of gaia_id or hip_id is provided
    if isnothing(gaia_id) && isnothing(hip_id)
        error("Either gaia_id or hip_id must be specified")
    end
    if !isnothing(gaia_id) && !isnothing(hip_id)
        error("Specify either gaia_id or hip_id, not both")
    end

    # allow passing in table directly
    if Tables.istable(catalog)
        # If hip_id provided, look up gaia_id
        if !isnothing(hip_id)
            hip_matches = findall(==(hip_id), catalog.hip_id)
            if isempty(hip_matches)
                error("The requested Hipparcos ID $hip_id was not found in the catalog.")
            end
            gaia_id = catalog.gaia_source_id[hip_matches[1]]
            @info "Resolved HIP $hip_id to Gaia DR3 source ID $gaia_id"
        end
        idx = findfirst(==(gaia_id), catalog.gaia_source_id)
        catalog = NamedTuple(catalog[idx,:])
    else
        # Load the catalog row for this system
        catalog = FITS(catalog, "r") do fits
            t = Table(fits[2])
            # If hip_id provided, look up gaia_id
            if !isnothing(hip_id)
                hip_matches = findall(==(hip_id), t.hip_id)
                if isempty(hip_matches)
                    error("The requested Hipparcos ID $hip_id was not found in the catalog file $catalog.")
                end
                gaia_id = t.gaia_source_id[hip_matches[1]]
                @info "Resolved HIP $hip_id to Gaia DR3 source ID $gaia_id"
            end
            idx = findfirst(==(gaia_id), t.gaia_source_id)
            if isnothing(idx)
                error("The requested gaia source ID $gaia_id was not found in the catalog file $catalog.")
            end
            return NamedTuple(t[idx])
        end
    end

    # Convert measurement epochs to MJD.
    # Careful: these are Julian years, not decimal years (T. Brant., private communications)
    J2000_mjd = 51544.5 # year J2000 in MJD
    catalog = (;
        catalog...,
        epoch_ra_hip_mjd=(catalog.epoch_ra_hip - 2000) * julian_year + J2000_mjd,
        epoch_dec_hip_mjd=(catalog.epoch_dec_hip - 2000) * julian_year + J2000_mjd,
        epoch_ra_dr2_mjd=(catalog.epoch_ra_dr2 - 2000) * julian_year + J2000_mjd,
        epoch_dec_dr2_mjd=(catalog.epoch_dec_dr2 - 2000) * julian_year + J2000_mjd,
        epoch_ra_dr3_mjd=(catalog.epoch_ra_dr3 - 2000) * julian_year + J2000_mjd,
        epoch_dec_dr3_mjd=(catalog.epoch_dec_dr3 - 2000) * julian_year + J2000_mjd,
    )


    if !hasproperty(catalog, :astrometric_chi2_al_dr3) || !hasproperty(catalog, :rv_nb_transits) 
        @warn "Column missing from catalog, querying Gaia DR3 TAP server (or using cached value)"

        dr3 = Octofitter._query_gaia_dr3(;gaia_id)
        catalog = (;
            catalog...,
            astrometric_chi2_al_dr3=dr3.astrometric_chi2_al,
            parallax_error=dr3.parallax_error,
            rv_nb_transits=dr3.rv_nb_transits,
            radial_velocity_error=dr3.radial_velocity_error,
        )
    end

    if isnan(catalog.hip_id)
        @warn "No Hipparcos data found; will skip HGCA and IAD modelling"
        hip_like = nothing
        hip_sol = nothing
        dist_hip = nothing
        dist_hg  = nothing
        hip_table = Table(
            @NamedTuple{iorb::Int64, epoch_yrs::Float64, parf::Float64, cosϕ::Float64, sinϕ::Float64, res::Float64, sres::Float64, reject::Bool, sres_renorm::Float64, epoch::Float64, x::Float64, y::Float64, z::Float64, vx::Float64, vy::Float64, vz::Float64, rv_kms::Float64, Δα✱::Float64, Δδ::Float64, plx_vs_time::Float64, α✱ₐ::Float64, δₐ::Float64, α✱ₘ::SubArray{Float64, 1, Matrix{Float64}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}, δₘ::SubArray{Float64, 1, Matrix{Float64}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}, scanAngle_rad::Float64, parallaxFactorAlongScan::Float64}[]
        )
        A_prepared_5_hip = fill(0.0, 0,0)
    else
        # Load the Hipparcos IAD data for epochs and scan angles
        hip_like = HipparcosIADLikelihood(;
            catalog.hip_id,
            ref_epoch_ra=catalog.epoch_ra_hip_mjd,
            ref_epoch_dec=catalog.epoch_dec_hip_mjd,
            variables=(Priors(),Derived())
        )
        A_prepared_5_hip = hip_like.A_prepared_5
        hip_table = hip_like.table
        hip_sol = hip_like.hip_sol

        # Following "Statistical properties of Hipparcos 2, caveats on its use, and a recalibration of the intermediate astrometric data"
        # by G Mirek Brandt,  Daniel Michalik,  Timothy D Brandt; 
        # we add 0.140 mas to the residuals and 2.25 mas additional dispersion to the unceratinties
        # this mitigates overfitting.
        hip_table.res .+= 0.140
        hip_table.sres_renorm .= hypot.(hip_table.sres_renorm, 2.25)

        # Precompute MvNormal distributions for correlation between ra and dec
        # Hipparcos epoch
        c = catalog.pmra_pmdec_hip[1] * catalog.pmra_hip_error[1] * catalog.pmdec_hip_error[1]
        dist_hip = MvNormal(
            @SVector([catalog.pmra_hip, catalog.pmdec_hip]),
            @SArray[
                catalog.pmra_hip_error[1]^2 c
                c catalog.pmdec_hip_error[1]^2
            ]
        )

        # Hipparcos - GAIA epoch
        c = catalog.pmra_pmdec_hg[1] * catalog.pmra_hg_error[1] * catalog.pmdec_hg_error[1]
        dist_hg = MvNormal(
            @SVector([catalog.pmra_hg, catalog.pmdec_hg]),
            @SArray [
                catalog.pmra_hg_error[1]^2 c
                c catalog.pmdec_hg_error[1]^2
            ]
        )
    end

    # Load the Gaia scanlaw etc
    # gaia_like = GaiaCatalogFitLikelihood(; gaia_id_dr3=gaia_id)

    # Besides epoch and catalog, I'm not sure we will really use this data table
    # except maybe for plotting
    # table = Table(;
    #     epoch=[hipparcos_catalog_epoch_mjd, meta_gaia_DR2.ref_epoch_mjd],
    #     catalog=[:hipparcos, :gaia],
    #     ra=[hip_like.hip_sol.radeg, gaia_like.gaia_sol.ra],
    #     dec=[hip_like.hip_sol.dedeg, gaia_like.gaia_sol.dec],
    #     plx=[hip_like.hip_sol.plx, gaia_like.gaia_sol.parallax],
    #     pmra=[hip_like.hip_sol.pm_ra, gaia_like.gaia_sol.pmra],
    #     pmdec=[hip_like.hip_sol.pm_de, gaia_like.gaia_sol.pmdec],
    # )

   

    # GAIA DR2 epoch
    c = catalog.pmra_pmdec_dr2[1] * catalog.pmra_dr2_error[1] * catalog.pmdec_dr2_error[1]
    dist_dr2 = MvNormal(
        @SVector([catalog.pmra_dr2, catalog.pmdec_dr2]),
        @SArray [
            catalog.pmra_dr2_error[1]^2 c
            c catalog.pmdec_dr2_error[1]^2
        ]
    )


    # GAIA DR3-DR2 epoch
    c = catalog.pmra_pmdec_dr32[1] * catalog.pmra_dr32_error[1] * catalog.pmdec_dr32_error[1]
    dist_dr32 = MvNormal(
        @SVector([catalog.pmra_dr32, catalog.pmdec_dr32]),
        @SArray [
            catalog.pmra_dr32_error[1]^2 c
            c catalog.pmdec_dr32_error[1]^2
        ]
    )

    # GAIA DR3 epoch
    c = catalog.pmra_pmdec_dr3[1] * catalog.pmra_dr3_error[1] * catalog.pmdec_dr3_error[1]
    dist_dr3 = MvNormal(
        @SVector([catalog.pmra_dr3, catalog.pmdec_dr3]),
        @SArray [
            catalog.pmra_dr3_error[1]^2 c
            c catalog.pmdec_dr3_error[1]^2
        ]
    )

    catalog = (; catalog..., dist_hip, dist_hg, dist_dr2, dist_dr32, dist_dr3)

    if isnothing(scanlaw_table)
        # @warn "No scan law table provided. We will fetch an approximate solution from the GOST webservice, but for best results please use the `scanninglaw` python package, installable via pip, to query the RA and Dec of this target and supply it as `scanlaw_table`. Run: `import astropy.coordinates, scanninglaw, pandas; o = astropy.coordinates.SkyCoord(158.30707896392835, 40.42555422701387,unit='deg');t = scanninglaw.times.Times(version='dr3_nominal'); t.query(o,return_angles=True)`"
        # Get predicted GAIA scan epochs and angles
        forecast_table = FlexTable(GOST_forecast(catalog.ra, catalog.dec))
        forecast_table.epoch = jd2mjd.(forecast_table.ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_)
        forecast_table.scanAngle_rad = forecast_table.scanAngle_rad_
    else
        @info "Scanlaw table from the `scanninglaw` python package was provided, will not use GOST."
        forecast_table = FlexTable(scanlaw_table)
        forecast_table.epoch = tcb_at_gaia_2mjd.(forecast_table.times)
        forecast_table.scanAngle_rad = deg2rad.(forecast_table.angles)

        earth_pos_vel = FlexTable(geocentre_position_query.(forecast_table.epoch))

        f = @. earth_pos_vel.x * sind(catalog.ra)-earth_pos_vel.y*cosd(catalog.ra)
        g = @. earth_pos_vel.x * cosd(catalog.ra) * sind(catalog.dec) + 
            earth_pos_vel.y * sind(catalog.ra) * sind(catalog.dec) -
            earth_pos_vel.z * cosd(catalog.dec)
        forecast_table.parallaxFactorAlongScan = @. f*sin(forecast_table.scanAngle_rad) + g*cos(forecast_table.scanAngle_rad)

    end
    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    forecast_table.cosϕ = cos.(π/2 .+ forecast_table.scanAngle_rad)
    forecast_table.sinϕ = sin.(π/2 .+ forecast_table.scanAngle_rad)

    # Get the Earth's position at those epochs
    earth_pos_vel = geocentre_position_query.(forecast_table.epoch)

    # merge the Gaia scan prediction and geocentre position results into one table
    gaia_table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)



    # Now remove any known gaps -- data sourced from HTOF.py; authors G.M. Brandt et al
    gaps_dr2 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiadr2_08252020.csv"), FlexTable)
    gaps_edr23 = CSV.read(joinpath(@__DIR__, "astrometric_gaps_gaiaedr3_12232020.csv"), FlexTable)
    gaps = Table(
        start_mjd=obmt2mjd.(vcat(gaps_dr2.start,gaps_edr23.start)),
        stop_mjd=obmt2mjd.(vcat(gaps_dr2.end,gaps_edr23.end)),
        note=[gaps_dr2.comment; gaps_edr23.description]
    )
    gaia_table = filter(eachrow(gaia_table)) do row
        row = row[]
        for gap in eachrow(gaps)
            gap = gap[]
            if gap.start_mjd <= row.epoch <= gap.stop_mjd
                @info "Detected known gap in Gaia scans; skipping." window=row.epoch note=gap.note
                return false
            end
        end
        return true
    end
    gaia_table = Table(map(dat->dat[], gaia_table))


    # Determine fraction of epochs in DR2 that overlap with DR3

    istart_dr2 = findfirst(>=(meta_gaia_DR2.start_mjd), vec(gaia_table.epoch))
    iend_dr2 = findlast(<=(meta_gaia_DR2.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr2)
        istart_dr2 = 1
    end
    if isnothing(iend_dr2)
        iend_dr2 = length(gaia_table.epoch)
    end

    istart_dr3 = findfirst(>=(meta_gaia_DR3.start_mjd), vec(gaia_table.epoch))
    iend_dr3 = findlast(<=(meta_gaia_DR3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(gaia_table.epoch)
    end

    min_epoch = +Inf
    max_epoch = -Inf
    min_epoch = min(min_epoch,meta_gaia_DR2.start_mjd)
    max_epoch = max(max_epoch,meta_gaia_DR2.stop_mjd)
    min_epoch = min(min_epoch,meta_gaia_DR3.start_mjd)
    max_epoch = max(max_epoch,meta_gaia_DR3.stop_mjd)
    gaia_table = Table(gaia_table[min_epoch .<= gaia_table.epoch .<= max_epoch,:])

    # DR2
    A_prepared_5_dr2 = prepare_A_5param(gaia_table, catalog.epoch_ra_dr2_mjd,  catalog.epoch_dec_dr2_mjd)
    
    # DR3
    A_prepared_5_dr3 = prepare_A_5param(gaia_table, catalog.epoch_ra_dr3_mjd,  catalog.epoch_dec_dr3_mjd)

    # This table serves only for plotting and keeping track of subsetting for cross-validation.
    # The actual likelihood calculations happen against the prepared linear system matrices above.
    table = Table(
        epoch=[
            isnothing(hip_like) ? NaN : mean(hip_like.table.epoch),
            isnothing(hip_like) ? NaN : years2mjd(catalog.epoch_ra_hip),
            isnothing(hip_like) ? NaN : years2mjd(catalog.epoch_dec_hip),
            isnothing(hip_like) ? NaN : years2mjd(catalog.epoch_ra_hg),
            isnothing(hip_like) ? NaN : years2mjd(catalog.epoch_dec_hg),
            years2mjd(catalog.epoch_ra_dr2),
            years2mjd(catalog.epoch_dec_dr2),
            years2mjd(catalog.epoch_ra_dr32),
            years2mjd(catalog.epoch_dec_dr32),
            years2mjd(catalog.epoch_ra_dr3),
            years2mjd(catalog.epoch_dec_dr3),
            years2mjd((catalog.epoch_dec_dr3+catalog.epoch_dec_dr2)/2),
        ],
        start_epoch=[
            isnothing(hip_like) ? 0 : minimum(hip_table.epoch),
            isnothing(hip_like) ? 0 : minimum(hip_table.epoch),
            isnothing(hip_like) ? 0 : minimum(hip_table.epoch),
            isnothing(hip_like) ? 0 : years2mjd(catalog.epoch_ra_hip),
            isnothing(hip_like) ? 0 : years2mjd(catalog.epoch_dec_hip),
            first(gaia_table.epoch[gaia_table.epoch.>=(meta_gaia_DR2.start_mjd)]),
            first(gaia_table.epoch[gaia_table.epoch.>=(meta_gaia_DR2.start_mjd)]),
            years2mjd(catalog.epoch_ra_dr2),
            years2mjd(catalog.epoch_dec_dr2),
            first(gaia_table.epoch[gaia_table.epoch.>=(meta_gaia_DR3.start_mjd)]),
            first(gaia_table.epoch[gaia_table.epoch.>=(meta_gaia_DR3.start_mjd)]),
            first(gaia_table.epoch[gaia_table.epoch.>=(meta_gaia_DR3.start_mjd)]),
        ],
        stop_epoch=[
            isnothing(hip_like) ? 0 : maximum(hip_table.epoch),
            isnothing(hip_like) ? 0 : maximum(hip_table.epoch),
            isnothing(hip_like) ? 0 : maximum(hip_table.epoch),
            isnothing(hip_like) ? 0 : years2mjd(catalog.epoch_ra_dr3),
            isnothing(hip_like) ? 0 : years2mjd(catalog.epoch_dec_dr3),
            last(gaia_table.epoch[gaia_table.epoch.<=(meta_gaia_DR2.stop_mjd)]),
            last(gaia_table.epoch[gaia_table.epoch.<=(meta_gaia_DR2.stop_mjd)]),
            years2mjd(catalog.epoch_ra_dr3),
            years2mjd(catalog.epoch_dec_dr3),
            last(gaia_table.epoch[gaia_table.epoch.<=(meta_gaia_DR3.stop_mjd)]),
            last(gaia_table.epoch[gaia_table.epoch.<=(meta_gaia_DR3.stop_mjd)]),
            last(gaia_table.epoch[gaia_table.epoch.<=(meta_gaia_DR3.stop_mjd)]),
        ],
        pm = [
            NaN,
            catalog.pmra_hip,
            catalog.pmdec_hip,
            catalog.pmra_hg,
            catalog.pmdec_hg, 
            catalog.pmra_dr2,
            catalog.pmdec_dr2,
            catalog.pmra_dr32,
            catalog.pmdec_dr32,
            catalog.pmra_dr3,
            catalog.pmdec_dr3,
            ueva_mode == :RUWE ? catalog.ruwe_dr3 : catalog.astrometric_excess_noise_dr3
        ],
        σ_pm = [
            NaN,
            catalog.pmra_hip_error,
            catalog.pmdec_hip_error,
            catalog.pmra_hg_error,
            catalog.pmdec_hg_error,
            catalog.pmra_dr2_error,
            catalog.pmdec_dr2_error,
            catalog.pmra_dr32_error,
            catalog.pmdec_dr32_error,
            catalog.pmra_dr3_error,
            catalog.pmdec_dr3_error,
            NaN
        ],
        kind=[
            :iad_hip,
            :ra_hip,
            :dec_hip,
            :ra_hg,
            :dec_hg,
            :ra_dr2,
            :dec_dr2,
            :ra_dr32,
            :dec_dr32,
            :ra_dr3,
            :dec_dr3,
            :ueva_dr3,
        ],

    )


    has_rv = include_rv && (
        hasproperty(catalog, :rv_ln_uncert_dr3) && !ismissing(catalog.rv_ln_uncert_dr3) && !ismissing(catalog.rv_ln_uncert_err_dr3) &&
        isfinite(catalog.rv_ln_uncert_dr3) && isfinite(catalog.rv_ln_uncert_err_dr3)
    )
    if has_rv
        push!(table.epoch, mean(gaia_table.epoch))  # or use RV-specific epoch if available
        push!(table.start_epoch, first(gaia_table.epoch))
        push!(table.stop_epoch, last(gaia_table.epoch))
        push!(table.pm, NaN)  # RV doesn't have a "pm" equivalent
        push!(table.σ_pm, NaN)
        push!(table.kind, :rv_dr3)
    end
    if isempty(hip_table)
        splice!(table, 1:5)
    end



    if isnothing(variables)

        len_epochs = length(gaia_table.epoch)
        astrometric_matched_transits_dr3 = catalog.astrometric_matched_transits_dr3
        missed_transits = Int(len_epochs - astrometric_matched_transits_dr3)
        dec = catalog.dec
        ra = catalog.ra
        if missed_transits < 0
            @warn "Transits missing from GOST (more matched transits than available options from GOST)"
            missed_transits = 0
        end
        @info "Count of missed or rejected transits:"  dr3=missed_transits

        variables = @variables begin
            σ_AL ~ truncated(Normal(catalog.sig_AL, catalog.sig_AL_sigma), lower=eps(), upper=10.0)
            σ_att ~ truncated(Normal(catalog.sig_att_radec, catalog.sig_att_radec_sigma), lower=eps(), upper=10.0)
            σ_calib ~ truncated(Normal(catalog.sig_cal, catalog.sig_cal_sigma), lower=eps(), upper=10.0)
            fluxratio = hasproperty(sys, :fluxratio) ? sys.fluxratio : 0.0
        end

        if(len_epochs) < astrometric_matched_transits_dr3
            @warn "Fewer epochs in GOST forecast than `astrometric_matched_transits` reported by Gaia. Results by be innaccurate."
            variables = vcat(variables, @variables begin
                transits = $(1:len_epochs)
            end)
        else
            # This is an optional approximation that can massievly speed up sampling -- 
            # sample the epochs randomly once and fix them for all remaining sampling
            if freeze_epochs
                transit_priorities = (randn(len_epochs)...,)
                transits = partialsortperm(SVector(transit_priorities), 1:astrometric_matched_transits_dr3, rev=true)
                variables = vcat(variables, @variables begin
                    transits = $transits
                end)
            else
                # Full model:
                # include the epochs of the Gaia observations as variables 
                # Our goal is to sample the *indices* of the subset of possible observation epochs
                # reported by GOST. We assume that the same subset of epochs used by DR2 are also used
                # in DR3. We assume that the RV epochs used by DR3 are a subset of the astrometry epochs
                # ie any rejected or skipped astrometry epochs are also skipped for RV. This is just to 
                # make the problem tractable.
                # 
                # We sample the discrete set of epochs used by a having a full set of continuous variables
                # corresponding to each possible observing epoch. At any given point, the epochs used in the model
                # are the values with the highest values.
                variables = vcat(variables, @variables begin
                    transit_priorities ~ MvNormal(zeros(len_epochs), I)
                    transits = partialsortperm(SVector(transit_priorities), 1:$astrometric_matched_transits_dr3, rev=true)
                end)
            end
        end

        if !isnothing(hip_like)
            variables_iad = @variables begin
                hip_iad_jitter ~ LogUniform(0.001, 100)
                iad_Δra     ~ Uniform(-1000, 1000)
                iad_Δdec    ~ Uniform(-1000, 1000)
                iad_Δplx    ~ Uniform(-10, 10)
                iad_Δpmra   ~ Uniform(-1000, 1000) 
                iad_Δpmdec  ~ Uniform(-1000, 1000)
                iad_pmra = $(hip_sol.pm_ra) + iad_Δpmra
                iad_pmdec = $(hip_sol.pm_de) + iad_Δpmdec
            end
            variables = vcat(variables, variables_iad)
        end

        if has_rv

            len_epochs = length(gaia_table.epoch)
            variables_rv = @variables begin
                σ_rv_per_transit ~ truncated(Normal(exp(catalog.rv_ln_uncert_dr3), exp(catalog.rv_ln_uncert_err_dr3)), lower=eps())
            end
            variables = vcat(variables, variables_rv)

            transits_rv = Int(len_epochs - catalog.rv_nb_transits)
            @info "Count of RV transits:"  transits_rv total_transits=len_epochs
           
            if transits_rv > 0
                missed_vars = @variables begin
                    transits_rv = partialsortperm(SVector(transit_priorities), 1:$transits_rv, rev=true)
                end
                variables = vcat(variables, missed_vars)
            end

        end

        @info "Added the following observation variables:"
        display(variables[1])
        display(variables[2])
    end
    (priors,derived)=variables


    return GaiaHipparcosUEVAJointObs{
        typeof(table),
        typeof(hip_table),
        typeof(gaia_table),
        typeof(catalog),
        typeof(hip_sol),
    }(
        table,
        priors,
        derived,
        hip_table,
        gaia_table,
        catalog,
        hip_sol,
        A_prepared_5_hip,
        A_prepared_5_dr2,
        A_prepared_5_dr3,
        include_iad,
        ueva_mode,
    )

end

function Octofitter.likeobj_from_epoch_subset(like::GaiaHipparcosUEVAJointObs, obs_inds)
    (;  table,
        priors,
        derived,
        hip_table,
        gaia_table,
        catalog,
        hip_sol,
        A_prepared_5_hip,
        A_prepared_5_dr2,
        A_prepared_5_dr3,
        include_iad,
        ueva_mode ) = like

    table = table[obs_inds,:]
    if  (
            :iad_hip ∉ like.table.kind &&
            :ra_hip ∉ like.table.kind &&
            :dec_hip ∉ like.table.kind &&
            :ra_hg ∉ like.table.kind &&
            :dec_hg ∉ like.table.kind
        )
        catalog = (;catalog..., dist_hip = nothing, dist_hg=nothing)
    end
    return GaiaHipparcosUEVAJointObs{
        typeof(table),
        typeof(hip_table),
        typeof(gaia_table),
        typeof(catalog),
        typeof(hip_sol),
    }(
        table,
        priors,
        derived,
        hip_table,
        gaia_table,
        catalog,
        hip_sol,
        A_prepared_5_hip,
        A_prepared_5_dr2,
        A_prepared_5_dr3,
        include_iad,
        ueva_mode,
    )
end

function ln_like(like::GaiaHipparcosUEVAJointObs, ctx::SystemObservationContext)
    (; θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = ctx 

    T = _system_number_type(θ_system)
    ll = zero(T)

    # TODO: optimize this, we only need to grab the epochs here -- it'll be faster
    if hasproperty(θ_obs, :transits)
        (;transits) = θ_obs 
        if eltype(transits) <: AbstractFloat
            transits = Int.(transits)
        end
        # The list of missed transits must be unique
        if length(unique(transits)) < length(transits)
            return nothing
        end
        ii = transits
        gaia_table = like.gaia_table[ii,:]
    else
        gaia_table = like.gaia_table
    end
    
    istart_dr2 = findfirst(>=(meta_gaia_DR2.start_mjd), vec(gaia_table.epoch))
    iend_dr2 = findlast(<=(meta_gaia_DR2.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr2)
        istart_dr2 = 1
    end
    if isnothing(iend_dr2)
        iend_dr2 = length(gaia_table.epoch)
    end

    istart_dr3 = findfirst(>=(meta_gaia_DR3.start_mjd), vec(gaia_table.epoch))
    iend_dr3 = findlast(<=(meta_gaia_DR3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(gaia_table.epoch)
    end

    @no_escape begin


        iad_resid  = @alloc(T, size(like.hip_table,1)); fill!(iad_resid, 0)
        Δα_mas_hip = @alloc(T, size(like.hip_table,1)); fill!(Δα_mas_hip, 0)
        Δδ_mas_hip = @alloc(T, size(like.hip_table,1)); fill!(Δδ_mas_hip, 0)
        Δα_mas_dr2 = @alloc(T, iend_dr2-istart_dr2+1); fill!(Δα_mas_dr2, 0)
        Δδ_mas_dr2 = @alloc(T, iend_dr2-istart_dr2+1); fill!(Δδ_mas_dr2, 0)
        Δα_mas_dr3 = @alloc(T, iend_dr3-istart_dr3+1); fill!(Δα_mas_dr3, 0)
        Δδ_mas_dr3 = @alloc(T, iend_dr3-istart_dr3+1); fill!(Δδ_mas_dr3, 0)
        buffers = (;iad_resid, Δα_mas_hip, Δδ_mas_hip, Δα_mas_dr2, Δδ_mas_dr2, Δα_mas_dr3, Δδ_mas_dr3)

        sim = simulate!(buffers, like, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

        if isnothing(sim)
            ll = convert(T,-Inf)
        else

            (; μ_h, μ_hg, μ_dr2, μ_dr32, μ_dr3, UEVA_model, UEVA_unc, μ_1_3, n_dr3, n_dr2) = sim       
            (;deflation_factor_dr3) = sim
            # Check if we have absolute orbits
            absolute_orbits = false
            for orbit in orbits
                absolute_orbits |= orbit isa AbsoluteVisual
            end

            # Get distribution objects from catalog
            dist_hip = like.catalog.dist_hip
            dist_hg = like.catalog.dist_hg
            dist_dr2 = like.catalog.dist_dr2
            dist_dr32 = like.catalog.dist_dr32
            dist_dr3 = like.catalog.dist_dr3

            # If the user is fitting jitters, we have to regenerate these distributions 
            if !isnothing(dist_hip) && (hasproperty(θ_obs, :σ_hip_pmra) || hasproperty(θ_obs, :σ_hip_pmdec))
                c = like.catalog.pmra_pmdec_hip[1] * like.catalog.pmra_hip_error[1] * like.catalog.pmdec_hip_error[1]
                pmra_var = like.catalog.pmra_hip_error[1]^2 + (hasproperty(θ_obs, :σ_hip_pmra) ? θ_obs.σ_hip_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_hip_error[1]^2 + (hasproperty(θ_obs, :σ_hip_pmdec) ? θ_obs.σ_hip_pmdec : zero(T))^2
                dist_hip = MvNormal(
                    @SVector([like.catalog.pmra_hip, like.catalog.pmdec_hip]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end

            if !isnothing(dist_hg) && (hasproperty(θ_obs, :σ_hg_pmra) || hasproperty(θ_obs, :σ_hg_pmdec))
                c = like.catalog.pmra_pmdec_hg[1] * like.catalog.pmra_hg_error[1] * like.catalog.pmdec_hg_error[1]
                pmra_var = like.catalog.pmra_hg_error[1]^2 + (hasproperty(θ_obs, :σ_hg_pmra) ? θ_obs.σ_hg_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_hg_error[1]^2 + (hasproperty(θ_obs, :σ_hg_pmdec) ? θ_obs.σ_hg_pmdec : zero(T))^2
                dist_hg = MvNormal(
                    @SVector([like.catalog.pmra_hg, like.catalog.pmdec_hg]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end

            if !isnothing(dist_dr2) && (hasproperty(θ_obs, :σ_dr2_pmra) || hasproperty(θ_obs, :σ_dr2_pmdec))
                c = like.catalog.pmra_pmdec_dr2[1] * like.catalog.pmra_dr2_error[1] * like.catalog.pmdec_dr2_error[1]
                pmra_var = like.catalog.pmra_dr2_error[1]^2 + (hasproperty(θ_obs, :σ_dr2_pmra) ? θ_obs.σ_dr2_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_dr2_error[1]^2 + (hasproperty(θ_obs, :σ_dr2_pmdec) ? θ_obs.σ_dr2_pmdec : zero(T))^2
                dist_dr2 = MvNormal(
                    @SVector([like.catalog.pmra_dr2, like.catalog.pmdec_dr2]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end

            if !isnothing(dist_dr32) && (hasproperty(θ_obs, :σ_dr32_pmra) || hasproperty(θ_obs, :σ_dr32_pmdec))
                c = like.catalog.pmra_pmdec_dr32[1] * like.catalog.pmra_dr32_error[1] * like.catalog.pmdec_dr32_error[1]
                pmra_var = like.catalog.pmra_dr32_error[1]^2 + (hasproperty(θ_obs, :σ_dr32_pmra) ? θ_obs.σ_dr32_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_dr32_error[1]^2 + (hasproperty(θ_obs, :σ_dr32_pmdec) ? θ_obs.σ_dr32_pmdec : zero(T))^2
                dist_dr32 = MvNormal(
                    @SVector([like.catalog.pmra_dr32, like.catalog.pmdec_dr32]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end

            if !isnothing(dist_dr3) && (hasproperty(θ_obs, :σ_dr3_pmra) || hasproperty(θ_obs, :σ_dr3_pmdec))
                c = like.catalog.pmra_pmdec_dr3[1] * like.catalog.pmra_dr3_error[1] * like.catalog.pmdec_dr3_error[1]
                pmra_var = like.catalog.pmra_dr3_error[1]^2 + (hasproperty(θ_obs, :σ_dr3_pmra) ? θ_obs.σ_dr3_pmra : zero(T))^2
                pmdec_var = like.catalog.pmdec_dr3_error[1]^2 + (hasproperty(θ_obs, :σ_dr3_pmdec) ? θ_obs.σ_dr3_pmdec : zero(T))^2
                dist_dr3 = MvNormal(
                    @SVector([like.catalog.pmra_dr2, like.catalog.pmdec_dr2]),
                    @SArray[
                        pmra_var c
                        c pmdec_var
                    ]
                )
            end


            # Apply nonlinear correction for absolute orbits
            if absolute_orbits && !isnothing(dist_hip)
                # Add nonlinear corrections to model values
                μ_hg += @SVector [
                    like.catalog.nonlinear_dpmra,
                    like.catalog.nonlinear_dpmdec,
                ]

                # Remove HGCA's nonlinear correction from Hipparcos epoch
                # Factor of two needed since dpmra is defined to the HG epoch
                μ_h += @SVector [
                    2like.catalog.nonlinear_dpmra,
                    2like.catalog.nonlinear_dpmdec,
                ]
            end

            # Define the order of all components in our unified system
            # Order: hip_iad, ra_hip, dec_hip, ra_hg, dec_hg, ra_dr2, dec_dr2, ra_dr32, dec_dr32, ra_dr3, dec_dr3, ueva_dr3
            component_flags = @SVector [
                :ra_hip, :dec_hip,
                :ra_hg, :dec_hg,
                :ra_dr2, :dec_dr2,
                :ra_dr32, :dec_dr32,
                :ra_dr3, :dec_dr3,
                :ueva_dr3
            ]

            # Build index mask for which components are present
            mask = [flag ∈ like.table.kind for flag in component_flags]
            indices = findall(mask)
            n_components = length(indices)

            # We handle the iad separately from the big covariance matrix below, since it's not correlated
            # and we don't want to factorize a massively bigger matrix than necessary
            if :iad_hip ∈ like.table.kind
                 for i in eachindex(iad_resid)
                    if like.hip_table.reject[i]
                        continue
                    end
                    (;hip_iad_jitter) = θ_obs
                    ll += logpdf(Normal(0, hypot(like.hip_table.sres_renorm[i], hip_iad_jitter) ), iad_resid[i])
                end
            end

            # # We handle the RV unceratinties separatley for the same reasons
            # if :rv_dr3 ∈ like.table.kind
            #     # The likelihood is based on chi-squared distribution
            #     # From the paper: ξ² follows χ²(N-1) distribution
            #     rv_chi2_dist = Chisq(sim.rv_dof)
                
            #     # Calculate log-likelihood
            #     # We want P(observing the catalog sample variance | our model)
            #     # The catalog reports error on median, we need to convert back to sample variance
            #     ε_catalog = like.catalog.radial_velocity_error
            #     N_rv = like.catalog.rv_nb_transits
            #     s_catalog_squared = (2 * N_rv / π) * (ε_catalog^2 - 0.113^2)  # Equation 4 from paper
                
            #     ξ_catalog_squared = (N_rv - 1) * s_catalog_squared / θ_obs.σ_rv_per_transit^2 # these are all in km/s

            #     @show s_catalog_squared ξ_catalog_squared rv_chi2_dist
                
            #     # Compare model to catalog
            #     ll += @show logpdf(rv_chi2_dist, ξ_catalog_squared)
            # end
            if :rv_dr3 ∈ like.table.kind
                ε_catalog = like.catalog.radial_velocity_error
                N_rv = like.catalog.rv_nb_transits
                σ_rv_per_transit = θ_obs.σ_rv_per_transit  # per-transit uncertainty in km/s
                
                # Convert catalog error to sample variance
                s_catalog_squared = (2 * N_rv / π) * (ε_catalog^2 - 0.113^2)

                # non centrality parameter
                ncp = (N_rv - 1) * sim.sample_variance / σ_rv_per_transit^2
                
                # Use NON-CENTRAL chi-squared with signal included
                rv_chi2_dist = NoncentralChisq(N_rv - 1, ncp)
                
                # Catalog's chi-squared statistic
                ξ_catalog_squared = (N_rv - 1) * s_catalog_squared / σ_rv_per_transit^2
                
                # Now this correctly penalizes high-amplitude models
                ll += logpdf(rv_chi2_dist, ξ_catalog_squared)

            end


            μ_dr3_cat, Σ_dr3 = params(dist_dr3)

            #############################################
            # Account for the UEVA-based potential uncertainty delfation of Gaia DR3 positions
            # we know the astrometric excess noise they applied and we are assuming a-priori that the
            # EAN is well-accounted for by a planet

            # DR3 position covariance at central epoch (already inflated by Gaia)
            σ_ra_dr3 = like.catalog.ra_error_central_dr3
            σ_dec_dr3 = like.catalog.dec_error_central_dr3
            ρ_radec_dr3 = like.catalog.ra_dec_corr_central_dr3

            Σ_pos_dr3 = @SMatrix [
                σ_ra_dr3^2                    ρ_radec_dr3*σ_ra_dr3*σ_dec_dr3
                ρ_radec_dr3*σ_ra_dr3*σ_dec_dr3    σ_dec_dr3^2
            ]

            # DR2 position covariance at central epoch
            σ_ra_dr2 = like.catalog.ra_error_central_dr2
            σ_dec_dr2 = like.catalog.dec_error_central_dr2
            ρ_radec_dr2 = like.catalog.ra_dec_corr_central_dr2

            Σ_pos_dr2 = @SMatrix [
                σ_ra_dr2^2                    ρ_radec_dr2*σ_ra_dr2*σ_dec_dr2
                ρ_radec_dr2*σ_ra_dr2*σ_dec_dr2    σ_dec_dr2^2
            ]
            
            # ρ_23 = like.catalog.rho_dr2_dr3
            ρ_dr3_dr2 = √(min(n_dr2, n_dr3) / max(n_dr2, n_dr3))
            # ρ_dr3_dr2 = θ_obs.ρ_dr3_dr2
            
            Σ_cross = @SMatrix [
                ρ_dr3_dr2*σ_ra_dr3*σ_ra_dr2                        ρ_dr3_dr2*ρ_radec_dr3*σ_ra_dr3*σ_dec_dr2
                ρ_dr3_dr2*ρ_radec_dr2*σ_dec_dr3*σ_ra_dr2          ρ_dr3_dr2*σ_dec_dr3*σ_dec_dr2
            ]
            
            # Time baselines for DR3-DR2 scaled position difference
            Δt_ra = (like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_dr2_mjd) / julian_year
            Δt_dec = (like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_dr2_mjd) / julian_year
            
            # Deflation adjustment for DR32 proper motions
            # Only the DR3-contributed terms get deflated
            # deflation_factor_dr3 = 1.0
            d = deflation_factor_dr3

            # Position covariance adjustment (in mas²)
            # ΔΣ_pos = (d^2 - 1) * Σ_pos_dr3 - 2 * (d - 1) * Σ_cross
            ΔΣ_pos = (d^2 - 1) * Σ_pos_dr3 - (d - 1) * (Σ_cross + Σ_cross')

            # Transform to proper motion covariance (mas²/yr²)
            # Different time baselines for RA and Dec
            Tr = @SMatrix [
                1/Δt_ra    0.0
                0.0        1/Δt_dec
            ]

            ΔΣ_dr32 = Tr * ΔΣ_pos * Tr'

            # Extract catalog parameters
            μ_h_cat, Σ_h = isnothing(dist_hip) ? (@SVector[0.,0.], @SMatrix zeros(2,2)) : params(dist_hip) 
            μ_hg_cat, Σ_hg = isnothing(dist_hg) ? (@SVector[0.,0.], @SMatrix zeros(2,2)) : params(dist_hg) 
            μ_dr2_cat, Σ_dr2 = params(dist_dr2)
            μ_dr32_cat, Σ_dr32 = params(dist_dr32)
            T = promote_type(
                eltype(Σ_dr32),
                eltype(Σ_h),
                eltype(Σ_hg),
                eltype(Σ_dr2),
                eltype(Σ_dr32),
                eltype(ΔΣ_dr32),
                eltype(Σ_dr3),
                typeof(deflation_factor_dr3)
            )
            Σ_h = SMatrix{2, 2, T, 4}(Σ_h)
            Σ_hg = SMatrix{2, 2, T, 4}(Σ_hg)
            Σ_dr2 = SMatrix{2, 2, T, 4}(Σ_dr2)
            # Apply deflation adjustment to DR32 covariance
            Σ_dr32 = SMatrix{2, 2, T, 4}(Σ_dr32 .+ ΔΣ_dr32)
            # Σ_dr32 = SMatrix{2, 2, T, 4}(Σ_dr32)
            # Σ_dr32 =SMatrix{2, 2, T, 4}( [
            #     Σ_dr32[1,1] 0
            #     0           Σ_dr32[1,1]
            # ] .+ ΔΣ_dr32)
            # Apply deflation adjustment to DR3 proper motions
            Σ_dr3 = SMatrix{2, 2, T, 4}(Σ_dr3 .* deflation_factor_dr3^2)

            

            μ_catalog_full = @SVector [
                μ_h_cat[1], μ_h_cat[2],     # ra_hip, dec_hip
                μ_hg_cat[1], μ_hg_cat[2],   # ra_hg, dec_hg
                μ_dr2_cat[1], μ_dr2_cat[2], # ra_dr2, dec_dr2
                μ_dr32_cat[1], μ_dr32_cat[2], # ra_dr32, dec_dr32
                μ_dr3_cat[1], μ_dr3_cat[2], # ra_dr3, dec_dr3
                μ_1_3                       # ueva_dr3 catalog value
            ]

            μ_model_full = @SVector [
                μ_h[1], μ_h[2],           # ra_hip, dec_hip
                μ_hg[1], μ_hg[2],         # ra_hg, dec_hg
                μ_dr2[1], μ_dr2[2],       # ra_dr2, dec_dr2
                μ_dr32[1], μ_dr32[2],     # ra_dr32, dec_dr32
                μ_dr3[1], μ_dr3[2],       # ra_dr3, dec_dr3
                UEVA_model                # ueva_dr3 model value
            ]

            # Extract selected components
            μ_catalog_selected = μ_catalog_full[indices]
            μ_model_selected = μ_model_full[indices]

            # Build full covariance matrix (11x11)
            # Initialize with zeros
            Σ_full = zeros(T, 11, 11)

            # Fill in the block diagonal elements (within-epoch covariances)
            Σ_full[1:2, 1:2] .= Σ_h     # Hipparcos
            Σ_full[3:4, 3:4] .= Σ_hg    # HGCA
            Σ_full[5:6, 5:6] .= Σ_dr2   # DR2
            Σ_full[7:8, 7:8] .= Σ_dr32  # DR3-DR2
            Σ_full[9:10, 9:10] .= Σ_dr3 # DR3

            # UEVA variance
            Σ_full[11, 11] = UEVA_unc^2

            # Add cross-epoch correlations between DR2 and DR3
            
            # Compute the cross-correlation matrix K
            K = ρ_dr3_dr2 * sqrt(Σ_dr2) * sqrt(Σ_dr3)'
            
            # Fill in the cross-correlation blocks
            Σ_full[5:6, 9:10] .= K      # DR2 -> DR3
            Σ_full[9:10, 5:6] .= K'     # DR3 -> DR2 (transpose)

            # Extract the submatrix for selected components
            Σ_selected = Σ_full[indices, indices]

            # @show (μ_catalog_selected .- μ_model_selected ) ./ sqrt.(diag(Σ_selected))
            # @show μ_catalog_selected μ_model_selected

            # Compute likelihood
            if n_components == 1
                # Single component - use univariate normal
                ll += logpdf(Normal(μ_catalog_selected[1], sqrt(Σ_selected[1,1])), μ_model_selected[1])
            else
                # Multiple components - use multivariate normal
                local dist_selected
                try
                    dist_selected = MvNormal(μ_catalog_selected, Hermitian(Σ_selected))
                    ll += logpdf(dist_selected, μ_model_selected)
                catch err
                    if err isa PosDefException
                        ll = convert(T, -Inf)
                    else
                        rethrow(err)
                    end
                end
            end
        end
    end

    if isnan(ll)
        return convert(T, -Inf)
    end

    return ll
end

function simulate(like::GaiaHipparcosUEVAJointObs, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    # TODO: optimize this, we only need to grab the epochs here -- it'll be faster
    if hasproperty(θ_obs, :transits)
        (;transits) = θ_obs 
        if eltype(transits) <: AbstractFloat
            transits = Int.(transits)
        end
        # The list of missed transits must be unique
        if length(unique(transits)) < length(transits)
            return nothing
        end
        ii = transits
        gaia_table = like.gaia_table[ii,:]
    else
        gaia_table = like.gaia_table
    end

    istart_dr2 = findfirst(>=(meta_gaia_DR2.start_mjd), vec(gaia_table.epoch))
    iend_dr2 = findlast(<=(meta_gaia_DR2.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr2)
        istart_dr2 = 1
    end
    if isnothing(iend_dr2)
        iend_dr2 = length(gaia_table.epoch)
    end

    istart_dr3 = findfirst(>=(meta_gaia_DR3.start_mjd), vec(gaia_table.epoch))
    iend_dr3 = findlast(<=(meta_gaia_DR3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(gaia_table.epoch)
    end

    iad_resid  = zeros(size(like.hip_table,1)); fill!(iad_resid, 0)
    Δα_mas_hip = zeros(size(like.hip_table,1)); fill!(Δα_mas_hip, 0)
    Δδ_mas_hip = zeros(size(like.hip_table,1)); fill!(Δδ_mas_hip, 0)
    Δα_mas_dr2 = zeros(iend_dr2-istart_dr2+1); fill!(Δα_mas_dr2, 0)
    Δδ_mas_dr2 = zeros(iend_dr2-istart_dr2+1); fill!(Δδ_mas_dr2, 0)
    Δα_mas_dr3 = zeros(iend_dr3-istart_dr3+1); fill!(Δα_mas_dr3, 0)
    Δδ_mas_dr3 = zeros(iend_dr3-istart_dr3+1); fill!(Δδ_mas_dr3, 0)

    buffers = (;iad_resid, Δα_mas_hip, Δδ_mas_hip, Δα_mas_dr2, Δδ_mas_dr2, Δα_mas_dr3, Δδ_mas_dr3)

    out = simulate!(buffers, like, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    return out 
end

function simulate!(buffers, like::GaiaHipparcosUEVAJointObs, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) 

    (;Δα_mas_hip, Δδ_mas_hip, Δα_mas_dr2, Δδ_mas_dr2, Δα_mas_dr3, Δδ_mas_dr3, iad_resid, ) = buffers


    T = _system_number_type(θ_system)

    # Generate simulated observations from this sample draw
    # Get Gaia noise parameters from observation variables
    (;σ_att, σ_AL, σ_calib,) = θ_obs
    σ_formal = sqrt(σ_att^2 + σ_AL^2)

    gaia_n_dof = like.catalog.astrometric_params_solved_dr3 == 31 ? 5 : 6

    # The gaia_table and A_prepared_5_dr3/A_prepared_5_dr2 include all available
    # visibility windows, not filtered to specifically be DR2 or DR3. 
    # Here we may further reject some more to marginalize over
    # unknown missed/rejected transits.
    # In theory these could be different between DR2 and DR3 but we assume they aren't.
    if hasproperty(θ_obs, :transits)
        (;transits) = θ_obs 
        if eltype(transits) <: AbstractFloat
            transits = Int.(transits)
        end
        # The list of missed transits must be unique
        if length(unique(transits)) < length(transits)
            return nothing
        end
        ii = transits
        gaia_table = like.gaia_table[ii,:]
        A_prepared_5_dr3 = view(like.A_prepared_5_dr3, ii,:)
        A_prepared_5_dr2 = view(like.A_prepared_5_dr2, ii,:)
    else
        gaia_table = like.gaia_table
        A_prepared_5_dr3 = like.A_prepared_5_dr3
        A_prepared_5_dr2 = like.A_prepared_5_dr2
    end

    if hasproperty(θ_obs, :transits_rv)
        (;transits_rv) = θ_obs 
        if eltype(transits_rv) <: AbstractFloat
            transits_rv = Int.(transits_rv)
        end
        # The list of missed transits must be unique
        if length(unique(transits_rv)) < length(transits_rv)
            return nothing
        end
        jj = collect(transits_rv)#Int.(sort(setdiff(1:length(like.gaia_table.epoch), transits_rv)))
        gaia_table_rv = like.gaia_table[jj,:]
    else
        gaia_table_rv = like.gaia_table
    end

    # Now we fit a no-planet (zero mass planet) sky path model to this data.
    # These should be fit using the appropriate catalog reference epoch so 
    # that they can be compared correctly.

    absolute_orbits = false
    for orbit in orbits
        absolute_orbits |= orbit isa AbsoluteVisual
        # TODO: could check in a more user-friendly way
        # that we don't have a mismatch of different orbit types
        # for different planets?
    end



    # Helper functions to either get the static pmra from the orbital elements,
    # or, if using an AbsoluteVisualOrbit, get the propagated pmra at the
    # current epoch accounting for barycentric motion.
    function propagate_astrom(orbits::NTuple{N,<:PlanetOrbits.AbsoluteVisualOrbit} where N, epoch_ra, epoch_dec)
        o = first(orbits)
        sol_ra = orbitsolve(o, epoch_ra)
        cmp_ra = sol_ra.compensated
        sol_dec = orbitsolve(o, epoch_dec)
        cmp_dec = sol_dec.compensated
        # Account for the instantaneous differential light travel time apparent acceleration.
        # Treat as linear for the duration of Gaia or Hipparcos
        t1 = max(epoch_ra, epoch_dec)
        Δt = 100
        t2 = t1 + Δt
        sol = epoch_ra >= epoch_dec ? sol_ra : sol_dec
        sol′ = orbitsolve(o,t2)
        # This isn't right! This is double counting the proper motion which already goes into ra/dec
        # Take change in delta_time and multiply it by pmra/pmdec
        diff_lt_app_pmra = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmra2
        diff_lt_app_pmdec = (sol′.compensated.t_em_days - sol.compensated.t_em_days - Δt)/Δt*sol.compensated.pmdec2
        return cmp_ra.ra2, cmp_dec.dec2, cmp_ra.pmra2+diff_lt_app_pmra, cmp_dec.pmdec2+diff_lt_app_pmdec
        # return (
        #     cmp_ra.ra2 - Δα_dr3/60/60/1000/cosd(cmp_dec.dec2),
        #     cmp_dec.dec2 - Δδ_dr3/60/60/1000,
        #     cmp_ra.pmra2+diff_lt_app_pmra - Δpmra_dr3,
        #     cmp_dec.pmdec2+diff_lt_app_pmdec - Δpmdec_dr3
        # )
    end
    function propagate_astrom(orbits::Tuple{}, _, _)
        return 0.0, 0.0, θ_system.pmra, θ_system.pmdec
    end
    function propagate_astrom(orbits::Any, _, _)
        return 0.0, 0.0, θ_system.pmra, θ_system.pmdec
    end




    ################################
    # DR3
    istart_dr3 = findfirst(>=(meta_gaia_DR3.start_mjd), vec(gaia_table.epoch))
    iend_dr3 = findlast(<=(meta_gaia_DR3.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr3)
        istart_dr3 = 1
    end
    if isnothing(iend_dr3)
        iend_dr3 = length(gaia_table.epoch)
    end
    gaia_table_dr3 = @views gaia_table[istart_dr3:iend_dr3]
    # gaia_table_dr3.epoch .+= Δepoch_dr3_days
    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        if planet_mass_msol == 0.0
            continue
        end
        if hasproperty(θ_obs, :fluxratio)
            if θ_obs.fluxratio isa Number
                fluxratio = θ_obs.fluxratio
            else
                fluxratio = θ_obs.fluxratio[i_planet]
            end
        else
            fluxratio = 0.0
        end
        _simulate_skypath_perturbations!(
            Δα_mas_dr3, Δδ_mas_dr3,
            gaia_table_dr3, orbit,
            planet_mass_msol, fluxratio,
            orbit_solutions[i_planet],
            -1, T,
        )
    end

    out_dr3 = fit_5param_prepared(
        view(A_prepared_5_dr3, istart_dr3:iend_dr3,:),
        view(gaia_table, istart_dr3:iend_dr3),
        Δα_mas_dr3, Δδ_mas_dr3, 0.0, σ_formal;
        include_chi2=Val(true)
    )
    Δα_dr3, Δδ_dr3, Δpmra_dr3, Δpmdec_dr3 = out_dr3.parameters
    # Rigorously propagate the linear proper motion component in spherical coordinates
    # Account for within-gaia differential light travel time 
    α_dr3₀, δ_dr3₀, pmra_dr3₀, pmdec_dr3₀ = propagate_astrom(orbits, like.catalog.epoch_ra_dr3_mjd, like.catalog.epoch_dec_dr3_mjd)
    μ_dr3 = @SVector [pmra_dr3₀ + Δpmra_dr3, pmdec_dr3₀ + Δpmdec_dr3]

    # Note: we shift the entire reference frame so that the proper motion is defined on the primary star
    # all proper motions derived below are shifted the perturbation in DR3 
    # This vastly improves sampling efficiency.
    # Leave Δpmdec_dr3 - Δpmdec_dr3 above as an explicit reminder about this ^

    # TODO: efficiency: since we assume all DR2 epochs are a subset of DR3, we
    # could re-use part of the _simulate_skypath_perturbations! done for DR3

    ################################
    # DR2
    istart_dr2 = findfirst(>=(meta_gaia_DR2.start_mjd), vec(gaia_table.epoch))
    iend_dr2 = findlast(<=(meta_gaia_DR2.stop_mjd), vec(gaia_table.epoch))
    if isnothing(istart_dr2)
        istart_dr2 = 1
    end
    if isnothing(iend_dr2)
        iend_dr2 = length(gaia_table.epoch)
    end
    gaia_table_dr2 = @views gaia_table[istart_dr2:iend_dr2]
    # gaia_table_dr2.epoch .+= Δepoch_dr2_days
    for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
        planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
        if planet_mass_msol == 0.0
            continue
        end
        if hasproperty(θ_obs, :fluxratio)
            if θ_obs.fluxratio isa Number
                fluxratio = θ_obs.fluxratio
            else
                fluxratio = θ_obs.fluxratio[i_planet]
            end
        else
            fluxratio = 0.0
        end
        _simulate_skypath_perturbations!(
            Δα_mas_dr2, Δδ_mas_dr2,
            gaia_table_dr2, orbit,
            planet_mass_msol, fluxratio,
            orbit_solutions[i_planet],
            -1, T
        )
    end

    out = fit_5param_prepared(view(A_prepared_5_dr2, istart_dr2:iend_dr2,:), view(gaia_table, istart_dr2:iend_dr2), Δα_mas_dr2, Δδ_mas_dr2)
    # out = fit_4param_prepared(hgca_like.gaialike.A_prepared_4, gaia_table, Δα_mas_dr2, Δδ_mas_dr2)
    Δα_dr2, Δδ_dr2, Δpmra_dr2, Δpmdec_dr2 = out.parameters
    # Rigorously propagate the linear proper motion component in spherical coordinates
    # Account for within-gaia differential light travel time 
    α_dr2₀, δ_dr2₀, pmra_dr2₀, pmdec_dr2₀ = propagate_astrom(orbits, like.catalog.epoch_ra_dr2_mjd, like.catalog.epoch_dec_dr2_mjd)
    μ_dr2 = @SVector [pmra_dr2₀ + Δpmra_dr2, pmdec_dr2₀ + Δpmdec_dr2]

        

    ################################
    # Hipparcos
    if isnothing(like.catalog.dist_hip)
        # type stable since dist_hip is part of the likelihood type parameters
        # ie. we statically know which of these branches will be taken.
        μ_h = @SVector [zero(T), zero(T)]
    else


        for (i_planet,(orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
            planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
            if planet_mass_msol == 0.0
                continue
            end
            if hasproperty(θ_obs, :fluxratio)
                if θ_obs.fluxratio isa Number
                    fluxratio = θ_obs.fluxratio
                else
                    fluxratio = θ_obs.fluxratio[i_planet]
                end
            else
                fluxratio = 0.0
            end
            _simulate_skypath_perturbations!(
                Δα_mas_hip, Δδ_mas_hip,
                like.hip_table, orbit,
                planet_mass_msol, fluxratio,
                orbit_solutions[i_planet],
                -1, T
            )
        end
        if like.include_iad
            out = fit_5param_prepared(like.A_prepared_5_hip, like.hip_table, Δα_mas_hip, Δδ_mas_hip, like.hip_table.res, like.hip_table.sres)
        else
            out = fit_5param_prepared(like.A_prepared_5_hip, like.hip_table, Δα_mas_hip, Δδ_mas_hip)
        end
        Δα_h, Δδ_h, Δpmra_h, Δpmdec_h = out.parameters
        α_h₀, δ_h₀, pmra_h₀, pmdec_h₀ = propagate_astrom(orbits, like.catalog.epoch_ra_hip_mjd, like.catalog.epoch_dec_hip_mjd)
        μ_h = @SVector [pmra_h₀ + Δpmra_h, pmdec_h₀ + Δpmdec_h]


        ################################
        # Hipparcos IAD
            

        # We can include the Hipparcos IAD in the following way
        #=
            * The Hipaprcos IAD from the Java tool (what we have) is not necessarily consistent with the composite catalog point used in the HGCA, which is more accurate
            * We want to keep fitting the composite point with well-calibrated uncertainties compared to Gaia; however, we can still additionally fit the Hipparcos IAD
            * using a flexible offset and linear trend, just like RV data.
            * We need the following user-supplied parameters:
                * iad_Δra ("zero point offset for pmra [mas]")
                * iad_Δdec ("zero point offset for pmdec [mas]")
                * iad_pmra ("slope in pmra [mas/yr] ")
                * iad_pmdec ("slope in pmra [mas/yr] ")
                * hip_iad_jitter (astrometric excess jitter term)
            * This includes information from the IAD in a way that's still consistent with the calibrated composite Hipparcos catalog values for ra, dec, pmra, pmdec.
            * This doesn't count any information twice, because we here we are only fitting curvature, while above we are only fitting positions and proper motions.

            Probably, one could put informed priors on most of these properies but we leave them wide open.
        =#
        # if like.include_iad
            (;iad_Δra,
                iad_Δdec,
                iad_pmra,
                iad_pmdec,
                iad_Δplx) = θ_obs


            # like.hip_table.res, like.hip_table.sres
            # Δα_mas_hip, Δδ_mas_hip

            α✱_models=[]
            δ_models=[]
            for i_epoch in eachindex(like.hip_table.epoch, Δα_mas_hip, Δδ_mas_hip)
                delta_time_julian_year = (like.hip_table.epoch[i_epoch] - hipparcos_catalog_epoch_mjd) / julian_year
                plx_at_epoch = like.hip_sol.plx + iad_Δplx
                # TODO: for very very nearby or high RV objects with specific IDs, our main Hipparcos code correctly
                # accounts for changing parallax vs time. This does not.
                α✱_model = iad_Δra - Δα_h + plx_at_epoch * (
                            like.hip_table.x[i_epoch] * sind(like.hip_sol.radeg) -
                            like.hip_table.y[i_epoch] * cosd(like.hip_sol.radeg)
                ) + delta_time_julian_year * (iad_pmra - Δpmra_h)
                δ_model = iad_Δdec - Δδ_h + plx_at_epoch * (
                            like.hip_table.x[i_epoch] * cosd(like.hip_sol.radeg) * sind(like.hip_sol.dedeg) +
                            like.hip_table.y[i_epoch] * sind(like.hip_sol.radeg) * sind(like.hip_sol.dedeg) -
                            like.hip_table.z[i_epoch] * cosd(like.hip_sol.dedeg)
                ) + delta_time_julian_year * (iad_pmdec - Δpmdec_h)
                
                α✱_model_with_perturbation = α✱_model + Δα_mas_hip[i_epoch]
                δ_model_with_perturbation = δ_model + Δδ_mas_hip[i_epoch]

                # We've hit the same problem again. As the perturbation gets huge, we
                # start to need to adjust these parameters to bring the measurements back over to near
                # the Hipparcos measurements.
                # But in this case, do we really care? We are only fitting curvature.
                # Can we subtract the average or something?




                # push!(α✱_models,α✱_model_with_perturbation)
                # push!(δ_models, δ_model_with_perturbation)


                # Calculate differences in milliarcseconds by mapping into a local tangent plane,
                # with the local coordinate defined by the Hipparcos solution at this *epoch* (not 
                # at the best solution)
                point = @SVector [
                    α✱_model_with_perturbation,
                    δ_model_with_perturbation,
                ]

                # Two points defining a line along which the star's position was measured
                line_point_1 = @SVector [
                    like.hip_table.α✱ₘ[i_epoch][1],
                    like.hip_table.δₘ[i_epoch][1],
                ]
                line_point_2 = @SVector [
                    like.hip_table.α✱ₘ[i_epoch][2],
                    like.hip_table.δₘ[i_epoch][2],
                ]
                # Distance from model star Delta position to this line where it was measured 
                # by the satellite
                resid = distance_point_to_line(point, line_point_1, line_point_2) # degrees

                iad_resid[i_epoch] = resid
            end

            # @show θ_system.ra  θ_system.dec iad_Δra iad_Δdec iad_pmra iad_pmdec
            # f,a,p=Main.scatterlines(
            #     α✱_models, δ_models
            # )
            # Main.stem(f[2,1], α✱_models, iad_resid)
            # f|>display
        # end

    end


    ################################
    # H-G (if H is available) and DR3-DR2

    # Simple linear approximation: don't deal with curvature & secular acceleration directly
    if absolute_orbits

        if isnothing(like.catalog.dist_hip)
            pmra_hg_model = zero(T)
            pmdec_hg_model = zero(T)
        else
            # HG
            Δα_hg_prop = (α_dr3₀ - α_h₀)*60*60*1000*cosd((δ_dr3₀ + δ_h₀)/2)
            Δδ_hg_prop = (δ_dr3₀ - δ_h₀)*60*60*1000
            pmra_hg_model = (Δα_dr3 - Δα_h + Δα_hg_prop) / (
                like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_hip_mjd
            )*julian_year
            pmdec_hg_model = (Δδ_dr3 - Δδ_h + Δδ_hg_prop) / (
                like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_hip_mjd
            )*julian_year
        end


        # DR3-DR2
        Δα_dr32_prop = (α_dr3₀ - α_dr2₀)*60*60*1000*cosd((δ_dr3₀ + δ_dr2₀)/2)
        Δδ_dr32_prop = (δ_dr3₀ - δ_dr2₀)*60*60*1000
        pmra_dr32_model = (Δα_dr3 - Δα_dr2 + Δα_dr32_prop) / (
            like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_dr2_mjd
        )*julian_year
        pmdec_dr32_model = (Δδ_dr3 - Δδ_dr2 + Δδ_dr32_prop) / (
            like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_dr2_mjd
        )*julian_year

    else
        if isnothing(like.catalog.dist_hip)
            pmra_hg_model = zero(T)
            pmdec_hg_model = zero(T)
        else
            pmra_hg_model = (Δα_dr3 - Δα_h) / (
                    like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_hip_mjd
            )*julian_year + θ_system.pmra
            pmdec_hg_model = (Δδ_dr3 - Δδ_h) / (
                like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_hip_mjd
            )*julian_year + θ_system.pmdec
        end

        pmra_dr32_model = (Δα_dr3 - Δα_dr2) / (
                like.catalog.epoch_ra_dr3_mjd - like.catalog.epoch_ra_dr2_mjd
        )*julian_year + θ_system.pmra
        pmdec_dr32_model = (Δδ_dr3 - Δδ_dr2) / (
            like.catalog.epoch_dec_dr3_mjd - like.catalog.epoch_dec_dr2_mjd
        )*julian_year + θ_system.pmdec


    end


    # μ_hg = @SVector [pmra_hg_model - Δpmra_dr3, pmdec_hg_model - Δpmdec_dr3]
    # μ_dr32 = @SVector [pmra_dr32_model - Δpmra_dr3, pmdec_dr32_model - Δpmdec_dr3]
    μ_hg = @SVector [pmra_hg_model, pmdec_hg_model]
    μ_dr32 = @SVector [pmra_dr32_model, pmdec_dr32_model]

    ##############################
    # DR3 UEVA calculation and uncertainty deflation
    # From Gaia catalog:
    (;
        astrometric_chi2_al_dr3,         
        astrometric_n_good_obs_al_dr3,   
        astrometric_matched_transits_dr3,
        astrometric_excess_noise_dr3,
        ruwe_dr3,
    ) = like.catalog

    N = astrometric_n_good_obs_al_dr3
    N_FoV = astrometric_matched_transits_dr3
    N_AL = N / N_FoV

    # Calculate Gaia's published UEVA (what they measured)
    if like.ueva_mode == :EAN
        UEVA_Gaia = astrometric_excess_noise_dr3^2 + σ_att^2 + σ_AL^2
    elseif like.ueva_mode == :RUWE
        u0 = 1/ruwe_dr3 * sqrt(astrometric_chi2_al_dr3/(N - gaia_n_dof))
        UEVA_Gaia = (ruwe_dr3 * u0)^2 * σ_formal^2
    else
        error("Unsupported mode (should be :EAN or :RUWE, was $(like.ueva_mode)")
    end

    # Calculate expected UEVA for a single star (Eq. D.8 from paper)
    μ_UEVA_single = (N_AL / (N - gaia_n_dof)) * 
                ((N_FoV - gaia_n_dof) * σ_calib^2 + N_FoV * σ_AL^2)

    # And its standard deviation (Eq. D.9)
    σ_UEVA_single = sqrt(
        2 * N_AL / (N - gaia_n_dof)^2 * (
            N_AL * (N_FoV - gaia_n_dof) * σ_calib^4 + 
            N_FoV * σ_AL^4 + 
            2 * N_FoV * σ_AL^2 * σ_calib^2
        )
    )

    μ_1_3 = UEVA_Gaia^(1/3)
    UEVA_unc = σ_UEVA_single * μ_UEVA_single^(-2/3) / 3 # divide by 3 due to cube root transformation

    # Calculate model-predicted UEVA from our fit
    chi2_astro_scaled = out_dr3.chi_squared_astro * N_AL
    UEVA_model_raw = (chi2_astro_scaled * σ_formal^2) / (N - gaia_n_dof)

    # For the UEVA likelihood, use cube-root transformation (Eq. 27, Sect 5.1.1)
    UEVA_model_1 = (chi2_astro_scaled * σ_formal^2) / (N_AL * N_FoV - gaia_n_dof)
    UEVA_model = cbrt(UEVA_model_1 + μ_UEVA_single)

    # Calculate the "deflation factor" -- the amount of Gaia's inflated uncertainties
    # that come from our now-explained companion model

    # What a 5-param fit would measure with this companion model
    UEVA_predicted = UEVA_model_raw + μ_UEVA_single

    # Deflation factor
    deflation_factor_raw = sqrt(μ_UEVA_single / UEVA_Gaia)
    # equivalent to :
    # deflation_factor_raw = sqrt(1 - UEVA_orbital_perturb / UEVA_Gaia) 


    # @show deflation_factor_raw 

    # Clamp to valid range
    deflation_factor_dr3 = deflation_factor_raw > 1.0 ? 1.0 : deflation_factor_raw

    # # for data simulation purposes, here is an estimate of what these parameters would produce for RUWE
    # # given everything we know about the gaia uncertainties etc.
    # # Given: UEVA, u0, σ_formal, calculate ruwe
    # # u0 = 1/ruwe_dr3*sqrt(astrometric_chi2_al_dr3/(astrometric_n_good_obs_al_dr3-gaia_n_dof))
    # # UEVA = ruwe_dr3^2 * u0^2 * σ_formal^2
    # # UEVA/( u0^2 * σ_formal^2) = ruwe_dr3^2
    # ruwe_dr3 = sqrt(UEVA_model/( u0^2 * σ_formal^2))

    # Forward-model the Gaia RV uncertainty using approach from the "paired" tool:
    # https://arxiv.org/pdf/2206.11275
    if :rv_dr3 ∈ like.table.kind
        # Get the per-transit RV uncertainty from catalog
        σ_rv_per_transit = exp(θ_obs.σ_rv_per_transit)  # Convert from log space to km/s
        
        # Simulate RV measurements at Gaia epochs
        rv_model = zeros(T, length(gaia_table_rv.epoch))
        
        for (i_planet, (orbit, θ_planet)) in enumerate(zip(orbits, θ_system.planets))
            planet_mass_msol = θ_planet.mass*Octofitter.mjup2msol
            if planet_mass_msol == 0.0
                continue
            end
            # Calculate RV at each epoch
            for (i, epoch) in enumerate(gaia_table_rv.epoch)
                sol = orbitsolve(orbit, epoch)
                # Add the RV component from this planet
                rv_model[i] = radvel(sol, planet_mass_msol)/1e3 # barycentric rv in km/s
            end
        end
        
        # Calculate sample variance
        rv_mean = mean(rv_model)
        sample_variance = sum((rv_model .- rv_mean).^2) / (length(rv_model) - 1)
        
        # The chi-squared statistic (equation 6)
        N_rv = like.catalog.rv_nb_transits
        ξ_squared = (N_rv - 1) * sample_variance / σ_rv_per_transit^2

        # Store for likelihood calculation
        # rv_chi2_stat = ξ_squared
        rv_dof = N_rv - 1

        ε_catalog = like.catalog.radial_velocity_error

        # Convert catalog error to sample variance
        s_catalog_squared = (2 * N_rv / π) * (ε_catalog^2 - 0.113^2)
    else
        # rv_chi2_stat = convert(T, NaN)
        rv_dof = convert(T, NaN)
        s_catalog_squared  = convert(T, NaN)
        rv_dof = convert(T, NaN)
        rv_mean = convert(T, NaN)
        sample_variance = convert(T, NaN)
        s_catalog_squared = convert(T, NaN)
    end

    # Adjust the reference frame such that, effectively, the pmra/pmdec system variables are referring to the primary
    # instead of the barycentre.
    # Specifically, the primary's proper motion at this epoch:
    μ_h    = μ_h     .- @SVector [Δpmra_dr3, Δpmdec_dr3,]
    μ_hg   = μ_hg    .- @SVector [Δpmra_dr3, Δpmdec_dr3,]
    μ_dr2  = μ_dr2   .- @SVector [Δpmra_dr3, Δpmdec_dr3,]
    μ_dr32 = μ_dr32  .- @SVector [Δpmra_dr3, Δpmdec_dr3,]
    μ_dr3  = μ_dr3   .- @SVector [Δpmra_dr3, Δpmdec_dr3,]


    return (;

        # UEVA: EAN/RUWE
        UEVA_model,
        UEVA_unc,
        μ_1_3,

        # Packaged up nicely
        μ_h,
        μ_hg,
        μ_dr2,
        μ_dr32,
        μ_dr3,
        μ = (@SVector [μ_h[1],μ_h[2],μ_hg[1],μ_hg[2],μ_dr2[1],μ_dr2[2],μ_dr32[1],μ_dr32[2],μ_dr3[1],μ_dr3[2],UEVA_model,sample_variance]),

        n_dr3 = iend_dr3 - istart_dr3 + 1,
        n_dr2 = iend_dr2 - istart_dr2 + 1,

        # rv_chi2_stat,
        rv_dof,
        rv_mean,
        sample_variance,
        s_catalog_squared,

        deflation_factor_dr3,

        # Individual
        # TODO: get rid of these
        pmra_hip_model=μ_h[1],
        pmdec_hip_model=μ_h[2],
        pmra_hg_model=μ_hg[1],
        pmdec_hg_model=μ_hg[2],
        pmra_dr2_model=μ_dr2[1],
        pmdec_dr2_model=μ_dr2[2],
        pmra_dr32_model=μ_dr32[1],
        pmdec_dr32_model=μ_dr32[2],
        pmra_dr3_model=μ_dr3[1],
        pmdec_dr3_model=μ_dr3[2],

        Δα_dr3, Δδ_dr3, Δpmra_dr3, Δpmdec_dr3


    )
end



# Generate new astrometry observations
function Octofitter.generate_from_params(like::GaiaHipparcosUEVAJointObs, ctx::SystemObservationContext; add_noise)
    (; θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start) = ctx

    sim = simulate(like, θ_system, θ_obs, orbits, orbit_solutions, orbit_solutions_i_epoch_start)
    (; μ_h, μ_hg, μ_dr2, μ_dr32, μ_dr3, UEVA_model, UEVA_unc, μ_1_3) = sim

    error("Simulating GaiaHipparcosUEVAJointObs not implemented yet")
    return like
end

# Backwards compatibility alias
const GaiaHipparcosUEVAJointLikelihood = GaiaHipparcosUEVAJointObs

# User-friendly alias
const G23HObs = GaiaHipparcosUEVAJointObs

export GaiaHipparcosUEVAJointObs, GaiaHipparcosUEVAJointLikelihood, G23HObs
