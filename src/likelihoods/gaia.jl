
const info_gaia_DR2 = (;
    start_mjd=mjd("2014-07-25"), # 10:30 UTC
    stop_mjd =mjd("2016-05-23"), # 11:35 UTC
)

const info_gaia_DR3 = (;
    start_mjd=mjd("2014-07-25"), # 10:30 UTC
    stop_mjd =mjd("2017-05-28"), # 11:35 UTC
)


abstract type GaiaCatalogLikelihood_v2 <: AbstractLikelihood end
struct GaiaDR3_v2{TCat,TTable,TDist} <: GaiaCatalogLikelihood_v2
    gaia_id::Int
    dr3::TCat
    table::TTable
    catalog_distribution::TDist
end
function GaiaDR3_v2(;gaia_id)

    # Query Gaia archive for DR3 solution
    dr3 = Octofitter._query_gaia_dr3(;gaia_id)

    # Get predicted GAIA scan epochs and angles
    forecast_table = FlexTable(GHOST_forecast(;gaia_id))
    forecast_table.epoch = jd2mjd.(forecast_table.var" ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]")

    # Calculate the scan angle using the same convention that Hipparcos uses,
    # namely psi = π/2 + scanAngle
    forecast_table.cosϕ = cos.(π/2 .+ forecast_table.var" scanAngle[rad]")
    forecast_table.sinϕ = sin.(π/2 .+ forecast_table.var" scanAngle[rad]")

    # Get the Earth's position at those epochs
    earth_pos_vel = geocentre_position_query.(forecast_table.epoch)

    # merge the Gaia scan prediction and geocentre position results into one table
    table = FlexTable(eachcol(forecast_table)..., eachcol(earth_pos_vel)...)

    # This will be hoisted
    μ = [
        dr3.parallax,
        dr3.ra,# deg
        dr3.dec,# deg
        dr3.pmra, 
        dr3.pmdec,
    ]
    σ = [
        dr3.parallax_error,
        dr3.ra_error/ 60/60/1000,
        dr3.dec_error / 60/60/1000,
        dr3.pmra_error ,
        dr3.pmdec_error,
    ]
    C = [
        # plx                   ra                      dec                     pmra                    pmdec
        1                       dr3.ra_parallax_corr    dr3.dec_parallax_corr   dr3.parallax_pmra_corr  dr3.parallax_pmdec_corr
        dr3.ra_parallax_corr    1                       dr3.ra_dec_corr         dr3.ra_pmra_corr        dr3.ra_pmdec_corr
        dr3.dec_parallax_corr   dr3.ra_dec_corr         1                       dr3.dec_pmra_corr       dr3.dec_pmdec_corr
        dr3.parallax_pmra_corr  dr3.ra_pmra_corr        dr3.dec_pmra_corr       1                       dr3.pmra_pmdec_corr
        dr3.parallax_pmdec_corr dr3.ra_pmdec_corr       dr3.dec_pmdec_corr      dr3.pmra_pmdec_corr     1
    ]
    Σ = Diagonal(σ) * C * Diagonal(σ)
    dist = MvNormal(μ,Hermitian(Σ))

    return GaiaDR3_v2(gaia_id, dr3, table, dist)
end




function ln_like(gaialike::GaiaCatalogLikelihood_v2, θ_system, orbits, num_epochs::Val{L}=Val(length(gaialike.table))) where L

    T = _system_number_type(θ_system)

    # TODO: expand to work with multiple, named planets
    planet_mass_msol = θ_system.planets.b.mass*Octofitter.mjup2msol
    # typical_along_scan_uncertainty_mas = θ_system.typical_along_scan_uncertainty_mas
    # along_scan_uncertainty_mas = sqrt(typical_along_scan_uncertainty_mas^2 + gaialike.dr3.astrometric_excess_noise)
    along_scan_uncertainty_mas = θ_system.typical_along_scan_uncertainty_mas


    sample_orbit = only(orbits)

    # Generate simulated observations from this sample draw
    (;α_model,δ_model,αₘ,δₘ) = _simulate_skypath_observations(gaialike, sample_orbit::AbsoluteVisual, planet_mass_msol, T)
    
    # Now we fit a no-planet (zero mass planet) sky path model to this data.
    # nested diff should be okay once we put this in a dedicated rule dedicated rule.
    func = OptimizationFunction((args...)->-_gaia_skypath_loglikelihood(args...), AutoForwardDiff())
    # func = OptimizationFunction((args...)->-_gaia_skypath_loglikelihood(args...), AutoFiniteDiff())
    guess = [
        sample_orbit.plx,
        sample_orbit.ra,
        sample_orbit.dec,
        # rv,
        sample_orbit.pmra,
        sample_orbit.pmdec,
    ]
    # Then iterate with qusi-Newton
    data = (;αₘ,δₘ,gaialike.table, sample_orbit.ref_epoch, along_scan_uncertainty_mas)
    prob = OptimizationProblem(func, guess, data)
    # could very well be a gradient problem! Better set up implicit diff sooner rather than later
    # sol = solve(prob,
    #     # LBFGS(
    #     #     m=4,
    #     #     # linesearch=OptimizationOptimJL.Optim.LineSearches.HagerZhang(linesearchmax=50),
    #     #     linesearch=OptimizationOptimJL.Optim.LineSearches.BackTracking(),
    #     #     alphaguess=OptimizationOptimJL.Optim.LineSearches.InitialHagerZhang()
    #     # ),
    #     ConjugateGradient(),
    #     x_abstol=NaN,
    #     x_reltol=NaN,
    #     f_abstol=NaN,
    #     f_reltol=NaN,
    #     g_abstol= 1e-1,
    #     g_reltol= 1e-1,
    #     # g_tol=1e-15, f_tol=NaN,  x_tol=NaN, x_reltol=NaN, 
    #     allow_f_increases=true,
    #     iterations=10000000
    # )
    # display(sol)
    # prob = OptimizationProblem(func, sol.u, data)
    sol = solve(prob,
        LBFGS(
            m=4,
            # linesearch=OptimizationOptimJL.Optim.LineSearches.HagerZhang(linesearchmax=50),
            linesearch=OptimizationOptimJL.Optim.LineSearches.BackTracking(),
            alphaguess=OptimizationOptimJL.Optim.LineSearches.InitialHagerZhang()
        ),
        # ConjugateGradient(),
        # x_abstol=NaN,
        # x_reltol=NaN,
        # f_abstol=NaN,
        # f_reltol=NaN,
        g_abstol= 1e-5,
        g_reltol= 1e-5,
        # g_tol=1e-15, f_tol=NaN,  x_tol=NaN, x_reltol=NaN, 
        allow_f_increases=true,
        iterations=10000
    )
    # sol = solve(prob, NelderMead, g_tol=1e-9)

    if sol.retcode != Optimization.ReturnCode.Success
        display(sol.original)
        @error("Optimization failed")
        return -Inf
    end

    # Version 1: Compare output means to catalog values and covariances
    # return logpdf(gaialike.catalog_distribution,sol.u)

    # Instead of comparing best fit to catalog mean, let's compute a distance metric between distrbutions.
    # For this, we will estimate the gaussian uncertainties using the hessian at the best fitting model
    
    H = ForwardDiff.hessian(u->_gaia_skypath_loglikelihood(u,data), sol.u)
    local gaia_post_covariance
    try
        gaia_post_covariance = -inv(H)
    catch 
        @warn "singular hessian encountered (maxlog=1)" maxlog=1
        return -Inf
    end
    # safesqrt(num) = num >=0 ? sqrt(num) : Inf
    # unc = safesqrt.(diag(Diagonal(-inv(H))))

    # We use the inverse Hessian from the gradient descent as the covariance matrix.
    # In case of numerical error (non-hermitian matrix) we regularize it by adding 
    # a value along the diagonals. 
    # P = MvNormal(sol.u, Hermitian(gaia_post_covariance))
    P = nothing
    for e in (0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6)#, 1e-5, 1e-4, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3, 1e4)
        try
            P = MvNormal(sol.u, Hermitian(gaia_post_covariance))
            break
        catch
            display(sol.original)
            @warn "Gaia model hessian matrix is not positive definite (maxlog=200)" maxlog=200 gaia_post_covariance
            for i in 1:length(sol.u)
                gaia_post_covariance[i] += e
            end
        end
    end
    if isnothing(P)
        return -Inf
    end
    Q = gaialike.catalog_distribution
    return -Distributions.kldivergence(P,Q)
    # return logpdf(Normal(0, 1), Distributions.kldivergence(P,Q)) # ???
end

function _simulate_skypath_observations(gaialike, sample_orbit::AbsoluteVisual, planet_mass_msol, T=Float64)
    # Compute the Ra and Dec sky path of the star, accounting for 3D barycentric
    # motion in spherical coordinates and perturbations from planets
    α_model = zeros(T, size(gaialike.table,1))
    δ_model = zeros(T, size(gaialike.table,1))
    for i in eachindex(gaialike.table.epoch)
        sol = orbitsolve(sample_orbit, gaialike.table.epoch[i])
        cmp = sol.compensated

        # Calculate the position of the star in cartesian coordinates
        # Get earth position in cartesian coordaintes (in units of AU)
        # Calculate apparent Ra and Dec.
        x_earth_pc = gaialike.table.x[i] / PlanetOrbits.pc2au
        y_earth_pc = gaialike.table.y[i] / PlanetOrbits.pc2au
        z_earth_pc = gaialike.table.z[i] / PlanetOrbits.pc2au
        x_diff_pc = cmp.x₂ - x_earth_pc
        y_diff_pc = cmp.y₂ - y_earth_pc
        z_diff_pc = cmp.z₂ - z_earth_pc
        distance_diff = sqrt(x_diff_pc^2 + y_diff_pc^2 + z_diff_pc^2)
        mydtor = π / 180
        ra_apparent_deg = ((atan(y_diff_pc,x_diff_pc)/mydtor + 360) % 360)
        arg = z_diff_pc / distance_diff
        arg = map(arg) do arg
            if 1.0 < arg < 1.0 + sqrt(eps(1.0))
                arg = 1.0
            end
            return arg
        end
        dec_apparent_deg = asin(arg) / mydtor

        # TODO: add perturbations from planets here
        for orb in (sample_orbit,)
            sol = orbitsolve(orb, gaialike.table.epoch[i])
            # Add perturbation from planet
            ra_apparent_deg += raoff(sol, planet_mass_msol)/60/60/1000/cos(sample_orbit.dec)
            dec_apparent_deg += decoff(sol, planet_mass_msol)/60/60/1000
        end

        α_model[i] = ra_apparent_deg
        δ_model[i] = dec_apparent_deg
    end

    # Now we need to generate a line (of any length) that passes through these points with an
    # angle given by psi. 

    # TODO: confirm we correctly accounted for a cos(delta) when generating these line segments
    αₘ = eachrow(@. T[-1, 1]'/60/60/1000/cos(sample_orbit.dec) * gaialike.table.sinϕ)
    for i in eachindex(α_model)
        αₘ[i] .+= α_model[i]
    end

    δₘ = eachrow(@. T[1, -1]'/60/60/1000 * gaialike.table.cosϕ)
    for i in eachindex(δ_model)
        δₘ[i] .+= δ_model[i]
    end

    return (;α_model,δ_model,αₘ,δₘ)

end


# This is an optimization objective function that gives the log-likelihood of Gaia skypath model
# matching some simulated data.
function _gaia_skypath_loglikelihood(args,(;αₘ,δₘ,table,ref_epoch,along_scan_uncertainty_mas))
    (
        plx,
        ra,
        dec,
        # rv,
        pmra,
        pmdec,
    ) = args
    inner_loop_skypath = AbsoluteVisual{KepOrbit}(;plx,ra,dec,rv=0,pmra,pmdec,
        ref_epoch=ref_epoch,
        M=1,a=1,e=0,i=0,ω=0,Ω=0,tp=0
    )
    ll = zero(eltype(args))
    for i in eachindex(table.epoch)
        sol = orbitsolve(inner_loop_skypath, table.epoch[i])
        cmp = sol.compensated
        x_earth_pc = table.x[i] / PlanetOrbits.pc2au
        y_earth_pc = table.y[i] / PlanetOrbits.pc2au
        z_earth_pc = table.z[i] / PlanetOrbits.pc2au
        x_diff_pc = cmp.x₂ - x_earth_pc
        y_diff_pc = cmp.y₂ - y_earth_pc
        z_diff_pc = cmp.z₂ - z_earth_pc
        distance_diff = sqrt(x_diff_pc^2 + y_diff_pc^2 + z_diff_pc^2)
        mydtor = π / 180
        ra_apparent_deg = ((atan(y_diff_pc,x_diff_pc)/mydtor + 360) % 360)
        arg = z_diff_pc / distance_diff
        arg = map(arg) do arg
            if 1.0 < arg < 1.0 + sqrt(eps(1.0))
                arg = 1.0
            end
            return arg
        end
        dec_apparent_deg = asin(arg) / mydtor
        # We don't add perturbations from planets here-- we are fitting Gaia's simplified
        # unperturbed model to our "true" perturbed data, to see if it reproduces the catalog
        # values
        point = @SVector [ra_apparent_deg, dec_apparent_deg]
        # Two points defining a line along which the star's position was measured
        line_point_1 = @SVector [αₘ[i][1],δₘ[i][1]]
        line_point_2 = @SVector [αₘ[i][2],δₘ[i][2]]
        # TODO: distance point to line should account for spherical coordinates!
        resid = distance_point_to_line(point, line_point_1, line_point_2) # degrees
        ll += logpdf(Normal(0,along_scan_uncertainty_mas/60/60/1000), resid)
        # if !isfinite(resid)
            # @error("non finite residual", resid, along_scan_uncertainty_mas, line_point_1[1], line_point_2[1])
            # error()
        # end
        # if !isfinite(ll)
            # @error("non finite ll", ll, along_scan_uncertainty_mas)
            # error()
        # end
    end
    return ll
end

# FOR testing at the moment
function _simulate_gaia_5_param_fit(args,(;αₘ,δₘ,table,ref_epoch,along_scan_uncertainty_mas))
    (
        plx,
        ra,
        dec,
        # rv,
        pmra,
        pmdec,
    ) = args
    inner_loop_skypath = AbsoluteVisual{KepOrbit}(;plx,ra,dec,rv=0,pmra,pmdec,
        ref_epoch=ref_epoch,
        M=1,a=1,e=0,i=0,ω=0,Ω=0,tp=0
    )
    display(inner_loop_skypath)
    ra_apparent_deg = zeros(length(table.epoch))
    dec_apparent_deg = zeros(length(table.epoch))
    resid = zeros(length(table.epoch))
    for i in eachindex(table.epoch)
        sol = orbitsolve(inner_loop_skypath, table.epoch[i])
        cmp = sol.compensated
        x_earth_pc = table.x[i] / PlanetOrbits.pc2au
        y_earth_pc = table.y[i] / PlanetOrbits.pc2au
        z_earth_pc = table.z[i] / PlanetOrbits.pc2au
        x_diff_pc = cmp.x₂ - x_earth_pc
        y_diff_pc = cmp.y₂ - y_earth_pc
        z_diff_pc = cmp.z₂ - z_earth_pc
        distance_diff = sqrt(x_diff_pc^2 + y_diff_pc^2 + z_diff_pc^2)
        mydtor = π / 180
        ra_apparent_deg[i] = ((atan(y_diff_pc,x_diff_pc)/mydtor + 360) % 360)
        arg = z_diff_pc / distance_diff
        arg = map(arg) do arg
            if 1.0 < arg < 1.0 + sqrt(eps(1.0))
                arg = 1.0
            end
            return arg
        end
        dec_apparent_deg[i] = asin(arg) / mydtor
        
        # We don't add perturbations from planets here-- we are fitting Gaia's simplified
        # unperturbed model to our "true" perturbed data, to see if it reproduces the catalog
        # values
        point = @SVector [ra_apparent_deg[i], dec_apparent_deg[i]]
        # Two points defining a line along which the star's position was measured
        line_point_1 = @SVector [αₘ[i][1],δₘ[i][1]]
        line_point_2 = @SVector [αₘ[i][2],δₘ[i][2]]
        # TODO: distance point to line should account for spherical coordinates!
        resid[i] = distance_point_to_line(point, line_point_1, line_point_2) # degrees

    end
    return ra_apparent_deg, dec_apparent_deg,resid
end


function _query_gaia_dr3(;gaia_id)
    fname = "_gaia_dr3/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiadr3.gaia_source WHERE source_id=$gaia_id"
            ],
            cookies=false,
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr3")
            mkdir("_gaia_dr3")
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end
    header_line, body_line = split(buf,"\n")
    headers = Symbol.(split(header_line,','))
    data = tryparse.(Float64, split(body_line,','))
    return namedtuple(headers, data)
end

function _query_gaia_dr2(;gaia_id)
    fname = "_gaia_dr2/source-$gaia_id.csv"
    if !isfile(fname)
        @info "Querying gea.esac.esa.int/tap-server" source_id=gaia_id
        resp = HTTP.get(
            "https://gea.esac.esa.int/tap-server/tap/sync",
            query=[
                "REQUEST"=>"doQuery",
                "LANG"=>"ADQL",
                "FORMAT"=>"CSV",
                "QUERY"=>"SELECT * FROM gaiadr2.gaia_source WHERE source_id=$gaia_id"
            ]
        )
        if resp.status != 200
            error("Error with GAIA query: $(resp.status)")
        end
        if !isdir("_gaia_dr2")
            mkdir("_gaia_dr2")
        end
        open(fname, write=true) do f
            write(f, resp.body)
        end
        buf = String(resp.body)
    else
        buf = read(fname, String)
    end
    header_line, body_line = split(buf,"\n")
    headers = Symbol.(split(header_line,','))
    data = tryparse.(Float64, split(body_line,','))
    return namedtuple(headers, data)
end

# TODO: can we add DR1 also? The challenge is that it has different source_id for the same object


"""
    geocentre_position_query(epoch_MJD)

Given a date+time in MJD format, return a named tuple of Earth position and velocity in AU 
on that date. The results are cached in a local `_geocentre_pos` directory for offline use.

Currently these are queried from the NASA HORIZONS online system, but this detail is subject
to change in future minor versions of the package.

The positions and velocities represent the Geocenter of the Earth relative to the solar sysytem
barycenter.
"""
function geocentre_position_query(epoch_MJD::Number)
    
    fname = "_geocentre_pos/HORIZONS-Geocenter-$(epoch_MJD)mjd.txt"
    if !isfile(fname)
        if !isdir("_geocentre_pos")
            mkdir("_geocentre_pos")
        end
        # Query the HORIZONS system for the Earth-Moon barycentre position at each specific epoch
        # of the dataset. 
        # This incurs a separate query for each epoch.
        t = DateTime(mjd2date(epoch_MJD))
        @info "Querying Earth Geocentre position from HORIZONS" epoch_MJD date=t
        HORIZONS.vec_tbl("Geocenter", t, t + Hour(1), Day(1); FILENAME=fname, CENTER="@ssb", REF_PLANE="FRAME", OUT_UNITS="AU-D", CSV_FORMAT=true, VEC_TABLE=2)
    end
        
    lines = readlines(fname)
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



"""
    forecast_table = GHOST_forecast(gaia_id=1234)

Given a Gaia ID, retreive a forecast of Gaia observations from the GHOST tool automatically.
See tool URL here: https://gaia.esac.esa.int/gost/

Please be aware that others  might be able to discover the target coordinates you searched for
(though not who performed the search) via information leaked to the external service.
"""
function GHOST_forecast(;gaia_id,  catalog=(datadep"HGCA_eDR3") * "/HGCA_vEDR3.fits")

    fname = "GHOST-$gaia_id.csv"
    if isfile(fname)
        forecast_table = CSV.read(fname, Table)
        return forecast_table
    end


    # Load the Hipparcos-GAIA catalog of accelerations as a table
    hgca = FITS(catalog, "r") do fits
        Table(fits[2])
    end
    idx = findfirst(==(gaia_id), hgca.gaia_source_id)

    #=
    curl \
    --form "srcname=eps indi test" \
    --form "inputmode=single" \
    --form "srcra=330.84022344234336" \
    --form "srdec=-56.785978554364995" \
    --form "from=2014-07-26T00:00:00" \
    --form "to=2017-05-28T00:00:00" \
    --cookie "JSESSIONID=73EE915E5F9FF8B5376D82FC116F765D" \
    'https://gaia.esac.esa.int/gost/GostServlet' 

    curl --cookie "JSESSIONID=73EE915E5F9FF8B5376D82FC116F765D" 'https://gaia.esac.esa.int/gost/export.jsp?id=73EE915E5F9FF8B5376D82FC116F765D/1722634273&format=csv'
    =#

    # Just pick a cookie ID to use. 
    # Might be better to let the service create one for us.

    cookiejar = HTTP.CookieJar()
    @info "Contacting the GAIA scan forecast tool GHOST: https://gaia.esac.esa.int/gost/"
    resp0 = HTTP.get(
        "https://gaia.esac.esa.int/gost/GostServlet",
        cookiejar=cookiejar,
    )
    if resp0.status != 200
        println(String(resp0.body))
        error("Could not contact the GAIA scan forecast tool GHOST https://gaia.esac.esa.int See above error message.")
    end
    formdata = Dict([
        "inputmode"=>"single",
        "srcra" => string(hgca.gaia_ra[idx]),
        "srdec" => string(hgca.gaia_dec[idx]),
        "from" => "2014-07-26T00:00:00",
        "to" => "2017-05-28T00:00:00",
    ])
    @info "Retrieving forecasted GAIA scans from GHOST: https://gaia.esac.esa.int/gost/"
    resp = HTTP.post(
        "https://gaia.esac.esa.int/gost/GostServlet",
        body=HTTP.Form(formdata),
        cookies=Dict([
            # "Cookie" => "$new_cookie",
            "JSESSIONID" => "73EE915E5F9FF8B5376D82FC116F765D",
            "noagreement" => "false"
        ]),
        cookiejar=HTTP.CookieJar()#cookiejar
    )
    if resp.status != 200 || contains(String(collect(resp.body)),"error")
        println(String(resp.body))
        error("Could not fetch GAIA scan forecast from GHOST. See above error message. Do you have an internet connection available?")
    end

    m = match(r"Submitted with id (\d+)", String(resp.body))
    response_id = m.captures[1]
    url = "https://gaia.esac.esa.int/gost/export.jsp?id=73EE915E5F9FF8B5376D82FC116F765D/$response_id&format=csv"
    resp_dat = HTTP.get(
        url,
        cookies=Dict([
            # "Cookie" => "$new_cookie",
            "JSESSIONID" => "73EE915E5F9FF8B5376D82FC116F765D",
            "noagreement" => "false"
        ]),
        cookiejar=HTTP.CookieJar()
    )
    @info "done"
    body = collect(resp_dat.body)
    
    # Save for offline use eg on clusters
    write(fname, body)

    io = IOBuffer(body)
    forecast_table = CSV.read(io, Table)

    return forecast_table
end




