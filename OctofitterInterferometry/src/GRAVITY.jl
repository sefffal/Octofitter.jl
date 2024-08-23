using LinearAlgebra
using Interpolations
using BlockArrays
using Distributions
using PDMats

#=
This file implements a variant of InterferometryLikelihood
that accounts for Fiber positioning throughput loss.
=#
const required_cols_grav_wide2 = (:epoch, :u, :v, :cps_data, :dcps, :index_cps1, :index_cps2, :index_cps3, :spectrum_var,)

include("GRAVITY-correlation.jl")

struct GRAVITYWideKPLikelihood{TTable<:Table,TInterp,spectrum_vars} <: AbstractInterferometryLikelihood
    table::TTable
    fiber_coupling_interpolator::TInterp
end
function GRAVITYWideKPLikelihood(
    observations...;
)
    input_table = Table(observations...)
    if :filename ∈ Tables.columnnames(input_table)
        rows = map(eachrow(input_table)) do row
            _prepare_input_row(row)
        end
        table = Table(rows)
    else
        table = input_table
    end

    if !issubset(required_cols_grav_wide2, Tables.columnnames(table))
        error("Expected columns $required_cols_grav_wide2, got $(Tables.columnnames(table))")
    end

    # Take some additional preparation steps per input-row.
    rows_with_kernel_phases = map(eachrow(table)) do row′
        row = row′[1]

        # Calculate the design matrix
        # TODO: replace hardcoded T with one calculated using cps_index1, 2, and 
        T = Int8[
            1 -1 0 1 0 0
            1 0 -1 0 1 0
            0 1 -1 0 0 1
            0 0 0 1 -1 1
        ]

        # We now generate a unified design matrix that works with all wavelengths.
        # Later, this will allow us to easily deal with correlations between both CPs & wavelengths.
        Λ = length(row.eff_wave)
        Tλ = zeros(Int8, Λ * size(T, 1), Λ * size(T, 2))

        # We now replicate our deign matrix T across wavelengths.
        # We will put each wavelength together, grouped by baseline.

        for baseline_i in axes(T, 1), baseline_j in axes(T, 2)
            for wavelength_i in 1:Λ
                Tλ[wavelength_i+(baseline_i-1)*Λ, wavelength_i+(baseline_j-1)*Λ] = T[baseline_i, baseline_j]
            end
        end

        # Determine a kernel phase basis using the cholesky factorization of the
        # design matrix * design matrixᵀ
        C, U=cholesky(Tλ*Tλ')
        P₁ = collect(C) ./ sqrt.(diag(C*C'))
        i_max = findfirst(<=(1e-5), diag(P₁))-1
        P₁ = P₁[:,1:i_max]'

        # We can pre-convert the CP uncertainties into KP uncertainties
        σ_kp =  P₁ * vec(row.dcps);

        return (; row..., Tλ, P₁, σ_kp)
    end
    table = Table(rows_with_kernel_phases)

    # Create an interpolator object that maps separation to fiber coupling efficiency
    @info "Pre-calculating fiber coupling efficiency over grid"
    sep_mas = 0:2:100
    λs = range(extrema(vec(table.eff_wave[1]))..., length=15)
    fiber_coupling = stack([
        fiber_coupling_fraction(sep_mas, λ)
        for λ in λs
    ])
    coupling_interp = LinearInterpolation((sep_mas, λs), fiber_coupling, extrapolation_bc=0.0)

    spectrum_vars = tuple(Symbol.(vec(table.spectrum_var))...)
    return GRAVITYWideKPLikelihood{typeof(table),typeof(coupling_interp),spectrum_vars}(table, coupling_interp)
end
GRAVITYWideKPLikelihood(observations::NamedTuple...) = GRAVITYWideKPLikelihood(observations)
export GRAVITYWideKPLikelihood


function _getparams(::GRAVITYWideKPLikelihood{TTable,TInterp,spectrum_vars}, θ_system) where
    {TTable,TInterp,spectrum_vars}
    spectrum_vals = tuple((
        getproperty(θ_system.planets[i_planet], spectrum_vars[i_planet]) for i_planet in eachindex(keys(θ_system.planets))
    )...)
    return (;spectrum_vals)
end

# theta = range(0, 100, length=250)
# Credit: W. Balmer, D. Bakely, and others.
function fiber_coupling_fraction(theta, lambda_w=2.2e-6)
    D = 8
    x = range(-D * 2, D * 2, length=500)#0)
    y = range(-D * 2, D * 2, length=500)#0)
    r = LinearAlgebra.norm.(x, y')
    m = r .< D / 2
    # arcseconds
    phase = reshape(x, :, 1, 1) ./ lambda_w .* reshape(theta, 1, 1, :) * 1e-3 / (180 / pi * 3600) * 2pi
    w_0 = 0.32D
    field_pup = @. m * exp(1im * phase)
    field_fiber = @. exp(-1 * r^2 / (2 * w_0^2))
    Inj = abs2.(sum(field_pup .* field_fiber, dims=(1, 2))) / abs(sum(m .* field_fiber))^2
    return Inj[1, 1, :]
end



"""
Visibliitiy modelling likelihood for point sources.
"""
function Octofitter.ln_like(vis::GRAVITYWideKPLikelihood, θ_system, orbits, num_epochs::Val{L}=Val(length(vis.table))) where {L}

    # Convoluted way to get either Float64 normally or a Dual{Float64} if using ForwardDiff
    T = Octofitter._system_number_type(θ_system)
    ll = zero(T)

    (;spectrum_vals) = _getparams(vis, θ_system)
    # contrasts = stack(spectrum_vals) # this could be more optimized

    # Assume that the fiber is positioned at the photocentre, averaged over wavelengths
    mean_constrasts_by_planet = zeros(T,length(spectrum_vals))
    for i_planet in 1:length(spectrum_vals)
        mean_constrasts_by_planet[i_planet] = mean(spectrum_vals[i_planet])
    end

    # Access the data here: 
    epochs = vis.table.epoch

    if length(epochs) > 0
        cps_model = zeros(T, size(vis.table.cps_data[1][:, 1]))
        cvis_model = zeros(complex(T), size(vis.table.u[1][:, 1]))
    end

    # Loop through epochs
    for i_epoch in eachindex(epochs)
        epoch = epochs[i_epoch]

        index_cps1 = vis.table.index_cps1[i_epoch]
        index_cps2 = vis.table.index_cps2[i_epoch]
        index_cps3 = vis.table.index_cps3[i_epoch]
        # use_vis2 = vis.table.use_vis2[i_epoch]

        # Re-use buffers between iterations if they are all the same shape (typical)
        if size(cps_model) == size(vis.table.cps_data[i_epoch][:, 1]) && size(cvis_model) == size(vis.table.u[i_epoch][:, 1])
            cps_model .= 0
            cvis_model .= 0
        else
            cps_model = zeros(T, size(vis.table.cps_data[i_epoch][:, 1]))
            cvis_model = zeros(complex(T), size(vis.table.u[i_epoch][:, 1]))
        end

        sols = [orbitsolve(orbits[i_planet], epoch) for i_planet in 1:min(length(θ_system.planets), length(orbits))]

        throughputs = zeros(T, length(sols), length(vis.table.eff_wave[i_epoch]))
        for i_planet in 1:length(sols)
            for i_wave in 1:length(vis.table.eff_wave[i_epoch])
                sol = sols[i_planet]
                flux_ratio = mean_constrasts_by_planet[i_planet]
                wavelength_m = vis.table.eff_wave[i_epoch][i_wave]
                # Model the fiber as placed at the photocentre of the two bodies
                secondary_offset_mas = projectedseparation(sol)
                # Now calculate throughput loss on the secondary due to it being offset wrt. the 
                # fiber (assumed to be at photocentre)
                fiber_offset_mas = (flux_ratio * secondary_offset_mas) / (1.0 + flux_ratio)
                coupling = vis.fiber_coupling_interpolator(fiber_offset_mas, wavelength_m)
                throughputs[i_planet,i_wave] = coupling
            end
        end

        Λ = length(vis.table.eff_wave[i_epoch])
        Len = Λ * size(vis.table.cps_data[i_epoch], 1)
        cp_resids = zeros(T, Len) # Fill this in a moment

        # Loop through wavelengths
        # The following is NOT threadsafe. DON'T multithread it!
        for i_wave in axes(vis.table.u[i_epoch], 2)
            u = @views vis.table.u[i_epoch][:, i_wave]
            v = @views vis.table.v[i_epoch][:, i_wave]
            cps_data = @views vis.table.cps_data[i_epoch][:, i_wave]
            σ_cp = @views vis.table.dcps[i_epoch][:, i_wave]
            # vis2_data = @views vis.table.vis2_data[i_epoch][:, i_wave]
            # dvis2 = @views vis.table.dvis2[i_epoch][:, i_wave]

            # to normalize complex visibilities 
            cvis_model .= 0
            cps_model .= 0
            norm_factor_model = zero(T)

            # Consider all planets
            for i_planet in eachindex(orbits)
                # All parameters relevant to this planet
                # Get model contrast parameter in this band (band provided as a symbol, e.g. :L along with data in table row.)
                contrast = spectrum_vals[i_planet][i_wave]
                throughput = throughputs[i_planet]
                Δra = raoff(sols[i_planet])  # in mas
                Δdec = decoff(sols[i_planet]) # in mas

                # add complex visibilities from all planets at a single epoch, for this wavelength
                cvis_bin!(cvis_model; Δdec, Δra, contrast=contrast * throughput, u, v)
                norm_factor_model += contrast
            end
            cvis_model .+= 1.0 # add contribution from the primary primary
            cvis_model .*= 1.0 / (1.0 + norm_factor_model)

            # Compute closure phases
            closurephase!(cps_model; vis=cvis_model, index_cps1, index_cps2, index_cps3)

            for i_T3 in eachindex(σ_cp, cps_data, cps_model)
                cp_resids[(i_T3-1)*Λ+i_wave] = cps_data[i_T3] - cps_model[i_T3]
                # if !isfinite(cp_resids[(i_T3-1)*Λ+i_wave] )
                #     @warn "non finite CP calculated"  cps_data[i_T3] cps_model[i_T3]
                # end
            end
        end
        # Done calculating the residuals for this epoch

        # Niave CP Only:
        # σ_cp = vec(vis.table[i_epoch].dcps) #sqrt.(σ_cp_jitter .^ 2 .+ vec(vis.table[i_epoch].dcps) .^ 2)
        # distribution = MvNormal(Diagonal(σ_cp))
        # ll += logpdf(distribution, cp_resids)

        # Kernel Phase likelihood with Jens's semi-analytic correlation matrix
        # This is 7x slower than naive CP modelling, but this implementation is 
        # still 5.5x faster than the straightforward approach.
        P₁ = vis.table.P₁[i_epoch]
        # σ_cp = vis.table[i_epoch].dcps

        ## Diagonalized covariance Kernphases
        if hasproperty(vis.table, :jitter)
            kp_jitter_name =  vis.table.jitter[i_epoch]
            kp_jitter = convert(T, getproperty(θ_system, kp_jitter_name))
            kp_jitter = max(eps(),kp_jitter)
        else
            kp_jitter = zero(T)
        end

        # Generate the semi-analytic correlation matrix from Jens
        # CT3_y = hasproperty(θ_system, :CT3_y) ? float(θ_system.CT3_y) : zero(T)
        if hasproperty(vis.table, :kp_Cy)
            kp_Cy =  vis.table.kp_Cy[i_epoch]
            kp_Cy = convert(T, getproperty(θ_system, kp_Cy))
            kp_Cy = max(eps(),kp_Cy)
        else
            kp_Cy = zero(T)
        end

        # We directly compute the kernel phase analog to that correlation matrix (digonal + spectral correlation within a given KP)
        
        # TODO: this might only be correct up to a scaling factor, verify
        C_kp = CKP(vis.table[i_epoch], kp_Cy)

        # Caculate KP uncertainties from CPs
        # σ_kp =  P₁ * vec(σ_cp); 
        # We have factored this out into the constructor to avoid doing it each iteration:
        σ_kp = vis.table[i_epoch].σ_kp

        # Calculate covariance from correlation
        Σ_kp = (Diagonal(σ_kp) * C_kp * Diagonal(σ_kp)' )
        # add KP jitter along the diagonal
        for i in axes(Σ_kp,1)
            Σ_kp[i,i] += kp_jitter^2
        end

        # Convert the CP residuals to KPs
        kernphase_resids = P₁ * vec(cp_resids)

        # We now exploit the block diagonal structure of the KP covariance matrix
        # (we have ensured this by construction). We can calculate the cholesky 
        # factorization for each block independently and then concatenate them together.
        # We know there are three KPs per wavelength with GRAVITY.
        N = length(vis.table.eff_wave)
        cholesky!(Hermitian(@view Σ_kp[1:N,1:N]));
        cholesky!(Hermitian(@view Σ_kp[N+1:2N,N+1:2N]));
        cholesky!(Hermitian(@view Σ_kp[2N+1:end,2N+1:end]));
        
        # Now we construct an MvNormal d1;istribution but avoid 
        # calling Cholesky on the full matrix
        Ch = Cholesky(UpperTriangular(Σ_kp));
        # need to pass matrix and cholesky factorization
        P = PDMat(Σ_kp, Ch);
        dist = MvNormal(P);
        ll += logpdf(dist,kernphase_resids)
    end


    return ll
end
