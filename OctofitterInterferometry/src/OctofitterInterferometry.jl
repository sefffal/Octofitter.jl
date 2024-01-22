module OctofitterInterferometry

using Octofitter
using PlanetOrbits
using Tables, TypedTables

using FITSIO, OIFITS

const required_cols = (:epoch, :u, :v, :cps_data, :dcps, :vis2_data, :dvis2, :index_cps1, :index_cps2, :index_cps3,:band,:use_vis2)
struct InterferometryLikelihood{TTable<:Table} <: Octofitter.AbstractLikelihood
    table::TTable
    function InterferometryLikelihood(observations...)
        input_table = Table(observations...) 
        if :filename ∈ Tables.columnnames(input_table)
            rows = map(eachrow(input_table)) do row
                row = only(row)
                f = FITS(row.filename, "r")     # open FITS file for reading
                wavs = read(OIDataBlock, f[2]) # idem
                vis2s = read(OIDataBlock, f[5]) # idem
                cps = read(OIDataBlock, f[6]) # idem
                #read data
                wav = wavs.eff_wave
                vis2 = vis2s.vis2data
                vis2_err = vis2s.vis2err
                ut = vis2s.ucoord
                vt = vis2s.vcoord
                vis2_index = vis2s.sta_index
                cp = cps.t3phi
                cp_err = cps.t3phierr
                cp_index = cps.sta_index
                #convert u,v to units of wavelength
                u = ut/wav
                v = vt/wav
                cp_inds1,cp_inds2,cp_inds3 = cp_indices(vis2_index=vis2_index,cp_index=cp_index)
                return (;
                    row.epoch,
                    row.band,
                    row.use_vis2,
                    u=u[:,1],
                    v=v[:,1],
                    cps_data=cp[1,:],
                    dcps=cp_err[1,:], 
                    vis2_data=vis2[1,:],
                    dvis2=vis2_err[1,:],
                    index_cps1=cp_inds1,
                    index_cps2=cp_inds2,
                    index_cps3=cp_inds3,
                )
            end
            table = Table(rows)
        else
            table = input_table
        end
        
        if !issubset(required_cols, Tables.columnnames(table))
            error("Expected columns $vis_cols")
        end
        return new{typeof(table)}(table)
    end
end
InterferometryLikelihood(observations::NamedTuple...) = InterferometryLikelihood(observations)
export InterferometryLikelihood

"""
Visibliitiy modelling likelihood for point sources.
"""
function Octofitter.ln_like(vis::InterferometryLikelihood, θ_system, orbits, num_epochs::Val{L}=Val(length(vis.table))) where L

    T = typeof(θ_system.M)
    ll = zero(T)

    # Access the data here: 
    epochs = vis.table.epoch
    band = vis.table.band
    u = vis.table.u
    v = vis.table.v
    cps_data = vis.table.cps_data
    dcps = vis.table.dcps
    vis2_data = vis.table.vis2_data
    dvis2 = vis.table.dvis2
    i_cps1 = vis.table.index_cps1
    i_cps2 = vis.table.index_cps2
    i_cps3 = vis.table.index_cps3
    use_vis2 = vis.table.use_vis2

    # Loop through planets
    for j in eachindex(epochs)
        
        epoch = epochs[j]
        this_band = band[j]
        
        cvis = 0.
        norm_factor = 0. #to normalize complex visibilities 
        for i in eachindex(orbits)
            # All parameters relevant to this planet
            θ_planet = θ_system.planets[i]
            
            # Get model contrast parameter in this band (band provided as a symbol, e.g. :L along with data in table row.)
            contrast = getproperty(θ_planet, this_band) #secondary/primary
            #contrast = contrast[i] 
            
            # orbit object pre-created from above parameters (shared between all likelihood functions)
            orbit = orbits[i]
            sol = orbitsolve(orbit, epoch)
            Δra  = raoff(sol)  # in mas
            Δdec = decoff(sol) # in mas
            
            #add complex visibilities from all planets at a single epoch
            cvis = cvis .+ cvis_bin(ddec=Δdec,dra=Δra,contrast=contrast,u=u[j],v=v[j])

            norm_factor += contrast
        end
        cvis = cvis .+ 1. #add contribution from the primary primary
        cvis *= 1. /(1. +norm_factor)
        #compute closure phase
        cps = closurephase(vis=cvis,index_cps1=i_cps1[j],index_cps2=i_cps2[j],index_cps3=i_cps3[j])
        lnlike_v2 = 0.
        if use_vis2[j]
            #compute squared visibilities
            vis2 = abs.(cvis).^2
            const_v2 = -sum(log.(2*pi*(dvis2[j].^2)))/2
            #calculate vis2 ln likelihood
            lnlike_v2 = lnlike_v2 .+ -0.5*sum((vis2_data[j].- vis2).^2 ./dvis2[j].^2) .+ const_v2
        end
        
        #calculate cp ln likelihood
        const_cp = -sum(log.(2*pi*(dcps[j].^2)))/2
        lnlike_cp = -0.5*sum((cps_data[j] .- (cps)).^2 ./dcps[j].^2) .+ const_cp

        # Accumulate into likelihood
        ll = ll .+ lnlike_cp .+ lnlike_v2
    end


    return ll
end

    
function cvis_bin(;ddec,dra,contrast,u,v)
    #u,v: baselines [wavelengths]
    #ddec: dec offset [mas]
    #dra: ra offset [mas]
    #contrast: secondary/primary  contrast 
    ##################################
    #returns complex visibilities of a point source at position ddec,dra
    
    l2 = contrast 
    # phase-factor
    phi= cos.(-2*pi*(u*dra + v*ddec)*pi/(180*3600*1000)) .+ sin.(-2*pi*(u*dra + v*ddec)*pi/(180*3600*1000))im
    cvis =  l2 * phi
    return cvis
end


function closurephase(;vis::AbstractVector,index_cps1::AbstractVector,index_cps2::AbstractVector,index_cps3::AbstractVector)
    #vis: complex visibilities
    #i_cps1,i_cps2,i_cps3: closure triangle indices 
    ##################################
    #returns closure phases [degrees]
    
    realt = real(vis)
    imagt = imag(vis)
    visphi = atan.(imagt,realt)*180/pi #convert to degrees
    visphi = mod.(visphi .+ 10980., 360.) .- 180. #phase in range (-180,180]
    cp = visphi[index_cps1] .+ visphi[index_cps2] .- visphi[index_cps3]
    return cp
end

function cp_indices(;vis2_index::Matrix{<:Int64}, cp_index::Matrix{<:Int64})
    """Extracts indices for calculating closure phases from visibility and closure phase station indices"""
    i_cps1 = zeros(Int64,size(cp_index)[2])
    i_cps2 = zeros(Int64,size(cp_index)[2])
    i_cps3 = zeros(Int64,size(cp_index)[2])
    
    nh = maximum(vis2_index) #number of stations making up your interferometer
    nb = Int64(nh*(nh-1)/2) #number of baselines
    ncp = Int64(nh*(nh-1)*(nh-2)/6) #total number of closure phases

    for i in range(1,size(cp_index)[2])
        for j in range(1,size(vis2_index)[2])
            if ((cp_index[1,i]==vis2_index[1,j])&&(cp_index[2,i]==vis2_index[2,j]))
                if (floor((j-1)/nb) == floor((i-1)/ncp)) 
                    i_cps1[i] = j
                end
            end
            if ((cp_index[2,i]==vis2_index[1,j])&&(cp_index[3,i]==vis2_index[2,j]))
                if (floor((j-1)/nb) == floor((i-1)/ncp))
                    i_cps2[i] = j
                end
            end
            if ((cp_index[1,i]==vis2_index[1,j])&&(cp_index[3,i]==vis2_index[2,j]))
                if (floor((j-1)/nb) == floor((i-1)/ncp))
                    i_cps3[i] = j
                end
           end
        end
    end
    return i_cps1,i_cps2,i_cps3
end

# Generate new observations for a system of possibly multiple planets
function Octofitter.generate_from_params(like::InterferometryLikelihood, θ_system, orbits::Vector{<:AbstractOrbit})

    # # Get epochs, uncertainties, and planet masses from observations and parameters
    epochs = like.table.epoch
    bands = like.table.band
    u = like.table.u
    v = like.table.v
    cps_data = like.table.cps_data
    dcps = like.table.dcps
    vis2_data = like.table.vis2_data
    dvis2 = like.table.dvis2
    i_cps1 = like.table.index_cps1
    i_cps2 = like.table.index_cps2
    i_cps3 = like.table.index_cps3
    use_vis2 = like.table.use_vis2
    cp_all = Any[];
    vis2_all = Any[];
    for j in eachindex(epochs)
        band = bands[j]
        epoch = epochs[j]
        cvis = 0.
        norm_factor = 0. #to normalize complex visibilities 
        for i in eachindex(orbits)
            # All parameters relevant to this planet
            θ_planet = θ_system.planets[i]

            # orbit object pre-created from above parameters (shared between all likelihood functions)
            orbit = orbits[i]
            contrast = getproperty(θ_planet, band) #secondary/primary
            sol = orbitsolve(orbit, epoch)
            Δra  = raoff(sol)  # in mas
            Δdec = decoff(sol) # in mas
            cvis = cvis .+ cvis_bin(ddec=Δdec,dra=Δra,contrast=contrast,u=u[j],v=v[j])
            norm_factor += contrast
        end
        cvis = cvis .+ 1. #add contribution from the primary 
        cvis *= 1. /(1. +norm_factor)
        #compute closure phase
        cp = closurephase(vis=cvis,index_cps1=i_cps1[j],index_cps2=i_cps2[j],index_cps3=i_cps3[j])
        #compute squared visibilities
        vis2 = abs.(cvis).^2
        
        cp_all = append!(cp_all,[cp])
        vis2_all = append!(vis2_all,[vis2])

    end
    new_vis_like_table = Table(epoch=epochs, u=u, v=v, cps_data=cp_all, dcps=dcps,
        vis2_data=vis2_all, dvis2=dvis2, index_cps1=i_cps1,
        index_cps2=i_cps2,index_cps3=i_cps3,band=bands,use_vis2=use_vis2)
    
    
    # return with same number of rows: band, epoch
    # position(s) of point sources according to orbits, θ_system
    return InterferometryLikelihood(new_vis_like_table)
end


end
