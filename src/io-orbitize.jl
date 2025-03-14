"""
This file contains functions for importing and exporting
chains in a format compatible with the python package Orbitize!.

It also contains functions for querying posteriors and datasets
from the public whereistheplanet.com website.
"""

using HDF5
using StringDistances
using MCMCChains: MCMCChains

# public loadhdf5, savehdf5, Whereistheplanet_search, Whereistheplanet_posterior, Whereistheplanet_astrom

"""
Whereistheplanet_search("targetname")

Search for an orbit posterior and/or astrometry hosted on whereistheplanet.com by a given target name.
If not found, a list of similar target names will be reported.
"""
function Whereistheplanet_search(target, catalog=datadep"Whereistheplanet")

    dirpath = joinpath(catalog, "whereistheplanet-master", "data")
    fnames = readdir(dirpath,join=true)
    avail_targets = map(fnames) do fname
        m = match(r"post_(.+)\.hdf5", fname)
        return !isnothing(m) ? m.captures[1] : ""
    end
    fname_matched_i = findfirst(==(target), avail_targets)

    if isnothing(fname_matched_i)
        avail_filt = avail_targets[avail_targets.!=""] 
        similarity = evaluate.(Ref(Levenshtein()), target, avail_filt)
        ii = sortperm(similarity)
        closest_3 = avail_filt[ii[1:min(3,end)]]
        @error "No results were found for the target $target."
        @info  "Here are a list of similar and available target names" closest_3
        error()
    end

    return fnames[fname_matched_i]

end

"""
    Whereistheplanet_astrom("targetname")

Load astrometry hosted on whereistheplanet.com by a given target name.
If not found, a list of similar target names will be reported.
"""
function Whereistheplanet_astrom(target, catalog=datadep"Whereistheplanet"; object=1)

    fname = Whereistheplanet_search(target, catalog)
    return h5open(fname, "r") do f
        records = read(f["data"])

        records = filter(row->row.object==object, records)

        # Group observations by type
        seppa = filter(row->row.quant_type=="seppa", records)
        radec = filter(row->row.quant_type=="radec", records)

        out = PlanetRelAstromLikelihood[]
        if length(seppa) > 0
            astrom_seppa = PlanetRelAstromLikelihood(map(seppa) do row
                cor=row.quant12_corr
                if !isfinite(cor)
                    cor = 0.0
                end
                (;row.epoch, sep=row.quant1, σ_sep=row.quant1_err, pa=deg2rad(row.quant2), σ_pa=deg2rad(row.quant2_err), cor)
            end)
            push!(out, astrom_seppa)
        end
        if length(radec) > 0
            astrom_radec = PlanetRelAstromLikelihood(map(radec) do row
                cor=row.quant12_corr
                if !isfinite(cor)
                    cor = 0.0
                end
                (;row.epoch, ra=row.quant1, σ_ra=row.quant1_err, dec=row.quant2, σ_dec=row.quant2_err, cor)
            end)
            push!(out, astrom_radec)
        end

        # return(out...)
    end
end

"""
Load an Orbitize! posterior from an HDF5 file and convert it into
an Octofitter-compatible chains format.
Both tools use the same orbit conventions so this is fairly straightforward.

If you pass `numchains` as a second argument, the array will be interpretted
as coming from multiple chains concatenated together, 
"""
function loadhdf5(fname_or_targetname, numchains=1; colnames=nothing)
    if !(occursin(".hdf5", fname_or_targetname) || occursin(".h5", fname_or_targetname))
        fname = Whereistheplanet_search(fname_or_targetname)
    else
        fname = fname_or_targetname
    end
    return h5open(fname, "r") do f
        # Standard orbitize basis assumed:
        # semi-major axis (sma), eccentricity (ecc), inclination (inc), argument of periastron (aop), position angle of the nodes (pan), epoch of periastron expressed as a fraction of the period past a reference epoch (tau), parallax (plx) and total system mass (mtot)
        arr = transpose(read(f["post"]))
        if numchains > 1
            lenchain = size(arr,1)÷numchains
            arr = stack(
                (
                    arr[(chain_no-1)*lenchain+1:chain_no*lenchain,:]
                    for chain_no in 1:numchains
                ),
            )
        end
        if isnothing(colnames)
            try
                colnames = read_attribute(f, "parameter_labels")
            catch
                @warn "`parameter_labels` not present, will have to fall back on a guess. You can also provide a vector of colname strings with the argument `colnames`"
                colnames = ["sma1","ecc1","inc1","aop1","pan1","tau1","plx","M"]
            end
        end
        num_planets = 0
        num_planets += "sma1" ∈ colnames
        num_planets += "sma2" ∈ colnames
        num_planets += "sma3" ∈ colnames
        num_planets += "sma4" ∈ colnames
        planet_keys = ["b", "c", "d", "e"]
        replace!(colnames,
            "sma1"=>"b_a",
            "ecc1"=>"b_e",
            "inc1"=>"b_i",
            "aop1"=>"b_ω",
            "pan1"=>"b_Ω",
            "tau1"=>"b_τ",
            "sma2"=>"c_a",
            "ecc2"=>"c_e",
            "inc2"=>"c_i",
            "aop2"=>"c_ω",
            "pan2"=>"c_Ω",
            "tau2"=>"c_τ",
            "sma3"=>"d_a",
            "ecc3"=>"d_e",
            "inc3"=>"d_i",
            "aop3"=>"d_ω",
            "pan3"=>"d_Ω",
            "tau3"=>"d_τ",
            "sma4"=>"e_a",
            "ecc4"=>"e_e",
            "inc4"=>"e_i",
            "aop4"=>"e_ω",
            "pan4"=>"e_Ω",
            "tau4"=>"e_τ",
            # "plx"
            "m1"=>"b_mass",
            "m2"=>"c_mass",
            "m3"=>"d_mass",
            "m4"=>"e_mass",
            "m0"=>"M_pri",
            "mtot"=>"M",
        )


        # We need to add total mass columns planet_M for each planet that
        # include the mass of each body and those interior to it
        for i_planet in 1:num_planets

            ind = findfirst(==("M_pri"), colnames)
            if isnothing(ind) 
                # We just have mtot
                continue
            end
            M_tot = arr[:,ind]
            for j_planet in 1:num_planets
                # Only add influce if sma is less than current
                sma_this = arr[:,findfirst(==(string(planet_keys[i_planet], "_a")),colnames)]
                sma_other = arr[:,findfirst(==(string(planet_keys[j_planet], "_a")),colnames)]
                mask_sma_lower =  sma_other .<= sma_this
                # mask_sma_lower = sma_other .== sma_this
                m_planet =  arr[:,findfirst(==(string(planet_keys[j_planet], "_mass")),colnames)]
                M_tot .+= mask_sma_lower .* m_planet
            end
            arr = hcat(arr, M_tot)
            push!(colnames, string(planet_keys[i_planet], "_M"))
        end

        for (col, colname) in zip(eachcol(arr), colnames)
            if endswith(colname, "_mass")
                col ./= mjup2msol
            end
        end


        chn =  MCMCChains.Chains(arr, Symbol.(colnames))
        # Calculate epoch of periastron passage from Orbitize tau variable:
        tau_ref_epoch = 58849
        if haskey(attrs(f),"tau_ref_epoch")
            tau_ref_epoch = attrs(f)["tau_ref_epoch"]
        end
        tps = map(1:num_planets) do i_planet
            a = chn[Symbol("$(planet_keys[i_planet])_a")]
            τ = chn[Symbol("$(planet_keys[i_planet])_τ")]
            k1 = :M
            k2 = Symbol("$(planet_keys[i_planet])_M")
            if haskey(chn,k1)
                M = chn[k1]
            else
                M = chn[k2]
            end
            period_days = @. √(a^3/M)*PlanetOrbits.kepler_year_to_julian_day_conversion_factor
            tp = @. τ * period_days + tau_ref_epoch
            return tp
        end
        tp_keys = map(1:num_planets) do i_planet
            Symbol("$(planet_keys[i_planet])_tp")
        end

        chn =  MCMCChains.Chains(hcat(arr, tps...), [Symbol.(colnames)..., tp_keys...])

        # Read additional attributes in and convert to named tuple
        metadata = [Symbol(k)=>v for (k,v) in attrs(f)]
        return MCMCChains.setinfo(chn, NamedTuple(metadata))
    end
end



"""
    savehdf5("filename.hdf5", chain)

For some limited cases, supports saving an Octofitter chain in a format used by Orbitize!
and used by whereistheplanet.com.

Currently only works with a single planet model. Does not currently export any raw data.
"""
function savehdf5(fname::AbstractString, model::Octofitter.LogDensityModel, chain::Chains, planet_key::Symbol=first(keys(model.system.planets)))
    return h5open(fname, "w") do f
        #  Just look up by planet key

        # Must calculate tau
        tau_ref_epoch = 58849
        tp = chain["$(planet_key)_tp"]
        a = chain[:b_a]
        M = chain[:M]
        period_days = @. √(a^3/M)*PlanetOrbits.kepler_year_to_julian_day_conversion_factor
        τ = @. mod((tp  - tau_ref_epoch)/period_days, 1)
        

        sma =  chain[:,"$(planet_key)_a",:]
        ecc =  chain[:,"$(planet_key)_e",:]
        inc =  chain[:,"$(planet_key)_i",:]
        aop =  chain[:,"$(planet_key)_ω",:]
        pan =  chain[:,"$(planet_key)_Ω",:]
        tau =  τ
        plx =  chain[:,"plx",:]
        mtot = chain[:,"M",:]

        dat = transpose(hcat(
            vec(collect(sma )),
            vec(collect(ecc )),
            vec(collect(inc )),
            vec(collect(aop )),
            vec(collect(pan )),
            vec(collect(tau )),
            vec(collect(plx )),
            vec(collect(mtot)),
        ))

        f["col_names"] = ["sma", "ecc", "inc", "aop", "pan", "tau", "plx", "mtot"]

        attrs(f)["tau_ref_epoch"] = tau_ref_epoch
        attrs(f)["sampler_name"] = "Octofitter"
        
        dset = create_dataset(f, "post", Float32, size(dat))
        write(dset, convert(Matrix{Float32}, dat))
        return
    end
end
