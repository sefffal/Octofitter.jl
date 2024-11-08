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
function loadhdf5(fname_or_targetname, numchains=1)
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
        # TODO: should read the `colnames` property instead of assuming:
        chn =  MCMCChains.Chains(
            arr,
            [
                :b_a,
                :b_e,
                :b_i,
                :b_ω,
                :b_Ω,
                :b_τ,
                :plx,
                :M
            ],
        )

        # Calculate epoch of periastron passage from Orbitize tau variable:
        tau_ref_epoch = 58849
        if haskey(attrs(f),"tau_ref_epoch")
            tau_ref_epoch = attrs(f)["tau_ref_epoch"]
        end

        a = chn[:b_a]
        τ = chn[:b_τ]
        M = chn[:M]
        period_days = @. √(a^3/M)*PlanetOrbits.kepler_year_to_julian_day_conversion_factor
            
        tp = @. τ * period_days + tau_ref_epoch

        chn =  MCMCChains.Chains(
            hcat(arr, tp),
            [
                :b_a,
                :b_e,
                :b_i,
                :b_ω,
                :b_Ω,
                :b_τ,
                :plx,
                :M,
                :b_tp
            ],
        )

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
        
        dset = create_dataset(f, "post", Float32, size(dat))
        write(dset, convert(Matrix{Float32}, dat))
        return
    end
end
