module OctofitterPairPlotsExt
using Octofitter
using PairPlots
using Tables
using NamedTupleTools
using MCMCChains

unitlengthprior_vars(::Octofitter.UnitLengthPrior{VX,VY}) where {VX, VY} = (VX, VY)

function Octofitter.octocorner(
    system::System,
    chains::Octofitter.MCMCChains.Chains...; 
    small=false
)
    labels = Dict{Symbol,Any}(
        :M => "total mass [M⊙]",
        :plx => "parallax [mas]",
    )
    prepared = map(chains) do chain
        chain_notinternal = MCMCChains.get_sections(chain, :parameters)
        chain_keys = string.(keys(chain_notinternal))
        table_cols = Pair{Symbol}[]
        if small
            push!(table_cols, :M => vec(chain_notinternal["M"]))
        else
            planetkeys = string.(keys(system.planets))
            ii = map(chain_keys) do k
                for pk in planetkeys
                    if startswith(k, pk*"_")
                        return false
                    end
                end
                return true
            end
            for chain_key in chain_keys[ii]
                push!(table_cols, Symbol(chain_key)=>vec(chain_notinternal[chain_key]))
            end
        end

        for (planetkey, planet) in pairs(system.planets)
            OrbitType = Octofitter.orbittype(planet)
            # els = Octofitter.construct_elements(chain,:b,:)
            
            pk_ = "$(planetkey)_"
            planet_var_keys = filter(startswith(pk_), chain_keys)
            planet_var_keys_chopped = chopprefix.(planet_var_keys, pk_)

            if small
                ii_splice = findall(map(planet_var_keys_chopped) do k
                    k ∉ ["a", "e", "i", "mass"]
                end)
            else
                # Remove x and y parameters used by UniformCircular
                # Also remove tp if tau is used
                ii_splice = findall(map(planet_var_keys_chopped) do k
                    k = Symbol(k)
                    if k == :tp && "τ" ∈ planet_var_keys_chopped
                        return true
                    end
                    for obs in planet.observations
                        if obs isa Octofitter.UnitLengthPrior
                            varx, vary = unitlengthprior_vars(obs)
                            if k == varx || k == vary
                                return true
                            end
                        end
                    end
                    return false
                end)
            end

            # remove variables that are held fixed.
            ii_splice_2 = findall(map(planet_var_keys) do k
                vals = vec(chain[k])
                all(==(first(vals)), vals)
            end)
            ii_splice = sort(union(ii_splice, ii_splice_2))

            splice!(planet_var_keys, ii_splice)
            splice!(planet_var_keys_chopped, ii_splice)

            labels[Symbol(pk_*"a")] = "$planetkey:\nsemi-major axis [au]"
            labels[Symbol(pk_*"i")] = "$planetkey:\ninclination [°]"
            labels[Symbol(pk_*"Ω")] = "$planetkey:\nlogintude of\nascending node [°]"
            labels[Symbol(pk_*"ω")] = "$planetkey:\nargument of\nperiapsis [°]"
            labels[Symbol(pk_*"e")] = "$planetkey:\neccentricity"
            labels[Symbol(pk_*"mass")] = "$planetkey:\nmass [Mⱼᵤₚ]"
            labels[Symbol(pk_*"tp")] = "$planetkey:\nperiastron [mjd]"
            labels[Symbol(pk_*"P")] = "$planetkey:\nperiod [yrs]"
            labels[Symbol(pk_*"τ")] = "$planetkey:\norbit fraction τ"

            for (k, pk) in zip(planet_var_keys_chopped, planet_var_keys)
                pks = Symbol(pk)

                # Add generic label if missing
                if pks ∉ keys(labels)
                    labels[pks] = "$planetkey:\n$k"
                end

                # Do any data conversions
                dat = vec(chain_notinternal[pk])
                if k == "i" || k == "Ω" || k == "ω"
                    dat = rad2deg.(rem2pi.(dat, RoundDown))
                end
                push!(table_cols, pks => dat)
            end

        end
        tbl = FlexTable(namedtuple(table_cols))
    end

    pairplot(prepared...;labels)
end


end