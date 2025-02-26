module OctofitterPairPlotsExt
using Octofitter
using PairPlots
using Tables
using NamedTupleTools
using MCMCChains

unitlengthprior_vars(::Octofitter.UnitLengthPrior{VX,VY}) where {VX, VY} = (VX, VY)

Octofitter.octocorner(
    model::Octofitter.LogDensityModel,
    chains::MCMCChains.Chains...; 
    kwargs...
) = octocorner(model.system, chains...; kwargs...)
function Octofitter.octocorner(
    system::System,
    chains::MCMCChains.Chains...; 
    small=false,
    labels=Dict{Symbol,Any}(),
    fname=small ? "$(system.name)-pairplot-small.png" : "$(system.name)-pairplot.png" ,
    viz=nothing,
    # Optionally include additional columns
    includecols=String[],
    excludecols=String[],
    bottomleft=true,
    topright=false,
    kwargs...
)
    labels_gen = Dict{Symbol,Any}(
        :M => "M [M⊙]\ntotal mass ",
        :plx => "plx [mas]\nparallax",
    )
    if length(chains) > 1
        if isnothing(viz)
            viz = PairPlots.multi_series_default_viz
        end
        colors = PairPlots.Makie.wong_colors()
    else
        if isnothing(viz)
            viz = PairPlots.single_series_default_viz
        end
        colors = [PairPlots.single_series_color]
    end
    colori = 1
    function preparechain((chain,viz)::Pair)
        prepped = _preparechain(chain)
        return prepped => viz
    end
    function preparechain(chain)
        prepped = _preparechain(chain)
        name = hasproperty(chain.info, :model_name) ? string(chain.info.model_name) : string(system.name)
        pair = PairPlots.Series(prepped,label=name,color=colors[mod1(colori,end)]; topright, bottomleft) => viz
        colori += 1
        return pair
    end
    function _preparechain(chain::Octofitter.MCMCChains.Chains)
        chain_notinternal = MCMCChains.get_sections(chain, :parameters)
        chain_keys = string.(keys(chain_notinternal))
        table_cols = Pair{Symbol}[]

        for colname in includecols
            if colname==:iter
                push!(table_cols, Symbol(:iter)=>repeat(1:size(chain,1),outer=size(chain,3)))
            else
                push!(table_cols, Symbol(colname)=>vec(chain[colname]))
            end
        end

        if small
            if haskey(chain_notinternal, :M)
                push!(table_cols, :M => vec(chain_notinternal["M"]))
            end
            planetkeys = string.(keys(system.planets))
            for obs in system.observations
                if hasproperty(obs,:table) && hasproperty(obs.table, :band)
                    bands = unique(obs.table.band)
                    for pk in planetkeys, band in bands
                        k = Symbol("$(pk)_$band")
                        if haskey(chain_notinternal, k)
                            push!(table_cols, band => vec(chain_notinternal[k]))
                        end
                    end
                end
            end
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
                keep_params_if_small = ["a", "e", "i", "mass", "A", "B", "F", "G"]
                for obs in planet.observations
                    if hasproperty(obs,:table) && hasproperty(obs.table, :band)
                        append!(keep_params_if_small, string.(unique(obs.table.band)))
                    end
                end
                ii_splice = findall(map(planet_var_keys_chopped) do k
                    k ∉ keep_params_if_small
                end)
            else
                # Remove x and y parameters used by UniformCircular
                # Also remove tp if tau is used
                ii_splice = findall(map(planet_var_keys_chopped) do k
                    k = Symbol(k)
                    if k == :tp && "τ" ∈ planet_var_keys_chopped
                        return true
                    end
                    if k == :tp && "θ" ∈ planet_var_keys_chopped
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


            for (k, pk) in zip(planet_var_keys_chopped, planet_var_keys)
                pks = Symbol(pk)

                if k=="a"
                    labels_gen[Symbol(pk_*"a")] = "$(planetkey)_a [au]\nsemi-major axis"
                end
                if k=="i"
                    labels_gen[Symbol(pk_*"i")] = "$(planetkey)_i [°]\ninclination"
                end
                if k=="Ω"
                    labels_gen[Symbol(pk_*"Ω")] = "$(planetkey)_Ω [°]\nlongitude of\nascending node"
                end
                if k=="ω"
                    labels_gen[Symbol(pk_*"ω")] = "$(planetkey)_ω [°]\nargument of\nperiapsis"
                end
                if k=="e"
                    labels_gen[Symbol(pk_*"e")] = "$(planetkey)_e\neccentricity"
                end
                if k=="mass"
                    labels_gen[Symbol(pk_*"mass")] = "$(planetkey)_mass [Mⱼᵤₚ]\nmass"
                end
                if k=="tp"
                    labels_gen[Symbol(pk_*"tp")] = "$(planetkey)_tp [mjd]\nepoch of periastron\npassage"
                end
                if k=="P"
                    labels_gen[Symbol(pk_*"P")] = "$(planetkey)_P [yrs]\nperiod"
                end
                if k=="τ"
                    labels_gen[Symbol(pk_*"τ")] = "$(planetkey)_τ\norbit fraction τ"
                end
                if k=="θ"
                    labels_gen[Symbol(pk_*"θ")] = "$(planetkey)_θ [°]\nposition angle\nat ref. epoch"
                end

                # Add generic label if missing
                if pks ∉ keys(labels_gen)
                    labels_gen[pks] = "$planetkey:\n$k"
                end

                # Do any data conversions
                dat = vec(chain_notinternal[pk])
                if k == "i" || k == "Ω" || k == "ω" || k == "θ" || k == "pa"
                    dat = rad2deg.(rem2pi.(dat, RoundDown))
                end
                push!(table_cols, pks => dat)
            end
        end
        # Remove duplicates (eg. if a column is added by the user to includecols
        # and would have been included anyways).
        # Also remove series where the std is 0.
        table_cols_unique = Pair{Symbol}[]
        for (k,v) in table_cols
            if k in excludecols
                continue
            end
            if length(unique(v)) == 1
                continue
            end
            found = false
            for (k2,v2) in table_cols_unique
                if k==k2
                    found = true
                    break
                end
            end
            if !found
                push!(table_cols_unique,(k=>v))
            end
        end
        tbl = FlexTable(namedtuple(table_cols_unique))
    end

    prepared = map(preparechain, chains)

    # Merge our generated labels with user labels
    labels = merge(labels_gen, labels)
    fig = pairplot(prepared...;labels, kwargs...)

    PairPlots.Makie.save(fname, fig, px_per_unit=3)

    return fig
end


end