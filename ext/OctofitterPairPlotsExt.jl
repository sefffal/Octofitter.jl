# ---------------------------------------------------
# OctofitterPairPlotsExt — corner plots for v2.
#
# The v1 extension kept a private suffix-matching label dictionary; here
# every label, unit, and radian→degree conversion comes from PlanetOrbits'
# `paraminfo` resolver table (the same one axis labels use), keyed by the
# flat `<owner>_<var>` chain naming. Columns are classified against the
# model — body names and observation names — not by string heuristics.
# ---------------------------------------------------
module OctofitterPairPlotsExt

using Octofitter
using Octofitter: LogDensityModel, normalizename, likelihoodname, UnitLengthPrior
using PlanetOrbits: paraminfo
using PairPlots
using MCMCChains

# Split a flat chain column `<owner>_<var>` against the model's owners
# (bodies and observations); system variables have no prefix.
function _colsplit(sys, col::Symbol)
    s = string(col)
    for owner in sys.bodynames
        pre = string(owner, '_')
        startswith(s, pre) && return (owner, Symbol(chopprefix(s, pre)), :body)
    end
    for obs in sys.observations
        owner = normalizename(likelihoodname(obs))
        pre = string(owner, '_')
        startswith(s, pre) && return (owner, Symbol(chopprefix(s, pre)), :obs)
    end
    return (:system, col, :system)
end

# The x/y helper pairs behind UniformCircular parametrizations: sampling
# machinery, not parameters anyone reads.
function _helper_columns(sys)
    out = Symbol[]
    for (owner, t) in sys.priorterms
        t isa UnitLengthPrior || continue
        for v in (t.varx, t.vary)
            push!(out, owner === :system ? v : Symbol(owner, :_, v))
        end
    end
    return out
end

function _label(owner, var, kind, col)
    info = paraminfo(var)
    if info === nothing
        return kind === :system ? string(col) : "$(col)\n($owner)"
    end
    head = isempty(info.unit) ? string(col) : "$(col) [$(info.unit)]"
    return "$head\n$(info.label)"
end

Octofitter.octocorner(model::LogDensityModel, chains::MCMCChains.Chains...; kwargs...) =
    _octocorner(model, chains...; kwargs...)

function _octocorner(model, chains...;
                     small=false,
                     truth=(),
                     labels=Dict{Symbol,Any}(),
                     viz=nothing,
                     includecols=Symbol[],
                     excludecols=Symbol[],
                     fname=nothing,
                     bottomleft=true, topright=false,
                     kwargs...)
    sys = model.system
    includecols = Symbol.(collect(includecols))
    excludecols = Symbol.(collect(excludecols))
    helpers = _helper_columns(sys)
    labels_gen = Dict{Symbol,Any}()

    function preparechain(chain::MCMCChains.Chains)
        params = MCMCChains.get_sections(chain, :parameters)
        cols = Pair{Symbol,Vector{Float64}}[]
        colnames = keys(params)
        # `tp` duplicates a `θ`/`M0` phase parametrization of the same body;
        # keep the sampled one.
        phase_owners = Set(o for c in colnames
                           for (o, v, k) in (_colsplit(sys, c),)
                           if k === :body && (v === :θ || v === :M0))
        for col in colnames
            dat = vec(params[col])
            owner, var, kind = _colsplit(sys, col)
            keep = col in includecols
            if !keep
                col in excludecols && continue
                col in helpers && continue
                length(unique(dat)) == 1 && continue
                kind === :body && var === :tp && owner in phase_owners && continue
                if small
                    (kind === :body && var in (:a, :e, :i, :mass)) || continue
                end
            end
            info = paraminfo(var)
            if info !== nothing && info.angle
                dat = rad2deg.(rem2pi.(dat, RoundDown))
            end
            labels_gen[col] = _label(owner, var, kind, col)
            push!(cols, col => collect(float.(dat)))
        end
        isempty(cols) && error("no plottable columns left after filtering")
        return (; cols...)
    end

    multi = length(chains) > 1
    vizn = viz === nothing ?
           (multi ? PairPlots.multi_series_default_viz : PairPlots.single_series_default_viz) :
           viz
    colors = multi ? PairPlots.Makie.wong_colors() : [PairPlots.single_series_color]
    prepared = map(enumerate(chains)) do (i, chain)
        name = hasproperty(chain.info, :model_name) ? string(chain.info.model_name) :
               string(sys.name)
        PairPlots.Series(preparechain(chain); label=name,
            color=colors[mod1(i, length(colors))], bottomleft, topright) => vizn
    end

    fig = pairplot(prepared..., truth...; labels=merge(labels_gen, labels), kwargs...)
    fname !== nothing && PairPlots.Makie.save(fname, fig, px_per_unit=3)
    return fig
end

end # module
