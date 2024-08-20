"""
    @planet [planet_name] [Orbit_Type] begin
        [prior_1] ~ [UnivariateDistribution]
        [prior_2] ~ [UnivariateDistribution]
        calculation_3 = [planet_name].[prior_1] + [planet_name].[prior_2] + [system.variable_x]
    end [likelihood_objects...]

Generate a Planet model named `planet_name`. A variable will be created with the name `[planet_name]`
in the current scope.
`Orbit_Type` specifies the orbit parameterization from PlanetOrbits.jl. You must provide all
input variables needed for the selected orbit type (see PlanetOrbits.jl documentation).
Following that is a block of variable assignments. Variables with a `~` will be free variables
with a prior distribution given by the right-hand-side (a UnivariateDistribution from Distributions.jl
or a `KDEDist`). 
Calculated quantities are also allowed. These may reference other variables using the planet name followed
by a dot and the name of the variable. Variables from other planets in a single system are not accessible.
You can access other variables in the current local scope, but these bindings are only guaranteed to be
resolved a single time. Note that using non-constant global variables in calculated expressions can lead 
to poor performance.
Finally, the planet model can be conditioned on data by supplying zero or more likelihood objects.
"""
macro planet(args...)
    name = args[1]
    orbit_type = args[2]
    variables_block_input = args[3]
    variables_block = filter(variables_block_input.args) do expr
        !(expr isa LineNumberNode)
    end
    likelihoods = args[4:end]
    variables = map(variables_block) do statement
        if statement.head == :call && statement.args[1] == :~
            varname = statement.args[2]
            expression = statement.args[3]
            return :($(esc(varname)) = (
                distribution = try
                    $(esc(expression))
                catch err
                    @error "There is an error in your prior specification for $($(Meta.quot(varname))). The right-hand side of the ~ must be a Distribution distribution. This is what we got:" expression=$(string(expression))
                    rethrow(err)
                end;
                if !(
                    distribution isa Distributions.Distribution ||
                    distribution isa Octofitter.Parameterization
                )
                    error("prior on variable $($(Meta.quot(varname))) must be a Distribution, Octofitter.Parameterization")
                end;
                distribution
            ))
        elseif statement.head == :(=)
            varname = statement.args[1]
            expression = statement.args[2]
            return :(
                $(esc(varname)) = ($(esc(:system)), $(esc(name))) -> $(esc(expression))
            )
            # TODO: can we constant-ify the arguments they passed in to avoid performance issues? Or wrap in a let block?
            # We would need to recurse through their expression and find all variabels, then put them in a let.
        else
            error("invalid statement encoutered $(statement.head)")
        end
    end
    return quote 
        $(esc(name)) = $Planet{$orbit_type}(
            Variables(;$(variables...))...,
            $((esc(l) for l in likelihoods)...);
            name=$(Meta.quot(name))
        )
    end
end
export @planet

"""
    @system [system_name] begin
        [prior_1] ~ [UnivariateDistribution]
        [prior_2] ~ [UnivariateDistribution]
        calculation_3 = system.[prior_1] + system.[prior_2]
    end [likelihood_objects...] [planet_models...]

Generate a System model named `system_name`. A variable will be created with the name `[system_name]`
in the current scope.
Following that is a block of variable assignments. Variables with a `~` will be free variables
with a prior distribution given by the right-hand-side (a UnivariateDistribution from Distributions.jl
or a `KDEDist`). 
Calculated quantities are also allowed. These may reference other variables using the planet name followed
by a dot and the name of the variable. Variables from other planets in a single system are not accessible.
You can access other variables in the current local scope, but these bindings are only guaranteed to be
resolved a single time. Note that using non-constant global variables in calculated expressions can lead 
to poor performance.
After the end of the variable block, the system model can be conditioned on data by supplying zero or more likelihood objects.
Finally, zero or more planet models can be attached to the system, potentially conditioned on likelihood objects of their own.
"""
macro system(args...)
    name = args[1]
    variables_block_input = args[2]
    variables_block = filter(variables_block_input.args) do expr
        !(expr isa LineNumberNode)
    end
    likelihoods = args[3:end]
    quote_vars = Symbol[]
    quote_vals = Expr[]
    variables = map(variables_block) do statement
        if statement.head == :call && statement.args[1] == :~
            varname = statement.args[2]
            expression = statement.args[3]
            return :($(esc(varname)) = (
                distribution = try
                    $(esc(expression))
                catch err
                    @error "There is an error in your prior specification for $($(Meta.quot(varname))). The right-hand side of the ~ must be a Distribution distribution. This is what we got:" expression=$(string(expression))
                    rethrow(err)
                end;
                if !(
                    distribution isa Distributions.Distribution ||
                    distribution isa Octofitter.Parameterization
                )
                    error("prior on variable $($(Meta.quot(varname))) must be a UnivariateDistribution, Octofitter.Parameterization")
                end;
                distribution
            ))
        elseif statement.head == :(=)
            varname = statement.args[1]
            expression = statement.args[2]
            expression = quasiquote!(expression, quote_vars, quote_vals)
            esc(:(
                $varname = system -> $expression
            ))
            # TODO: can we constant-ify the arguments they passed in to avoid performance issues? Or wrap in a let block?
            # We would need to recurse through their expression and find all variabels, then put them in a let.
        else
            error("invalid statement encoutered $(statement.head)")
        end
    end

    # Create constant bindings for all `$` interpolated variables
    # to ensure performance.
    # We just interpolate them into the global scope as constants,
    # evaluated immediately in the calling scope.
    quoted_var_pairs = map(quote_vars, quote_vals) do var, val
        return quote
            const $(esc(var)) = $(esc(val.args[1]))
        end
    end
    return quote
        $(quoted_var_pairs...)
        $(esc(name)) = $System(
            $Variables(;$(variables...))...,
            $((esc(o) for o in likelihoods)...);
            name=$(Meta.quot(name))
        )
    end
end
export @system

# Copied from BenchmarkTools
# We use this for $variable interpolation into models.
# Users can do this to avoid type-instabilities / global variable 
# access when using constants etc. in their models.
# This only applies to deterministic variables, sampled variables already
# resolve any variables at model-creation time.
"""
    quasiquote!(expr::Expr, vars::Vector{Symbol}, vals::Vector{Expr})

Replace every interpolated value in `expr` with a placeholder variable and
store the resulting variable / value pairings in `vars` and `vals`.
"""
quasiquote!(ex, _...) = ex
function quasiquote!(ex::Expr, vars::Vector{Symbol}, vals::Vector{Expr})
    if ex.head === :($)
        var = isa(ex.args[1], Symbol) ? gensym(ex.args[1]) : gensym()
        push!(vars, var)
        push!(vals, ex)
        return var
    elseif ex.head !== :quote
        for i in 1:length(ex.args)
            ex.args[i] = quasiquote!(ex.args[i], vars, vals)
        end
    end
    return ex
end

