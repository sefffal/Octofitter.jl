"""
    @variables begin
        [prior_1] ~ [UnivariateDistribution]
        [prior_2] ~ [UnivariateDistribution]
        calculation_3 = obs.[prior_1] + obs.[prior_2]
    end
"""
macro variables(variables_block_input)
    variables_block = filter(variables_block_input.args) do expr
        !(expr isa LineNumberNode)
    end
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
                    @error "There is an error in your prior specification for $($(Meta.quot(varname))). The right-hand side of the ~ must be a Distributions.jl distribution. Did you run `using Distributions`? This is what we got:" expression=$(string(expression))
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
                $varname = ((super, this) -> $expression
            )))
            # TODO: can we constant-ify the arguments they passed in to avoid performance issues? Or wrap in a let block?
            # We would need to recurse through their expression and find all variabels, then put them in a let.
        else
            error("invalid statement encoutered $(statement.head)")
        end
    end

    # We allow users to interpolate local variables into the model definition 
    # with `$`.
    # We wrap the model definition in a local anonymous function and use 
    # the arguments to pass in these local variables.
    quote_vars = map(quote_vals) do val
        return esc(val.args[1])
    end
    quoted_vals_escaped = map(quote_vals) do val
        return quote
            $(esc(val.args[1]))
        end
    end
    # Drop any duplicates (if the same variable is interpolated in twice)
    ii = []
    for j in eachindex(quote_vars)
        found = false
        for i in ii
            if quote_vars[i] == quote_vars[j]
                found = true
                break
            end
        end
        if !found
            push!(ii, j)
        end
    end
    quote_vars = quote_vars[ii]
    quoted_vals_escaped = quoted_vals_escaped[ii]
    return quote
        (function($(quote_vars...),)
            return $Variables(;$(variables...))
        end)($(quoted_vals_escaped...))
    end
end
export @variables

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
        var = ex.args[1]
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

