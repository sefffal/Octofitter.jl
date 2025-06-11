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
    
    # Global captured variables for all derived expressions
    all_quote_vars = Symbol[]
    all_quote_vals = Any[]
    
    priors = []
    derived_vars = OrderedDict{Symbol,Any}()
    
    for statement in variables_block
        if statement.head == :call && statement.args[1] == :~
            varname = statement.args[2]
            expression = statement.args[3]
            push!(priors, :($(esc(varname)) = ( 
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
            )))
        elseif statement.head == :(=)
            varname = statement.args[1]
            expression = statement.args[2]
            
            # Process the expression to extract interpolated variables
            local_quote_vars = Symbol[]
            local_quote_vals = Any[]
            processed_expr = quasiquote!(deepcopy(expression), local_quote_vars, local_quote_vals)
            
            # Add to global captured variables (avoiding duplicates)
            for (var, val) in zip(local_quote_vars, local_quote_vals)
                if !(var in all_quote_vars)
                    push!(all_quote_vars, var)
                    push!(all_quote_vals, val)
                end
            end
            
            # Store the processed expression directly (not wrapped in a function)
            derived_vars[varname] = processed_expr
        else
            error("invalid statement encountered $(statement.head)")
        end
    end

    # Create escaped versions of captured variable names and values
    escaped_vars = [esc(var) for var in all_quote_vars]
    # Properly escape the full interpolation expression to evaluate in caller's scope
    escaped_vals = [esc(val.args[1]) for val in all_quote_vals]
    
    # Build the Derived constructor call with OrderedDict
    derived_dict_expr = if isempty(derived_vars)
        :(OrderedDict{Symbol,Any}())
    else
        :(OrderedDict{Symbol,Any}($([:($(QuoteNode(k)) => $(Meta.quot(v))) for (k,v) in derived_vars]...)))
    end
    
    # Don't escape the variable names in the tuple - they should be symbols
    captured_names_tuple = Expr(:tuple, [QuoteNode(var) for var in all_quote_vars]...)
    
    return quote
        # Create the function with unescaped parameter names
        local captured_vals = tuple($([esc(val.args[1]) for val in all_quote_vals]...))
        local captured_names = $captured_names_tuple
        (
            Priors(;$(priors...)), 
            Derived(
                $derived_dict_expr,
                captured_names,
                captured_vals
            )
        )
    end
end

export @variables

# Adapted from BenchmarkTools
# We use this for $variable interpolation into models.
# Users can do this to avoid type-instabilities / global variable 
# access when using constants etc. in their models.
# This only applies to deterministic variables, sampled variables already
# resolve any variables at model-creation time.
function quasiquote!(ex, vars::Vector{Symbol}, vals::Vector)
    ex
end
function quasiquote!(ex::Symbol, vars::Vector{Symbol}, vals::Vector)
    ex
end
function quasiquote!(ex::Expr, vars::Vector{Symbol}, vals::Vector)
    if ex.head === :($)
        # Don't use gensym - keep the original variable name
        var = isa(ex.args[1], Symbol) ? ex.args[1] : gensym()
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