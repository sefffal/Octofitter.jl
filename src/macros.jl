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
    priors_varnames = []
    derived_vars = OrderedDict{Symbol,Any}()
    user_likelihoods = []
    
    # Track variable names for duplicate detection
    seen_prior_vars = Set{Symbol}()
    seen_derived_vars = Set{Symbol}()
    
    for statement in variables_block
        if statement.head == :call && statement.args[1] == :~
            # Check if LHS is a distribution (Distribution(...) ~ expression)
            if statement.args[2] isa Expr && statement.args[2].head == :call
                # This is a user likelihood: Distribution(...) ~ expression
                dist_expr = statement.args[2]
                rhs_expr = statement.args[3]
                
                # Generate unique symbol for the derived variable
                derived_sym = Symbol("rhs_",generate_userlike_name(rhs_expr))
                
                # Check for duplicate derived variable names (including generated ones)
                if derived_sym in seen_derived_vars
                    error("Generated derived variable name '$derived_sym' conflicts with existing variable. Please use a different expression or variable name.")
                end
                push!(seen_derived_vars, derived_sym)
                
                # Process the RHS expression for variable capture
                local_quote_vars = Symbol[]
                local_quote_vals = Any[]
                processed_expr = quasiquote!(deepcopy(rhs_expr), local_quote_vars, local_quote_vals)
                
                # Add to global captured variables
                for (var, val) in zip(local_quote_vars, local_quote_vals)
                    if !(var in all_quote_vars)
                        push!(all_quote_vars, var)
                        push!(all_quote_vals, val)
                    end
                end
                
                # Add to derived variables
                derived_vars[derived_sym] = processed_expr
                
                # Create UserLikelihood
                # Generate name from distribution type and expression
                like_name = generate_userlike_name(rhs_expr)
                
                push!(user_likelihoods, quote
                    distribution = try
                        $(esc(dist_expr))
                    catch err
                        @error "Error creating distribution for user likelihood" expression=$(string(dist_expr))
                        rethrow(err)
                    end
                    if !(distribution isa Distributions.Distribution)
                        error("Left-hand side of ~ must be a Distribution when used for user likelihood")
                    end
                    UserLikelihood(distribution, $(Meta.quot(derived_sym)), $(string(like_name)))
                end)
            else
                # Regular prior: varname ~ Distribution
                varname = statement.args[2]
                expression = statement.args[3]
                
                # Check for duplicate prior variable names
                if varname in seen_prior_vars
                    error("Duplicate prior variable '$varname'. Each variable can only be defined once with ~.")
                end
                if varname in seen_derived_vars
                    error("Variable '$varname' is already defined as a derived variable (=). Each variable can only be defined once.")
                end
                push!(seen_prior_vars, varname)
                
                push!(priors, :( 
                    distribution = try
                        $(esc(expression))
                    catch err
                        @error "There is an error in your prior specification for $($(Meta.quot(varname))). The right-hand side of the ~ must be a Distributions.jl distribution. Did you run `using Distributions`? This is what we got:" expression=$(string(expression))
                        rethrow(err)
                    end;
                    if distribution isa Tuple || distribution isa AbstractArray
                        distribution = Product([distribution...])
                    end;
                    if !(
                        distribution isa Distributions.Distribution ||
                        distribution isa Octofitter.Parameterization
                    )
                        error("prior on variable $($(Meta.quot(varname))) must be a UnivariateDistribution, Octofitter.Parameterization")
                    end;
                    distribution
                ))
                push!(priors_varnames, varname)
            end
        elseif statement.head == :(=)
            varname = statement.args[1]
            expression = statement.args[2]
            
            # Check for duplicate derived variable names
            if varname in seen_derived_vars
                error("Duplicate derived variable '$varname'. Each variable can only be defined once with =.")
            end
            if varname in seen_prior_vars
                error("Variable '$varname' is already defined as a prior variable (~). Each variable can only be defined once.")
            end
            push!(seen_derived_vars, varname)
            
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
        priors_out = []
        derived_out = []
        likelihoods_out = AbstractLikelihood[]
        priors_evaled = [$(priors...)]
        for (varname, prior) in zip($priors_varnames, priors_evaled)
            out = expandparam(varname, prior)
            append!(priors_out, out.priors)
            append!(derived_out, out.derived)
            append!(likelihoods_out, out.likelihoods)
        end
        
        # Add derived_vars to derived_out - THIS IS THE FIX
        for (k, v) in $derived_dict_expr
            push!(derived_out, k => v)
        end
        
        # Add user likelihoods
        user_likes = [$(user_likelihoods...)]
        append!(likelihoods_out, user_likes)
        
        if isempty(likelihoods_out)
            (
                Priors(;[l=>r for (l,r) in priors_out]...), 
                Derived(
                    OrderedDict{Symbol,Any}([l=>r for (l,r) in derived_out]),
                    captured_names,
                    captured_vals
                ),
            )
        else
            (
                Priors(;[l=>r for (l,r) in priors_out]...), 
                Derived(
                    OrderedDict{Symbol,Any}([l=>r for (l,r) in derived_out]),
                    captured_names,
                    captured_vals
                ),
                likelihoods_out...
            )
        end
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

# Helper function to create human-readable names for user likelihoods
function generate_userlike_name(rhs_expr)
   
    # Simplify expression string
    expr_str = string(rhs_expr)
    
    # Replace common mathematical operators with words
    expr_str = replace(expr_str, "^" => "_pow_")
    expr_str = replace(expr_str, "*" => "_times_")
    expr_str = replace(expr_str, "/" => "_over_")
    expr_str = replace(expr_str, "+" => "_plus_")
    expr_str = replace(expr_str, "-" => "_minus_")
    expr_str = replace(expr_str, "(" => "")
    expr_str = replace(expr_str, ")" => "")
    expr_str = replace(expr_str, " " => "")
    
    # Truncate if too long
    if length(expr_str) > 20
        expr_str = expr_str[1:20] * "_etc"
    end
    
    # Create descriptive name
    
    return normalizename(expr_str)
end