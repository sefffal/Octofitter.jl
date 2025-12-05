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
            varname = statement.args[2]
            expression = statement.args[3]

            # Check if LHS is a distribution (Distribution(...) ~ expression)
            u = union(seen_prior_vars, seen_derived_vars)
            if varname in u || expression in u || (statement.args[2] isa Expr && statement.args[2].head == :call)
                # This is a user likelihood: Distribution(...) ~ expression
                lhs_expr = varname
                rhs_expr = expression
                
                # Generate unique symbol for the derived variable
                derived_sym_lhs = Symbol("lhs_",generate_userlike_name(rhs_expr))
                derived_sym_rhs = Symbol("rhs_",generate_userlike_name(rhs_expr))
                
                # Check for duplicate derived variable names (including generated ones)
                if derived_sym_rhs in seen_derived_vars
                    error("Generated derived variable name '$derived_sym_rhs' conflicts with existing variable. Please use a different expression or variable name.")
                end
                if derived_sym_lhs in seen_derived_vars
                    error("Generated derived variable name '$derived_sym_lhs' conflicts with existing variable. Please use a different expression or variable name.")
                end
                push!(seen_derived_vars, derived_sym_lhs)
                push!(seen_derived_vars, derived_sym_rhs)
                
                # Process the RHS expression for variable capture
                local_quote_vars = Symbol[]
                local_quote_vals = Any[]
                processed_expr_lhs = quasiquote!(deepcopy(lhs_expr), local_quote_vars, local_quote_vals)
                processed_expr_rhs = quasiquote!(deepcopy(rhs_expr), local_quote_vars, local_quote_vals)
                
                # Add to global captured variables
                for (var, val) in zip(local_quote_vars, local_quote_vals)
                    if !(var in all_quote_vars)
                        push!(all_quote_vars, var)
                        push!(all_quote_vals, val)
                    end
                end
                
                # Add to derived variables
                derived_vars[derived_sym_rhs] = processed_expr_rhs
                derived_vars[derived_sym_lhs] = processed_expr_lhs
                
                # Create UserLikelihood
                # Generate name from distribution type and expression
                like_name = generate_userlike_name(rhs_expr)
                
                push!(user_likelihoods, UserLikelihood(derived_sym_lhs, derived_sym_rhs, like_name))
                # push!(user_likelihoods, quote
                #     distribution = try
                #         $(esc(dist_expr))
                #     catch err
                #         @error "Error creating distribution for user likelihood" expression=$(string(lhs_expr))
                #         rethrow(err)
                #     end
                #     if !(distribution isa Distributions.Distribution)
                #         error("Left-hand side of ~ must be a Distribution when used for user likelihood")
                #     end
                #     UserLikelihood(distribution, $(Meta.quot(derived_sym)), $(string(like_name)))
                # end)
            else
                # Regular prior: varname ~ Distribution
                
                # Check for duplicate prior variable names
                if varname in seen_prior_vars
                    error("Duplicate prior variable '$varname'. Each variable can only be defined once with ~.")
                end
                # if varname in seen_derived_vars
                #     error("Variable '$varname' is already defined as a derived variable (=). Each variable can only be defined once.")
                # end
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


"""
    Base.vcat for Variables blocks
    
Concatenates two or more @variables blocks together, combining their priors,
derived variables, and likelihoods. Throws an error if any variable names
are duplicated across the blocks.

Example:
```julia
vars1 = @variables begin
    a ~ Normal()
    b = 2a
end

vars2 = @variables begin
    d ~ Normal()
    e = 2d^2
end

vars = vcat(vars1, vars2)
```
"""
function Base.vcat(vars1::Tuple, vars2::Tuple, vars_rest::Tuple...)
    # Check that these are actually @variables outputs
    # They should have at least Priors and Derived as first two elements
    if length(vars1) < 2 || !isa(vars1[1], Priors) || !isa(vars1[2], Derived)
        throw(ArgumentError("First argument does not appear to be from @variables macro"))
    end
    if length(vars2) < 2 || !isa(vars2[1], Priors) || !isa(vars2[2], Derived)
        throw(ArgumentError("Second argument does not appear to be from @variables macro"))
    end
    
    # Start with vars1 and vars2
    result = _vcat_two_variables(vars1, vars2)
    
    # If there are more variable blocks, concatenate them iteratively
    for vars in vars_rest
        if length(vars) < 2 || !isa(vars[1], Priors) || !isa(vars[2], Derived)
            throw(ArgumentError("Additional argument does not appear to be from @variables macro"))
        end
        result = _vcat_two_variables(result, vars)
    end
    
    return result
end

# Helper function to concatenate exactly two variable blocks
function _vcat_two_variables(vars1::Tuple, vars2::Tuple)
    priors1 = vars1[1]
    derived1 = vars1[2]
    likelihoods1 = length(vars1) > 2 ? vars1[3:end] : ()
    
    priors2 = vars2[1]
    derived2 = vars2[2]
    likelihoods2 = length(vars2) > 2 ? vars2[3:end] : ()
    
    # Check for duplicate variable names
    prior_names1 = Set(keys(priors1.priors))
    prior_names2 = Set(keys(priors2.priors))
    derived_names1 = Set(keys(derived1.variables))
    derived_names2 = Set(keys(derived2.variables))
    
    # Check for duplicates between priors
    duplicate_priors = intersect(prior_names1, prior_names2)
    if !isempty(duplicate_priors)
        error("Duplicate prior variable(s) found when concatenating @variables blocks: $(join(duplicate_priors, ", "))")
    end
    
    # Check for duplicates between derived variables
    duplicate_derived = intersect(derived_names1, derived_names2)
    if !isempty(duplicate_derived)
        error("Duplicate derived variable(s) found when concatenating @variables blocks: $(join(duplicate_derived, ", "))")
    end
    
    # Check for cross-duplicates (prior in one, derived in another)
    cross_duplicates1 = intersect(prior_names1, derived_names2)
    if !isempty(cross_duplicates1)
        error("Variable(s) defined as prior in first block but derived in second block: $(join(cross_duplicates1, ", "))")
    end
    
    cross_duplicates2 = intersect(prior_names2, derived_names1)
    if !isempty(cross_duplicates2)
        error("Variable(s) defined as prior in second block but derived in first block: $(join(cross_duplicates2, ", "))")
    end
    
    # Merge priors
    merged_priors_dict = OrderedDict{Symbol,Distribution}()
    for (k, v) in priors1.priors
        merged_priors_dict[k] = v
    end
    for (k, v) in priors2.priors
        merged_priors_dict[k] = v
    end
    merged_priors = Priors(merged_priors_dict)
    
    # Merge derived variables
    merged_derived_dict = OrderedDict{Symbol,Any}()
    for (k, v) in derived1.variables
        merged_derived_dict[k] = v
    end
    for (k, v) in derived2.variables
        merged_derived_dict[k] = v
    end
    
    # Merge captured variables
    # We need to combine the captured names and values from both blocks
    captured_names = (derived1.captured_names..., derived2.captured_names...)
    captured_vals = (derived1.captured_vals..., derived2.captured_vals...)
    
    # Remove duplicates in captured variables (keeping first occurrence)
    unique_captures = OrderedDict{Symbol,Any}()
    for (name, val) in zip(captured_names, captured_vals)
        if !haskey(unique_captures, name)
            unique_captures[name] = val
        end
    end
    
    merged_captured_names = tuple(keys(unique_captures)...)
    merged_captured_vals = tuple(values(unique_captures)...)
    
    merged_derived = Derived(merged_derived_dict, merged_captured_names, merged_captured_vals)
    
    # Merge likelihoods
    merged_likelihoods = (likelihoods1..., likelihoods2...)
    
    # Return in the same format as @variables macro
    if isempty(merged_likelihoods)
        return (merged_priors, merged_derived)
    else
        return (merged_priors, merged_derived, merged_likelihoods...)
    end
end

# Also support concatenating more than 2 blocks at once using varargs
Base.vcat(vars::Tuple...) = Base.vcat(vars[1], vars[2], vars[3:end]...)