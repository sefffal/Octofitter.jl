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
                    distribution isa Distributions.UnivariateDistribution ||
                    distribution isa Octofitter.Parameterization
                )
                    error("prior on variable $($(Meta.quot(varname))) must be a UnivariateDistribution, Octofitter.Parameterization")
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
                    distribution isa Distributions.UnivariateDistribution ||
                    distribution isa Octofitter.Parameterization
                )
                    error("prior on variable $($(Meta.quot(varname))) must be a UnivariateDistribution, Octofitter.Parameterization")
                end;
                distribution
            ))
        elseif statement.head == :(=)
            varname = statement.args[1]
            expression = statement.args[2]
            esc(:(
                $varname = system -> $expression
            ))
            # TODO: can we constant-ify the arguments they passed in to avoid performance issues? Or wrap in a let block?
            # We would need to recurse through their expression and find all variabels, then put them in a let.
        else
            error("invalid statement encoutered $(statement.head)")
        end
    end
    return quote 
        $(esc(name)) = $System(
            Variables(;$(variables...))...,
            $((esc(o) for o in likelihoods)...);
            name=$(Meta.quot(name))
        )
    end
end
export @system



# Copied from ModellingToolkit.

"""
    @named x = f(...)

For variable assignments like the above, this is equivalent to:
`x = f(..., name=:x)`.

This shorthand copied from ModellingToolkit is a handy way to
cut down on typing and make sure that names of objects are kept
in sync with variable names.
"""
macro named(expr)
    name, call = split_assign(expr)
    if Meta.isexpr(name, :ref)
        name, idxs = name.args
        check_name(name)
        esc(_named_idxs(name, idxs, :($(gensym()) -> $call)))
    else
        check_name(name)
        esc(:($name = $(_named(name, call))))
    end
end
export @named 

macro named(name::Symbol, idxs, call)
    esc(_named_idxs(name, idxs, call))
end

function _named(name, call, runtime=false)
    has_kw = false
    call isa Expr || throw(Meta.ParseError("The rhs must be an Expr. Got $call."))
    if length(call.args) >= 2 && call.args[2] isa Expr
        # canonicalize to use `:parameters`
        if call.args[2].head === :kw
            call.args[2] = Expr(:parameters, Expr(:kw, call.args[2].args...))
            has_kw = true
        elseif call.args[2].head === :parameters
            has_kw = true
        end
    end

    if !has_kw
        param = Expr(:parameters)
        if length(call.args) == 1
            push!(call.args, param)
        else
            insert!(call.args, 2, param)
        end
    end

    kws = call.args[2].args

    if !any(kw->(kw isa Symbol ? kw : kw.args[1]) == :name, kws) # don't overwrite `name` kwarg
        pushfirst!(kws, Expr(:kw, :name, runtime ? name : Meta.quot(name)))
    end
    call
end

function _named_idxs(name::Symbol, idxs, call)
    if call.head !== :->
        throw(ArgumentError("Not an anonymous function"))
    end
    if !isa(call.args[1], Symbol)
        throw(ArgumentError("not a single-argument anonymous function"))
    end
    sym, ex = call.args
    ex = Base.Cartesian.poplinenum(ex)
    ex = _named(:(Symbol($(Meta.quot(name)), :_, $sym)), ex, true)
    ex = Base.Cartesian.poplinenum(ex)
    :($name = $map($sym->$ex, $idxs))
end

check_name(name) = name isa Symbol || throw(Meta.ParseError("The lhs must be a symbol (a) or a ref (a[1:10]). Got $name."))


function split_assign(expr)
    if !(expr isa Expr && expr.head === :(=) && expr.args[2].head === :call)
        throw(ArgumentError("expression should be of the form `sys = foo(a, b)`"))
    end
    name, call = expr.args
end
