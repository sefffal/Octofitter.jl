"""

    astrom = CSV.read("octofitter/input.csv", AstrometryLikelihood)

    vars = Octofitter.@planet b VisualOrbit astrom begin
        a  ~ LogUniform(2.5, 25)
        a2 = b.a^2
        i ~ Sine()
        e ~ truncated(
            Normal(0.2, 0.2),
            lower=0
        )
        Ω ~ UniformCircular()
        ω ~ UniformCircular()
        τ ~ UniformCircular(1.0)
        mass  ~ Uniform(0, 50)
        H_band = cooling_track_flux(system.age_Myr, b.mass)
    end

Generate a Planet model 
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
            distribution = eval(expression)
            if !(distribution isa Distributions.UnivariateDistribution || distribution isa Octofitter.Parameterization)
                error("prior on variable $varname must be a UnivariateDistribution")
            end
            return :( $(esc(varname)) = $distribution)
        elseif statement.head == :(=)
            varname = statement.args[1]
            expression = statement.args[2]
            return esc(:(
                $varname = (system, $name) -> $expression
            ))
            # TODO: can we constant-ify the arguments they passed in to avoid performance issues? Or wrap in a let block?
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
            distribution = statement.args[3]
            # distribution = eval(expression)
            # if !(distribution isa Distributions.UnivariateDistribution || distribution isa Octofitter.Parameterization)
            #     error("prior on variable $varname must be a UnivariateDistribution")
            # end
            return :( $(esc(varname)) = $(esc(distribution)))
        elseif statement.head == :(=)
            varname = statement.args[1]
            expression = statement.args[2]
            return :(
                $(esc(varname)) = (system) -> $(esc(expression))
            )
        else
            error("invalid statement encoutered $(statement.head)")
        end
    end
    return quote 
        $(esc(name)) = $System(
            Variables(;$(variables...)),
            $((esc(o) for o in likelihoods)...);
            name=$(Meta.quot(name))
        )
    end
end
export @system



# Copied from ModellingToolkit.

"""
    @named x = f(...)

For variable assiments like the above, this is equivalent to:
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
