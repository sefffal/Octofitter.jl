# ---------------------------------------------------
# Variable namespaces
#
# `Priors` and `Derived` are what a `@variables` block compiles to: an
# ordered set of random variables, and an ordered set of *unevaluated*
# expressions computed from them. Keeping the expressions unevaluated is what
# lets `System` decide statically which system lines must wait until after
# the body blocks, and lets the code generator inline them.
# ---------------------------------------------------

"""
    Priors(key=dist, ...)

Zero or more named prior distributions, in declaration order.
"""
struct Priors
    priors::OrderedDict{Symbol,Distribution}
end
Priors(; priors...) = Priors(OrderedDict(priors))

function Base.show(io::IO, mime::MIME"text/plain", priors::Priors)
    println(io, "Priors:")
    for (k, prior) in priors.priors
        print(io, @sprintf("%20s ~ ", k))
        show(io, mime, prior)
        print(io, "\n")
    end
end

"""
    Derived(key=expr, ...)

Zero or more named expressions computed from the priors (and from each
other, in order). Expressions must be pure and differentiable.
`captured_names`/`captured_vals` carry values interpolated into the block
with `\$` so they become compile-time constants rather than global lookups.
"""
struct Derived
    variables::OrderedDict{Symbol}
    captured_names::Tuple
    captured_vals::Tuple
end
function Derived(; _captured=nothing, variables...)
    captured_names = isnothing(_captured) ? () : _captured[1]
    captured_vals = isnothing(_captured) ? () : _captured[2]
    return Derived(OrderedDict(variables), captured_names, captured_vals)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize det::Derived)
    print(io, "Derived:\n  ")
    for (k, expr) in det.variables
        print(io, @sprintf("%20s = %s\n", k, string(expr)))
    end
    print(io, "\n")
end

# ---------------------------------------------------
# Parameterizations
# ---------------------------------------------------

abstract type Parameterization end

"""
    UniformCircular(domain=2π)

A variable on a continuous periodic domain. Introduces two normal variables
`<name>x`, `<name>y` and derives the angle as `atan(y, x)`, so the sampler
can wrap freely instead of hitting a hard boundary.

    @variables begin
        Ω ~ UniformCircular()
    end
"""
struct UniformCircular <: Parameterization
    domain::Float64
end
UniformCircular() = UniformCircular(2π)
export UniformCircular

# Each entry of a `@variables` block expands to some priors, some derived
# variables, and possibly some prior-shaped likelihood terms.
expandparam(var, d::Distribution) = (priors=[(var, d)], derived=[], likelihoods=[])
expandparam(var, n::Number) = (priors=[], derived=[(var, n)], likelihoods=[])
expandparam(var, f::Expr) = (priors=[], derived=[(var, f)], likelihoods=[])

function expandparam(var::Symbol, p::UniformCircular)
    varx = Symbol("$(var)x")
    vary = Symbol("$(var)y")
    return (
        priors=[(varx, Normal(0, 1)), (vary, Normal(0, 1))],
        derived=[(var, :(atan($vary, $varx) / 2π * $(p.domain)))],
        # Without a prior on the vector's length the pair is pinched at the
        # origin, where the angle is undefined.
        likelihoods=[UnitLengthPrior{varx,vary}(varx, vary)],
    )
end

# ---------------------------------------------------

# ---------------------------------------------------
# Identifier hygiene for names that become NamedTuple keys (from CSV.jl).
# ---------------------------------------------------

const RESERVED = Set(["local", "global", "export", "let",
    "for", "struct", "while", "const", "continue", "import",
    "function", "if", "else", "try", "begin", "break", "catch",
    "return", "using", "baremodule", "macro", "finally",
    "module", "elseif", "end", "quote", "do"])

function normalizename(name::String)::Symbol
    uname = strip(Base.Unicode.normalize(name))
    id = Base.isidentifier(uname) ? uname : map(c -> Base.is_id_char(c) ? c : '_', uname)
    cleansed = string((isempty(id) || !Base.is_id_start_char(id[1]) || id in RESERVED) ? "_" : "", id)
    return Symbol(replace(cleansed, r"(_)\1+" => "_"))
end
normalizename(name::Symbol) = normalizename(String(name))
