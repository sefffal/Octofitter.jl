# ---------------------------------------------------
# Model nodes: Body, Orbit, System
#
# A model is a flat list of nodes plus a list of observations. Each node owns
# a variable namespace (`@variables`); each observation names the references
# it measures. There is no per-companion observation scope and no `basis=`:
# what used to be `Visual{KepOrbit}` vs `AbsoluteVisual{KepOrbit}` vs
# `RadialVelocityOrbit` is now decided by which frame variables the system
# block supplies, and what used to be a hand-rolled parametrization is a
# constructor group on `PlanetOrbits.Orbit`.
#
# The 1:1 correspondence with PlanetOrbits is the teaching device:
#   Octofitter.Body  → one PlanetOrbits.Body (+ one PlanetOrbits.Orbit, if `about=`)
#   Octofitter.Orbit → one PlanetOrbits.Orbit
#   Octofitter.System → one PlanetOrbits.System, rebuilt every sample
# ---------------------------------------------------

"""
    Body(; name, about=nothing, variables)

One gravitating body in the model, and — when `about` is given — the orbit
that places it. `about` takes the same grammar `PlanetOrbits.Orbit` uses:

    A = Body(name="A", variables=@variables begin mass ~ … end)          # root
    b = Body(name="b", about=A,      variables=@variables begin … end)   # astrocentric
    c = Body(name="c", about=(A, b), variables=@variables begin … end)   # Jacobi

Exactly one body may omit `about`; it is the root of the hierarchy. The
convention is never inferred from semi-major axis — `about=A` and
`about=(A, b)` are materially different models under `KeplerianApprox`.

# Variables
Names in the block are read by position in three groups:

  - `mass` [M⊙] — the body's mass. Zero is allowed (a test particle).
  - `flux`, or `flux_<band>` — per-band flux, used to form photocentres.
    Setting the host's flux to 1.0 makes every other body's flux a contrast
    ratio.
  - orbital elements — any keyword `PlanetOrbits.Orbit` accepts (`a` or `P`;
    `e`/`ω` or `secosω`/`sesinω` or `ecosω`/`esinω`; `i`, `Ω`; `tp` or
    `M0`+`epoch` or `θ`+`epoch`; or a full Cartesian state). Supplying the
    wrong combination is a mechanical error from the constructor.

Every other name is an ordinary local, free to be used by later lines in the
same block (e.g. `mass_jup ~ LogUniform(…)` then `mass = mass_jup * mjup2msol`).
A body block sees `system.*`; it does not see its siblings.
"""
struct Body{TP<:Priors,TD<:Derived,TL<:Tuple}
    name::Symbol
    about::Union{Nothing,Tuple{Vararg{Symbol}}}
    priors::TP
    derived::TD
    likelihoods::TL
end

function Body(; name::Union{Symbol,AbstractString},
              about=nothing,
              variables::Tuple=(Priors(), Derived()))
    (priors, derived, extra...) = variables
    priors::Priors
    derived::Derived
    foreach(l -> l::AbstractObs, extra)
    return Body(Symbol(name), _aboutspec(about), priors, derived, extra)
end

"""
    Orbit(; name, exterior, about, variables)

One hierarchy row whose exterior is not a single body — the case
[`Body`](@ref)'s fused `about=` cannot express:

    Orbit(name="AB", exterior=(Ba, Bb), about=(Aa, Ab),
          variables=@variables begin P ~ …; e ~ … end)

`name` gives the row its own chain columns (`AB_P`, `AB_e`), which is how
binary-star people talk about that orbit. Its variables are read exactly as a
`Body`'s orbital elements are; it has no `mass` or `flux` of its own, since
every mass in the row already belongs to a body.
"""
struct Orbit{TP<:Priors,TD<:Derived,TL<:Tuple}
    name::Symbol
    exterior::Tuple{Vararg{Symbol}}
    about::Tuple{Vararg{Symbol}}
    priors::TP
    derived::TD
    likelihoods::TL
end

function Orbit(; name::Union{Symbol,AbstractString},
               exterior, about,
               variables::Tuple=(Priors(), Derived()))
    (priors, derived, extra...) = variables
    priors::Priors
    derived::Derived
    foreach(l -> l::AbstractObs, extra)
    ext = _aboutspec(exterior)
    int = _aboutspec(about)
    ext === nothing && error("Orbit $name: `exterior` is required")
    int === nothing && error("Orbit $name: `about` is required")
    return Orbit(Symbol(name), ext, int, priors, derived, extra)
end

export Body, Orbit

const ModelNode = Union{Body,Orbit}

# `about=` / `exterior=` accept a node, a Symbol, or a tuple of either.
_aboutspec(::Nothing) = nothing
_aboutspec(b::ModelNode) = (b.name,)
_aboutspec(s::Symbol) = (s,)
_aboutspec(s::AbstractString) = (Symbol(s),)
_aboutspec(t::Tuple) = map(_aboutspec1, t)
_aboutspec(v::AbstractVector) = _aboutspec(Tuple(v))
@noinline _aboutspec(@nospecialize x) = error(
    "an orbit endpoint must be a `Body` model node, a `Symbol` naming one, or a tuple " *
    "of those meaning their barycentre (e.g. `about=(A, b)`); got a value of type $(typeof(x))")
_aboutspec1(b::ModelNode) = b.name
_aboutspec1(s::Symbol) = s
_aboutspec1(s::AbstractString) = Symbol(s)
@noinline _aboutspec1(@nospecialize x) = error(
    "orbit endpoint members must be `Body` model nodes or `Symbol`s; got $(typeof(x))")

# Anywhere a reference is expected, a model node names its body.
refspec(b::Body) = BodyRefSpec{b.name}()

# --- variable-name bookkeeping ------------------------------------------------

varnames(n::ModelNode) = (collect(keys(n.priors.priors))...,
                          collect(keys(n.derived.variables))...)
priornames(n::ModelNode) = collect(keys(n.priors.priors))

# Element keywords forwarded verbatim to `PlanetOrbits.Orbit`. Anything else in
# a node's block is a local — users need intermediates, so unknown names are
# not an error here; the constructor's exactly-one-of-each-group validation is
# what catches a typo'd element (`aa=` gives "supply exactly one of `a` or `P`").
const ORBIT_ELEMENT_KEYWORDS = (
    :a, :P,                                     # size
    :e, :ω, :secosω, :sesinω, :ecosω, :esinω,   # shape
    :i, :Ω,                                     # orientation
    :tp, :M0, :θ, :epoch,                       # phase
    :x, :y, :z, :vx, :vy, :vz,                  # Cartesian state
    :M,                                         # compatibility override
)

# Frame variables read off the *system* block. Which ones are present chooses
# the frame — none, parallax-only, or a full absolute frame — so `basis=` has
# nothing left to say.
const FRAME_VARIABLES = (:plx, :ra, :dec, :pmra, :pmdec, :rv, :ref_epoch)

_element_vars(n::ModelNode) = filter(in(ORBIT_ELEMENT_KEYWORDS), varnames(n))
_flux_vars(n::Body) = filter(v -> v === :flux || startswith(string(v), "flux_"), varnames(n))
_flux_band(v::Symbol) = v === :flux ? :default : Symbol(string(v)[6:end])

# ---------------------------------------------------
# System
# ---------------------------------------------------

"""
    System(; name, bodies, observations=(), variables)

A model of one system: its bodies (and any explicit [`Orbit`](@ref) rows),
its observations, and its system-level variables.

    sys = System(
        name="HD 12345",
        bodies=[A, b],
        observations=[astrom, rvs],
        variables=@variables begin
            plx ~ Uniform(1, 100)
        end
    )

# Frame
Which of `plx`, `ra`, `dec`, `pmra`, `pmdec`, `rv`, `ref_epoch` the system
block defines chooses the frame: none gives physical-unit observables only,
`plx` alone gives angular observables in mas, and the full set gives a
rigorously propagated absolute frame. There is no `basis=` keyword.

# Deferred lines
A system line that mentions a body by name (`mut_inc = b.i - c.i`) is
evaluated *after* every body block, so couplings, period ratios and global
`LL +=` constraints need no new syntax. Bodies look up and outward; the
system block looks down and inward, after the fact; siblings never see each
other.

# Observations
All observations live here. Each one names its own references
(`target=`, `ref=`), so there is nothing left for per-companion attachment
to express.
"""
struct System{TN<:Tuple,TO<:Tuple,TT<:Tuple,TP<:Priors,TD<:Derived,TM}
    name::Symbol
    nodes::TN               # Body and Orbit nodes, in declaration order
    observations::TO        # real observations, in declaration order
    # Prior-shaped terms emitted by `@variables` (`lhs ~ dist`, `LL +=`,
    # `UniformCircular`), paired with the namespace they read from: a node's
    # name, or `:system`. They own no parameters of their own.
    priorterms::TT
    priors::TP
    derived::TD
    # --- solve configuration (see PlanetOrbits' `orbitsolve!`) ---
    method::TM
    observing_geometry::Bool
    barycentric_lighttime::Bool
    # --- derived at construction, all of it structural (no sample data) ---
    bodynames::Vector{Symbol}       # PlanetOrbits body order
    rows::Vector{Tuple{Symbol,Tuple{Vararg{Symbol}},Tuple{Vararg{Symbol}}}}
                                    # (owning node name, exterior names, interior names)
    deferred::Vector{Symbol}        # system derived vars evaluated after the body blocks
    framevars::Vector{Symbol}       # subset of FRAME_VARIABLES the system block defines
end

function System(; name::Union{Symbol,AbstractString},
                bodies,
                observations=(),
                variables::Tuple=(Priors(), Derived()),
                method=PlanetOrbits.KeplerianApprox(),
                observing_geometry::Bool=true,
                barycentric_lighttime::Bool=true)
    (priors, derived, extra...) = variables
    priors::Priors
    derived::Derived
    nodes = Tuple(bodies)
    for n in nodes
        n isa ModelNode || error(
            "`bodies` takes `Body` and `Orbit` model nodes; got a value of type $(typeof(n)). " *
            "(Observations go in `observations=`, not here.)")
    end
    obs = Tuple(observations)
    for o in obs
        o::AbstractObs
    end
    priorterms = reduce(nodes; init=map(l -> (:system => l), extra)) do acc, n
        (acc..., map(l -> (n.name => l), n.likelihoods)...)
    end

    bodynodes = filter(n -> n isa Body, nodes)
    bodynames = Symbol[n.name for n in bodynodes]
    allnames = Symbol[n.name for n in nodes]
    if length(unique(allnames)) != length(allnames)
        dups = [n for n in unique(allnames) if count(==(n), allnames) > 1]
        error("System $name: duplicate node name(s) $(join(dups, ", ")). Body and Orbit " *
              "names share one namespace and give the chain its column names, so they " *
              "must be unique.")
    end

    # Rows: one per `Body` with `about=`, plus every explicit `Orbit` node.
    rows = Tuple{Symbol,Tuple{Vararg{Symbol}},Tuple{Vararg{Symbol}}}[]
    roots = Symbol[]
    for n in nodes
        if n isa Body
            if n.about === nothing
                push!(roots, n.name)
            else
                push!(rows, (n.name, (n.name,), n.about))
            end
        else
            push!(rows, (n.name, n.exterior, n.about))
        end
    end
    isempty(roots) && error(
        "System $name: every body has an `about=`, so the hierarchy has no root. " *
        "Exactly one body — normally the host star — must omit it.")
    length(roots) == 1 || error(
        "System $name: bodies $(join(roots, ", ")) all omit `about=`, but a system has " *
        "exactly one root. Give the others an `about=`.")
    for (owner, ext, int) in rows, nm in (ext..., int...)
        nm in bodynames || error(
            "System $name: node $owner refers to :$nm, which is not a `Body` in this " *
            "system (its bodies are $(join(bodynames, ", "))).")
    end
    length(rows) == length(bodynames) - 1 || error(
        "System $name: $(length(bodynames)) bodies need exactly $(length(bodynames)-1) " *
        "orbits to determine every body's motion; got $(length(rows)). Add an explicit " *
        "`Orbit` node, or drop one.")

    # A system line mentioning a node by name is deferred until after the body
    # blocks (§8.4 tier 1). Detection is a static walk over the *unevaluated*
    # expressions the `@variables` macro already stored, and only qualified
    # access (`b.i`) counts — a bare `b` is far more likely to be a local.
    deferred = Symbol[]
    for (k, expr) in derived.variables
        _mentions_node(expr, allnames) && push!(deferred, k)
    end
    # A body block reading a deferred system variable is a cycle: the system
    # line needs the body, and the body needs the system line.
    for n in nodes, (k, expr) in n.derived.variables
        for d in _system_refs(expr)
            d in deferred && error(
                "System $name: circular definition — the system variable `$d` is " *
                "deferred because it refers to a body, but $(n.name)'s `$k` reads " *
                "`system.$d`. Break the cycle by hoisting `$d`'s inputs to the " *
                "system block.")
        end
    end

    sysnames = Symbol[collect(keys(priors.priors))...; collect(keys(derived.variables))...]
    framevars = Symbol[v for v in FRAME_VARIABLES if v in sysnames]
    if :plx in framevars && length(framevars) > 1 && length(framevars) < length(FRAME_VARIABLES)
        missing_ = setdiff(FRAME_VARIABLES, framevars)
        error("System $name: an absolute frame needs all of " *
              "$(join(FRAME_VARIABLES, ", ")); missing $(join(missing_, ", ")). " *
              "Supply them, or define `plx` alone for a parallax-only model.")
    end

    # Reference names used by observations must exist.
    for o in obs, s in refspecs(o), nm in _refnames(s)
        nm in bodynames || error(
            "System $name: observation \"$(likelihoodname(o))\" references :$nm, which " *
            "is not a `Body` in this system (its bodies are $(join(bodynames, ", "))).")
    end
    obsnames = String[likelihoodname(o) for o in obs]
    length(unique(obsnames)) == length(obsnames) || error(
        "System $name: duplicate observation name(s). Each observation needs a unique " *
        "`name=`, since it labels the observation's own variables in the chain.")

    return System(Symbol(name), nodes, obs, priorterms, priors, derived,
        method, observing_geometry, barycentric_lighttime,
        bodynames, rows, deferred, framevars)
end

export System

# --- static expression walks --------------------------------------------------

_mentions_node(x, names) = false
function _mentions_node(ex::Expr, names)
    if ex.head === :. && ex.args[1] isa Symbol && ex.args[1] in names
        return true
    end
    return any(a -> _mentions_node(a, names), ex.args)
end

"""Names `Y` such that the expression reads `system.Y`."""
_system_refs(x) = Symbol[]
function _system_refs(ex::Expr)
    out = Symbol[]
    if ex.head === :. && ex.args[1] === :system && ex.args[2] isa QuoteNode &&
       ex.args[2].value isa Symbol
        push!(out, ex.args[2].value)
    end
    for a in ex.args
        append!(out, _system_refs(a))
    end
    return out
end

# --- accessors ----------------------------------------------------------------

nbodies(sys::System) = length(sys.bodynames)
nrows(sys::System) = length(sys.rows)
bodynodes(sys::System) = filter(n -> n isa Body, sys.nodes)
_node(sys::System, name::Symbol) = sys.nodes[findfirst(n -> n.name === name, sys.nodes)]

"""
    list_parameter_names(system)

Flat vector of the model's free parameter names, in the order the sampler
sees them: system, then each body/orbit node, then each observation. Node and
observation variables are prefixed with their owner's name.
"""
function list_parameter_names(sys::System)
    out = Symbol[]
    append!(out, keys(sys.priors.priors))
    for n in sys.nodes, k in keys(n.priors.priors)
        push!(out, Symbol(n.name, :_, k))
    end
    for o in sys.observations
        hasproperty(o, :priors) || continue
        for k in keys(o.priors.priors)
            push!(out, Symbol(normalizename(likelihoodname(o)), :_, k))
        end
    end
    return out
end

function _count_epochs(sys::System)::Int
    n = 0
    for o in sys.observations
        n += hasproperty(o, :table) ? Tables.rowcount(o.table) : 1
    end
    return n
end

# --- display ------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", @nospecialize node::Body)
    println(io, "Body ", node.name,
        node.about === nothing ? "  (root)" : "  about " * _fmtmembers(node.about))
    _show_vars(io, node)
end
function Base.show(io::IO, ::MIME"text/plain", @nospecialize node::Orbit)
    println(io, "Orbit ", node.name, "  ", _fmtmembers(node.exterior),
        " about ", _fmtmembers(node.about))
    _show_vars(io, node)
end
_fmtmembers(t::Tuple{Symbol}) = string(t[1])
_fmtmembers(t::Tuple) = "barycentre(" * join(t, ", ") * ")"

function _show_vars(io::IO, node)
    for (k, d) in node.priors.priors
        println(io, @sprintf("  %14s ~ %s", k, d))
    end
    for (k, e) in node.derived.variables
        println(io, @sprintf("  %14s = %s", k, e))
    end
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize sys::System)
    println(io, "System model ", sys.name, " — ", nbodies(sys), " bodies, ",
        nrows(sys), " orbits, ", length(sys.observations), " observations")
    println(io, "  frame: ", isempty(sys.framevars) ? "none (physical units)" :
                sys.framevars == [:plx] ? "parallax" : "absolute")
    for (k, d) in sys.priors.priors
        println(io, @sprintf("  %14s ~ %s", k, d))
    end
    for (k, e) in sys.derived.variables
        println(io, @sprintf("  %14s = %s%s", k, e, k in sys.deferred ? "   [deferred]" : ""))
    end
    for n in sys.nodes
        print(io, "  ")
        show(io, MIME"text/plain"(), n)
    end
    for o in sys.observations
        println(io, "  ", typeof(o).name.name, " \"", likelihoodname(o), "\"",
            isempty(refspecs(o)) ? "" : "  " * join(map(_refstr, refspecs(o)), " vs "))
    end
    for (owner, t) in sys.priorterms
        println(io, "  ", typeof(t).name.name, " [", owner, "]")
    end
    return nothing
end
Base.show(io::IO, @nospecialize sys::System) =
    print(io, "System model ", sys.name, " (", nbodies(sys), " bodies, ",
        length(sys.observations), " observations)")
