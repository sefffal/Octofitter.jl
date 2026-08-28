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
_refname(b::Body) = b.name

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

# --- the interim system -------------------------------------------------------
#
# A deferred system line can reach a `PlanetOrbits.System` built from this
# model's own bodies *before* the frame those lines are computing exists. The
# enabling fact is that body-vs-barycentre kinematics need a distance, not a
# sky position: they are fully defined at `Parallax` level, so the interim
# solve never depends on the frame it is helping to construct.
#
# Built by codegen from the same body declarations as the final system, so the
# two cannot drift apart — which is the reason this lives in the code
# generator rather than being documented as a user-space pattern.

"""
Name a **deferred** system line uses to reach the frame-free interim system.
See the "Deferred lines" section of [`System`](@ref) and [`AnchoredFrame`](@ref).
"""
const INTERIM_SYSTEM_VAR = :system_interim

"""
Optional system variable designating the parallax the interim system is built
at. Needed only for *angular* interim observables (`raoff`, `pmra`, …); the
physical-unit ones (`posx`, `velx`, …) are frame-free.

Falls back to `plx` when that is itself non-deferred. In the anchored
parameterization `plx` is derived *from* the interim and so cannot be, which
is exactly why this second name exists: write `plx_interim = plx_A` beside the
anchor's sampled parallax.
"""
const INTERIM_PARALLAX_VAR = :plx_interim

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

A system may contain a single `Body` and no orbits: the lone body *is* the
system barycentre, and its motion is purely the frame's. That is the natural
model for fitting the absolute astrometry of an isolated star (e.g. the
first [`HipparcosIADObs`](@ref) tutorial) — no zero-mass companion needed.

# Deferred lines
A system line that mentions a body by name (`mut_inc = b.i - c.i`) is
evaluated *after* every body block, so couplings, period ratios and global
`LL +=` constraints need no new syntax. Bodies look up and outward; the
system block looks down and inward, after the fact; siblings never see each
other.

Deferred lines also see `system_interim`: a `PlanetOrbits.System` built from
this model's own bodies, which they may solve. It carries no absolute frame —
which is the point, since what deferred lines most usefully compute *is* the
frame:

    plx_interim = plx_A                                          # not deferred
    Δ    = anchor_offsets(system_interim, :A, 57388.0)           # deferred
    pmra = pmra_A - Δ.pmra                                       # …and so is this

Body–barycentre kinematics need a distance, not a sky position, so the interim
is complete at `Parallax` level and there is no circularity. Give it that
distance with a non-deferred `plx_interim` (or `plx`); without one it is
frame-free and only the physical-unit observables (`posx`, `velx`, …) work.
[`AnchoredFrame`](@ref) is the packaged form of this pattern. `system_interim`
is visible *only* in deferred system lines: a body block cannot see it (it is
built from the bodies), and an observation block must not, because
observations are evaluated against the final system rather than this one.
Reading it anywhere else is an error naming the block.

# Observations
All observations live here. Each one names its own references
(`target=`, `ref=`), so there is nothing left for per-companion attachment
to express.

# Corrections
`observing_geometry` and `barycentric_lighttime` gate PlanetOrbits' two
source-side precision corrections. Each takes `:auto` (the default), `:on` or
`:off` — `true`/`false` are accepted as aliases for the latter two.

`:auto` measures rather than guesses: at build, it draws ~300 parameter sets
from the priors, solves each with the correction on and off, and compares
**each observation's own model predictions** against that observation's own
tightest uncertainty. The correction is declined only if the resulting
*accumulated* bias — per-point impact × √(number of data points), since these
corrections are common-mode and so grow as √n rather than averaging down —
stays under 0.1σ for every observation. Needs are OR'd, so one µas dataset
keeps it on for the whole model. A model in which no observation type can
report comparable predictions resolves `:on` without taking a single draw.

The decision is taken once, here — never per draw, which would make the log
density discontinuous — and recorded on the built system as `sys.corrections`,
from where it travels into chain metadata. Tune with `correction_draws`,
`correction_seed`, `correction_threshold`, and silence the build log with
`verbosity=0`.

`einstein_rv` (`:on` by default, or `:off`) is the third and is *not*
tri-state: it chooses whether radial velocities are predicted from the
spectroscopic `radvel` or the kinematic `velz`. That is a counterfactual, for
measuring what the Einstein term does to a posterior — not something the data
can rule out, and not a provenance declaration, since no pipeline can have
removed the part of the term that varies with the orbit.

See the "How Octofitter Computes Orbits" manual page, and
[`recheck_corrections`](@ref) for the after-sampling check.
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
    # The *resolved* booleans, whatever the user wrote. `:auto` is decided
    # once, here, by the prior-predictive impact test in `corrections.jl`;
    # everything downstream (codegen's literal interpolation, cross-validation,
    # the plotting API) keeps reading a plain `Bool`.
    observing_geometry::Bool
    barycentric_lighttime::Bool
    # Whether radial velocities are predicted from the spectroscopic `radvel`
    # (the default) or the kinematic `velz`. A counterfactual physics switch,
    # not a provenance declaration: no pipeline can have removed the
    # orbit-varying part of the Einstein term, so there is nothing per-series
    # to declare and it belongs here beside the other two solve settings.
    einstein_rv::Bool
    # How those booleans came to be what they are, kept with the model so a
    # posterior is never separable from the corrections that produced it.
    corrections::CorrectionReport
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
                observing_geometry=:auto,
                barycentric_lighttime=:auto,
                einstein_rv=:on,
                correction_draws::Int=CORRECTION_DRAWS,
                correction_seed::UInt64=CORRECTION_SEED,
                correction_threshold::Float64=CORRECTION_BIAS_THRESHOLD,
                verbosity::Int=1)
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
    for n in nodes
        if n isa Body
            n.about === nothing || push!(rows, (n.name, (n.name,), n.about))
        else
            push!(rows, (n.name, n.exterior, n.about))
        end
    end
    # A body is *placed* by whichever row has it on the exterior side. That is
    # not the same as "has an `about=`": in a 2+2 quadruple, `Ba` omits
    # `about=` and is placed by the wide `Orbit` node's `exterior=(Ba, Bb)`.
    # Testing placement rather than the keyword is what makes set exteriors
    # work at all.
    placed = Symbol[nm for (_, ext, _) in rows for nm in ext]
    roots = setdiff(bodynames, placed)
    isempty(roots) && error(
        "System $name: every body is placed by some orbit, so the hierarchy has no root. " *
        "Exactly one body — normally the host star — must be left unplaced.")
    length(roots) == 1 || error(
        "System $name: bodies $(join(roots, ", ")) are not placed by any orbit, but a " *
        "system has exactly one root. Give the others an `about=`, or list them in an " *
        "`Orbit` node's `exterior=`.")
    # Note there is deliberately no "each body appears in exactly one
    # exterior" rule: in a 2+2 quadruple `Bb` appears both as the exterior of
    # its own tight row and as a member of the wide row's *set* exterior,
    # which places the pair's barycentre rather than `Bb` itself. Redundancy
    # that this cannot see shows up as a singular hierarchy matrix, which
    # PlanetOrbits reports with the offending rows listed.
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
    #
    # Deferral is **transitive**: a system line that reads a deferred system
    # variable must itself wait. Without the fixed point, `f = g ? 0 : 1`
    # where `g = u < b.a` is evaluated in step 1 while `g` only exists after
    # step 3, and the failure is an `UndefVarError` naming a variable the user
    # did write — at evaluation time, inside the sampler.
    #
    # Reading `system_interim` defers a line for the same reason: the interim
    # system is built from the bodies, so it does not exist until they do.
    deferred = Symbol[]
    let changed = true
        while changed
            changed = false
            for (k, expr) in derived.variables
                k in deferred && continue
                if _mentions_node(expr, allnames) || _mentions_symbol(expr, deferred) ||
                   _mentions_symbol(expr, (INTERIM_SYSTEM_VAR,))
                    push!(deferred, k)
                    changed = true
                end
            end
        end
    end
    # Back into declaration order: the generated code evaluates them in the
    # order of this vector, and a later line may read an earlier one.
    deferred = Symbol[k for k in keys(derived.variables) if k in deferred]
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
    # A node block sees only `system.*` — never a sibling. Writing a mass ratio
    # the natural way (`mass = q * A.mass` inside b's block) is the commonest
    # way to discover that, and without this check it surfaces as
    # `UndefVarError: A not defined in Octofitter`, thrown out of a
    # `RuntimeGeneratedFunction` at `LogDensityModel` construction time — which
    # names neither the block nor the rule.
    for n in nodes, (k, expr) in n.derived.variables
        siblings = Tuple(m for m in allnames if m !== n.name)
        _mentions_node(expr, siblings) && error(
            "System $name: $(n.name)'s `$k` refers to another node by name, but a node " *
            "block cannot see a sibling — each one is evaluated with only `system.*` in " *
            "scope, so that the evaluation order is well defined. Either hoist the " *
            "shared quantity into the system block and read it as `system.…`, or, if it " *
            "is a dynamical mass, express it through the hierarchy with `about=` rather " *
            "than by arithmetic. (Nodes here: $(join(allnames, ", ")).)")
    end

    # `system_interim` is a *deferred system line* scope only. It is built from
    # the bodies, so a body block asking for it is the same cycle the check
    # above catches; and an observation block asking for it would be reading a
    # system that has since been reframed, which is a different object from the
    # one the observation is actually evaluated against.
    for n in nodes, (k, expr) in n.derived.variables
        _mentions_symbol(expr, (INTERIM_SYSTEM_VAR,)) && error(
            "System $name: $(n.name)'s `$k` reads `$INTERIM_SYSTEM_VAR`, which only " *
            "exists in *deferred system lines*. The interim system is built from the " *
            "bodies, so it cannot be visible inside a body block. Move the calculation " *
            "to the system block, where it will be deferred automatically, and read it " *
            "back here as `system.…` — unless that is itself circular, in which case " *
            "what you want is a hierarchy row (`about=`), not arithmetic.")
    end
    for o in obs
        hasproperty(o, :derived) || continue
        for (k, expr) in o.derived.variables
            _mentions_symbol(expr, (INTERIM_SYSTEM_VAR,)) && error(
                "System $name: observation \"$(likelihoodname(o))\"'s `$k` reads " *
                "`$INTERIM_SYSTEM_VAR`, which only exists in deferred system lines. " *
                "Observations are evaluated against the *final* system — the one the " *
                "deferred lines finished building — so the interim is never the right " *
                "object here. Hoist the quantity into a system line and read it as " *
                "`system.…`.")
        end
    end
    if any(k -> _mentions_symbol(derived.variables[k], (INTERIM_SYSTEM_VAR,)), deferred)
        INTERIM_PARALLAX_VAR in deferred && error(
            "System $name: `$INTERIM_PARALLAX_VAR` is deferred (it mentions a body, or " *
            "reads `$INTERIM_SYSTEM_VAR`), but it is what the interim system is *built " *
            "at* — so it must be known before the bodies are. Define it from sampled " *
            "quantities only; in the anchored parameterization that is the anchor " *
            "source's own catalog parallax, e.g. `$INTERIM_PARALLAX_VAR = plx_A`.")
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
    # …and so must the bands a photocentre names. A typo'd band is otherwise
    # only found on the first likelihood evaluation, inside the sampler.
    declaredbands = unique(Symbol[_flux_band(v) for n in bodynodes for v in _flux_vars(n)])
    for o in obs, s in refspecs(o), band in _refbands(s)
        band in declaredbands || error(
            "System $name: observation \"$(likelihoodname(o))\" asks for the :$band " *
            "photocentre, but no body declares a `flux_$band` variable" *
            (isempty(declaredbands) ? " (no body declares any flux)." :
             " (the bands declared here are $(join(declaredbands, ", "))).") *
            " Photocentre weights come from the bodies' fluxes.")
    end
    # A photocentre with no band named picks the system's sole band, so it
    # needs exactly one to exist. `_refbands` returns `()` for that spec (it
    # names no band to check), which left the commonest spelling of all —
    # `GaiaDR4AstromObs`'s default `target=Photocentre` — falling through to a
    # PlanetOrbits error on the first likelihood evaluation, from inside the
    # `LogDensityModel` constructor, phrased in PlanetOrbits' own
    # `Body(…, flux=(G=1.0,))` spelling rather than the `@variables` one.
    for o in obs, s in refspecs(o)
        s isa PhotocentreSpec{nothing} || continue
        if isempty(declaredbands)
            error("System $name: observation \"$(likelihoodname(o))\" observes a " *
                  "`Photocentre`, but no body declares a flux. Photocentre weights " *
                  "come from the bodies' fluxes, so give at least one body a `flux` " *
                  "variable in its `@variables` block — e.g. `flux = 1.0` on the host " *
                  "and `flux = 0.0` on a dark companion, or `flux_G`/`flux_H` etc. for " *
                  "a named band.")
        elseif length(declaredbands) > 1
            error("System $name: observation \"$(likelihoodname(o))\" observes a " *
                  "`Photocentre` with no band, but this system declares several " *
                  "($(join(declaredbands, ", "))). Name one, e.g. " *
                  "`Photocentre(:$(first(declaredbands)))`.")
        end
    end
    # Settings that are only wrong in company — see `check_siblings`.
    let ctx = (; name=Symbol(name), framevars, sysnames)
        for o in obs
            check_siblings(o, obs, ctx)
        end
    end
    obsnames = String[likelihoodname(o) for o in obs]
    length(unique(obsnames)) == length(obsnames) || error(
        "System $name: duplicate observation name(s). Each observation needs a unique " *
        "`name=`, since it labels the observation's own variables in the chain.")

    og = _normalize_flag(observing_geometry, "observing_geometry")
    blt = _normalize_flag(barycentric_lighttime, "barycentric_lighttime")
    # No `:auto` here: the Einstein term is not a precision tier that data can
    # rule out — it is what `radvel` *means*. `:off` exists to measure what it
    # does to a posterior, which is a question only the user can be asking.
    ein = _normalize_flag(einstein_rv, "einstein_rv")
    ein === :auto && error(
        "`einstein_rv` takes `:on` (the default) or `:off`, not `:auto`. The " *
        "Einstein term cannot be ruled out by the data the way a geometry " *
        "correction can: no reduction pipeline can have removed the part of it " *
        "that varies with the orbit, and the constant part is absorbed by each " *
        "series' offset either way. `:off` is an experiment — it predicts from " *
        "the kinematic `velz` — not a declaration about your data. The " *
        "correction report prices the term for every RV series regardless.")

    # The impact test needs a fully-built `System` to draw from and solve, so
    # build one at the accurate setting first, measure, and rebuild with the
    # verdict. Both are cheap immutable values; nothing is mutated in place.
    proto = System(Symbol(name), nodes, obs, priorterms, priors, derived,
        method, true, true, ein === :on, CorrectionReport(),
        bodynames, rows, deferred, framevars)
    report = _impact_report(proto, og, blt;
        ndraws=correction_draws, minsuccess=CORRECTION_MIN_SUCCESS,
        threshold=correction_threshold, seed=correction_seed, verbosity)
    _log_corrections(name, report, verbosity)

    return System(Symbol(name), nodes, obs, priorterms, priors, derived,
        method, _resolved_flag(report, :observing_geometry),
        _resolved_flag(report, :barycentric_lighttime), ein === :on, report,
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

"""Whether the expression mentions any of `names` as a bare identifier."""
_mentions_symbol(x, names) = x isa Symbol && x in names
function _mentions_symbol(ex::Expr, names)
    # `a.b` is a field access, not a use of the local `b`.
    ex.head === :. && return _mentions_symbol(ex.args[1], names)
    return any(a -> _mentions_symbol(a, names), ex.args)
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

"""Whether any deferred system line reads the interim system."""
_uses_interim(sys::System) = any(sys.deferred) do k
    _mentions_symbol(sys.derived.variables[k], (INTERIM_SYSTEM_VAR,))
end

"""
Non-deferred system variable the interim system takes its parallax from, or
`nothing` for a frame-free (`NoFrame`) interim. See [`INTERIM_PARALLAX_VAR`](@ref).
"""
function _interim_plx(sys::System)
    # Priors are all consumed before any derived line, so a prior name is
    # never deferred and needs no check of its own.
    have = k -> (k in keys(sys.priors.priors)) ||
                (k in keys(sys.derived.variables) && !(k in sys.deferred))
    have(INTERIM_PARALLAX_VAR) && return INTERIM_PARALLAX_VAR
    have(:plx) && return :plx
    return nothing
end
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
        # `_refdesc`, not a flat "a vs b vs c": a type that overrides it says
        # something true about its references (`b < c` for an ordering prior,
        # `b.flux_H` for photometry), and the flat join says something false.
        desc = _refdesc(o)
        println(io, "  ", typeof(o).name.name, " \"", likelihoodname(o), "\"",
            isempty(desc) ? "" : "  " * desc)
    end
    for (owner, t) in sys.priorterms
        println(io, "  ", typeof(t).name.name, " [", owner, "]")
    end
    return nothing
end
Base.show(io::IO, @nospecialize sys::System) =
    print(io, "System model ", sys.name, " (", nbodies(sys), " bodies, ",
        length(sys.observations), " observations)")
