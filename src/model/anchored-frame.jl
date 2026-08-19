# ---------------------------------------------------
# The anchored frame
#
# A system's `AbsoluteFrame` describes its **barycentre**. Catalogues describe
# **sources**. For a single star with unseen companions the difference is a
# nuisance the wide priors absorb; for two catalogue sources belonging to one
# physical system it is the entire signal, because what binds the wide pair
# together is that both sources' astrometry must be predicted from one frame
# plus each source's own modelled motion about the barycentre.
#
# The barycentric parameterization ("sample `pmra`, `pmdec`, … directly") is
# correct but badly conditioned: the frame proper motion is degenerate with
# the reflex of whichever body dominates the light, so the sampler explores a
# ridge. `G23HObs(...; frame_shift=true)` fixes that for one source by
# redefining the frame proper motion to mean "the primary's, as DR3 measured
# it" — but it does so *inside one observation*, which is why it cannot be
# right for two (see `frame_shift`'s own docstring).
#
# The anchored parameterization does the same reconditioning in the *model*,
# where it belongs: sample the anchor source's observed catalogue quantities,
# and derive the barycentric ones by subtracting the model's own motion of the
# anchor body about the barycentre. Every source's prediction then comes from
# the one shared frame, so two sources differ by exactly the relative motion
# the model puts between them — the §4 bug is structurally impossible rather
# than merely avoided.
#
# The mechanism underneath is general: `system_interim` in codegen (a
# `PlanetOrbits.System` over this model's own bodies, exposed to deferred
# system lines). This file is sugar over that, and everything it emits can be
# written by hand.
# ---------------------------------------------------

# Solve the interim system at a handful of epochs into the caller's arena. Not
# `PlanetOrbits.orbitsolve(sys, t)`, which heap-allocates a whole `Trajectory`
# — this runs once per sample inside `arr2nt`, on the sampler's hot path.
#
# `f` receives the whole `Trajectory`, so a secant window can read all three of
# its epochs from one solve.
@inline function _interim_solve(f::F, posys::PlanetOrbits.System{NB,NR,T},
                                epochs::SVector{N};
                                method, observing_geometry,
                                barycentric_lighttime) where {F,N,NB,NR,T}
    buf = Bumper.default_buffer()
    S = promote_type(T, eltype(epochs))
    @no_escape buf begin
        traj = PlanetOrbits.Trajectory(BumpAlloc(buf), S, posys, epochs)
        PlanetOrbits.orbitsolve!(traj, posys;
            method, observing_geometry, barycentric_lighttime)
        # `f` must return only scalars: everything the arena holds dies here.
        f(traj)
    end
end

@inline _interim_solve(f::F, posys::PlanetOrbits.System, epoch::Real; kwargs...) where {F} =
    _interim_solve(traj -> f(traj[1]), posys, SVector(float(epoch)); kwargs...)

"""
    anchor_offsets(system_interim, anchor, epoch)

The **anchor body's motion about the system barycentre**, at `epoch` [MJD], in
the units the catalogue quantities are in. Subtracting these from an observed
catalogue solution for the anchor source gives the corresponding barycentric
frame quantities — see [`AnchoredFrame`](@ref), which emits exactly that.

`anchor` is a body name (a `Symbol`, or a `Body` model node). Returns

| field | units | meaning |
|---|---|---|
| `ra_cosdec` | mas | Δα⋅cos δ, great-circle, i.e. what `raoff` returns |
| `dec` | mas | Δδ |
| `pmra` | mas/yr | μα\\* of the anchor relative to the barycentre |
| `pmdec` | mas/yr | μδ |
| `rv` | m/s | *spectroscopic* line-of-sight velocity (`radvel`) |
| `dz` | AU | line-of-sight displacement, +receding — for the parallax |

The first four are angular and so need the interim system to carry a parallax:
give the model a `plx_interim` variable (or a non-deferred `plx`), or they
fail with a `MethodError` on `NoFrame`. `rv` and `dz` are frame-free.

This is the one-pass primitive. The anchored map itself is
[`anchored_frame`](@ref), which calls this twice — once at `Parallax` level
and once on the frame that produced — because two terms of the reflex are
invisible to a system with no space motion. See its docstring for both, with
their measured sizes.

`method`, `observing_geometry` and `barycentric_lighttime` are PlanetOrbits'
solve settings and default to the same values `orbitsolve` uses. Note the
default propagator is Keplerian regardless of what the *model* is configured
with: pass `method=` to match an N-body model, or accept that the frame is
derived from the Keplerian approximation of the anchor's reflex (which is a
correction to a correction, and below any catalogue's precision for the
hierarchies where the two propagators visibly differ).

Two notes on exactness, both immaterial at present accuracies but recorded so
nobody has to re-derive them:

  - `rv` is `radvel`, not `velz` — a catalogue radial velocity is the
    spectroscopic quantity, so the Einstein term belongs in the subtraction.
    On a `Parallax`-level system that term is built from orbital velocities
    alone; [`anchored_frame`](@ref)'s second pass is what supplies the rest.
  - the interim is solved with no observer position, so `ra_cosdec`/`dec` are
    offsets as seen from the solar-system barycentre. That is what a catalogue
    solution reports, having already removed the annual parallax, so this is
    the right quantity rather than an approximation to one.

# `window`

`window=(t1, t2)` [MJD] replaces the *instantaneous* reflex in `ra_cosdec`,
`dec`, `pmra` and `pmdec` with the **secant** through the reflex offset at `t1`
and `t2`, evaluated at `epoch`. `rv` and `dz` are unaffected: they are not part
of the degeneracy the window exists to break, and their instantaneous values at
`epoch` are what a catalogue solution means.

The reason to want it is what the subtraction is *for*. The quantity the data's
linear terms absorb — and therefore the thing that has to come out to
decorrelate the sampler — is the average drift across the observing window, not
the derivative at one instant. The two agree at long period. At short period
they could hardly differ more: a 10 M_jup companion at P = 10 d has α ≈ 2 mas
and an instantaneous reflex proper motion of ~475 mas/yr, so subtracting the
instantaneous value pushes the anchored `pmra_A` far outside any sane prior box
and clips the short-period, high-mass region the data actually allow. The
secant → 0 there, exactly as the data's sensitivity to the wiggle-as-a-drift
does.

Choose `t1`, `t2` to bracket every dataset in the model. There is no default,
deliberately: the right window is a statement about which observations the
shear is meant to decorrelate, and a guessed one would be silently wrong.

The map stays triangular with unit diagonal — the secant is still built from
the body variables and the anchored parallax alone — so the Jacobian argument
in [`AnchoredFrame`](@ref) is unchanged.

See also [`AnchoredFrame`](@ref), `PlanetOrbits.reframe`.
"""
function anchor_offsets(posys::PlanetOrbits.System, anchor, epoch;
                        window=nothing,
                        method=PlanetOrbits.KeplerianApprox(),
                        observing_geometry::Bool=true,
                        barycentric_lighttime::Bool=true)
    a = _anchorname(anchor)
    ref = PlanetOrbits.barycentre(posys)
    isnothing(window) || return _anchor_offsets_secant(posys, a, ref, epoch, window;
        method, observing_geometry, barycentric_lighttime)
    return _interim_solve(posys, epoch; method, observing_geometry,
                          barycentric_lighttime) do sol
        (; ra_cosdec=PlanetOrbits.raoff(sol, a, ref),
           dec=PlanetOrbits.decoff(sol, a, ref),
           pmra=PlanetOrbits.pmra(sol, a, ref),
           pmdec=PlanetOrbits.pmdec(sol, a, ref),
           rv=PlanetOrbits.radvel(sol, a, ref),
           dz=PlanetOrbits.posz(sol, a, ref))
    end
end

# The secant form. One solve at (t1, epoch, t2): the window ends give the slope,
# the middle epoch gives `rv`/`dz`, and the offsets are read off the secant LINE
# at `epoch` rather than the curve — a line has one slope everywhere, so taking
# the offset from the curve and the rate from the chord would describe no single
# straight motion, and the shear would stop being a pure shear.
@inline function _anchor_offsets_secant(posys, a::Symbol, ref, epoch, window;
                                        method, observing_geometry,
                                        barycentric_lighttime)
    t1, t2 = _check_window(window)
    Δt_yr = (t2 - t1) / julian_year
    return _interim_solve(posys, SVector(float(t1), float(epoch), float(t2));
                          method, observing_geometry,
                          barycentric_lighttime) do traj
        s1, s0, s2 = traj[1], traj[2], traj[3]
        r1 = PlanetOrbits.raoff(s1, a, ref)
        d1 = PlanetOrbits.decoff(s1, a, ref)
        μα = (PlanetOrbits.raoff(s2, a, ref) - r1) / Δt_yr
        μδ = (PlanetOrbits.decoff(s2, a, ref) - d1) / Δt_yr
        τ = (epoch - t1) / julian_year
        (; ra_cosdec=r1 + μα * τ,
           dec=d1 + μδ * τ,
           pmra=μα,
           pmdec=μδ,
           rv=PlanetOrbits.radvel(s0, a, ref),
           dz=PlanetOrbits.posz(s0, a, ref))
    end
end

function _check_window(window)
    length(window) == 2 || error(
        "the anchored frame's `window` is the two epochs [MJD] the secant is " *
        "taken through, as `(t1, t2)`; got $(length(window)) value(s).")
    t1, t2 = window
    (isfinite(t1) && isfinite(t2)) || error(
        "the anchored frame's `window` epochs must be finite; got ($t1, $t2).")
    t2 > t1 || error(
        "the anchored frame's `window` must satisfy t2 > t1; got ($t1, $t2). " *
        "It is a time span in MJD, and its width divides the secant.")
    return (t1, t2)
end

"""
    anchored_frame(system_interim, anchor, ref_epoch, correct; ra, dec, plx, pmra, pmdec, rv)

The **barycentric** frame quantities implied by an anchor source's catalogue
solution, as `(; ra, dec, plx, pmra, pmdec, rv)` in `System`'s units — the map
[`AnchoredFrame`](@ref) emits, and the thing to call if you are writing the
deferred lines yourself.

The six keywords are the *inputs*: the anchor source's observed values for
whichever quantities are anchored, and the already-barycentric values for the
rest. `correct` is a six-tuple of `Bool` in `(ra, dec, plx, pmra, pmdec, rv)`
order saying which is which; an uncorrected quantity is returned unchanged.
All six are required either way, because the second pass needs a complete
frame.

# Two passes, and why

The naive map subtracts [`anchor_offsets`](@ref) from the catalogue values.
Those offsets come from the interim system, which is `Parallax`-level — it has
a distance but no space motion — while the *final* system computes the same
reflex through a full `AbsoluteFrame`. Two terms differ, and neither is
negligible for a wide pair:

  - **proper motion.** An angular offset ρ from the barycentre shrinks as the
    system recedes, at ρ⋅(ṙ/d). The interim's distance is constant, so it
    misses exactly that. Measured on a 750 au pair at 13.5 pc receding at
    27 km/s: **4.0 × 10⁻³ mas/yr**, which is 10% of a DR3 proper-motion σ.
  - **radial velocity.** The interim's Einstein term is built from orbital
    velocities alone, missing ½|v_frame|²/c — **2.40 m/s** for a 38 km/s space
    velocity, and a constant, since the frame's own speed does not vary with
    the orbit.

So the offsets are recomputed on a `PlanetOrbits.reframe`d interim carrying
the first pass's `AbsoluteFrame` — the *same* frame the final system will be
built with, so the second pass computes the reflex on a system identical to
the one the likelihood uses. What is left is the second-order feedback of the
frame on its own offsets, ~10⁻⁴ of the first-pass residual: measured at
**3 × 10⁻⁷ mas/yr and 2 × 10⁻⁴ m/s** on the 750 au pair, i.e. 10⁻⁵ of a DR3 σ.

`reframe` rather than a second `System(...)` because nothing about the orbit
may be recomputed here — see its docstring.

`refine=false` stops after the first pass. The reason to want that is the
Jacobian: the one-pass map is *exactly* triangular with unit diagonal (the
correction it subtracts is built from the body variables and the anchored
parallax, and from nothing else), while the second pass makes each barycentric
quantity depend weakly on all six inputs. Measured, that perturbs the diagonal
by 10⁻⁹–10⁻⁴ and the determinant by 10⁻⁹ — far below anything a posterior can
resolve, and much smaller than the 0.1σ interpretation error it removes, which
is why refinement is the default. See [`AnchoredFrame`](@ref)'s "Priors and the
Jacobian".
"""
function anchored_frame(posys::PlanetOrbits.System, anchor, ref_epoch,
                        correct::NTuple{6,Bool}=ntuple(_ -> true, 6);
                        ra, dec, plx, pmra, pmdec, rv, refine::Bool=true,
                        window=nothing,
                        method=PlanetOrbits.KeplerianApprox(),
                        observing_geometry::Bool=true,
                        barycentric_lighttime::Bool=true)
    a = _anchorname(anchor)
    cat = (; ra, dec, plx, pmra, pmdec, rv)
    solve = (p) -> _anchor_map(p, a, float(ref_epoch), correct, cat;
        window, method, observing_geometry, barycentric_lighttime)
    f1 = solve(posys)
    refine || return f1
    f2 = PlanetOrbits.reframe(posys; f1.ra, f1.dec, f1.plx, f1.pmra, f1.pmdec, f1.rv,
        ref_epoch=float(ref_epoch))
    return solve(f2)
end

@inline function _anchor_map(posys, a::Symbol, epoch, correct, cat; kwargs...)
    posys.frame isa PlanetOrbits.NoFrame && _err_interim_frame()
    off = anchor_offsets(posys, a, epoch; kwargs...)
    return (;
        ra = correct[1] ? cat.ra - off.ra_cosdec / 3.6e6 / cosd(cat.dec) : cat.ra,
        dec = correct[2] ? cat.dec - off.dec / 3.6e6 : cat.dec,
        plx = correct[3] ? barycentre_parallax(cat.plx, off.dz) : cat.plx,
        pmra = correct[4] ? cat.pmra - off.pmra : cat.pmra,
        pmdec = correct[5] ? cat.pmdec - off.pmdec : cat.pmdec,
        rv = correct[6] ? cat.rv - off.rv : cat.rv,
    )
end

@noinline _err_interim_frame() = error(
    "the anchored frame needs the interim system to carry a parallax: the " *
    "anchor's offset from the barycentre is subtracted in mas and mas/yr, and " *
    "an AU has no angular size without a distance. Give the system block a " *
    "non-deferred `$INTERIM_PARALLAX_VAR` (the anchor source's own catalogue " *
    "parallax is the right value) or a non-deferred `plx`.")

_anchorname(s::Symbol) = s
_anchorname(b::Body) = b.name
_anchorname(s::AbstractString) = Symbol(s)
@noinline _anchorname(@nospecialize x) = error(
    "the frame anchor must be a `Body` model node or a `Symbol` naming one; " *
    "got a value of type $(typeof(x))")

"""
    barycentre_parallax(plx_anchor, dz_au)

Barycentre parallax [mas], given the anchor source's parallax [mas] and the
anchor's line-of-sight displacement from the barycentre [AU, +receding] —
i.e. `anchor_offsets(...).dz`.

Exact rather than a series: `d_bary = d_anchor − dz` in AU, inverted back to
mas. The correction is about `dz⋅plx²/2.06e8` mas — 1 part in 10⁵ for a body
an AU from the barycentre at 100 pc, and 1 part in 300 for a 750 au wide pair
at 13 pc, which is where it starts to matter.
"""
@inline barycentre_parallax(plx_anchor, dz_au) =
    inv(inv(plx_anchor) - dz_au / (1000 * PlanetOrbits.pc2au))

export anchor_offsets, anchored_frame, barycentre_parallax

"""
    AnchoredFrame(anchor; ref_epoch, variables, window=nothing, ra=true, dec=true, plx=true, pmra=true, pmdec=true, rv=true)

A `variables=` block for [`System`](@ref) in which the absolute frame is
parameterized by an **anchor source's observed catalogue solution** rather than
by the system barycentre's.

The frame still *means* the barycentre — nothing downstream changes — but what
the sampler moves is the anchor's own (ra, dec, plx, pmra, pmdec, rv), with the
model's motion of the anchor body about the barycentre subtracted to recover
the barycentric values. That subtraction reads the interim system, so it sees
the same bodies the final system is built from.

    A = Body(name="A", variables=@variables begin mass ~ truncated(Normal(1.29, 0.02), lower=0) end)
    b = Body(name="b", about=A, variables=@variables begin … end)
    B = Body(name="B", about=A, variables=@variables begin … end)   # the wide companion

    sys = System(
        name="ups And", bodies=[A, b, B], observations=[g23h_A, g23h_B],
        variables=AnchoredFrame(A; ref_epoch=mjd("2016-01-01"), variables=@variables begin
            ra_A     ~ Normal(24.19928, 1/3.6e6)
            dec_A    ~ Normal(41.40546, 1/3.6e6)
            plx_A    ~ Normal(74.19, 0.20)
            pmra_A   ~ Normal(-172.57, 0.05)
            pmdec_A  ~ Normal(-381.32, 0.04)
            rv_A     ~ Normal(-26900.0, 500.0)
        end))

The sampled names default to `<quantity>_<anchor>`; pass a `Symbol` to any of
the six keywords to rename one, or `false` to leave that quantity
**unanchored** — in which case your block must define the barycentric `ra` /
`dec` / … itself, in the ordinary way. Mixing is legal and is the point of the
per-quantity keywords: anchoring proper motion, where the degeneracy is, while
sampling the barycentric parallax directly is a perfectly reasonable model.

`ref_epoch` [MJD] is both the frame's reference epoch and the epoch the
anchor's motion is evaluated at, because that is what makes the two
consistent: the catalogue solution being anchored to is itself referred to
that epoch. Anchoring to *several* catalogue epochs is not a thing —
one frame has one reference epoch, and a source's motion between catalogue
epochs is what the observations are for.

`window=(t1, t2)` [MJD] switches the position and proper-motion corrections
from the anchor's *instantaneous* reflex at `ref_epoch` to the **secant**
through its reflex offset at `t1` and `t2`. Bracket every dataset in the model.
Prefer it whenever the prior admits periods shorter than the observing
baseline: the instantaneous reflex proper motion diverges as 2πα/P while the
data's sensitivity to a reflex-as-drift does not, so subtracting it clips the
short-period, high-mass corner against the anchored prior box. See
[`anchor_offsets`](@ref) for the full argument and the sizes.

# What this expands to

Ordinary system lines. Nothing here is privileged, `Base.show` on the system
prints them, and a model that wants something slightly different should write
them out rather than reach for another keyword:

    plx_interim = plx_A                                  # interim scale
    _frame      = anchored_frame(system_interim, :A, ref_epoch,   # deferred
                     (true, true, true, true, true, true);
                     ra=ra_A, dec=dec_A, plx=plx_A,
                     pmra=pmra_A, pmdec=pmdec_A, rv=rv_A)
                     # ... plus `window=(t1, t2)` when one is given
    ra    = _frame.ra
    dec   = _frame.dec
    plx   = _frame.plx
    pmra  = _frame.pmra
    pmdec = _frame.pmdec
    rv    = _frame.rv

`_frame` mentions `system_interim`, so it is deferred automatically, and so is
everything that reads it. `plx_interim` is *not*: it has to be known before
the bodies are, which is why the interim is built at the anchor's parallax
rather than the barycentre's — a distinction worth ~1 part in 10⁵ (see
[`barycentre_parallax`](@ref)).

A quantity passed `false` is left out of the tuple's `correct` mask and its
*barycentric* value, which your block must define, is passed through the map
untouched. All six are still needed, because the map's second pass builds a
frame out of them; see [`anchored_frame`](@ref).

# Priors and the Jacobian

The map (anchored) → (barycentric) is **triangular with unit diagonal** in ra,
dec, pmra, pmdec and rv: each barycentric quantity is its anchored counterpart
minus a correction built from the *body* variables and `plx_A` only. So a
prior declared on `pmra_A` transfers to `pmra` with volume factor 1, and needs
no `LL +=` correction.

That is exact for the one-pass map (`refine=false`). The default second pass —
which is what makes `pmra_A` mean the anchor's proper motion to better than
10⁻⁶ mas/yr rather than 4 × 10⁻³ — couples each output weakly to all six
inputs, perturbing the diagonal by 10⁻⁹–10⁻⁴ and the determinant by 10⁻⁹. That
is a slowly-varying tilt five orders below anything a posterior resolves, and
much smaller than the interpretation error it buys off; see
[`anchored_frame`](@ref).

Parallax is the one structural exception, and only because it composes as a
reciprocal:
`∂plx/∂plx_A = (plx/plx_A)²`, which differs from 1 by `2⋅dz⋅plx/2.06e8`, with
`dz` the anchor's line-of-sight offset from the barycentre — not the pair's
separation, which is larger by the inverse mass fraction. That is 1e-8 for a
planet-mass companion an AU out and 3e-5 for the 750 au stellar pair above. Add
`LL += 2*log(plx/plx_A)` if you want the declared prior to be exactly a prior
on the barycentric parallax; at these sizes it is a slowly-varying tilt across
the prior's own support, not something a posterior can notice.

A physical prior on a *derived* barycentric quantity is one line, and needs no
Jacobian either way — `~` on an expression is a likelihood term, not a change
of sampling coordinates:

    variables=AnchoredFrame(A; ref_epoch=…, variables=@variables begin
        …
        plx ~ Normal(74.05, 0.15)     # a prior on the *barycentre's* parallax
    end)

# Why not `frame_shift=true`

`frame_shift` reconditions the same degeneracy inside one `G23HObs`, by
subtracting that observation's own `Δpm` from every channel. With two sources
on one frame, each subtracts its *own* `Δpm` while redefining the *same*
`pmra`, so both predict one number where the catalogue has two — for ups And,
13σ and 55σ apart. `AnchoredFrame` gets the conditioning without the
redefinition, so it composes; `G23HObs` errors if `frame_shift=true` and more
than one of them shares the frame.

See also [`anchor_offsets`](@ref), [`System`](@ref).
"""
function AnchoredFrame(anchor;
                       ref_epoch,
                       variables::Tuple=(Priors(), Derived()),
                       ra=true, dec=true, plx=true, pmra=true, pmdec=true, rv=true,
                       refine::Bool=true, window=nothing)
    a = _anchorname(anchor)
    (priors, derived, extra...) = variables
    priors::Priors
    derived::Derived

    given = (; ra, dec, plx, pmra, pmdec, rv)
    anchored = NamedTuple{keys(given)}(
        map(k -> _anchorvar(given[k], k, a), keys(given)))

    declared = Set{Symbol}((keys(priors.priors)..., keys(derived.variables)...))
    for (q, v) in pairs(anchored)
        isnothing(v) && continue
        v in declared || error(
            "AnchoredFrame: `$q` is anchored to the sampled variable `$v`, but the " *
            "`variables=` block defines no `$v`. Declare it (e.g. `$v ~ Normal(…)` " *
            "with the catalogue value and uncertainty for source $a), rename it with " *
            "`$q=:<name>`, or pass `$q=false` to leave $q barycentric and define it " *
            "yourself.")
        q in declared && error(
            "AnchoredFrame: the block defines `$q` directly *and* anchors it to `$v`. " *
            "Those are the two alternatives, not a pair: pass `$q=false` to keep your " *
            "own definition, or drop the `$q =` line and let the anchoring derive it.")
    end
    all(isnothing, anchored) && error(
        "AnchoredFrame: every quantity was passed as `false`, so nothing is anchored " *
        "and this call would do nothing. Use a plain `@variables` block instead.")
    # The map needs a complete frame, anchored or not: its second pass builds
    # one, and the un-anchored quantities are passed through it untouched. So
    # anything not anchored has to be defined in the block — otherwise the
    # generated line fails as an `UndefVarError` inside the sampler, naming a
    # variable the user never wrote.
    for (q, v) in pairs(anchored)
        (isnothing(v) && !(q in declared)) && error(
            "AnchoredFrame: `$q` is not anchored (`$q=false`), so the barycentric " *
            "`$q` has to come from the block — but nothing defines it. The anchored " *
            "map needs all six of ra, dec, plx, pmra, pmdec, rv either way, because " *
            "it refines itself against the frame they make. Define `$q`, or drop " *
            "`$q=false` and anchor it.")
    end

    # Appended in evaluation order; every line reads the ones above it.
    out = Pair{Symbol,Any}[]
    if !isnothing(anchored.plx) && !(INTERIM_PARALLAX_VAR in declared)
        push!(out, INTERIM_PARALLAX_VAR => anchored.plx)
    end
    # One call, shared by all six — and it is the whole map, not a set of
    # offsets, because closing it takes two passes (see `anchored_frame`).
    # Each quantity's input is the anchored variable if it is anchored, and
    # the user's own barycentric one if it is not; `correct` says which.
    #
    # A NamedTuple-valued derived variable is fine here: chain flattening
    # keeps numbers and arrays and skips everything else, so this never
    # becomes a chain column, and it is recomputed from the flat parameter
    # vector on the way back in.
    QS = (:ra, :dec, :plx, :pmra, :pmdec, :rv)
    correct = Expr(:tuple, (!isnothing(anchored[q]) for q in QS)...)
    inputs = Expr(:parameters,
        (Expr(:kw, q, something(anchored[q], q)) for q in QS)...)
    push!(inputs.args, Expr(:kw, :refine, refine))
    # Validate here, at model-build time, rather than letting a bad window
    # surface from inside `arr2nt` on the sampler's first draw.
    isnothing(window) ||
        push!(inputs.args, Expr(:kw, :window, map(float, Tuple(_check_window(window)))))
    push!(out, :_frame => Expr(:call, :anchored_frame, inputs,
        INTERIM_SYSTEM_VAR, QuoteNode(a), float(ref_epoch), correct))
    for q in QS
        isnothing(anchored[q]) || push!(out, q => :(_frame.$q))
    end
    :ref_epoch in declared || push!(out, :ref_epoch => float(ref_epoch))

    vars = OrderedDict{Symbol,Any}(derived.variables)
    for (k, v) in out
        vars[k] = v
    end
    return (priors, Derived(vars, derived.captured_names, derived.captured_vals), extra...)
end

AnchoredFrame(; anchor, kwargs...) = AnchoredFrame(anchor; kwargs...)

export AnchoredFrame

_anchorvar(v::Bool, q::Symbol, a::Symbol) = v ? Symbol(q, :_, a) : nothing
_anchorvar(v::Symbol, ::Symbol, ::Symbol) = v
_anchorvar(::Nothing, ::Symbol, ::Symbol) = nothing
@noinline _anchorvar(@nospecialize(v), q::Symbol, ::Symbol) = error(
    "AnchoredFrame: `$q` takes `true` (anchor, with the default name), `false` " *
    "(leave $q barycentric), or a `Symbol` naming the sampled variable to anchor " *
    "to; got a value of type $(typeof(v)).")
