# ---------------------------------------------------
# Correction flags: `:auto | :on | :off`, decided once at build
#
# PlanetOrbits' two precision opt-outs (`observing_geometry`,
# `barycentric_lighttime`) are assertions about the *data*, which is exactly
# the thing Octofitter can see and PlanetOrbits cannot. So here they are
# tri-state and default to `:auto`, and `:auto` is decided by measurement
# rather than by a rule of thumb.
#
# Why not a closed-form rule. The ρ-scaling formulas on PlanetOrbits'
# "Precision opt-outs" page were derived for star-plus-planet topologies. v2's
# whole point is hierarchies — moons, 2+2 quadruples, N-body — where "the" ρ is
# not well defined, and the fourth correction (line-of-sight projection) does
# not scale with ρ at all. What *is* always well defined is the manual's own
# "Checking the difference" recipe: solve both ways and compare against your
# own uncertainties. `:auto` mechanizes that.
#
# Why once, at build, and never per draw. Toggling a correction inside the
# likelihood would make the log density and its gradient discontinuous across
# parameter space, which HMC pays for in divergences; and the flags are
# interpolated into the generated code as literals (see `codegen.jl`), which a
# per-draw decision would defeat. The *resolved booleans* are what the built
# `System` stores, so everything downstream — codegen, cross-validation,
# plotting — keeps reading a plain `Bool` and needs no changes.
#
# Why the decision is recorded rather than merely made. A posterior is not
# separable from the corrections that produced it. The verdict, the measured
# impacts, the margin, the draw count and the seed all travel on the built
# model and into the chain's metadata, so two models under comparison show
# their (possibly different) verdicts side by side.
# ---------------------------------------------------

"""
    Octofitter.ObsImpact

What one observation had to say about one correction flag.

`impact` is the per-point figure — the 99th percentile over draws of
`max_epochs |Δ prediction| / σ`. `bias` is `impact · √n`, which is the
quantity the decision is actually made on: a *common-mode* per-point shift of
`b·σ` across `n` points moves a posterior mean by about `b·√n·σ`, so the same
per-point impact matters far more to a 10⁵-transit dataset than to a
20-epoch one. `NaN` means the observation could not say.
"""
struct ObsImpact
    name::String
    impact::Float64
    n::Int
    bias::Float64
end
ObsImpact(name, impact, n) = ObsImpact(name, impact, n, impact * sqrt(n))

"""
    Octofitter.CorrectionDecision

What `:auto` decided for one correction flag, and the evidence for it.

Fields: `flag`, `resolved` (the `Bool` actually used), `source` (`:user` when
given explicitly, `:auto` when measured), `impacts` (a worst-first vector of
[`ObsImpact`](@ref Octofitter.ObsImpact)), `ndraws`/`nfailed`, `threshold`
(the accumulated-bias limit, in σ of the posterior mean), `seed`, and a
human-readable `note`.

An impact of `NaN` means that observation type does not report predictions the
test can compare, which is read as "keep the correction on".
"""
struct CorrectionDecision
    flag::Symbol
    resolved::Bool
    source::Symbol
    impacts::Vector{ObsImpact}
    ndraws::Int
    nfailed::Int
    threshold::Float64
    seed::UInt64
    note::String
end

"""
    Octofitter.CorrectionReport

The [`CorrectionDecision`](@ref)s for a built model, plus informational
advisories that carry no decision — the size of the `radvel` Einstein term and
of the secular-acceleration drift against each radial-velocity series' own
uncertainty. Stored on the `System` and propagated into chain metadata.
"""
struct CorrectionReport
    decisions::Vector{CorrectionDecision}
    advisories::Vector{String}
end
CorrectionReport() = CorrectionReport(CorrectionDecision[], String[])

const CORRECTION_FLAGS = (:observing_geometry, :barycentric_lighttime)

# The decision criterion, in σ of the *posterior mean*.
#
# A flat per-point margin ("the correction must be 10³ below your error bar")
# is the wrong shape, because these corrections are common-mode: they push
# every epoch the same way, so their effect on a fitted parameter accumulates
# as √n rather than averaging down. A 0.001σ per-point bias is invisible in 10
# points and a 0.3σ shift in 10⁵ Gaia transits. So the test is on the
# accumulated bias, `impact · √n`, against one threshold — which makes the
# effective per-point margin ~100 for a 100-point dataset and ~300 for a
# DR4-scale one, automatically, instead of by a number somebody picked.
#
# 0.1σ is the "you would not notice this in a corner plot" line.
const CORRECTION_BIAS_THRESHOLD = 0.1
const CORRECTION_DRAWS = 300
const CORRECTION_MIN_SUCCESS = 100
const CORRECTION_SEED = 0x0c70f177e5000001

# The boolean a report resolved for one flag. Every path through
# `_impact_report` emits a decision for both flags, so the `true` fallback is
# a belt-and-braces "accurate by default", not a live branch.
function _resolved_flag(r::CorrectionReport, flag::Symbol)
    for d in r.decisions
        d.flag === flag && return d.resolved
    end
    return true
end

# One line, ASCII, and **under 68 characters** — this is what `savechain`
# hands to a FITS header card, and a longer string value is silently dropped
# rather than truncated. The two flags at their longest spelling come to 63.
# (The resolved booleans are also written to the header under their own names,
# so what this string adds is *how* they were decided.)
function Base.show(io::IO, r::CorrectionReport)
    isempty(r.decisions) && return print(io, "CorrectionReport()")
    join(io, ("$(d.flag)=$(d.resolved)($(d.source))" for d in r.decisions), " ")
end

function Base.show(io::IO, ::MIME"text/plain", r::CorrectionReport)
    println(io, "CorrectionReport:")
    for d in r.decisions
        println(io, "  ", _correction_logline(d))
    end
    for a in r.advisories
        println(io, "  · ", a)
    end
end

# ---------------------------------------------------
# Flag normalization
# ---------------------------------------------------

_normalize_flag(v::Symbol, name) =
    v in (:auto, :on, :off) ? v : _err_flag(v, name)
_normalize_flag(v::Bool, ::Any) = v ? :on : :off
_normalize_flag(v, name) = _err_flag(v, name)

@noinline _err_flag(v, name) = error(
    "`$name` takes `:auto`, `:on` or `:off` (and `true`/`false` as aliases for " *
    "the latter two); got $(repr(v)). `:auto` — the default — measures whether " *
    "the correction moves this model's predictions by anything comparable to " *
    "your uncertainties, and decides once, at build.")

# ---------------------------------------------------
# What an observation can say about a correction
#
# One generic function, with a conservative default. An observation type that
# does not implement it forces every correction on, which is the setting it
# had before `:auto` existed — `:auto` can only ever make a model cheaper than
# the old default, never less accurate.
# ---------------------------------------------------

"""
    correction_impact(obs, ctx_a, ctx_b) -> (; delta, sigma, n)

What `obs` has to say about a correction: `delta`, the largest change in its
own model predictions between two solves of the *same* parameters; `sigma`,
the tightest per-point uncertainty that change should be judged against; and
`n`, how many numbers were compared.

Implement this for an observation type to let `observing_geometry=:auto` and
`barycentric_lighttime=:auto` decide whether a correction is worth its cost
for data like yours, and declare [`has_correction_impact`](@ref) `true` for
it. Reuse `simulate` rather than re-deriving which observables matter:
whatever the likelihood actually fits is what should be compared.

**`sigma` comes from you, not from the caller.** An earlier version had the
resolver guess it from a list of plausible σ column names, which meant a
renamed column silently became "cannot say" and quietly turned the correction
back on — a failure mode invisible to everybody. Read it from your own table.

The default returns `delta = sigma = NaN`, meaning "cannot say", which
resolves the flag **on**.
"""
correction_impact(::AbstractObs, ::ObsContext, ::ObsContext) =
    (; delta=NaN, sigma=NaN, n=0)

"""
    has_correction_impact(::Type{<:AbstractObs}) -> Bool

Whether this observation type implements [`correction_impact`](@ref).

Consulted *before* any draws are taken: a model in which nothing can answer
resolves both flags `:on` immediately, rather than solving three hundred prior
draws to reach a conclusion that was fixed before it started. Declare it
beside your `correction_impact` method.
"""
has_correction_impact(::Type{<:AbstractObs}) = false
has_correction_impact(obs::AbstractObs) = has_correction_impact(typeof(obs))

"""
    correction_advisories(obs, ctx) -> Tuple of (label, value, σ)

Magnitudes worth reporting to the user about `obs` at one draw, with no
decision attached — e.g. how large the `radvel` Einstein term is against this
series' uncertainty. The `:auto` machinery medians them over draws and prints
them; a model with explicit flags computes none of this.
"""
correction_advisories(::AbstractObs, ::ObsContext) = ()

# Whether every observation that was asked returned "cannot say". `all` over an
# empty collection is `true`, which is the wrong answer here, so this tracks
# whether anything was asked at all.
function _all_cannot_say(ratios)
    asked = false
    for d in values(ratios), v in values(d)
        asked = true
        all(isnan, v) || return false
    end
    return asked
end

# Largest absolute difference across the named prediction vectors of two
# `simulate` results, and how many numbers were compared. The `σ` comes from
# the caller — see `correction_impact`.
function _simulate_impact(a, b, keys, σ)
    m = 0.0
    n = 0
    for k in keys
        va = getproperty(a, k)
        vb = getproperty(b, k)
        for i in eachindex(va)
            d = abs(float(va[i] - vb[i]))
            isfinite(d) || return (; delta=NaN, sigma=float(σ), n=0)
            m = max(m, d)
            n += 1
        end
    end
    return (; delta=m, sigma=float(σ), n)
end

# Tightest finite, positive value of one uncertainty column. `NaN`, not `Inf`,
# when there is none: "I do not know" and "infinitely tolerant" are opposite
# answers, and only the first is honest.
function _tightest(vals...)
    best = NaN
    for v in vals, x in v
        y = float(x)
        (isfinite(y) && y > 0) || continue
        best = isnan(best) ? y : min(best, y)
    end
    return best
end

# ---------------------------------------------------
# The prior-predictive impact test
# ---------------------------------------------------

# Three solves per draw: the accurate baseline, and each candidate flag turned
# off in turn. ~10 µs each — negligible next to anything else that happens at
# build.
#
# `sys` is an `Octofitter.System`, deliberately not annotated: this file is
# included before `nodes.jl` so that the `System` struct can carry a
# `CorrectionReport` field.
function _impact_report(sys, og::Symbol, blt::Symbol;
                        ndraws::Int, minsuccess::Int, threshold::Float64,
                        seed::UInt64, verbosity::Int)
    flags = Symbol[]
    og === :auto && push!(flags, :observing_geometry)
    blt === :auto && push!(flags, :barycentric_lighttime)

    decisions = CorrectionDecision[]
    advisories = String[]

    # Explicit settings are recorded, not measured.
    for (flag, spec) in ((:observing_geometry, og), (:barycentric_lighttime, blt))
        spec === :auto && continue
        push!(decisions, CorrectionDecision(flag, spec === :on, :user,
            ObsImpact[], 0, 0, NaN, seed, "set explicitly"))
    end
    isempty(flags) && return CorrectionReport(decisions, advisories)

    _resolve_on(note) = CorrectionReport(
        vcat(decisions, [CorrectionDecision(f, true, :auto, ObsImpact[], 0, 0,
                                            threshold, seed, note) for f in flags]),
        advisories)

    dataobs = [o for o in sys.observations if !_isprior(o)]
    # Nothing to test against, so accuracy wins: a simulation-only model is
    # usually being used to *generate* the truth.
    isempty(dataobs) && return _resolve_on("no observations to test against")

    # Static capability check, before a single solve. An observation type that
    # implements no `correction_impact` can only ever answer "cannot say",
    # which resolves the flag on — so if nothing in the model can answer, the
    # conclusion is fixed before the test starts and running it is pure waste.
    # This is the common case for the most expensive models there are (G23H,
    # images, interferometry), which is exactly where the waste would hurt.
    if !any(has_correction_impact, dataobs)
        kinds = join(unique(string(nameof(typeof(o))) for o in dataobs), ", ")
        return _resolve_on("no observation can report comparable predictions " *
                           "($kinds); nothing to measure, so keeping the correction on")
    end

    unique_ep, maps = epoch_plan(sys)
    arr2nt = make_arr2nt(sys)
    sampler = make_prior_sampler(sys)
    build = @RuntimeGeneratedFunction(:(function (θ)
        T = $_system_number_type(θ)
        $(_build_system_expr(sys)...)
    end))
    rng = Random.Xoshiro(seed)

    # ratios[flag][obs name] = per-point |Δ|/σ over draws; npoints is fixed by
    # the data, so it is recorded once.
    ratios = Dict(f => Dict{String,Vector{Float64}}() for f in flags)
    npoints = Dict{String,Int}()
    advice = Dict{Tuple{String,String},Vector{Float64}}()
    advice_σ = Dict{Tuple{String,String},Float64}()
    nfailed = 0
    nok = 0

    for _ in 1:ndraws
        local θ, posys
        try
            θ = arr2nt(collect(sampler(rng)))
            posys = build(θ)
        catch err
            err isa InterruptException && rethrow()
            nfailed += 1
            continue
        end
        try
            base = PlanetOrbits.orbitsolve(posys, unique_ep; method=sys.method,
                observing_geometry=true, barycentric_lighttime=true)
            alt = Dict{Symbol,Any}()
            for f in flags
                alt[f] = PlanetOrbits.orbitsolve(posys, unique_ep; method=sys.method,
                    observing_geometry=(f !== :observing_geometry),
                    barycentric_lighttime=(f !== :barycentric_lighttime))
            end
            for o in dataobs
                nm = likelihoodname(o)
                key = normalizename(nm)
                θ_obs = hasproperty(θ.observations, key) ?
                        getproperty(θ.observations, key) : (;)
                ctx_base = _obsctx(sys, θ, θ_obs, posys, base, maps[o])
                for f in flags
                    ctx_alt = _obsctx(sys, θ, θ_obs, posys, alt[f], maps[o])
                    r = correction_impact(o, ctx_base, ctx_alt)
                    push!(get!(Vector{Float64}, ratios[f], nm),
                          float(r.delta) / float(r.sigma))
                    npoints[nm] = max(get(npoints, nm, 0), r.n)
                end
                for (label, value, σ) in correction_advisories(o, ctx_base)
                    push!(get!(Vector{Float64}, advice, (nm, label)), float(value))
                    advice_σ[(nm, label)] = float(σ)
                end
            end
            nok += 1
        catch err
            err isa InterruptException && rethrow()
            nfailed += 1
        end
        # Runtime backstop behind the static check: a type that *claims* to
        # answer but returns NaN anyway (a σ column of zeros, say) would
        # otherwise buy three hundred draws of the same non-answer.
        nok == 1 && _all_cannot_say(ratios) && break
    end

    # A third of the requested draws, capped at `minsuccess`. Scaling with
    # `ndraws` matters: a caller who asks for 40 draws has already accepted a
    # coarser measurement, and should not silently get the fallback instead.
    need = max(1, min(minsuccess, cld(ndraws, 3)))
    if nok < need && !_all_cannot_say(ratios)
        verbosity >= 1 && @warn """
        Correction flags left on: too few prior draws could be solved ($nok of $ndraws).
        Usually this means the priors admit parameter combinations the orbit
        constructor rejects. That costs accuracy, never correctness — but it is
        worth knowing about.
        """
        return _resolve_on("only $nok of $ndraws prior draws could be solved " *
                           "(needed $need); keeping the correction on")
    end

    for flag in flags
        per_obs = ObsImpact[]
        for (nm, vals) in ratios[flag]
            push!(per_obs, ObsImpact(nm, _p99(vals), get(npoints, nm, 0)))
        end
        # Worst-first, with "cannot say" sorting to the top: it is the strongest
        # statement in the set, not the weakest.
        sort!(per_obs; by=x -> isnan(x.bias) ? Inf : x.bias, rev=true)
        worst = isempty(per_obs) ? ObsImpact("(none)", NaN, 0) : first(per_obs)
        # OR across observations: one observation that needs the correction
        # needs it for the whole model.
        needed = isempty(per_obs) || isnan(worst.bias) || worst.bias >= threshold
        note = if isnan(worst.bias)
            "\"$(worst.name)\" does not report comparable predictions"
        elseif needed
            "would bias \"$(worst.name)\" by $(_g(worst.bias))σ " *
            "($(_g(worst.impact))σ per point × √$(worst.n))"
        elseif iszero(worst.bias)
            "changes no prediction at all"
        else
            "worst accumulated bias $(_g(worst.bias))σ, " *
            "$(_g(threshold / worst.bias))× inside the $(_g(threshold))σ limit"
        end
        push!(decisions, CorrectionDecision(flag, needed, :auto, per_obs,
            nok, nfailed, threshold, seed, note))
    end

    for ((nm, label), vals) in advice
        σ = get(advice_σ, (nm, label), NaN)
        med = _median(vals)
        push!(advisories, "\"$nm\": $label ≈ $(_g(med))" *
            (isfinite(σ) && σ > 0 ? " ($(_g(med / σ))× its tightest σ = $(_g(σ)))" : ""))
    end
    sort!(advisories)

    return CorrectionReport(decisions, advisories)
end

# 99th percentile, NaN-preserving: a single "cannot say" draw must not be
# averaged away.
function _p99(v::Vector{Float64})
    isempty(v) && return NaN
    any(isnan, v) && return NaN
    s = sort(v)
    return s[clamp(ceil(Int, 0.99 * length(s)), 1, length(s))]
end

function _median(v::Vector{Float64})
    isempty(v) && return NaN
    s = sort(filter(isfinite, v))
    isempty(s) && return NaN
    n = length(s)
    return isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

_g(x) = isfinite(x) ? string(round(x; sigdigits=3)) : string(x)

function _correction_logline(d::CorrectionDecision)
    tag = d.source === :user ? "user" : "auto"
    base = "$(d.flag) = $(d.resolved) ($tag)"
    d.source === :user && return base
    draws = d.ndraws == 0 ? "no draws needed" :
            "$(d.ndraws) prior draws" * (d.nfailed > 0 ? " ($(d.nfailed) rejected)" : "")
    return base * ": " * d.note * " — over " * draws *
           " (seed 0x$(string(d.seed, base=16)))"
end

function _log_corrections(sysname, r::CorrectionReport, verbosity::Int)
    verbosity >= 1 || return
    for d in r.decisions
        d.source === :user && continue
        @info "[$sysname] " * _correction_logline(d)
    end
    for a in r.advisories
        @info "[$sysname] " * a
    end
end

# ---------------------------------------------------
# Post-sampling re-check
#
# The impact test is one function, runnable at up to three stages: prior draws
# at build (today's behaviour, and the only one that *decides*), variational
# draws at `initialize!` (not implemented), and posterior draws after sampling
# (below, which only *reports*).
#
# If the prior-based decision proves noisy in practice — the expected failure
# mode is wide priors resolving a flag *on* that the posterior never needs,
# i.e. wasted compute rather than a wrong answer — the targeted remedy is a
# re-decision at `initialize!` from its mixture-of-Gaussians variational
# approximation. That stage has a safety asymmetry worth writing down: because
# build resolved conservatively, an `initialize!` re-decision could only ever
# **de-escalate**, which is safe by construction (nothing has sampled yet with
# a correction missing), and it would trigger a model rebuild — which the
# codegen literals force anyway. It would keep the same 10³ margin: an MoG
# approximation under-covers multimodal and heavy-tailed posteriors, so the
# variational stage earns no friendlier threshold. Ship prior-based first.
# ---------------------------------------------------

"""
    recheck_corrections(model, chain; ndraws=300, seed=…, verbosity=1)

Re-run the build-time impact test on draws from the sampled posterior, and
report whether the corrections the model was built with are still the right
ones. Returns a [`CorrectionReport`](@ref).

Two directions, and they mean different things:

- **Escalation (warning).** A flag resolved *off* under the priors would
  resolve *on* under the posterior — the posterior concentrated somewhere the
  correction matters, so the results may be biased and the fit should be
  re-run with the flag on. Expected to be rare.
- **De-escalation (hint).** A flag resolved *on* — usually because broad
  priors covered parameter space the posterior never visits — is comfortably
  unneeded under the posterior. Nothing is wrong; the next fit of the same
  data can set it `:off` and sample faster. This is the common case.

Draw from the **target** chain only. Tempered chains sit closer to the prior
and would re-trigger exactly the broad-prior inclusions the posterior has
ruled out.
"""
function recheck_corrections(model, chain; ndraws::Int=CORRECTION_DRAWS,
                             seed::UInt64=CORRECTION_SEED, verbosity::Int=1,
                             threshold::Float64=CORRECTION_BIAS_THRESHOLD)
    sys = model.system
    dataobs = [o for o in sys.observations if !_isprior(o)]
    isempty(dataobs) && return CorrectionReport()
    any(has_correction_impact, dataobs) || return CorrectionReport()
    unique_ep, maps = epoch_plan(sys)
    build = @RuntimeGeneratedFunction(:(function (θ)
        T = $_system_number_type(θ)
        $(_build_system_expr(sys)...)
    end))

    n = size(chain, 1) * size(chain, 3)
    rng = Random.Xoshiro(seed)
    idx = n <= ndraws ? collect(1:n) : rand(rng, 1:n, ndraws)
    θs = mcmcchain2result(model, chain, idx)

    ratios = Dict(f => Dict{String,Vector{Float64}}() for f in CORRECTION_FLAGS)
    npoints = Dict{String,Int}()
    nok = 0
    nfailed = 0
    for θ in θs
        local posys
        try
            posys = build(θ)
        catch err
            err isa InterruptException && rethrow()
            nfailed += 1
            continue
        end
        try
            base = PlanetOrbits.orbitsolve(posys, unique_ep; method=sys.method,
                observing_geometry=true, barycentric_lighttime=true)
            for f in CORRECTION_FLAGS
                alt = PlanetOrbits.orbitsolve(posys, unique_ep; method=sys.method,
                    observing_geometry=(f !== :observing_geometry),
                    barycentric_lighttime=(f !== :barycentric_lighttime))
                for o in dataobs
                    nm = likelihoodname(o)
                    key = normalizename(nm)
                    θ_obs = hasproperty(θ.observations, key) ?
                            getproperty(θ.observations, key) : (;)
                    r = correction_impact(o,
                        _obsctx(sys, θ, θ_obs, posys, base, maps[o]),
                        _obsctx(sys, θ, θ_obs, posys, alt, maps[o]))
                    push!(get!(Vector{Float64}, ratios[f], nm),
                          float(r.delta) / float(r.sigma))
                    npoints[nm] = max(get(npoints, nm, 0), r.n)
                end
            end
            nok += 1
        catch err
            err isa InterruptException && rethrow()
            nfailed += 1
        end
    end

    decisions = CorrectionDecision[]
    for (f, built) in ((:observing_geometry, sys.observing_geometry),
                       (:barycentric_lighttime, sys.barycentric_lighttime))
        per_obs = [ObsImpact(nm, _p99(v), get(npoints, nm, 0)) for (nm, v) in ratios[f]]
        sort!(per_obs; by=x -> isnan(x.bias) ? Inf : x.bias, rev=true)
        worst = isempty(per_obs) ? ObsImpact("(none)", NaN, 0) : first(per_obs)
        # The identical criterion the build-time test used, so the two verdicts
        # are comparable rather than merely similar.
        needed = isempty(per_obs) || isnan(worst.bias) || worst.bias >= threshold
        note = needed == built ? "unchanged under the posterior" :
               needed ? "**escalation**: needed under the posterior but the fit ran without it" :
               "de-escalation: not needed under the posterior"
        push!(decisions, CorrectionDecision(f, needed, :posterior, per_obs,
            nok, nfailed, threshold, seed, note))
        verbosity >= 1 || continue
        if needed && !built
            @warn """
            $f resolved `off` under the priors but is needed under the posterior.
            The posterior concentrated where this correction matters — it would bias
            "$(worst.name)" by $(_g(worst.bias))σ. These results may be biased; re-run
            with `$f=:on`.
            """
        elseif !needed && built
            @info """
            $f was on and changed nothing the posterior can see (worst accumulated
            bias $(_g(worst.bias))σ, against a $(_g(threshold))σ limit). Setting
            `$f=:off` will sample faster next time.
            """
        end
    end
    return CorrectionReport(decisions, String[])
end

export recheck_corrections
