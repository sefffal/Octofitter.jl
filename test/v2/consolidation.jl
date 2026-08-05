# Regressions for the consolidation pass — the API problems the documentation
# authors found by trying to *use* the ported surface, and the fixes for them.
#
# Two kinds of thing live here:
#
#   * **Error-message contracts.** A name retired by v2, a declaration spec used
#     as a query, a sibling reference in a body block. Each of these failed as an
#     `UndefVarError` or a raw `MethodError` before, which is accurate and
#     useless. The tests assert the *message*, because the message is the fix.
#   * **Silent-wrongness fixes.** `generate_from_params` dropping the `:cor`
#     column and the O'Neil wrapper; `extract_fixed_params` locating a parameter
#     by matching values rather than by layout. None of these threw — they
#     returned something plausible and wrong.
#
# Self-contained: no helper from another test file.

using Test
using Octofitter
using Octofitter: Octofitter, Body, System
using Octofitter.TypedTables: Table
using Distributions
using PlanetOrbits
using LinearAlgebra
using Random
using Statistics

# ---------------------------------------------------------------------------

const CONS_EPOCHS = collect(range(56000.0, 58500.0, length=6))

"""Host + one companion, everything sampled that needs to be."""
function cons_bodies(; hostmass_prior=truncated(Normal(1.0, 0.05), lower=0.2))
    A = Body(name="A", variables=@variables begin
        mass ~ hostmass_prior
        flux_H = 1.0
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass = 2mjup
        a ~ Uniform(1.0, 20.0)
        e ~ Uniform(0.0, 0.6)
        i ~ Sine()
        ω ~ UniformCircular()
        Ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 56000.0
        flux_H = 1e-3
    end)
    return A, b
end

function cons_astrom_table(; cor=nothing)
    n = length(CONS_EPOCHS)
    cols = (; epoch=CONS_EPOCHS,
              ra=collect(range(100.0, 60.0, length=n)),
              dec=collect(range(50.0, -30.0, length=n)),
              σ_ra=fill(5.0, n), σ_dec=fill(7.0, n))
    return isnothing(cor) ? Table(cols) : Table(merge(cols, (; cor=fill(cor, n))))
end

function cons_system(; obs, sysvars=@variables begin plx ~ truncated(Normal(25.0, 0.5), lower=1.0) end)
    A, b = cons_bodies()
    return A, b, System(name="cons", bodies=[A, b], observations=obs, variables=sysvars)
end

# ---------------------------------------------------------------------------
@testset "retired names error, naming the replacement" begin
    # `Planet` and `θ_at_epoch_to_tperi` are the two an old script hits first.
    # Left undefined they are bare `UndefVarError`s; the point of the stub is
    # the message. Every other retired name is simply gone — the name is free
    # again — so there is nothing to assert about it here.
    @test_throws r"`Body`" Planet(name="b")
    @test_throws r"`θ` \(position" θ_at_epoch_to_tperi(0.1, 58849.0; M=1.0, e=0.1, a=5.0)
    @test !isdefined(Octofitter, :PlanetRelAstromObs)
    @test !isdefined(Octofitter, :StarAbsoluteRVObs)
    @test !isdefined(Octofitter, :PlanetRelativeRVObs)
    @test !isdefined(Octofitter, :MarginalizedStarAbsoluteRVObs)
    @test !isdefined(Octofitter, :masspostplot)
    # The retired *aliases* are gone too, so those names are free for reuse.
    @test !isdefined(Octofitter, :PhotometryLikelihood)
    @test !isdefined(Octofitter, :PlanetOrderPrior)
    @test !isdefined(Octofitter, :ObsPriorAstromONeil2019)
    # `rvpostplot`/`rvpostplot_animated` are *not* retired — so without a Makie
    # backend loaded they give the "load a backend" error every plotting entry
    # point gives.
    @test_throws r"requires a Makie backend" rvpostplot(nothing, nothing)
    @test_throws r"requires a Makie backend" rvpostplot_animated(nothing, nothing)
    for f in (dotplot, gaiastarplot, gaiatimeplot, skytrackplot, hipparcosplot)
        @test_throws r"requires a Makie backend" f(nothing, nothing)
    end
    # The HGCA pair predates this file but belongs to the same contract.
    @test_throws r"HGCAObs" HGCAInstantaneousObs(gaia_id=1)
    @test_throws r"G23HObs" GaiaCatalogFitObs()

end

# ---------------------------------------------------------------------------
@testset "declaration specs rejected as query arguments" begin
    A = PlanetOrbits.Body(mass=1.0, flux=(; H=1.0), name=:A)
    b = PlanetOrbits.Body(mass=1e-3, flux=(; H=1e-3), name=:b)
    posys = PlanetOrbits.System((A, b),
        (PlanetOrbits.Orbit(b, about=A; a=3.0, e=0.1, i=0.5, ω=0.2, Ω=0.3, tp=55000.0),);
        plx=25.0)
    sol = PlanetOrbits.orbitsolve(posys, 56000.0)

    # `radvel(sol, :A, Barycentre)` is the obvious thing to write after reading
    # any tutorial, and used to be a MethodError with a page of candidates.
    for f in (radvel, raoff, decoff, posx, projectedseparation)
        @test_throws r"declaration spec, not a resolved reference" f(sol, :A, Barycentre)
        @test_throws r"declaration spec, not a resolved reference" f(sol, Photocentre(:H), :A)
    end
    msg = try; radvel(sol, :A, Barycentre); catch e; sprint(showerror, e); end
    @test occursin("resolveref(posys, Barycentre)", msg)
    @test occursin("PlanetOrbits.barycentre(posys)", msg)

    # …and the spelling the message recommends works.
    @test radvel(sol, :A, Octofitter.resolveref(posys, Barycentre)) isa Real
    @test raoff(sol, Octofitter.resolveref(posys, Photocentre(:H)), PlanetOrbits.BodyRef(1)) isa Real
end

# ---------------------------------------------------------------------------
@testset "`@variables` rejects its two syntactic traps" begin
    # `$` on a `~` line. Parsed from a string, since a literal `$` in this file
    # would be interpolated by the test file's own reader.
    ex = Meta.parse("@variables begin\n  mass ~ LogUniform(0.01, \$host_mass)\nend")
    err = try; @eval(Main, $ex); nothing; catch e; sprint(showerror, e); end
    @test !isnothing(err)
    @test occursin("not supported on a `~` line", err)
    # The message shows the spelling that works.
    @test occursin("mass ~ LogUniform(0.01, host_mass)", err)

    # An unparenthesized `@variables` swallowing a following keyword. Both the
    # in-call spelling (which reaches the macro as one `:(=)` fragment) and the
    # parenthesized-but-still-wrong spelling (two arguments) get the message.
    for src in ("f(variables=@variables begin\n a ~ Uniform(1,2)\nend, method=:AHL21)",
                "@variables(begin\n a ~ Uniform(1,2)\nend, method=:AHL21)")
        ex2 = Meta.parse(src)
        err2 = try; @eval(Main, $ex2); nothing; catch e; sprint(showerror, e); end
        @test !isnothing(err2)
        @test occursin("single `begin ... end` block", err2)
    end

    # `$` on an `=` line is still supported — that asymmetry is the whole point.
    local scale = 3.0
    v = @variables begin
        a ~ Uniform(1.0, 2.0)
        b = a * $scale
    end
    @test v isa Tuple{Octofitter.Priors,Octofitter.Derived}
end

# ---------------------------------------------------------------------------
@testset "a node block cannot read a sibling" begin
    A = Body(name="A", variables=@variables begin mass ~ Uniform(0.5, 2.0) end)
    b = Body(name="b", about=A, variables=@variables begin
        q ~ Uniform(0.0, 1.0)
        mass = q * A.mass          # the natural way to write a mass ratio
        a = 3.0; e = 0.0; i = 0.0; ω = 0.0; Ω = 0.0; tp = 55000.0
    end)
    err = try
        System(name="sib", bodies=[A, b], observations=(),
               variables=@variables begin plx = 25.0 end)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err)
    @test occursin("cannot see a sibling", err)
    @test occursin("about=", err)          # names the fix, not just the rule
end

# ---------------------------------------------------------------------------
@testset "an unnamed `Photocentre` needs exactly one declared band" begin
    # No flux anywhere: this is what every v1 Gaia DR4 model migrated as-is hits,
    # since `target=Photocentre` is that observation's default.
    A = Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Body(name="b", about=A, variables=@variables begin
        mass = 1e-3; a = 3.0; e = 0.0; i = 0.0; ω = 0.0; Ω = 0.0; tp = 55000.0
    end)
    o = RelAstromObs(cons_astrom_table(); target=Photocentre, ref=:A, name="x")
    err = try
        System(name="noflux", bodies=[A, b], observations=[o],
               variables=@variables begin plx = 25.0 end)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err)
    @test occursin("no body declares a flux", err)
    @test occursin("@variables", err)      # in the user's spelling, not PlanetOrbits'

    # Several bands and no band named: also a build-time error rather than a
    # runtime one, and it lists the choices.
    A2 = Body(name="A", variables=@variables begin mass = 1.0; flux_H = 1.0; flux_K = 1.0 end)
    b2 = Body(name="b", about=A2, variables=@variables begin
        mass = 1e-3; a = 3.0; e = 0.0; i = 0.0; ω = 0.0; Ω = 0.0; tp = 55000.0
        flux_H = 1e-3; flux_K = 2e-3
    end)
    o2 = RelAstromObs(cons_astrom_table(); target=Photocentre, ref=:A, name="x")
    err2 = try
        System(name="twoband", bodies=[A2, b2], observations=[o2],
               variables=@variables begin plx = 25.0 end)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err2)
    @test occursin("declares several", err2)

    # One band declared: builds, and the photocentre resolves.
    A3 = Body(name="A", variables=@variables begin mass = 1.0; flux_H = 1.0 end)
    b3 = Body(name="b", about=A3, variables=@variables begin
        mass = 1e-3; a = 3.0; e = 0.0; i = 0.0; ω = 0.0; Ω = 0.0; tp = 55000.0
        flux_H = 1e-3
    end)
    o3 = RelAstromObs(cons_astrom_table(); target=Photocentre, ref=:A, name="x")
    sys3 = System(name="oneband", bodies=[A3, b3], observations=[o3],
                  variables=@variables begin plx ~ truncated(Normal(25.0, 0.5), lower=1.0) end)
    m3 = Octofitter.LogDensityModel(sys3; verbosity=0)
    @test isfinite(m3.ℓπcallback(m3.link(m3.sample_priors(Random.Xoshiro(1)))))
end

# ---------------------------------------------------------------------------
@testset "fixed-parameter lookup is by layout, not by matching values" begin
    A = Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.2)
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass = 0.0
        a ~ Uniform(1.0, 20.0)
        e ~ Uniform(0.0, 0.5)
        i = 0.5; ω = 0.0; Ω = 0.0; tp = 55000.0
    end)
    # A vector-valued prior on the observation. This is the case sentinel
    # matching *cannot* express at all: it searched the flat vector for one
    # value and returned one index, so only the first element of `w` was ever
    # addressable.
    obs = RelAstromObs(cons_astrom_table(); target=b, ref=A, name="GPI cam",
        variables=@variables begin
            w ~ Product([Uniform(0.0, 1.0), Uniform(0.0, 1.0), Uniform(0.0, 1.0)])
        end)
    sys = System(name="slots", bodies=[A, b], observations=[obs],
                 variables=@variables begin plx ~ truncated(Normal(25.0, 0.5), lower=1.0) end)
    model = Octofitter.LogDensityModel(sys; verbosity=0)

    slots = Octofitter._prior_slots(sys)
    @test [s[3] for s in slots] == [:plx, :mass, :a, :e, :w]
    @test [s[1] for s in slots] == [:system, :bodies, :bodies, :bodies, :observations]
    @test [s[4] for s in slots] == [1:1, 2:2, 3:3, 4:4, 5:7]

    # Scalar lookups land on the slot the layout says, not on whichever entry
    # happened to hold the same number.
    _, idx = Octofitter.extract_fixed_params(model, (; bodies=(; b=(; a=7.0))))
    @test idx == [3]
    _, idx2 = Octofitter.extract_fixed_params(model, (; plx=25.0, bodies=(; A=(; mass=1.0))))
    @test sort(idx2) == [1, 2]

    # Every element of the vector prior is addressable, and lands in order.
    vals3, idx3 = Octofitter.extract_fixed_params(model,
        (; observations=(; var"GPI cam"=(; w=[0.1, 0.2, 0.3]))))
    @test idx3 == [5, 6, 7]
    @test vals3 == [0.1, 0.2, 0.3]
    # …under the normalized name too, which is how it appears in a chain.
    @test last(Octofitter.extract_fixed_params(model,
        (; observations=(; GPI_cam=(; w=[0.1, 0.2, 0.3]))))) == [5, 6, 7]

    @test isempty(last(Octofitter.extract_fixed_params(model, NamedTuple())))

    # `planets=` names the v1 spelling and says what replaced it.
    @test_throws r"nested under `bodies`" Octofitter.extract_fixed_params(
        model, (; planets=(; b=(; a=1.0))))
    # A derived variable cannot be fixed, and the error lists what can be.
    @test_throws r"Could not find free parameter" Octofitter.extract_fixed_params(
        model, (; bodies=(; b=(; i=0.5))))
    # A scalar for a vector-valued prior is rejected rather than silently
    # writing one slot.
    @test_throws r"vector-valued prior" Octofitter.extract_fixed_params(
        model, (; observations=(; GPI_cam=(; w=0.5))))

    # Flat names come from the same layout.
    @test Octofitter._flat_param_names(model) ==
        ["plx", "bodies.A.mass", "bodies.b.a", "bodies.b.e",
         "observations.GPI_cam.w[1]", "observations.GPI_cam.w[2]", "observations.GPI_cam.w[3]"]
end

# ---------------------------------------------------------------------------
@testset "generate_from_params keeps the noise model it replicates" begin
    Random.seed!(1234)

    # --- the `:cor` column survives, and the draw uses it -------------------
    tab = cons_astrom_table(cor=0.8)
    A, b, sys = cons_system(obs=[RelAstromObs(tab; target=:b, ref=:A, name="astrom")])
    θ = drawfrompriors(sys)
    sim = generate_from_params(sys, θ; add_noise=true)
    newtab = sim.observations[1].table
    @test hasproperty(newtab, :cor)
    @test all(newtab.cor .== 0.8)

    # The residuals of a correlated draw are themselves correlated. With ρ=0.95
    # over many replicates the sample correlation of (ra−model, dec−model) has
    # to come out near ρ; independent `randn()` per component would give ~0.
    tab95 = cons_astrom_table(cor=0.95)
    A2, b2, sys95 = cons_system(obs=[RelAstromObs(tab95; target=:b, ref=:A, name="astrom")])
    θ95 = drawfrompriors(sys95)
    truth = generate_from_params(sys95, θ95; add_noise=false).observations[1].table
    dr = Float64[]; dd = Float64[]
    for _ in 1:400
        t = generate_from_params(sys95, θ95; add_noise=true).observations[1].table
        append!(dr, (t.ra .- truth.ra) ./ tab95.σ_ra)
        append!(dd, (t.dec .- truth.dec) ./ tab95.σ_dec)
    end
    @test cor(dr, dd) > 0.85

    # …and a table with no `:cor` column does not grow one.
    A3, b3, sysnc = cons_system(obs=[RelAstromObs(cons_astrom_table(); target=:b, ref=:A, name="astrom")])
    θnc = drawfrompriors(sysnc)
    @test !hasproperty(generate_from_params(sysnc, θnc; add_noise=true).observations[1].table, :cor)

    # --- the O'Neil wrapper survives ---------------------------------------
    # Unwrapping returns a *different model*: no Jacobian term, and the
    # likelihood's name changes from `obspri_astrom` to `astrom`, so a chain fit
    # to the replicate does not line up with one fit to the original.
    A4, b4 = cons_bodies()
    inner = RelAstromObs(cons_astrom_table(); target=b4, ref=A4, name="astrom")
    wrapped = ObsPriorONeil2019(inner)
    sysw = System(name="oneil", bodies=[A4, b4], observations=[wrapped],
                  variables=@variables begin plx ~ truncated(Normal(25.0, 0.5), lower=1.0) end)
    θw = drawfrompriors(sysw)
    simw = generate_from_params(sysw, θw; add_noise=false)
    @test simw.observations[1] isa ObsPriorONeil2019
    @test Octofitter.likelihoodname(simw.observations[1]) == Octofitter.likelihoodname(wrapped)
    @test simw.observations[1].orbit === wrapped.orbit
    # The replicate is a model in its own right, with the same parameter shape.
    mw = Octofitter.LogDensityModel(sysw; verbosity=0)
    ms = Octofitter.LogDensityModel(simw; verbosity=0)
    @test mw.D == ms.D
end

# ---------------------------------------------------------------------------
@testset "GOST forecasts with no scans are rejected, not cached" begin
    empty_tbl = Table(; ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_=Float64[])
    # The failure mode this replaces was a `DimensionMismatch` thrown from deep
    # inside the G23H forecast pool, with nothing naming the cache file.
    @test_throws r"no scans for ra=" Octofitter._check_gost_nonempty(
        empty_tbl, 53.2, -9.4, :dr3, "the GOST service")
    @test_throws r"Delete that file" Octofitter._check_gost_nonempty(
        empty_tbl, 53.2, -9.4, :dr3, "cached file /tmp/x.csv"; cached=true)
    good = Table(; ObservationTimeAtBarycentre_BarycentricJulianDateInTCB_=[2.457e6, 2.4571e6])
    @test isnothing(Octofitter._check_gost_nonempty(good, 53.2, -9.4, :dr3, "s"))
end

# ---------------------------------------------------------------------------
@testset "orbitize! export recovers `tp` when the model did not sample it" begin
    # The recommended v2 phase spelling is `θ` + `epoch`, which produces neither
    # a `tp` prior nor a `b_tp` chain column — but `tp` is determined by the
    # elements, so the export derives it instead of refusing.
    A, b, sys = cons_system(obs=[RelAstromObs(cons_astrom_table(); target=:b, ref=:A, name="astrom")])
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    chain = Octofitter.result2mcmcchain(
        [model.arr2nt(model.sample_priors(Random.Xoshiro(k))) for k in 1:5])
    @test !haskey(chain, :b_tp)

    tp = Octofitter._tp_from_chain(model, chain, :b)
    @test length(tp) == 5
    @test all(isfinite, tp)
    # Cross-check against the reconstructed systems directly.
    @test tp ≈ [PlanetOrbits.periastron(s, 1) for s in construct_system(model, chain, :)]

    mktempdir() do dir
        f = joinpath(dir, "out.hdf5")
        @test Octofitter.savehdf5(f, model, chain) == f
        @test isfile(f)
    end
end
