# The inference-analysis machinery over a fitted model: prior-only models,
# data subsetting for cross-validation, pointwise likelihoods, forward
# simulation, and injection-recovery bookkeeping.
#
# Everything here is deliberately tiny — two bodies, a handful of epochs, no
# sampling unless the integration gate is on. The point is the plumbing, not
# the statistics.

const AM_EPOCHS_1 = [56000.0, 56500.0, 57000.0, 57500.0]
const AM_EPOCHS_2 = [58000.0, 58500.0]

# A fixed truth to generate exact "data" from, so simulate-then-fit round
# trips are checkable to machine precision.
function am_truth(; plx=25.0)
    A = PlanetOrbits.Body(mass=1.0 - 5mjup, name=:A)
    b = PlanetOrbits.Body(mass=5mjup, name=:b)
    return PlanetOrbits.System((A, b),
        (PlanetOrbits.Orbit(b, about=A; a=3.0, e=0.2, i=0.9, ω=1.2, Ω=2.0, tp=55000.0),);
        plx)
end

function am_astrom(epochs, name; σ=1.0,
                   variables=(Octofitter.Priors(), Octofitter.Derived()))
    tr = orbitsolve(am_truth(), epochs)
    ra = [raoff(tr[k], :b, :A) for k in eachindex(epochs)]
    dec = [decoff(tr[k], :b, :A) for k in eachindex(epochs)]
    return RelAstromObs((epoch=epochs, ra=ra, dec=dec,
            σ_ra=fill(σ, length(epochs)), σ_dec=fill(σ, length(epochs)));
        target=:b, ref=:A, name, variables)
end

# The model under test: fixed elements except `a`, so it has exactly one free
# body parameter plus one per-observation jitter.
function am_system(; obs=nothing)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0 - 5mjup
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 5mjup
        a ~ Uniform(1, 8)
        e = 0.2
        i = 0.9
        ω = 1.2
        Ω = 2.0
        tp = 55000.0
    end)
    observations = isnothing(obs) ? (
        am_astrom(AM_EPOCHS_1, "astrom1"),
        am_astrom(AM_EPOCHS_2, "astrom2";
            variables=@variables begin jitter ~ LogUniform(0.1, 10.0) end),
    ) : obs
    return Octofitter.System(name="am", bodies=[A, b], observations=observations,
        variables=@variables begin plx = 25.0 end)
end

# A `UniformCircular` variant: the case where the model carries a
# `@variables`-emitted prior term (a `UnitLengthPrior`) alongside its data.
function am_system_circular()
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 - 5mjup end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 5mjup
        a = 3.0
        e = 0.2
        i = 0.9
        ω ~ UniformCircular()
        Ω = 2.0
        tp = 55000.0
    end)
    return Octofitter.System(name="amu", bodies=[A, b],
        observations=[am_astrom(AM_EPOCHS_1, "astrom1")],
        variables=@variables begin plx = 25.0 end)
end

# An observation type that carries data but implements neither subsetting nor
# forward simulation — the shape of every not-yet-ported likelihood, and of
# `MarginalizedRVObs`, which structurally cannot implement the first.
struct AMUnsubsettable{T} <: Octofitter.AbstractObs
    table::T
    name::String
end
Octofitter.ln_like(::AMUnsubsettable, ctx::Octofitter.ObsContext) =
    zero(Octofitter._system_number_type(ctx.θ_system))

# ---------------------------------------------------

@testset "prior_only_model" begin
    sys = am_system()
    prior = Octofitter.prior_only_model(sys)

    # Same parameter vector: the blanked observations keep their variables, so
    # a prior-only chain lines up column for column with a full one.
    @test Octofitter.list_parameter_names(prior) == Octofitter.list_parameter_names(sys)
    @test length(Octofitter._list_priors(prior)) == length(Octofitter._list_priors(sys))
    @test all(o -> o isa Octofitter.BlankLikelihood, prior.observations)
    @test [Octofitter.likelihoodname(o) for o in prior.observations] ==
          [Octofitter.likelihoodname(o) for o in sys.observations]

    # And the likelihood really is gone.
    θ = Octofitter.make_arr2nt(prior)(Octofitter.sample_priors(Random.Xoshiro(3), prior))
    @test Octofitter.make_ln_like(prior)(prior, θ) == 0.0

    # `exclude_all` additionally drops the prior-shaped terms — the point of
    # the flag, since those are exactly what a `log_Z0` reference must not
    # contain — while leaving the free variables, and hence the parameter
    # space, untouched.
    sysu = am_system_circular()
    @test length(sysu.priorterms) == 1
    @test length(Octofitter.prior_only_model(sysu).priorterms) == 1
    bare = Octofitter.prior_only_model(sysu; exclude_all=true)
    @test isempty(bare.priorterms)
    @test Octofitter.list_parameter_names(bare) == Octofitter.list_parameter_names(sysu)
end

@testset "likelihood-object folds" begin
    sys = am_system()
    @test Octofitter._count_likeobj(sys) == 2

    kfold = Octofitter.generate_kfold_systems(sys)
    @test length(kfold) == 2
    @test [Octofitter.likelihoodname(o) for o in kfold[1].observations] == ["astrom2"]
    @test [Octofitter.likelihoodname(o) for o in kfold[2].observations] == ["astrom1"]

    per = Octofitter.generate_systems_per_like(sys)
    @test length(per) == 2
    @test [Octofitter.likelihoodname(o) for o in per[1].observations] == ["astrom1"]
    @test [Octofitter.likelihoodname(o) for o in per[2].observations] == ["astrom2"]

    filt = Octofitter.generate_system_filtered_like(
        o -> Octofitter.likelihoodname(o) == "astrom2", sys)
    @test [Octofitter.likelihoodname(o) for o in filt.observations] == ["astrom2"]

    # Prior-shaped terms carry no data, so a fold never drops one.
    sysu = am_system_circular()
    @test Octofitter._count_likeobj(sysu) == 1
    @test length(Octofitter.generate_kfold_systems(sysu)) == 1
    @test length(Octofitter.generate_kfold_systems(sysu)[1].priorterms) == 1
end

@testset "epoch groups keep the rows they name" begin
    sys = am_system()
    total = length(AM_EPOCHS_1) + length(AM_EPOCHS_2)

    # This is the v1 inversion regression: v1 passed the *complement* of the
    # wanted rows to `likeobj_from_epoch_subset`, which keeps (not drops) what
    # it is given, so asking for row `i` produced every row except `i`.
    systems, eps = Octofitter.generate_system_per_epoch(sys)
    @test length(systems) == total
    @test eps == vcat(AM_EPOCHS_1, AM_EPOCHS_2)
    for (k, s) in enumerate(systems)
        @test Octofitter._count_epochs(s) == 1
        @test only(Octofitter.epoch_plan(s)[1]) == eps[k]
    end
    # Row 1 lives in astrom1, the last row in astrom2 — and each system carries
    # only the observation that owns its row.
    @test [Octofitter.likelihoodname(o) for o in systems[1].observations] == ["astrom1"]
    @test [Octofitter.likelihoodname(o) for o in systems[end].observations] == ["astrom2"]

    cum, cum_eps = Octofitter.generate_cumulative_system_per_epoch(sys)
    @test length(cum) == total
    @test [Octofitter._count_epochs(s) for s in cum] == collect(1:total)
    @test cum_eps[end] == vcat(AM_EPOCHS_1, AM_EPOCHS_2)

    # Grouping across observation boundaries picks rows from both.
    grouped, gep = Octofitter.generate_systems_with_epoch_groups(
        sys, [[1, 5]], i -> "_g$i")
    @test length(grouped) == 1
    @test [Octofitter.likelihoodname(o) for o in grouped[1].observations] ==
          ["astrom1", "astrom2"]
    @test gep[1] == [AM_EPOCHS_1[1], AM_EPOCHS_2[1]]
end

@testset "pointwise_like" begin
    sys = am_system()
    model = Octofitter.LogDensityModel(sys; verbosity=0)

    # A two-sample "chain", built by hand so no sampler is needed.
    θs = [model.arr2nt(Octofitter.sample_priors(Random.Xoshiro(s), sys)) for s in (1, 2)]
    chain = Octofitter.result2mcmcchain(θs)

    LL, eps = Octofitter.pointwise_like(model, chain; verbosity=0)
    @test size(LL) == (2, length(AM_EPOCHS_1) + length(AM_EPOCHS_2))
    @test eps == vcat(AM_EPOCHS_1, AM_EPOCHS_2)
    @test all(isfinite, LL)

    # The defining property PSIS-LOO needs: the columns are a decomposition of
    # the model's own log-likelihood, not overlapping pieces of it.
    lnlike = Octofitter.make_ln_like(sys)
    for (i, θ) in enumerate(θs)
        @test sum(LL[i, :]) ≈ lnlike(sys, θ) rtol = 1e-10
    end

    # And it agrees with the (much slower) build-a-system-per-epoch route.
    per_epoch, _ = Octofitter.generate_system_per_epoch(sys)
    for (j, s) in enumerate(per_epoch)
        f = Octofitter.make_ln_like(s)
        for (i, θ) in enumerate(θs)
            @test LL[i, j] ≈ f(s, θ) rtol = 1e-10
        end
    end
end

@testset "pointwise_like excludes prior-shaped terms" begin
    sysu = am_system_circular()
    model = Octofitter.LogDensityModel(sysu; verbosity=0)
    θs = [model.arr2nt(Octofitter.sample_priors(Random.Xoshiro(s), sysu)) for s in (1, 2)]
    chain = Octofitter.result2mcmcchain(θs)

    LL, _ = Octofitter.pointwise_like(model, chain; verbosity=0)
    # One column per data row and no more: the `UnitLengthPrior` behind
    # `UniformCircular` is not a data point. v1 gave it a column of its own and
    # added it to every other column as well.
    @test size(LL, 2) == length(AM_EPOCHS_1)

    lnlike = Octofitter.make_ln_like(sysu)
    for (i, θ) in enumerate(θs)
        ulp = Octofitter.ln_like(sysu.priorterms[1][2],
            Octofitter.ObsContext(θ, θ.bodies.b, nothing, nothing, Int[]))
        @test sum(LL[i, :]) + ulp ≈ lnlike(sysu, θ) rtol = 1e-10
    end
end

@testset "subsetting failures name the observation" begin
    o = AMUnsubsettable(Table(epoch=[57000.0], val=[1.0]), "stubborn")
    err = try
        Octofitter._subset(o, [1]; what="cross-validation")
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("stubborn", err.msg)
    @test occursin("cross-validation", err.msg)
    @test occursin("AMUnsubsettable", err.msg)
end

@testset "generate_from_params round-trips" begin
    sys = am_system()
    θ = Octofitter.make_arr2nt(sys)(Octofitter.sample_priors(Random.Xoshiro(11), sys))

    sim = Octofitter.generate_from_params(sys, θ; add_noise=false)
    @test sim isa Octofitter.System
    @test [Octofitter.likelihoodname(o) for o in sim.observations] ==
          [Octofitter.likelihoodname(o) for o in sys.observations]
    @test Octofitter.list_parameter_names(sim) == Octofitter.list_parameter_names(sys)
    @test sim.observations[1].table.epoch == AM_EPOCHS_1

    # Noiseless data are the model's own prediction — check that directly
    # against PlanetOrbits rather than trusting the likelihood.
    posys = Octofitter._build_po_system(sys, θ)
    tr = orbitsolve(posys, AM_EPOCHS_1; Octofitter._solvekw(sys)...)
    @test sim.observations[1].table.ra ≈ [raoff(tr[k], :b, :A) for k in eachindex(AM_EPOCHS_1)]
    @test sim.observations[1].table.dec ≈ [decoff(tr[k], :b, :A) for k in eachindex(AM_EPOCHS_1)]

    # With noise the data move, but stay within a few σ.
    Random.seed!(4)
    noisy = Octofitter.generate_from_params(sys, θ; add_noise=true)
    @test noisy.observations[1].table.ra != sim.observations[1].table.ra
    @test all(abs.(noisy.observations[1].table.ra .- sim.observations[1].table.ra) .<
              8 .* sim.observations[1].table.σ_ra)

    # An observation that carries data but cannot simulate is an error, not a
    # silent pass-through of the real measurements.
    ctx = Octofitter.ObsContext(θ, (;), posys, nothing, Int[])
    @test_throws r"does not implement" Octofitter.generate_from_params(
        AMUnsubsettable(Table(epoch=[57000.0]), "stubborn"), ctx; add_noise=false)
    # ...while a prior-shaped term passes straight through.
    blank = Octofitter.BlankLikelihood()
    @test Octofitter.generate_from_params(blank, ctx; add_noise=false) === blank
end

@testset "seeded noise, and the caller's stream left alone" begin
    sys = am_system()
    θ = Octofitter.make_arr2nt(sys)(Octofitter.sample_priors(Random.Xoshiro(11), sys))

    # `generate_from_params` draws its noise from the default RNG (the per-type
    # methods take no `rng`), so a trial is only reproducible if something
    # seeds around it.
    sim(seed) = Octofitter._with_seed(seed) do
        Octofitter.generate_from_params(sys, θ; add_noise=true)
    end
    @test sim(UInt64(7)).observations[1].table.ra == sim(UInt64(7)).observations[1].table.ra
    @test sim(UInt64(7)).observations[1].table.ra != sim(UInt64(8)).observations[1].table.ra

    # ...and a script that seeded before calling still gets its own stream.
    Random.seed!(123)
    expected = randn()
    Random.seed!(123)
    sim(UInt64(7))
    @test randn() == expected
end

@testset "completeness bookkeeping" begin
    masses = [1.0, 10.0]
    seps = [3.0, 9.0, 27.0]
    jobs = completeness_jobs(; masses, separations=seps, n_trials=2)
    @test length(jobs) == 2 * 3 * 2
    @test jobs[1].i_mass == 1 && jobs[1].i_sep == 1 && jobs[1].i_trial == 1
    @test jobs[1].mass == 1.0 && jobs[1].separation == 3.0
    # Seeds are a pure function of the grid position, so re-running one array
    # task reproduces exactly the same trial.
    @test length(unique(j.seed for j in jobs)) == length(jobs)
    @test completeness_jobs(; masses, separations=seps, n_trials=2)[4].seed == jobs[4].seed

    # Assembly is pure counting, and tolerates a grid cell with no results at
    # all (a cluster job that timed out).
    results = [CompletenessResult(jobs[k], (; a=1.0), (;)) for k in 1:3]
    cmap = assemble_completeness(results, (c, θ) -> true; masses, separations=seps)
    @test size(cmap.completeness) == (2, 3)
    @test cmap.n_total[1, 1] == 2 && cmap.n_total[1, 2] == 1
    @test cmap.completeness[1, 1] == 1.0
    @test cmap.n_total[2, 3] == 0
    @test cmap.completeness[2, 3] == 0.0

    half = assemble_completeness(results, (c, θ) -> false; masses, separations=seps)
    @test all(iszero, half.completeness)
end

@testset "parameter overrides address the right slot" begin
    sys = am_system()
    slots = Octofitter._prior_slots(sys)
    # No system priors here (`plx` is derived), then node priors, then
    # observation priors — the order `make_arr2nt` consumes them in.
    @test slots == [(:bodies, :b, :a, 1:1), (:observations, :astrom2, :jitter, 2:2)]

    θ = Octofitter.sample_priors(Random.Xoshiro(5), sys)
    Octofitter._apply_overrides!(θ, sys, (; bodies=(; b=(; a=4.25))))
    @test θ[1] == 4.25
    Octofitter._apply_overrides!(θ, sys, (; observations=(; astrom2=(; jitter=2.5))))
    @test θ[2] == 2.5
    @test Octofitter.make_arr2nt(sys)(θ).bodies.b.a == 4.25

    # A derived variable is not a slot, and saying so beats writing a value
    # nothing will ever read.
    @test_throws r"not a free \(prior\) variable" Octofitter._apply_overrides!(
        θ, sys, (; bodies=(; b=(; e=0.5))))
    # The v1 spelling gets a pointer rather than a missing-key error.
    @test_throws r"nested under `bodies`" Octofitter._apply_overrides!(
        θ, sys, (; planets=(; b=(; a=4.0))))
end

@testset "starting points at the truth" begin
    sys = am_system()
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θ = Octofitter.sample_priors(Random.Xoshiro(9), sys)
    Octofitter._initialize_at_truth!(model, θ; ndraws=4)
    @test length(model.starting_points) == 4
    @test all(p -> p == model.starting_points[1], model.starting_points)
    @test collect(model.invlink(model.starting_points[1])) ≈ θ
end

# Anything that actually samples is slow enough to belong behind the
# integration gate.
if get(ENV, "OCTOFITTER_TEST_MODE", "all") in ("all", "integration")
    @testset "one SBC trial" begin
        sys = am_system()
        θ = Octofitter.sample_priors(Random.Xoshiro(2), sys)
        priorsample, ranks, chain = Octofitter.calibrationhmc(
            sys; rng=Random.Xoshiro(2), θ, verbosity=0,
            adaptation=200, iterations=200)
        @test haskey(ranks, "b_a")
        @test priorsample["b_a"] ≈ θ[1]
        @test 0 <= ranks["b_a"] <= 100
        @test size(chain, 1) == 200
    end

    @testset "one completeness trial" begin
        sys = am_system()
        job = only(completeness_jobs(masses=[5.0], separations=[3.0], n_trials=1))
        res = run_completeness_trial(job, sys,
            m -> octofit(m; adaptation=200, iterations=200, verbosity=0);
            inject=(mass, sep) -> (; bodies=(; b=(; a=sep))), verbosity=0)
        @test res.job === job
        @test res.θ_true.bodies.b.a == 3.0
        @test size(res.chain, 1) == 200

        cmap = assemble_completeness([res],
            (c, θ) -> abs(median(vec(c["b_a"])) - θ.bodies.b.a) < 1.0;
            masses=[5.0], separations=[3.0])
        @test cmap.n_total == fill(1, 1, 1)
    end

    # `sbctrial` is what the SBC tutorial's per-trial cluster script actually
    # calls, and everything it does beyond `calibrationhmc` is serialization.
    # That is exactly where it can break without any signature changing: the
    # sampler settings a user passes may hold a distribution or a function,
    # which TOML cannot represent, so they are stringified on the way out.
    @testset "sbctrial writes its four files" begin
        sys = am_system()
        θ = Octofitter.sample_priors(Random.Xoshiro(2), sys)
        mktempdir() do dir
            saveas = joinpath(dir, "trial-1")
            Octofitter.sbctrial(sys, (;
                    θ, verbosity=0, target_accept=0.85,
                    adaptation=200, iterations=200), saveas)

            for suffix in ("_sampler_parameters.toml", "_parameters.toml",
                           "_rank_stats.toml", "_chains.fits")
                @test isfile(saveas * suffix)
            end

            params = Octofitter.TOML.parsefile(saveas * "_parameters.toml")
            @test params["b_a"] ≈ θ[1]
            ranks = Octofitter.TOML.parsefile(saveas * "_rank_stats.toml")
            @test 0 <= ranks["b_a"] <= 100

            # `target_accept` is consumed as a positional argument, so it is
            # not among the forwarded settings; the rest round-trip as strings.
            settings = Octofitter.TOML.parsefile(saveas * "_sampler_parameters.toml")
            @test settings["iterations"] == "200"
            @test !haskey(settings, "target_accept")

            @test size(Octofitter.loadchain(saveas * "_chains.fits"), 1) == 200
        end
    end
end
