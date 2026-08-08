# OctofitterImages — v2 observation types.
#
# The two things worth proving here are (a) that the multi-target `ImageObs`
# still computes the Ruffio/Mawet matched-filter term v1 computed, one target
# at a time, and (b) that nothing about it privileges a host any more: the
# reference is a ref like any other, and every source's flux is its own body's
# variable.

using Test
using OctofitterImages
using Octofitter
using PlanetOrbits
using AstroImages
using Distributions
using ForwardDiff
using Random
using Statistics
using TypedTables

const EPOCHS = [56000.0, 56800.0, 57600.0]

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# A radially varying 1σ contrast, as a callable — what the `contrast` column
# holds. A struct rather than a closure so the column stays concretely typed.
struct RadialContrast{T}
    σ0::T
    slope::T
end
(c::RadialContrast)(r_px) = c.σ0 + c.slope * r_px

# A smooth, deterministic stand-in for a reduced image: a couple of Gaussian
# blobs plus ripple, so interpolated values differ from point to point.
function fixture_image(; halfwidth=30, amp=1.0)
    ax = float(-halfwidth):float(halfwidth)
    data = [amp * (exp(-((x - 8.0)^2 + (y + 5.0)^2) / 18) +
                   0.4exp(-((x + 11.0)^2 + (y - 3.0)^2) / 40) +
                   0.02cos(x / 3) * sin(y / 4))
            for x in ax, y in ax]
    return AstroImages.recenter(AstroImage(data))
end

# A logL surface peaked at (x₀, y₀) px, falling off with slope 1/s along each
# axis. Deliberately *not* a Gaussian: with the peak on a pixel centre this
# surface is a sum of two piecewise-linear functions whose kinks sit on grid
# nodes, so bilinear interpolation reproduces it exactly and the expected
# log-likelihood is analytic rather than a tolerance.
fixture_map_value(x, y, x0, y0, s) = -(abs(x - x0) + abs(y - y0)) / s
function fixture_map(x0, y0; s=8.0, halfwidth=30)
    ax = float(-halfwidth):float(halfwidth)
    return AstroImages.recenter(AstroImage(
        [fixture_map_value(x, y, x0, y0, s) for x in ax, y in ax]))
end

function image_table(; platescale=10.0, contrast=RadialContrast(0.05, 0.002), epochs=EPOCHS)
    return Table(
        epoch=collect(epochs),
        image=[fixture_image() for _ in epochs],
        platescale=fill(platescale, length(epochs)),
        contrast=fill(contrast, length(epochs)),
    )
end

# Model system: host `A` plus companions `b` (inner, massive) and `c` (outer).
# `b` being massive is deliberate — it is the configuration whose apparent
# positions v1 tried to patch up with its superposition loop.
function model_system(; obs, flux_b=3.8, flux_c=1.1, flux_A=1.0)
    A = Octofitter.Body(name="A", variables=Octofitter.@variables begin
        mass = 1.0
        flux_H = $flux_A
    end)
    b = Octofitter.Body(name="b", about=A, variables=Octofitter.@variables begin
        mass = 30mjup
        flux_H = $flux_b
        a = 3.0
        e = 0.1
        i = 0.6
        ω = 1.2
        Ω = 2.0
        tp = 55000.0
    end)
    c = Octofitter.Body(name="c", about=A, variables=Octofitter.@variables begin
        mass = 8mjup
        flux_H = $flux_c
        a = 9.0
        e = 0.05
        i = 1.1
        ω = 0.3
        Ω = 1.0
        tp = 56000.0
    end)
    return Octofitter.System(name="img", bodies=[A, b, c], observations=obs,
        variables=Octofitter.@variables begin
            plx = 25.0
        end)
end

# The same system with the bodies' fluxes declared as a bare `flux`, i.e. in
# the single unnamed band.
function model_system_oneband(; obs, flux_b=3.8, flux_c=1.1, flux_A=1.0)
    A = Octofitter.Body(name="A", variables=Octofitter.@variables begin
        mass = 1.0
        flux = $flux_A
    end)
    b = Octofitter.Body(name="b", about=A, variables=Octofitter.@variables begin
        mass = 30mjup
        flux = $flux_b
        a = 3.0; e = 0.1; i = 0.6; ω = 1.2; Ω = 2.0; tp = 55000.0
    end)
    c = Octofitter.Body(name="c", about=A, variables=Octofitter.@variables begin
        mass = 8mjup
        flux = $flux_c
        a = 9.0; e = 0.05; i = 1.1; ω = 0.3; Ω = 1.0; tp = 56000.0
    end)
    return Octofitter.System(name="img", bodies=[A, b, c], observations=obs,
        variables=Octofitter.@variables begin
            plx = 25.0
        end)
end

# Everything a likelihood needs, for calling `ln_like`/`simulate` directly.
function obs_context(sys, obs; θ=Float64[])
    nt = Octofitter.make_arr2nt(sys)(θ)
    lnl = Octofitter.make_ln_like(sys)
    posys = lnl.build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep)
    # An observation that declares no variables of its own gets no namespace.
    sym = Symbol(Octofitter.likelihoodname(obs))
    θ_obs = hasproperty(nt.observations, sym) ? getproperty(nt.observations, sym) : (;)
    return Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs]), posys, traj, maps[obs]
end

system_ln_like(sys) = Octofitter.make_ln_like(sys)(sys, Octofitter.make_arr2nt(sys)(Float64[]))

# A literal transcription of the v1 `ln_like(::ImageObs, ::PlanetObservationContext)`
# inner loop, minus the epicyclic superposition — which is identically zero for
# a single-planet system, and which `raoff(sol, target, ref)` replaces exactly
# for every other one. Used to pin the ported likelihood to v1's arithmetic.
function v1_image_ln_like(obs, traj, idx; target, ref, flux, platescale_mult=1.0, northangle=0.0)
    ll = 0.0
    for i_epoch in eachindex(obs.table.epoch)
        sol = traj[idx[i_epoch]]
        ra_raw = raoff(sol, target, ref)
        dec_raw = decoff(sol, target, ref)
        cos_θ = cos(northangle)
        sin_θ = sin(northangle)
        ra_rotated = ra_raw * cos_θ - dec_raw * sin_θ
        dec_rotated = ra_raw * sin_θ + dec_raw * cos_θ
        x = -ra_rotated
        y = +dec_rotated
        platescale = obs.table.platescale[i_epoch] * platescale_mult
        f̃ₓ = obs.table.imageinterp[i_epoch](x / platescale, y / platescale)
        r = √(x^2 + y^2)
        σₓ = obs.table.contrast[i_epoch](r / platescale)
        if !isfinite(f̃ₓ)
            f̃ₓ = zero(typeof(f̃ₓ))
        end
        if !isfinite(σₓ) || iszero(σₓ)
            return -Inf
        end
        ll += -1 / (2 * σₓ^2) * (flux^2 - 2 * flux * f̃ₓ)
    end
    return ll
end

lnlike_allocations(obs, ctx) = @allocated Octofitter.ln_like(obs, ctx)

# ---------------------------------------------------------------------------

@testset "OctofitterImages" begin

    @testset "construction, references, and display" begin
        tab = image_table()
        obs = ImageObs(tab; targets=(:b, :c), ref=:A, band=:H, name="SPHERE")
        @test Octofitter.refspecs(obs) === (Octofitter.BodyRefSpec{:b}(),
            Octofitter.BodyRefSpec{:c}(), Octofitter.BodyRefSpec{:A}())
        @test Octofitter.epochs(obs) == EPOCHS
        @test occursin("(b, c) vs A", sprint(show, MIME"text/plain"(), obs))

        # A single target need not be wrapped in a tuple.
        one = ImageObs(tab; targets=:b, ref=:A, band=:H, name="one")
        @test Octofitter.refspecs(one) === (Octofitter.BodyRefSpec{:b}(),
            Octofitter.BodyRefSpec{:A}())
        @test occursin("b vs A", sprint(show, MIME"text/plain"(), one))
        # …and a vector of them is accepted, but stored as a tuple: a Vector of
        # specs would allocate and lose inference on every evaluation.
        @test ImageObs(tab; targets=[:b, :c], ref=:A, band=:H, name="v").targets isa Tuple

        @test_throws r"`targets` is empty" ImageObs(tab; targets=(), ref=:A, name="e")
        @test_throws r"Expected columns" ImageObs(
            Table(epoch=EPOCHS, platescale=fill(10.0, 3)); targets=:b, ref=:A, name="e")

        # A barycentre emits nothing, so it cannot be a source in an image.
        @test_throws r"barycentre has no" ImageObs(tab; targets=Barycentre, ref=:A, name="e")
        # …but it is a perfectly good reference point.
        @test ImageObs(tab; targets=:b, ref=Barycentre, band=:H, name="e") isa ImageObs

        # An unresolved pair adopts the image's band rather than silently
        # weighting its position by a different one.
        ph = ImageObs(tab; targets=Photocentre((:b, :c)), ref=:A, band=:H, name="p")
        @test ph.targets === (Photocentre(:H, (:b, :c)),)
        @test_throws r"weighted by band" ImageObs(tab; targets=Photocentre(:K, (:b, :c)),
            ref=:A, band=:H, name="p")

        # v1 put the flux on the observation. That cannot describe an image
        # with two companions in it, so it is an error rather than a silent
        # no-op.
        @test_throws r"per-body variables in v2" ImageObs(tab; targets=:b, ref=:A,
            band=:H, name="legacy", variables=Octofitter.@variables begin
                flux ~ Normal(3.8, 0.5)
            end)
    end

    @testset "one target reproduces the v1 matched-filter likelihood" begin
        # The claim under test: with one target, the (A, b) quadratic form
        # collapses to v1's per-epoch −1/(2σₓ²)(f² − 2f f̃ₓ), exactly.
        tab = image_table()
        obs = ImageObs(tab; targets=(:b,), ref=:A, band=:H, name="SPHERE")
        sys = model_system(obs=[obs], flux_b=3.8)
        ctx, _, traj, idx = obs_context(sys, obs)

        expected = v1_image_ln_like(obs, traj, idx; target=:b, ref=:A, flux=3.8)
        @test isfinite(expected)
        @test Octofitter.ln_like(obs, ctx) ≈ expected rtol = 1e-12
        # …and it is the whole system likelihood, since nothing else is in it.
        @test system_ln_like(sys) ≈ expected rtol = 1e-12

        # Same, with the instrument calibration engaged, in both directions.
        obs2 = ImageObs(tab; targets=(:b,), ref=:A, band=:H, name="SPHERE",
            variables=Octofitter.@variables begin
                platescale = 1.04
                northangle = 0.07
            end)
        sys2 = model_system(obs=[obs2], flux_b=3.8)
        ctx2, _, traj2, idx2 = obs_context(sys2, obs2)
        expected2 = v1_image_ln_like(obs2, traj2, idx2; target=:b, ref=:A,
            flux=3.8, platescale_mult=1.04, northangle=0.07)
        @test Octofitter.ln_like(obs2, ctx2) ≈ expected2 rtol = 1e-12
        @test !(expected2 ≈ expected)
    end

    @testset "the flux comes from the body, and only from the body" begin
        tab = image_table()
        obs = ImageObs(tab; targets=(:b,), ref=:A, band=:H, name="SPHERE")
        sys = model_system(obs=[obs], flux_b=2.0)
        ctx, _, traj, idx = obs_context(sys, obs)
        @test Octofitter.ln_like(obs, ctx) ≈
              v1_image_ln_like(obs, traj, idx; target=:b, ref=:A, flux=2.0) rtol = 1e-12

        # Changing `c`'s flux cannot move a likelihood that does not model `c`…
        @test system_ln_like(model_system(obs=[obs], flux_b=2.0, flux_c=1.0)) ==
              system_ln_like(model_system(obs=[obs], flux_b=2.0, flux_c=9.0))
        # …but changing `b`'s does.
        @test system_ln_like(model_system(obs=[obs], flux_b=2.5)) !=
              system_ln_like(model_system(obs=[obs], flux_b=2.0))

        # The bare `flux` variable is the single unnamed band, and `band=`
        # may then be omitted.
        obsnb = ImageObs(tab; targets=(:b,), ref=:A, name="SPHERE")
        sysnb = model_system_oneband(obs=[obsnb], flux_b=2.0)
        ctxnb, _, trajnb, idxnb = obs_context(sysnb, obsnb)
        @test Octofitter.ln_like(obsnb, ctxnb) ≈
              v1_image_ln_like(obsnb, trajnb, idxnb; target=:b, ref=:A, flux=2.0) rtol = 1e-12
    end

    @testset "several companions in one image" begin
        # The headline change: one image, one likelihood, N sources — each with
        # its own body flux. Under the diagonal-A reduction that must equal the
        # sum of the individual terms, which is exactly the statement that this
        # reduction treats the companions as independent.
        tab = image_table()
        both = ImageObs(tab; targets=(:b, :c), ref=:A, band=:H, name="SPHERE")
        sysb = model_system(obs=[both], flux_b=3.8, flux_c=1.1)
        ctxb, _, traj, idx = obs_context(sysb, both)

        ll_b = v1_image_ln_like(both, traj, idx; target=:b, ref=:A, flux=3.8)
        ll_c = v1_image_ln_like(both, traj, idx; target=:c, ref=:A, flux=1.1)
        @test Octofitter.ln_like(both, ctxb) ≈ ll_b + ll_c rtol = 1e-12

        # A dark companion contributes nothing at all: f_j = 0 zeroes both the
        # fᵀAf and the fᵀb term.
        @test system_ln_like(model_system(obs=[both], flux_b=3.8, flux_c=0.0)) ≈
              ll_b rtol = 1e-12

        # And `targets` is a *model* statement, not a summary of who is bright:
        # dropping `c` from the list is not the same model as `c` being in it.
        onlyb = ImageObs(tab; targets=(:b,), ref=:A, band=:H, name="SPHERE")
        @test system_ln_like(model_system(obs=[onlyb], flux_b=3.8, flux_c=1.1)) ≈
              ll_b rtol = 1e-12
        @test !(system_ln_like(model_system(obs=[both], flux_b=3.8, flux_c=1.1)) ≈ ll_b)
    end

    @testset "no privileged primary" begin
        tab = image_table()
        # A companion measured against *another companion* — a query v1's
        # per-planet contexts could not express at all.
        cb = ImageObs(tab; targets=(:c,), ref=:b, band=:H, name="SPHERE")
        sys = model_system(obs=[cb], flux_c=1.1)
        ctx, _, traj, idx = obs_context(sys, cb)
        @test Octofitter.ln_like(cb, ctx) ≈
              v1_image_ln_like(cb, traj, idx; target=:c, ref=:b, flux=1.1) rtol = 1e-12

        # The reference genuinely changes the modelled position: measured from
        # the host, from the barycentre and from `b` are three different tracks.
        ca = ImageObs(tab; targets=(:c,), ref=:A, band=:H, name="SPHERE")
        cbary = ImageObs(tab; targets=(:c,), ref=Barycentre, band=:H, name="SPHERE")
        xs = map((cb, ca, cbary)) do o
            c, = obs_context(model_system(obs=[o]), o)
            Octofitter.simulate(o, c).x
        end
        @test xs[1] != xs[2]
        @test xs[2] != xs[3]

        # An unresolved pair is one source at the pair's photocentre, whose
        # flux is the sum of its members'.
        ph = ImageObs(tab; targets=(Photocentre(:H, (:b, :c)),), ref=:A, band=:H, name="SPHERE")
        sysp = model_system(obs=[ph], flux_b=3.0, flux_c=1.0)
        ctxp, posysp, trajp, idxp = obs_context(sysp, ph)
        @test Octofitter.ln_like(ph, ctxp) ≈ v1_image_ln_like(ph, trajp, idxp;
            target=PlanetOrbits.photocentre(posysp, :b, :c; band=:H), ref=:A,
            flux=4.0) rtol = 1e-12
    end

    @testset "off the detector, and unusable contrast" begin
        # A source predicted outside the image reads zero flux there, leaving
        # −f²/2σ² — evidence against a bright companion, not a missing term and
        # not a NaN. (This is v1's behaviour, kept deliberately.)
        tiny = Table(epoch=[57000.0], image=[fixture_image(halfwidth=3)],
            platescale=[0.05], contrast=[RadialContrast(0.05, 0.0)])
        obs = ImageObs(tiny; targets=(:b,), ref=:A, band=:H, name="tiny")
        sys = model_system(obs=[obs], flux_b=3.0)
        ctx, _, _, _ = obs_context(sys, obs)
        sim = Octofitter.simulate(obs, ctx)
        @test !isfinite(obs.table.imageinterp[1](sim.x[1], sim.y[1]))
        ll = Octofitter.ln_like(obs, ctx)
        @test isfinite(ll)
        @test ll ≈ -3.0^2 / (2 * 0.05^2) rtol = 1e-12

        # A σ the data cannot supply is a different thing, and kills the sample.
        zeroσ = Table(epoch=[57000.0], image=[fixture_image()], platescale=[10.0],
            contrast=[RadialContrast(0.0, 0.0)])
        obsz = ImageObs(zeroσ; targets=(:b,), ref=:A, band=:H, name="z")
        ctxz, _, _, _ = obs_context(model_system(obs=[obsz]), obsz)
        @test Octofitter.ln_like(obsz, ctxz) == -Inf
    end

    @testset "contrast maps and measured contrast" begin
        # A 2-D contrast map is read at the source position rather than at its
        # separation.
        cmapdata = AstroImages.recenter(AstroImage(
            [0.04 + 0.001 * hypot(x, y) for x in -30.0:30.0, y in -30.0:30.0]))
        tab = Table(epoch=[57000.0], image=[fixture_image()], platescale=[10.0],
            contrastmap=[cmapdata])
        obs = ImageObs(tab; targets=(:b,), ref=:A, band=:H, name="cm")
        @test hasproperty(obs.table, :contrastmapinterp)
        ctx, _, _, _ = obs_context(model_system(obs=[obs]), obs)
        @test isfinite(Octofitter.ln_like(obs, ctx))

        # With neither column, the contrast is measured from the image itself.
        bare = Table(epoch=[57000.0], image=[fixture_image()], platescale=[10.0])
        obsb = ImageObs(bare; targets=(:b,), ref=:A, band=:H, name="auto")
        @test hasproperty(obsb.table, :contrast)
        ctxb, _, _, _ = obs_context(model_system(obs=[obsb]), obsb)
        @test isfinite(Octofitter.ln_like(obsb, ctxb))

        # `contrast`/`contrast_interp` themselves still work on a real image.
        c = OctofitterImages.contrast(fixture_image())
        @test length(c.separation) == length(c.contrast)
        @test all(>=(0), filter(isfinite, c.contrast))
        @test OctofitterImages.contrast_interp(fixture_image())(5.0) > 0
        @test size(OctofitterImages.imgsep(fixture_image())) == (61, 61)
    end

    @testset "bands" begin
        tab = image_table()
        A = Octofitter.Body(name="A", variables=Octofitter.@variables begin
            mass = 1.0
            flux_H = 1.0
            flux_K = 2.0
        end)
        b = Octofitter.Body(name="b", about=A, variables=Octofitter.@variables begin
            mass = 30mjup
            flux_H = 3.8
            flux_K = 1.9
            a = 3.0; e = 0.1; i = 0.6; ω = 1.2; Ω = 2.0; tp = 55000.0
        end)
        mk(band) = Octofitter.System(name="two", bodies=[A, b],
            observations=[ImageObs(tab; targets=(:b,), ref=:A, band, name="SPHERE")],
            variables=Octofitter.@variables begin
                plx = 25.0
            end)
        @test_throws r"several bands are defined" system_ln_like(mk(nothing))
        llH = system_ln_like(mk(:H))
        llK = system_ln_like(mk(:K))
        @test isfinite(llH) && isfinite(llK) && llH != llK
        @test_throws r"no band :J" system_ln_like(mk(:J))
    end

    @testset "hot-loop hygiene and ForwardDiff" begin
        tab = image_table()
        obs = ImageObs(tab; targets=(:b, :c), ref=:A, band=:H, name="SPHERE")
        sys = model_system(obs=[obs])
        ctx, _, _, _ = obs_context(sys, obs)

        # Refs resolve to a tuple outside the epoch loop and cost nothing.
        Octofitter.resolverefs(ctx, obs.targets)
        @test (@allocated Octofitter.resolverefs(ctx, obs.targets)) == 0
        @test Octofitter.resolverefs(ctx, obs.targets) isa Tuple

        # …and the epoch loop itself never reaches the heap.
        lnlike_allocations(obs, ctx)
        @test lnlike_allocations(obs, ctx) == 0

        # Differentiable through the fluxes and through the calibration.
        obsd = ImageObs(tab; targets=(:b, :c), ref=:A, band=:H, name="SPHERE",
            variables=Octofitter.@variables begin
                platescale ~ truncated(Normal(1, 0.01), lower=0)
                northangle ~ Normal(0, 0.01)
            end)
        A = Octofitter.Body(name="A", variables=Octofitter.@variables begin
            mass = 1.0
            flux_H = 1.0
        end)
        b = Octofitter.Body(name="b", about=A, variables=Octofitter.@variables begin
            mass = 30mjup
            flux_H ~ Normal(3.8, 0.5)
            a ~ truncated(Normal(3.0, 0.2), lower=0)
            e = 0.1; i = 0.6; ω = 1.2; Ω = 2.0; tp = 55000.0
        end)
        c = Octofitter.Body(name="c", about=A, variables=Octofitter.@variables begin
            mass = 8mjup
            flux_H ~ Normal(1.1, 0.5)
            a ~ truncated(Normal(9.0, 0.2), lower=0)
            e = 0.05; i = 1.1; ω = 0.3; Ω = 1.0; tp = 56000.0
        end)
        sysd = Octofitter.System(name="ad", bodies=[A, b, c], observations=[obsd],
            variables=Octofitter.@variables begin
                plx = 25.0
            end)
        arr2nt = Octofitter.make_arr2nt(sysd)
        lnl = Octofitter.make_ln_like(sysd)
        θv = Octofitter.sample_priors(Random.Xoshiro(3), sysd)
        f(θ) = lnl(sysd, arr2nt(θ))
        @test isfinite(f(θv))
        g = ForwardDiff.gradient(f, θv)
        @test all(isfinite, g)
        @test all(!iszero, g)
        # A finite difference in each coordinate, to be sure the AD path and
        # the primal path see the same model.
        for k in eachindex(θv)
            h = 1e-6 * max(abs(θv[k]), 1)
            δ = zeros(length(θv)); δ[k] = h
            @test (f(θv .+ δ) - f(θv .- δ)) / 2h ≈ g[k] rtol = 1e-4
        end
    end

    @testset "epoch subsetting" begin
        tab = image_table()
        obs = ImageObs(tab; targets=(:b, :c), ref=:A, band=:H, name="SPHERE")
        sub = Octofitter.likeobj_from_epoch_subset(obs, [1, 3])
        @test sub isa ImageObs
        @test Octofitter.epochs(sub) == EPOCHS[[1, 3]]
        @test sub.targets === obs.targets
        @test sub.ref === obs.ref
        @test sub.band === obs.band

        # The parts add up: two disjoint subsets recover the whole.
        whole = system_ln_like(model_system(obs=[obs]))
        s2 = Octofitter.likeobj_from_epoch_subset(obs, [2])
        @test system_ln_like(model_system(obs=[sub])) +
              system_ln_like(model_system(obs=[s2])) ≈ whole rtol = 1e-12
    end

    @testset "generating replicate images" begin
        tab = image_table()
        obs = ImageObs(tab; targets=(:b, :c), ref=:A, band=:H, name="SPHERE")
        sys = model_system(obs=[obs], flux_b=3.8, flux_c=1.1)
        ctx, _, _, _ = obs_context(sys, obs)

        rep = Octofitter.generate_from_params(obs, ctx; add_noise=false)
        @test rep isa ImageObs
        @test rep.targets === obs.targets && rep.ref === obs.ref
        @test Octofitter.epochs(rep) == EPOCHS

        # A noiseless replicate puts each source's flux exactly at its own
        # modelled position, so re-evaluating there recovers f̃ₓ = f and the
        # likelihood is +½ Σ f²/σ² — the peak of the profiled form.
        sim = Octofitter.simulate(obs, ctx)
        ctxr, _, _, _ = obs_context(model_system(obs=[rep], flux_b=3.8, flux_c=1.1), rep)
        expected = 0.0
        for i in eachindex(obs.table.epoch), j in eachindex(sim.flux)
            σ = obs.table.contrast[i](hypot(sim.x[j, i], sim.y[j, i]))
            expected += sim.flux[j]^2 / (2σ^2)
        end
        @test Octofitter.ln_like(rep, ctxr) ≈ expected rtol = 1e-9

        # With noise the images differ, but the geometry does not.
        noisy = Octofitter.generate_from_params(obs, ctx; add_noise=true)
        @test collect(noisy.table.image[1]) != collect(rep.table.image[1])
        @test std(collect(noisy.table.image[1])) > 0
    end

    @testset "a real reduced image stack" begin
        # Everything above runs on synthetic arrays. This one walks the path
        # the docs walk a user through — a multi-extension FITS off disk,
        # `AstroImages.recenter`, and a contrast measured from the data — so
        # that a change in AstroImages' array or dimension types shows up here
        # rather than in a docs build. Guarded on the file, since the
        # subpackage is testable on its own.
        fitspath = joinpath(@__DIR__, "..", "..", "image-examples-1.fits")
        if !isfile(fitspath)
            @warn "image-examples-1.fits not found; skipping the real-image tests" fitspath
        else
            imgs = AstroImages.load(fitspath, :)
            dat = Table(
                epoch=[1238.6, 1584.7, 3220.0, 7495.9, 7610.4],
                image=[AstroImages.recenter(img) for img in imgs],
                platescale=fill(10.0, length(imgs)),
            )
            obs = ImageObs(dat; targets=(:b, :c), ref=:A, band=:H, name="SPHERE")
            sys = model_system(obs=[obs], flux_b=3.8, flux_c=1.1)
            ctx, _, traj, idx = obs_context(sys, obs)

            # Both companions land on the detector at every epoch, so this
            # exercises real matched-filter readouts. Some of them fall in the
            # masked region the reduction left as NaN — that is the branch v1
            # read as "no flux measured here", and it must not poison the sum.
            sim = Octofitter.simulate(obs, ctx)
            @test all(isfinite, sim.x) && all(isfinite, sim.y)
            xs, ys = parent.(dims(obs.table.image[1]))
            @test all(first(xs) .<= sim.x .<= last(xs)) &&
                  all(first(ys) .<= sim.y .<= last(ys))
            @test any(isfinite(obs.table.imageinterp[i](sim.x[j, i], sim.y[j, i]))
                      for i in eachindex(obs.table.epoch), j in 1:2)

            ll = Octofitter.ln_like(obs, ctx)
            @test isfinite(ll)
            # Additivity on real data: the diagonal-A form is still the sum of
            # v1's per-companion terms.
            @test ll ≈ v1_image_ln_like(obs, traj, idx; target=:b, ref=:A, flux=3.8) +
                       v1_image_ln_like(obs, traj, idx; target=:c, ref=:A, flux=1.1) rtol = 1e-12
            lnlike_allocations(obs, ctx)
            @test lnlike_allocations(obs, ctx) == 0
        end
    end

    # -----------------------------------------------------------------------
    # Log-likelihood maps
    # -----------------------------------------------------------------------

    @testset "LogLikelihoodMapObs" begin
        # Build the maps so each peaks exactly where `b` is modelled.
        probe = ImageObs(image_table(); targets=(:b,), ref=:A, band=:H, name="probe")
        pctx, _, _, _ = obs_context(model_system(obs=[probe]), probe)
        psim = Octofitter.simulate(probe, pctx)

        # Peaks on pixel centres, so the surface is bilinear-exact and the
        # expected log-likelihood is analytic.
        xpk = [round(psim.x[1, i]) for i in eachindex(EPOCHS)]
        ypk = [round(psim.y[1, i]) for i in eachindex(EPOCHS)]
        mapdat = Table(epoch=collect(EPOCHS), platescale=fill(10.0, length(EPOCHS)),
            map=[fixture_map(xpk[i], ypk[i]) for i in eachindex(EPOCHS)])
        obs = LogLikelihoodMapObs(mapdat; target=:b, ref=:A, name="GRAVITY")
        @test Octofitter.refspecs(obs) === (Octofitter.BodyRefSpec{:b}(),
            Octofitter.BodyRefSpec{:A}())
        @test occursin("b vs A", sprint(show, MIME"text/plain"(), obs))
        @test Octofitter.epochs(obs) == EPOCHS

        sys = model_system(obs=[obs])
        ctx, _, _, _ = obs_context(sys, obs)
        sim = Octofitter.simulate(obs, ctx)

        # Exactly the surfaces, read at the modelled positions — the plumbing
        # claim, free of interpolation error.
        exact = sum(obs.table.mapinterp[i](sim.x[i], sim.y[i]) for i in eachindex(EPOCHS))
        @test Octofitter.ln_like(obs, ctx) == exact
        @test system_ln_like(sys) ≈ exact rtol = 1e-12
        # …and the value is the analytic one: each surface read at the position
        # the model puts the companion, offset from that map's own peak.
        analytic = sum(fixture_map_value(psim.x[1, i], psim.y[1, i], xpk[i], ypk[i], 8.0)
                       for i in eachindex(EPOCHS))
        @test Octofitter.ln_like(obs, ctx) ≈ analytic rtol = 1e-12
        @test analytic > -0.2       # sub-pixel: the maps do peak on the model

        # Displace the maps by 2 px and the likelihood drops by the slope times
        # the displacement, exactly.
        shifted = Table(epoch=collect(EPOCHS), platescale=fill(10.0, length(EPOCHS)),
            map=[fixture_map(xpk[i] + 2.0, ypk[i]) for i in eachindex(EPOCHS)])
        obs2 = LogLikelihoodMapObs(shifted; target=:b, ref=:A, name="GRAVITY")
        ctx2, _, _, _ = obs_context(model_system(obs=[obs2]), obs2)
        @test Octofitter.ln_like(obs2, ctx2) ≈
              sum(fixture_map_value(psim.x[1, i], psim.y[1, i], xpk[i] + 2, ypk[i], 8.0)
                  for i in eachindex(EPOCHS)) rtol = 1e-12
        @test Octofitter.ln_like(obs2, ctx2) < Octofitter.ln_like(obs, ctx)

        # …and it is the *reference* that defines where zero is: measured from
        # `b` instead of `A`, the same maps say something else entirely.
        obs3 = LogLikelihoodMapObs(mapdat; target=:c, ref=:b, name="GRAVITY")
        ctx3, _, _, _ = obs_context(model_system(obs=[obs3]), obs3)
        @test Octofitter.ln_like(obs3, ctx3) < -1.0

        # Off the map, the fill value applies rather than a NaN.
        far = Table(epoch=[57000.0], platescale=[1e-4],
            map=[fixture_map(0.0, 0.0, halfwidth=3)])
        obsf = LogLikelihoodMapObs(far; target=:b, ref=:A, name="far")
        ctxf, _, _, _ = obs_context(model_system(obs=[obsf]), obsf)
        @test Octofitter.ln_like(obsf, ctxf) ≈ obsf.table.fillvalue[1] rtol = 1e-12

        # Subsetting, for cross-validation.
        sub = Octofitter.likeobj_from_epoch_subset(obs, [2, 3])
        @test sub isa LogLikelihoodMapObs
        @test Octofitter.epochs(sub) == EPOCHS[[2, 3]]
        @test sub.target === obs.target && sub.ref === obs.ref

        # Replicates: the measured surface, recentred on the model. The peak of
        # each replicate must land on the pixel nearest the modelled position,
        # and the replicate must therefore fit better than the maps it came
        # from. (Its *value* is not analytic — translating by a fractional
        # offset resamples the surface through its own interpolator.)
        rep = Octofitter.generate_from_params(obs2, ctx2; add_noise=false)
        @test rep isa LogLikelihoodMapObs
        ctxr, _, _, _ = obs_context(model_system(obs=[rep]), rep)
        for i in eachindex(EPOCHS)
            m = rep.table.map[i]
            ipk = argmax(collect(m))
            xs, ys = parent.(dims(m))
            @test hypot(xs[ipk[1]] - psim.x[1, i], ys[ipk[2]] - psim.y[1, i]) <= 1.0
        end
        @test Octofitter.ln_like(rep, ctxr) > Octofitter.ln_like(obs2, ctx2)
        @test_throws r"no sampling distribution" Octofitter.generate_from_params(
            obs2, ctx2; add_noise=true)
    end
end

# ---------------------------------------------------------------------------
# Plotting protocol
#
# The 1-D reduction of each type: for an image, the matched-filter readout at
# each modelled source position against that source's fitted flux; for a
# likelihood map, the shortfall below the map's own peak.
# ---------------------------------------------------------------------------
@testset "plotting protocol" begin
    @testset "ImageObs exposes one flux channel per target" begin
        obs = ImageObs(image_table(); targets=(:b, :c), ref=:A, band=:H, name="SPHERE")
        sys = model_system(obs=[obs])
        @test isempty(Octofitter.unplottable_observations(sys))
        chs = Octofitter.plotchannels(obs)
        @test [ch.name for ch in chs] == [:flux_b, :flux_c]
        # Flux is epoch-independent and no PlanetOrbits observable produces it.
        @test all(ch -> ch.query === nothing, chs)

        ctx, _, _, _ = obs_context(sys, obs)
        r = Octofitter.residuals(obs, ctx)
        @test keys(r) == (:flux_b, :flux_c)
        @test length(r.flux_b.epoch) == length(EPOCHS)
        @test r.flux_b.model == fill(3.8, length(EPOCHS))   # b's `flux_H`
        @test r.flux_c.model == fill(1.1, length(EPOCHS))
        @test r.flux_b.resid ≈ r.flux_b.data .- r.flux_b.model
        # The implied χ² is the estimator's, up to the additive f̃²/2σ² that
        # `ln_like` drops as a constant.
        sim = Octofitter.simulate(obs, ctx)
        dropped = sum(r.flux_b.data .^ 2 ./ (2 .* r.flux_b.σ .^ 2)) +
                  sum(r.flux_c.data .^ 2 ./ (2 .* r.flux_c.σ .^ 2))
        χ² = sum(r.flux_b.resid .^ 2 ./ (2 .* r.flux_b.σ .^ 2)) +
             sum(r.flux_c.resid .^ 2 ./ (2 .* r.flux_c.σ .^ 2))
        @test Octofitter.ln_like(obs, ctx) ≈ -χ² + dropped rtol = 1e-10
        @test sim.flux[1] == 3.8
    end

    @testset "LogLikelihoodMapObs reports its shortfall below the map peak" begin
        obs = LogLikelihoodMapObs(
            Table(epoch=EPOCHS, map=[fixture_map(3.0, -4.0) for _ in EPOCHS],
                platescale=fill(10.0, length(EPOCHS)));
            target=:b, ref=:A, name="likemap")
        sys = model_system(obs=[obs])
        @test isempty(Octofitter.unplottable_observations(sys))
        @test [ch.name for ch in Octofitter.plotchannels(obs)] == [:lnlike]
        # A logL surface has no (data, model, σ) triple, so it declares its own
        # panel rather than being drawn as a residual series with an invented σ.
        @test !isempty(Octofitter.defaultpanels(obs))

        ctx, _, _, _ = obs_context(sys, obs)
        r = Octofitter.residuals(obs, ctx).lnlike
        # `data` is the map's own maximum: the best score any orbit could get.
        @test r.data == fill(fixture_map_value(3.0, -4.0, 3.0, -4.0, 8.0), length(EPOCHS))
        # `model` is what this sample scores, and they sum to `ln_like`.
        @test sum(r.model) ≈ Octofitter.ln_like(obs, ctx) rtol = 1e-10
        # The shortfall is non-negative, and zero only at the peak.
        @test all(>=(0), r.resid)
        @test r.σ == ones(length(EPOCHS))
    end
end
