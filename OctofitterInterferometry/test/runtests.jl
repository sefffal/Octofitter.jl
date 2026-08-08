using Test
using Octofitter
using OctofitterInterferometry
using PlanetOrbits
using Distributions
using LinearAlgebra
using TypedTables
using ForwardDiff
using Random

const OI = OctofitterInterferometry

# ---------------------------------------------------------------------------
# Fixtures
#
# A four-telescope array: six baselines ordered (1,2) (1,3) (1,4) (2,3) (2,4)
# (3,4), and the four closure triangles over them. The index vectors say which
# baseline phases combine as +1 +2 −3, and they match `GRAVITY_T3_DESIGN` row
# for row — which is what lets the same fixture drive both the closure-phase
# and the kernel-phase paths.
# ---------------------------------------------------------------------------

const IDX_CPS1 = [1, 1, 2, 4]
const IDX_CPS2 = [4, 5, 6, 6]
const IDX_CPS3 = [2, 3, 3, 5]

const EPOCHS = [59000.0, 59200.0, 59450.0]

"Baselines [m] of a compact four-telescope array, as (u, v) pairs."
const BASELINES = [
    (46.6, 12.3), (-9.1, 54.8), (33.7, -41.2),
    (-55.7, 42.5), (-12.9, -53.5), (42.8, -96.0),
]

function uv_tables(eff_wave; jitter_seed=0.0)
    Λ = length(eff_wave)
    u = zeros(length(BASELINES), Λ)
    v = zeros(length(BASELINES), Λ)
    for (b, (ub, vb)) in enumerate(BASELINES), l in 1:Λ
        u[b, l] = (ub + jitter_seed * b) / eff_wave[l]
        v[b, l] = (vb - jitter_seed * b) / eff_wave[l]
    end
    return u, v
end

"""
One synthetic exposure. The closure phases are zeros — the value of the *data*
never matters to the tests below, only that the model and the reference
implementation see the same numbers.
"""
function exposure(epoch; Λ=3, use_vis2=false, jitter_seed=0.0, extra...)
    eff_wave = collect(range(2.0e-6, 2.4e-6, length=Λ))
    u, v = uv_tables(eff_wave; jitter_seed)
    n_cp = length(IDX_CPS1)
    n_bl = length(BASELINES)
    return (;
        epoch,
        eff_wave,
        u, v,
        cps_data=fill(3.0, n_cp, Λ),
        dcps=fill(2.5, n_cp, Λ),
        vis2_data=fill(0.9, n_bl, Λ),
        dvis2=fill(0.05, n_bl, Λ),
        index_cps1=IDX_CPS1, index_cps2=IDX_CPS2, index_cps3=IDX_CPS3,
        use_vis2,
        extra...
    )
end

synthetic_table(; kwargs...) = Table([exposure(e; kwargs...) for e in EPOCHS])

function model_system(; obs, flux_A=1.0, flux_b=0.3, flux_c=nothing, a_b=2.0, sysvars=nothing)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.2
        flux_K = $flux_A
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 20mjup
        flux_K = $flux_b
        a = $a_b
        e = 0.12
        i = 0.6
        ω = 1.0
        Ω = 0.4
        tp = 58800.0
    end)
    bodies = Any[A, b]
    if !isnothing(flux_c)
        c = Octofitter.Body(name="c", about=A, variables=@variables begin
            mass = 8mjup
            flux_K = $flux_c
            a = 5.5
            e = 0.03
            i = 0.9
            ω = 0.2
            Ω = 1.7
            tp = 57900.0
        end)
        push!(bodies, c)
    end
    vars = isnothing(sysvars) ? (@variables begin
        plx = 100.0
    end) : sysvars
    return Octofitter.System(name="sys", bodies=bodies, observations=obs, variables=vars)
end

"Build an `ObsContext` for `obs` inside `sys`, at the parameters' point values."
function ctx_for(sys, obs)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep)
    # `normalizename`, not `Symbol`: codegen sanitizes the observation's name
    # into a NamedTuple key, so "GRAVITY-WIDE" is stored as `GRAVITY_WIDE`.
    sym = Octofitter.normalizename(Octofitter.likelihoodname(obs))
    θ_obs = hasproperty(nt.observations, sym) ? getproperty(nt.observations, sym) : (;)
    return Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs]), posys, traj
end

# ---------------------------------------------------------------------------
# The v1 implementation, transcribed
#
# This is v1's `ln_like(::InterferometryObs, ::SystemObservationContext)` inner
# loop, copied line for line from
# `packages/Octofitter-v1-ref/OctofitterInterferometry/src/OctofitterInterferometry.jl`
# and given a v2 trajectory to read positions from. The epicyclic superposition
# loop is dropped because it is a no-op for a single companion — the case this
# check covers, and the only case it was ever right for.
#
# Everything else is kept exactly as v1 wrote it: `cvis_model .+= 1.0` for the
# host, the `1/(1 + Σ contrast)` normalization spelled as a multiplication by
# the reciprocal, and the two-pass accumulation of the closure-phase term.
# ---------------------------------------------------------------------------

function v1_ln_like(table, traj, epoch_index, contrasts; σ_cp_jitter=0.0)
    ll = 0.0
    for i_epoch in eachindex(table.epoch)
        sol = traj[epoch_index[i_epoch]]
        cps_model = zeros(size(table.cps_data[i_epoch][:, 1]))
        cvis_model = zeros(ComplexF64, size(table.u[i_epoch][:, 1]))
        for i_wave in axes(table.u[i_epoch], 2)
            u = @views table.u[i_epoch][:, i_wave]
            v = @views table.v[i_epoch][:, i_wave]
            cps_data = @views table.cps_data[i_epoch][:, i_wave]
            σ_cp = @views table.dcps[i_epoch][:, i_wave]

            cvis_model .= 0
            cps_model .= 0
            norm_factor_model = 0.0
            for (i_planet, contrast) in enumerate(contrasts)
                # v1: `raoff(sol)` of the planet's own per-planet solution,
                # i.e. the planet relative to the host.
                Δra = raoff(sol, PlanetOrbits.BodyRef(i_planet + 1), PlanetOrbits.BodyRef(1))
                Δdec = decoff(sol, PlanetOrbits.BodyRef(i_planet + 1), PlanetOrbits.BodyRef(1))
                OI.cvis_bin!(cvis_model; Δdec, Δra, contrast, u, v)
                norm_factor_model += contrast
            end
            cvis_model .+= 1.0
            cvis_model .*= 1.0 / (1.0 + norm_factor_model)
            OI.closurephase!(cps_model; vis=cvis_model,
                index_cps1=table.index_cps1[i_epoch],
                index_cps2=table.index_cps2[i_epoch],
                index_cps3=table.index_cps3[i_epoch])

            const_cp = 0.0
            for I in eachindex(σ_cp)
                const_cp -= log(2π * (σ_cp[I]^2 + σ_cp_jitter^2))
            end
            const_cp /= 2
            lnlike_cp = 0.0
            for I in eachindex(σ_cp, cps_data, cps_model)
                σ² = σ_cp_jitter^2 + σ_cp[I]^2
                lnlike_cp -= 0.5 * (cps_data[I] - cps_model[I])^2 / σ²
            end
            lnlike_cp += const_cp
            ll += lnlike_cp
        end
    end
    return ll
end

# ---------------------------------------------------------------------------

@testset "OctofitterInterferometry" begin

@testset "construction, references and display" begin
    tab = synthetic_table()
    obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="NIRISS")
    @test Octofitter.refspecs(obs) ==
          (Octofitter.refspec(:A), Octofitter.refspec(:b), Octofitter.refspec(:A))
    @test Octofitter.epochs(obs) == EPOCHS
    @test sprint(show, MIME"text/plain"(), obs) |> s -> occursin("(A, b) vs A", s)

    # The v1 `flux` vector on the observation is gone; targets are bodies and
    # carry their own band fluxes.
    @test_throws r"must name bodies" InterferometryObs(
        tab; targets=(Photocentre(:K),), ref=:A, band=:K, name="x")
    @test_throws r"at least one source" InterferometryObs(
        tab; targets=(), ref=:A, band=:K, name="x")

    # Rows are sorted by epoch, like every other v2 observation.
    shuffled = Table([exposure(e) for e in reverse(EPOCHS)])
    @test InterferometryObs(shuffled; targets=(:A, :b), ref=:A, band=:K, name="x").table.epoch == EPOCHS

    # `Body` model nodes name themselves, which is the spelling the docs use.
    Abody = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.2
        flux_K = 1.0
    end)
    @test Octofitter.refspecs(
        InterferometryObs(tab; targets=(Abody, :b), ref=Abody, band=:K, name="x")) ==
          (Octofitter.refspec(:A), Octofitter.refspec(:b), Octofitter.refspec(:A))

    # `fiber_pointing` without `fiber_coupling` would silently do nothing.
    @test_throws r"no effect" InterferometryObs(
        tab; targets=(:A, :b), ref=:A, band=:K, name="x", fiber_pointing=:A)
end

@testset "band selection" begin
    tab = synthetic_table()
    obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, name="vis")   # no band
    sys = model_system(obs=[obs])
    ctx, _, _ = ctx_for(sys, obs)
    # One band declared, so `band=nothing` picks it.
    @test isfinite(Octofitter.ln_like(obs, ctx))

    # A band nothing declares is caught by PlanetOrbits, naming the bands there are.
    obs2 = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:H, name="vis")
    sys2 = model_system(obs=[obs2])
    ctx2, _, _ = ctx_for(sys2, obs2)
    @test_throws r"no band :H" Octofitter.ln_like(obs2, ctx2)

    # A band every target is dark in is 0/0, and says so rather than
    # returning NaN out of the middle of a sampler.
    obs3 = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="vis")
    sys3 = model_system(obs=[obs3], flux_A=0.0, flux_b=0.0)
    ctx3, _, _ = ctx_for(sys3, obs3)
    @test_throws r"zero flux" Octofitter.ln_like(obs3, ctx3)
end

@testset "the host's flux_K = 1.0 reproduces v1 bit-for-bit" begin
    # The single most valuable check in this port. v1 hard-coded the host into
    # the visibility sum as `cvis_model .+= 1.0` and normalized by
    # `1/(1 + Σ contrast)`; v2 puts the host in `targets` with an ordinary
    # `flux_K`. Setting that flux to 1.0 must recover v1 *exactly* — not
    # approximately — or the port has moved the physics.
    for jitter in (0.0, 4.0), contrast in (0.3, 0.02)
        variables = @variables begin
            σ_cp_jitter = $jitter
        end
        tab = synthetic_table()
        obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K,
                                name="NIRISS", variables)
        sys = model_system(obs=[obs], flux_A=1.0, flux_b=contrast)
        ctx, _, traj = ctx_for(sys, obs)

        expected = v1_ln_like(obs.table, traj, ctx.epoch_index, (contrast,);
                              σ_cp_jitter=jitter)
        @test Octofitter.ln_like(obs, ctx) === expected
    end
end

@testset "two companions agree with v1 to roundoff" begin
    # With more than one companion the *order* of the visibility sum differs —
    # v1 accumulated the companions first and added the host last — so the
    # agreement is to roundoff rather than bit-for-bit. (The epicyclic loop v1
    # would also have run here is gone; it is a no-op for these coplanar
    # astrocentric orbits only because both companions orbit the same body,
    # which is exactly the assumption it encoded.)
    tab = synthetic_table()
    obs = InterferometryObs(tab; targets=(:A, :b, :c), ref=:A, band=:K, name="NIRISS")
    sys = model_system(obs=[obs], flux_A=1.0, flux_b=0.3, flux_c=0.11)
    ctx, _, traj = ctx_for(sys, obs)
    @test Octofitter.ln_like(obs, ctx) ≈ v1_ln_like(obs.table, traj, ctx.epoch_index, (0.3, 0.11)) rtol = 1e-10
end

@testset "no body is privileged" begin
    tab = synthetic_table()
    # Closure phases are invariant to the phase centre: shifting it multiplies
    # every complex visibility by one global phase, which cancels in the
    # triangle sum. So the same three sources measured against A, against b or
    # against the barycentre give the same closure phases — the property that
    # makes `ref` a free choice, and the reason no "primary" is needed.
    #
    # Modulo 360°, because `closurephase!` folds each baseline phase into
    # (−180°, 180°] before summing and does not fold the sum. That is v1's
    # behaviour, kept; see the warning in `InterferometryObs`'s docstring.
    models = map((:A, :b, :c, Barycentre, Photocentre(:K))) do r
        obs = InterferometryObs(tab; targets=(:A, :b, :c), ref=r, band=:K, name="vis")
        sys = model_system(obs=[obs], flux_A=0.7, flux_b=1.0, flux_c=0.4)
        ctx, _, _ = ctx_for(sys, obs)
        [s.cps_model for s in Octofitter.simulate(obs, ctx)]
    end
    for m in models, i in eachindex(m)
        @test all(≈(0.0, atol=1e-8), mod.(m[i] .- models[1][i] .+ 180, 360) .- 180)
    end

    # …and the visibility is normalized, so a global rescaling of every band
    # flux changes nothing either. v1 could not express this: the host's flux
    # was the literal constant 1.
    scaled = map((1.0, 3.7, 0.02)) do k
        obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="vis")
        sys = model_system(obs=[obs], flux_A=1.0k, flux_b=0.3k)
        ctx, _, _ = ctx_for(sys, obs)
        Octofitter.ln_like(obs, ctx)
    end
    @test all(≈(scaled[1], rtol=1e-9), scaled)
end

@testset "squared visibilities" begin
    tab = synthetic_table(use_vis2=true)
    obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="vis")
    sys = model_system(obs=[obs])
    ctx, _, traj = ctx_for(sys, obs)
    ll = Octofitter.ln_like(obs, ctx)
    # v1's vis² branch read an uninitialized variable and could only throw, so
    # there is nothing to reproduce; check instead that it is the ordinary
    # Gaussian term on top of the closure-phase one.
    @test isfinite(ll)
    ll_cp_only = Octofitter.ln_like(
        InterferometryObs(synthetic_table(use_vis2=false); targets=(:A, :b), ref=:A,
                          band=:K, name="vis"), ctx)
    sim = Octofitter.simulate(obs, ctx)
    expect_v2 = sum(eachindex(EPOCHS)) do i
        r = tab.vis2_data[i] .- sim[i].vis2_model
        sum(@. -0.5 * (r^2 / tab.dvis2[i]^2 + log(2π * tab.dvis2[i]^2)))
    end
    @test ll - ll_cp_only ≈ expect_v2 rtol = 1e-10
end

@testset "epoch subsetting and simulation" begin
    tab = synthetic_table()
    obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="vis")
    sys = model_system(obs=[obs])
    ctx, _, _ = ctx_for(sys, obs)

    sub = Octofitter.likeobj_from_epoch_subset(obs, [2])
    @test sub.table.epoch == [EPOCHS[2]]
    @test sub.targets === obs.targets && sub.ref === obs.ref && sub.band === obs.band
    # Subsetting must not re-read the input files or rebuild the design matrix;
    # it rebuilds the struct around a sliced table.
    subsys = model_system(obs=[sub])
    subctx, _, _ = ctx_for(subsys, sub)
    per_epoch = map(1:3) do i
        o = Octofitter.likeobj_from_epoch_subset(obs, [i])
        s = model_system(obs=[o])
        c, _, _ = ctx_for(s, o)
        Octofitter.ln_like(o, c)
    end
    @test sum(per_epoch) ≈ Octofitter.ln_like(obs, ctx) rtol = 1e-12

    # Simulated data reproduce themselves: fit the generated closure phases and
    # the residuals are zero, so the likelihood is exactly its normalization.
    gen = Octofitter.generate_from_params(obs, ctx; add_noise=false)
    gensys = model_system(obs=[gen])
    genctx, _, _ = ctx_for(gensys, gen)
    n_cp = length(IDX_CPS1)
    Λ = length(tab.eff_wave[1])
    expected = -0.5 * length(EPOCHS) * Λ * n_cp * log(2π * 2.5^2)
    @test Octofitter.ln_like(gen, genctx) ≈ expected rtol = 1e-10

    # …and `add_noise` is honoured, which v1 accepted and ignored.
    noisy = Octofitter.generate_from_params(obs, ctx; add_noise=true)
    @test noisy.table.cps_data != gen.table.cps_data
end

@testset "ForwardDiff and the hot loop" begin
    tab = synthetic_table()
    obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="vis",
        variables=@variables begin
            σ_cp_jitter ~ Uniform(0.1, 10)
            platescale = 1.0
            northangle = 0.0
        end)
    sysvars = @variables begin
        plx = 100.0
    end
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.2
        flux_K = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 20mjup
        flux_K ~ Uniform(0, 1)
        a ~ Uniform(1, 4)
        e = 0.12
        i = 0.6
        ω = 1.0
        Ω = 0.4
        tp = 58800.0
    end)
    sys = Octofitter.System(name="s", bodies=[A, b], observations=[obs], variables=sysvars)
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θt = model.link(Octofitter.sample_priors(Random.Xoshiro(3), sys))
    val, g = model.∇ℓπcallback(θt)
    @test isfinite(val)
    @test all(isfinite, g)
    @test any(!iszero, g)
end

@testset "ForwardDiff through the kernel-phase path" begin
    # The kernel-phase branch takes a different route to its number: the
    # covariance is assembled into an owned matrix, factorized in place, and
    # wrapped in a `PDMat`, with a `try`/`catch` around the whole thing. None of
    # that may type-assert to `Float64` or the gradient silently dies — and the
    # `-Inf` fallback must be reachable without throwing out of a dual number.
    tab = synthetic_table(Λ=3, jitter=:kp_jit, kp_Cy=:kp_cy)
    obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="grav",
        kernel_phases=true,
        variables=@variables begin
            kp_jit ~ Uniform(0.1, 10)
            kp_cy = 0.04
        end)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.2
        flux_K = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 20mjup
        flux_K ~ Uniform(0, 1)
        a ~ Uniform(1, 4)
        e = 0.12
        i = 0.6
        ω = 1.0
        Ω = 0.4
        tp = 58800.0
    end)
    sys = Octofitter.System(name="s", bodies=[A, b], observations=[obs],
        variables=@variables begin
            plx = 100.0
        end)
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θt = model.link(Octofitter.sample_priors(Random.Xoshiro(5), sys))
    val, g = model.∇ℓπcallback(θt)
    @test isfinite(val)
    @test all(isfinite, g)
    @test any(!iszero, g)
end

@testset "platescale is a divisor now, and northangle rotates by −na" begin
    # v1's interferometry likelihood *multiplied* its offsets by `platescale`,
    # the reciprocal of the convention every other likelihood uses. The shared
    # `sky_offset` front-end uses the majority convention, so this pins the new
    # sense — it is a user-visible change, not an implementation detail.
    tab = synthetic_table()
    mk(ps, na) = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="vis",
        variables=@variables begin
            platescale = $ps
            northangle = $na
        end)
    base = mk(1.0, 0.0)
    ctx_base, _, _ = ctx_for(model_system(obs=[base]), base)
    front = OI._front_end(base, ctx_base, Float64)
    off0, _ = OI._epoch_offsets(ctx_base, 1, front)

    scaled = mk(2.0, 0.0)
    ctx_s, _, _ = ctx_for(model_system(obs=[scaled]), scaled)
    off_s, _ = OI._epoch_offsets(ctx_s, 1, OI._front_end(scaled, ctx_s, Float64))
    for j in eachindex(off0)
        @test off_s[j][1] ≈ off0[j][1] / 2 rtol = 1e-12
        @test off_s[j][2] ≈ off0[j][2] / 2 rtol = 1e-12
    end

    rotated = mk(1.0, 0.3)
    ctx_r, _, _ = ctx_for(model_system(obs=[rotated]), rotated)
    off_r, _ = OI._epoch_offsets(ctx_r, 1, OI._front_end(rotated, ctx_r, Float64))
    for j in eachindex(off0)
        # Position angle (N through E) rotated by −northangle: the separation is
        # preserved and atan(Δα✱, Δδ) decreases by `northangle`.
        @test hypot(off_r[j]...) ≈ hypot(off0[j]...) rtol = 1e-12
        iszero(hypot(off0[j]...)) && continue
        Δpa = atan(off_r[j][1], off_r[j][2]) - atan(off0[j][1], off0[j][2])
        @test mod(Δpa + π, 2π) - π ≈ -0.3 rtol = 1e-10
    end
end

@testset "no allocations in the closure-phase hot loop" begin
    tab = synthetic_table()
    obs = InterferometryObs(tab; targets=(:A, :b, :c), ref=:A, band=:K, name="vis")
    sys = model_system(obs=[obs], flux_c=0.1)
    ctx, _, _ = ctx_for(sys, obs)
    f(obs, ctx) = Octofitter.ln_like(obs, ctx)
    f(obs, ctx)
    # References resolve outside the epoch loop and every temporary comes from
    # `ctx.buf`; nothing here may reach the heap.
    @test (@allocated f(obs, ctx)) == 0
    @test Core.Compiler.return_type(
        Octofitter.ln_like, Tuple{typeof(obs),typeof(ctx)}) === Float64
end

@testset "the GRAVITY-WIDE preset" begin
    tab = synthetic_table(Λ=4, jitter=:kp_jit)
    obs = GRAVITYWideKPObs(tab; targets=(:A, :b), ref=:A, band=:K,
        variables=@variables begin
            kp_jit = 1.5
        end)
    @test obs isa InterferometryObs
    @test obs.phases isa OI.KernelPhases
    @test obs.fiber_pointing === Photocentre(:K)
    @test Octofitter.likelihoodname(obs) == "GRAVITY-WIDE"
    sys = model_system(obs=[obs], a_b=0.35)
    ctx, _, _ = ctx_for(sys, obs)
    @test isfinite(Octofitter.ln_like(obs, ctx))
    @test Core.Compiler.return_type(
        Octofitter.ln_like, Tuple{typeof(obs),typeof(ctx)}) === Float64
end

# ---------------------------------------------------------------------------
# Kernel phases
# ---------------------------------------------------------------------------

@testset "kernel-phase basis" begin
    Λ = 4
    Tλ, P₁ = OI.kernel_phase_basis(OI.GRAVITY_T3_DESIGN, Λ)
    @test size(Tλ) == (4Λ, 6Λ)
    @test size(P₁) == (3Λ, 4Λ)
    # The four closure phases over six baselines span only three independent
    # directions per channel, and `P₁`'s rows are a basis of exactly that span.
    # (Its rows are orthogonal but not unit — v1 normalizes by the Cholesky
    # factor's row norms, which leaves ‖·‖² = 4/3 — so this checks the span,
    # not orthonormality.)
    @test rank(Tλ) == 3Λ
    @test rank(P₁) == 3Λ
    @test P₁ * P₁' ≈ (4 / 3) * I(3Λ)
    # …and that span is the range of Tλ Tλ', so stacking adds no rank.
    @test rank([Tλ * Tλ' P₁']) == 3Λ
end

@testset "the block Cholesky agrees with a dense MvNormal" begin
    # The whole point of the block factorization is that it is *the same
    # answer* as factorizing the full covariance. v1 sized its blocks by
    # `length(table.eff_wave)` — the number of exposures in the table, not the
    # number of spectral channels in this one — so the `Cholesky` it handed to
    # `PDMat` was not a factorization of Σ except by coincidence. This is the
    # check that pins it.
    Λ = 5
    n_kp = 3Λ
    σ = collect(range(1.0, 2.0, length=n_kp))
    C = zeros(n_kp, n_kp)
    for k in 1:3
        blk = (k-1)*Λ+1:k*Λ
        C[blk, blk] .= 0.3
        for i in blk
            C[i, i] = 1.0
        end
    end
    Σ = OI._kp_covariance!(zeros(n_kp, n_kp), C, σ, 0.7)
    resid = collect(range(-2.0, 3.0, length=n_kp))
    dense = logpdf(MvNormal(Symmetric(copy(Σ))), resid)
    @test OI._kp_logpdf(copy(Σ), resid, 3) ≈ dense rtol = 1e-10
    # A non-positive-definite draw is a `-Inf`, not a thrown exception.
    bad = copy(Σ); bad[1, 1] = -1.0
    @test OI._kp_logpdf(bad, resid, 3) == -Inf
end

@testset "kernel-phase likelihood" begin
    tab = synthetic_table(Λ=4, jitter=:kp_jit, kp_Cy=:kp_cy)
    obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="grav",
        kernel_phases=true,
        variables=@variables begin
            kp_jit = 1.5
            kp_cy = 0.04
        end)
    @test obs.phases isa OI.KernelPhases
    sys = model_system(obs=[obs])
    ctx, _, _ = ctx_for(sys, obs)
    ll = Octofitter.ln_like(obs, ctx)
    @test isfinite(ll)

    # Reconstruct the same number by hand from the residuals, the basis and a
    # plain dense MvNormal — no block structure assumed anywhere.
    sim = Octofitter.simulate(obs, ctx)
    expect = 0.0
    for i in eachindex(EPOCHS)
        Λ = length(tab.eff_wave[i])
        n_T3 = size(tab.cps_data[i], 1)
        r = zeros(n_T3 * Λ)
        for i_T3 in 1:n_T3, i_wave in 1:Λ
            r[(i_T3-1)*Λ+i_wave] = tab.cps_data[i][i_T3, i_wave] - sim[i].cps_model[i_T3, i_wave]
        end
        P₁ = obs.table.P₁[i]
        σ_kp = obs.table.σ_kp[i]
        C = OI.CKP(obs.table[i], 0.04)
        Σ = OI._kp_covariance!(zeros(size(C)), C, σ_kp, 1.5)
        expect += logpdf(MvNormal(Symmetric(Σ)), P₁ * r)
    end
    @test ll ≈ expect rtol = 1e-10

    # `correlated=false` drops the spectral correlation but keeps the basis.
    obs0 = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="grav",
        kernel_phases=true, kp_correlation=false,
        variables=@variables begin
            kp_jit = 1.5
            kp_cy = 0.04
        end)
    sys0 = model_system(obs=[obs0])
    ctx0, _, _ = ctx_for(sys0, obs0)
    @test isfinite(Octofitter.ln_like(obs0, ctx0))
    @test Octofitter.ln_like(obs0, ctx0) != ll

    # Kernel phases need GRAVITY's geometry; say so rather than mis-blocking.
    weird = Table([(; exposure(EPOCHS[1])...,
        cps_data=fill(3.0, 3, 3), dcps=fill(2.5, 3, 3),
        index_cps1=[1, 1, 2], index_cps2=[4, 5, 6], index_cps3=[2, 3, 3])])
    @test_throws r"6 baselines forming 4 triangles" InterferometryObs(
        weird; targets=(:A, :b), ref=:A, band=:K, name="x", kernel_phases=true)
end

# ---------------------------------------------------------------------------
# Fibre coupling
# ---------------------------------------------------------------------------

@testset "fibre coupling is a reference, not a two-body formula" begin
    tab = synthetic_table(Λ=3)
    mk(pointing) = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K,
        name="grav", fiber_coupling=true, fiber_pointing=pointing)

    # Pointing at the host: the host sits on the fibre axis and keeps its full
    # throughput, while the companion is attenuated by its own separation. v1
    # had this exactly backwards — it applied the *host's* throughput to the
    # companion and left the host at 1.0.
    obs = mk(:A)
    @test Octofitter.refspecs(obs) ==
          (Octofitter.refspec(:A), Octofitter.refspec(:b), Octofitter.refspec(:A),
           Octofitter.refspec(:A))
    @test occursin("fiber at A", sprint(show, MIME"text/plain"(), obs))
    # 0.35 AU at 100 mas parallax is ~35 mas, inside the coupling grid and out
    # far enough that the throughput is well below 1.
    sys = model_system(obs=[obs], a_b=0.35)
    ctx, posys, traj = ctx_for(sys, obs)

    front = OI._front_end(obs, ctx, Float64)
    offsets, fiberoff = OI._epoch_offsets(ctx, 1, front)
    w = OI._effective_weights(obs, front.fluxes, offsets, fiberoff, tab.eff_wave[1][1])
    @test w[1] == front.fluxes[1] * obs.fiber_coupling(0.0, tab.eff_wave[1][1])
    sep_b = hypot(offsets[2]...)
    @test w[2] ≈ front.fluxes[2] * obs.fiber_coupling(sep_b, tab.eff_wave[1][1])
    # Real attenuation, not a rounding effect: this companion is tens of mas out.
    @test 0 < w[2] / front.fluxes[2] < 0.9

    # The default pointing is the band photocentre, which for three bodies is
    # defined where v1's `f·ρ/(1+f)` was not.
    obs3 = InterferometryObs(tab; targets=(:A, :b, :c), ref=:A, band=:K,
                             name="grav", fiber_coupling=true)
    @test obs3.fiber_pointing === Photocentre(:K)
    sys3 = model_system(obs=[obs3], flux_c=0.2, a_b=0.35)
    ctx3, _, _ = ctx_for(sys3, obs3)
    @test isfinite(Octofitter.ln_like(obs3, ctx3))

    # Reweighting the sources changes the visibility, so it must change the
    # answer: the coupling is not a no-op the way v1's `f·ρ/(1+f)` nearly was
    # for a faint companion.
    plain = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="grav")
    @test Octofitter.ln_like(obs, ctx) != Octofitter.ln_like(plain, ctx)
end

# ---------------------------------------------------------------------------
# OI-FITS loading
# ---------------------------------------------------------------------------

AMI_DIR = joinpath(@__DIR__, "..", "..", "examples", "AMI_data")

if isdir(AMI_DIR)
    @testset "OI-FITS input" begin
        data = Table([
            (; filename=joinpath(AMI_DIR, "Sim_data_2023_1_.oifits"),
               epoch=mjd("2023-06-01"), use_vis2=false),
            (; filename=joinpath(AMI_DIR, "Sim_data_2024_1_.oifits"),
               epoch=mjd("2024-06-01"), use_vis2=false),
        ])
        obs = InterferometryObs(data; targets=(:A, :b), ref=:A, band=:K,
                                name="NIRISS-AMI")
        @test length(obs.table.epoch) == 2
        @test size(obs.table.u[1], 1) == size(obs.table.v[1], 1)
        @test size(obs.table.cps_data[1]) == size(obs.table.dcps[1])
        @test all(>(0), obs.table.eff_wave[1])

        sys = model_system(obs=[obs])
        ctx, _, traj = ctx_for(sys, obs)
        @test isfinite(Octofitter.ln_like(obs, ctx))
        # Real data, real u/v: still bit-for-bit v1.
        @test Octofitter.ln_like(obs, ctx) ===
              v1_ln_like(obs.table, traj, ctx.epoch_index, (0.3,))
    end
end

end

# ---------------------------------------------------------------------------
# Plotting protocol
#
# One exposure is an (n_triangles × n_channels) block, so the channel flattens
# it and repeats the epoch: every measured closure phase is one point.
# ---------------------------------------------------------------------------
@testset "plotting protocol" begin
    tab = synthetic_table()
    obs = InterferometryObs(tab; targets=(:A, :b), ref=:A, band=:K, name="vis",
        variables=Octofitter.@variables begin
            σ_cp_jitter = 0.2
        end)
    sys = model_system(obs=[obs])
    @test isempty(Octofitter.unplottable_observations(sys))
    chs = Octofitter.plotchannels(obs)
    @test :closure_phase in [ch.name for ch in chs]
    # No PlanetOrbits observable produces a closure phase.
    @test all(ch -> ch.query === nothing, chs)

    ctx, _, _ = ctx_for(sys, obs)
    r = Octofitter.residuals(obs, ctx)
    sims = Octofitter.simulate(obs, ctx)
    ntot = sum(length(tab.cps_data[i]) for i in eachindex(tab.epoch))
    @test length(r.closure_phase.data) == ntot
    @test length(unique(r.closure_phase.epoch)) == length(tab.epoch)
    @test r.closure_phase.resid ≈ r.closure_phase.data .- r.closure_phase.model
    @test r.closure_phase.model[1:length(sims[1].cps_model)] ≈ vec(sims[1].cps_model)
    # σ_eff carries the fitted closure-phase jitter, in degrees.
    @test r.closure_phase.σ_eff ≈ hypot.(r.closure_phase.σ, 0.2)
    # Squared visibilities appear only when the table asks for them.
    if any(tab.use_vis2 .=== true)
        @test haskey(r, :vis2)
    else
        @test !haskey(r, :vis2)
    end
end
