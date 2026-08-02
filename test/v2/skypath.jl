# The shared sky-path / LSQ helpers migrated out of `src/legacy/`.
#
# The five-parameter refit is portable linear algebra, so the gate is
# agreement with the legacy implementation on random inputs — the functions
# below are transcribed verbatim from `src/legacy/gaia-utils-fitting.jl`,
# with Bumper's `@no_escape`/`@alloc` replaced by ordinary arrays (an
# allocation strategy, not arithmetic) and the `table.cosϕ`/`table.sinϕ`
# column reads replaced by the vectors themselves.

const JY = Octofitter.julian_year

function legacy_prepare_A_5param(table, ref_ra, ref_dec)
    n_obs = size(table, 1)
    A = zeros(n_obs, 5)
    for i in 1:n_obs
        A[i, 1] = table.cosϕ[i]
        A[i, 2] = table.sinϕ[i]
        A[i, 3] = -table.parallaxFactorAlongScan[i]
        A[i, 4] = table.cosϕ[i] * (table.epoch[i] - ref_ra) / JY
        A[i, 5] = table.sinϕ[i] * (table.epoch[i] - ref_dec) / JY
    end
    return A
end

function legacy_prepare_A_4param(table, ref_ra, ref_dec)
    n_obs = size(table, 1)
    A = zeros(n_obs, 4)
    for i in 1:n_obs
        A[i, 1] = table.cosϕ[i]
        A[i, 2] = table.sinϕ[i]
        A[i, 3] = table.cosϕ[i] * (table.epoch[i] - ref_ra) / JY
        A[i, 4] = table.sinϕ[i] * (table.epoch[i] - ref_dec) / JY
    end
    return A
end

function legacy_fit_5param_prepared(A_prepared, table, Δα_mas, Δδ_mas,
                                    residuals=0.0, σ_formal=0.0; include_chi2=Val(false))
    n_obs = size(A_prepared, 1)
    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))
    b = zeros(T, n_obs)
    A_dense = zeros(T, n_obs, size(A_prepared, 2))
    @. b = Δα_mas * table.cosϕ + Δδ_mas * table.sinϕ + residuals
    if σ_formal isa Number
        @. A_dense = A_prepared
    else
        @. A_dense = A_prepared / σ_formal
        @. b = b / σ_formal
    end
    x = A_dense \ b
    parameters = @SVector [x[1], x[2], x[4], x[5], x[3]]
    include_chi2 == Val(true) || return (; parameters)
    σ_formal == 0 && error("Asked for `include_chi2=true` but `σ_formal==0`")
    resid = b .- A_dense * x
    chi_squared_astro = σ_formal isa Number ?
                        dot(resid, resid) / (σ_formal * σ_formal) : dot(resid, resid)
    dof = n_obs - 5
    return (; parameters, chi_squared_astro, chi2_reduced=chi_squared_astro / dof, dof)
end

function legacy_fit_5param_pinv(pinv_A, table, Δα_mas, Δδ_mas, residuals=0.0)
    n_obs = size(pinv_A, 2)
    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))
    b = zeros(T, n_obs)
    @. b = Δα_mas * table.cosϕ + Δδ_mas * table.sinϕ + residuals
    x_buf = pinv_A * b
    return (; parameters=@SVector [x_buf[1], x_buf[2], x_buf[4], x_buf[5], x_buf[3]])
end

function legacy_fit_4param_prepared(A_factored, table, Δα_mas, Δδ_mas, σ_formal=0.0)
    n_obs = size(table, 1)
    T = promote_type(eltype(Δα_mas), eltype(Δδ_mas))
    b = zeros(T, n_obs)
    for i in 1:n_obs
        b[i] = Δα_mas[i] * table.cosϕ[i] + Δδ_mas[i] * table.sinϕ[i]
    end
    if σ_formal != 0.0
        @. b *= 1 / σ_formal
    end
    x = A_factored \ b
    return (; parameters=@SVector [x[1], x[2], x[3], x[4]])
end

@testset "five-parameter refit matches the legacy implementation" begin
    rng = Random.Xoshiro(20260802)
    for trial in 1:8
        n = 40 + trial
        ϕ = 2π .* rand(rng, n)
        tab = Table(cosϕ=cos.(ϕ), sinϕ=sin.(ϕ),
            epoch=sort(56000 .+ 1200 .* rand(rng, n)),
            parallaxFactorAlongScan=0.8 .* randn(rng, n))
        ref_ra, ref_dec = 57388.5, 57390.0
        A_new = Octofitter.prepare_A_5param(tab.cosϕ, tab.sinϕ, tab.epoch,
            tab.parallaxFactorAlongScan, ref_ra, ref_dec)
        A_old = legacy_prepare_A_5param(tab, ref_ra, ref_dec)
        @test A_new == A_old

        A4_new = Octofitter.prepare_A_4param(tab.cosϕ, tab.sinϕ, tab.epoch, ref_ra, ref_dec)
        @test A4_new == legacy_prepare_A_4param(tab, ref_ra, ref_dec)

        Δα = 3 .* randn(rng, n)
        Δδ = 3 .* randn(rng, n)
        res = 0.5 .* randn(rng, n)
        σvec = 0.4 .+ rand(rng, n)
        σscal = 0.7

        # unweighted, scalar σ, and per-epoch σ; with and without residuals
        for (r, σ) in ((0.0, 0.0), (res, 0.0), (res, σscal), (res, σvec), (0.0, σvec))
            new = Octofitter.fit_5param_prepared(A_new, tab.cosϕ, tab.sinϕ, Δα, Δδ, r, σ)
            old = legacy_fit_5param_prepared(A_old, tab, Δα, Δδ, r, σ)
            @test new.parameters ≈ old.parameters rtol = 1e-10
        end

        # χ² channel
        newc = Octofitter.fit_5param_prepared(A_new, tab.cosϕ, tab.sinϕ, Δα, Δδ, res, σvec;
            include_chi2=Val(true))
        oldc = legacy_fit_5param_prepared(A_old, tab, Δα, Δδ, res, σvec; include_chi2=Val(true))
        @test newc.parameters ≈ oldc.parameters rtol = 1e-10
        @test newc.chi_squared_astro ≈ oldc.chi_squared_astro rtol = 1e-10
        @test newc.chi2_reduced ≈ oldc.chi2_reduced rtol = 1e-10
        @test newc.dof == oldc.dof == n - 5
        @test_throws r"σ_formal==0" Octofitter.fit_5param_prepared(
            A_new, tab.cosϕ, tab.sinϕ, Δα, Δδ, res, 0.0; include_chi2=Val(true))

        # the cached pseudo-inverse route is the same weighted answer
        Q_new = Octofitter.prepare_pinv_5param(A_new, σvec)
        Q_old = pinv(A_old ./ σvec) ./ permutedims(σvec)
        @test Q_new ≈ Q_old rtol = 1e-10
        pn = Octofitter.fit_5param_pinv(Q_new, tab.cosϕ, tab.sinϕ, Δα, Δδ, res)
        po = legacy_fit_5param_pinv(Q_old, tab, Δα, Δδ, res)
        @test pn.parameters ≈ po.parameters rtol = 1e-10
        # …and it agrees with the factorizing weighted solve it replaces
        @test pn.parameters ≈ Octofitter.fit_5param_prepared(
            A_new, tab.cosϕ, tab.sinϕ, Δα, Δδ, res, σvec).parameters rtol = 1e-8

        # four-parameter counterpart
        f4n = Octofitter.fit_4param_prepared(A4_new, tab.cosϕ, tab.sinϕ, Δα, Δδ, σscal)
        f4o = legacy_fit_4param_prepared(A4_new, tab, Δα, Δδ, σscal)
        @test f4n.parameters ≈ f4o.parameters rtol = 1e-10
    end

    # Rejected scans must be given σ = Inf explicitly rather than sneaking
    # back in at weight 1/|σ|.
    A = Octofitter.prepare_A_5param(cos.(1:10), sin.(1:10), collect(56000.0:10:56090.0),
        zeros(10), 57388.5, 57388.5)
    @test_throws r"every σ must be positive" Octofitter.prepare_pinv_5param(A, [-1.0; ones(9)])
    @test Octofitter.prepare_pinv_5param(A, [Inf; ones(9)]) isa Matrix

    # The refit is differentiable: it sits inside a likelihood.
    #
    # NB the parallax column must not be collinear with the scan-angle
    # columns. A "plausible" fixture with `parallax_factor = c·sinϕ` makes A
    # rank-4 with cond ~1e16, and then AD reports gradients of order 1e15
    # while finite differences report O(1) — which reads as an AD bug in the
    # solver and is a degenerate design matrix. (The legacy implementation
    # produces the identical garbage, which is how this was diagnosed.)
    let n = 30
        ϕ = range(0, 2π, length=n)
        cosϕ = cos.(ϕ); sinϕ = sin.(ϕ)
        ep = collect(range(56000.0, 57200.0, length=n))
        plxfac = 0.5 .* sin.(2π .* (ep .- 56000.0) ./ 365.25)   # annual, not per-scan
        A = Octofitter.prepare_A_5param(cosϕ, sinϕ, ep, plxfac, 56600.0, 56600.0)
        @test rank(A) == 5
        f(p) = sum(Octofitter.fit_5param_prepared(A, cosϕ, sinϕ,
            p[1] .* cosϕ .+ p[2], p[3] .* sinϕ, 0.0, 0.0).parameters)
        g = Octofitter.ForwardDiff.gradient(f, [2.0, 1.0, -3.0])
        gfd = FiniteDiff.finite_difference_gradient(f, [2.0, 1.0, -3.0])
        @test norm(g .- gfd) / norm(g) < 1e-6
    end
end

@testset "accumulate_offsets! is the refs-based sky-path accumulator" begin
    # Build a real context so the helper is tested against the same thing a
    # likelihood sees, not a mock.
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
        flux_G = 1.0
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 30mjup
        flux_G = 0.3
        a = 4.0; e = 0.15; i = 0.8; ω = 0.4; Ω = 1.7; tp = 55500.0
    end)
    ep = collect(range(56000.0, 58400.0, length=11))
    obs = RelAstromObs((epoch=ep, ra=zeros(length(ep)), dec=zeros(length(ep)),
        σ_ra=ones(length(ep)), σ_dec=ones(length(ep))); target=:b, ref=:A, name="astrom")
    sys = Octofitter.System(name="acc", bodies=[A, b], observations=[obs],
        variables=@variables begin plx = 30.0 end)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    eps, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, eps)
    # this observation declares no variables of its own
    ctx = Octofitter.ObsContext(nt, (;), posys, traj, maps[obs])

    tgt = Octofitter.ref(ctx, Photocentre(:G, (:A, :b)))
    rf = Octofitter.ref(ctx, Barycentre)
    expect_ra = [raoff(Octofitter.solutionat(ctx, i), tgt, rf) for i in eachindex(ep)]
    expect_dec = [decoff(Octofitter.solutionat(ctx, i), tgt, rf) for i in eachindex(ep)]

    # whole table
    ra = zeros(length(ep)); dec = zeros(length(ep))
    Octofitter.accumulate_offsets!(ra, dec, ctx, tgt, rf)
    @test ra == expect_ra
    @test dec == expect_dec
    # it *accumulates*, so a second pass doubles it
    Octofitter.accumulate_offsets!(ra, dec, ctx, tgt, rf)
    @test ra ≈ 2 .* expect_ra

    # a row range, with the buffers indexed by position within the range —
    # the shape a likelihood whose table concatenates instrument channels needs
    rows = 4:9
    ra2 = zeros(length(rows)); dec2 = zeros(length(rows))
    Octofitter.accumulate_offsets!(ra2, dec2, ctx, tgt, rf, rows)
    @test ra2 == expect_ra[rows]
    @test dec2 == expect_dec[rows]
    @test_throws r"indexed by position within" Octofitter.accumulate_offsets!(
        zeros(length(ep)), zeros(length(ep)), ctx, tgt, rf, rows)

    # a hand-built WeightedPoint goes in the same slot as a resolved spec —
    # this is the tier-2/tier-3 entry point
    f = PlanetOrbits.fluxes(posys, :G)
    wp = PlanetOrbits.photocentre(f .* SVector(1.0, 1.0))
    ra3 = zeros(length(ep)); dec3 = zeros(length(ep))
    Octofitter.accumulate_offsets!(ra3, dec3, ctx, wp, rf)
    @test ra3 ≈ expect_ra
    # …and a bare body ref works too (single-member degradation)
    ra4 = zeros(length(ep)); dec4 = zeros(length(ep))
    Octofitter.accumulate_offsets!(ra4, dec4, ctx, Octofitter.ref(ctx, Octofitter.refspec(:A)), rf)
    @test ra4 == [raoff(Octofitter.solutionat(ctx, i), :A, rf) for i in eachindex(ep)]
end
