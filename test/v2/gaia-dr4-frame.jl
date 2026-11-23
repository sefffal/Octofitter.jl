# `GaiaDR4AstromObs` on a system with an absolute frame: the reference point
# is the barycentre's propagated direction (PlanetOrbits' `frame_ra`/
# `frame_dec`), projected about the data's tangent point, with the parallax at
# the epoch. Checked against independent closed forms rather than against the
# implementation's own arithmetic:
#
#   - Barnard-like kinematics: switching the radial velocity on must add
#     exactly the first-order perspective term to every abscissa;
#   - a slow, distant star with rv = 0: the absolute-frame prediction must
#     reproduce the linear parallax-only model to the gnomonic terms;
#   - the build-time contract (no observation-level reference point under an
#     absolute frame, a tangent point from the data, `ref=Barycentre`);
#   - parallax-only models are untouched by the tangent-point columns.
#
# Offline and self-contained.

using Octofitter
using Octofitter.TypedTables: Table
using PlanetOrbits
using Distributions, LinearAlgebra, Random, Test
using FiniteDiff

# A DR4-like epoch set: 140 transits over the 66-month DR4 window
# (2014-07-25 .. 2020-01-20), irregularly spaced, with scan angles that wrap
# the circle and parallax factors of the right size. The values need not be a
# real scan law: the closed forms below take the same ψ and f_al as data.
const DR4F_N = 140
const DR4F_EPOCHS = [56863.0 + (k - 1) * 14.3 + 3.0 * sin(1.7k) for k in 1:DR4F_N]
const DR4F_PSI = [mod(rad2deg(0.37k + 1.1), 360.0) for k in 1:DR4F_N]   # degrees
const DR4F_FAL = [0.9 * sin(2π * (DR4F_EPOCHS[k] - 56863.0) / 365.25 + 0.4) for k in 1:DR4F_N]
# The DR4 reference epoch, J2017.5. Also the frame's `ref_epoch` throughout.
const DR4F_REF = 57936.375

# Barnard's star, roughly (Gaia DR3; rv from the literature).
const DR4F_BARNARD = (ra=269.4485, dec=4.7394, plx=546.98,
                 pmra=-801.55, pmdec=10362.39, rv=-110.1e3)   # rv in m/s
const DR4F_AV = 4.740470446                                         # km/s per (au/yr)

function dr4f_table(; ra0=nothing, dec0=nothing)
    cols = (epoch=copy(DR4F_EPOCHS), scan_pos_angle=copy(DR4F_PSI),
            parallax_factor_al=copy(DR4F_FAL), centroid_pos_al=zeros(DR4F_N),
            centroid_pos_error_al=fill(0.08, DR4F_N), outlier_flag=zeros(Int, DR4F_N))
    ra0 === nothing && return Table(cols)
    return Table(merge(cols, (ra0=fill(Float64(ra0), DR4F_N), dec0=fill(Float64(dec0), DR4F_N))))
end

# A single star, so the photocentre *is* the barycentre and the prediction is
# the reference-point model alone.
dr4f_star() = Octofitter.Body(name="A", variables=@variables begin
    mass = 0.16
    flux = 1.0
end)

function dr4f_absolute_system(obs, k; lighttime=:off, name="absframe")
    return Octofitter.System(name=name, bodies=[dr4f_star()], observations=[obs],
        barycentric_lighttime=lighttime, observing_geometry=:off,
        variables=@variables begin
            plx = $(k.plx)
            ra = $(k.ra)
            dec = $(k.dec)
            pmra = $(k.pmra)
            pmdec = $(k.pmdec)
            rv = $(k.rv)
            ref_epoch = $DR4F_REF
        end)
end

function dr4f_parallax_system(obs, plx; name="plxonly")
    return Octofitter.System(name=name, bodies=[dr4f_star()], observations=[obs],
        barycentric_lighttime=:off, observing_geometry=:off,
        variables=@variables begin
            plx = $plx
        end)
end

# The modelled abscissae of a system's one observation, in table order, solved
# with the system's own solve configuration.
function dr4f_simulate(sys)
    obs = only(sys.observations)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    ep, maps = Octofitter.epoch_plan(sys)
    traj = orbitsolve(posys, ep; Octofitter._solvekw(sys)...)
    key = Symbol(Octofitter.likelihoodname(obs))
    # An observation with no variables has no entry in θ.
    obsns = hasproperty(nt, :observations) ? nt.observations : (;)
    θ_obs = hasproperty(obsns, key) ? getproperty(obsns, key) : (;)
    ctx = Octofitter.ObsContext(nt, θ_obs, posys, traj, maps[obs])
    return Octofitter.simulate(obs, ctx)
end

dr4f_lnlike(sys) = Octofitter.make_ln_like(sys)(sys, Octofitter.make_arr2nt(sys)(Float64[]))

@testset "Barnard-like perspective curvature reaches the abscissae" begin
    # Tangent point at the frame's own position at `ref_epoch`, as DR4's
    # (ra0, dec0) is (near) the source's catalogue position.
    obs = GaiaDR4AstromObs(dr4f_table(ra0=DR4F_BARNARD.ra, dec0=DR4F_BARNARD.dec);
        target=Photocentre, ref=Barycentre, name="GaiaDR4")
    w_rv = dr4f_simulate(dr4f_absolute_system(obs, DR4F_BARNARD)).along_scan
    w_0 = dr4f_simulate(dr4f_absolute_system(obs, merge(DR4F_BARNARD, (; rv=0.0)))).along_scan
    Δw = w_rv .- w_0

    # Independent first-order formula: with the barycentre receding/approaching
    # at v_r the local-plane motion is μ t / (1 + μ_r t) and the parallax
    # ϖ / (1 + μ_r t), with μ_r = ϖ[rad] v_r / A_v [1/yr]. To first order in
    # μ_r t,
    #     Δw = −μ_r t (ϖ f_al + (μα* sin ψ + μδ cos ψ) t).
    μ_r = deg2rad(DR4F_BARNARD.plx / 3.6e6) * (DR4F_BARNARD.rv / 1e3) / DR4F_AV
    t = (DR4F_EPOCHS .- DR4F_REF) ./ PlanetOrbits.year2day_julian
    s, c = sind.(DR4F_PSI), cosd.(DR4F_PSI)
    Δw_1 = @. -μ_r * t * (DR4F_BARNARD.plx * DR4F_FAL + (DR4F_BARNARD.pmra * s + DR4F_BARNARD.pmdec * c) * t)

    # The effect is several mas — 50–100 times the per-transit error of a
    # bright star — so it is not something a linear reference point can be
    # excused for omitting.
    @test maximum(abs, Δw) > 3.0
    # Agreement to second order: the next term in μ_r t is μ μ_r² t³, which is
    # ≲ 1 µas over this window (and the gnomonic and parallax-curvature terms
    # are smaller still).
    resid1 = maximum(abs, Δw .- Δw_1)
    Δw_2 = @. Δw_1 + μ_r^2 * t^3 * (DR4F_BARNARD.pmra * s + DR4F_BARNARD.pmdec * c)
    resid2 = maximum(abs, Δw .- Δw_2)
    @info "Barnard-like perspective curvature" max_abs_Δw = maximum(abs, Δw) first_order_residual = resid1 second_order_residual = resid2
    @test resid1 < 2e-3
    @test resid1 < 1e-3 * maximum(abs, Δw)
    # …and the residual is the expected second-order term, not noise: adding
    # it takes the residual down by more than an order of magnitude.
    @test resid2 < resid1 / 10

    # The parallax the abscissae use is the epoch's: 1000/|r0 + v t| for the
    # light-time-free worldline, computed here from the catalogue values
    # independently of PlanetOrbits.
    sys = dr4f_absolute_system(obs, DR4F_BARNARD)
    nt = Octofitter.make_arr2nt(sys)(Float64[])
    posys = Octofitter.make_ln_like(sys).build(nt)
    traj = orbitsolve(posys, DR4F_EPOCHS; Octofitter._solvekw(sys)...)
    ϖ_t = [Octofitter._dr4_frame_plx(traj[k]) for k in eachindex(DR4F_EPOCHS)]
    au_per_pc = 206264.80624709636
    d0 = 1000 / DR4F_BARNARD.plx                                          # pc
    vr = DR4F_BARNARD.rv / 1e3 / DR4F_AV / au_per_pc                          # pc / yr
    vt = hypot(DR4F_BARNARD.pmra, DR4F_BARNARD.pmdec) / 1000 / au_per_pc * d0  # pc / yr
    ϖ_ref = @. 1000 / sqrt((d0 + vr * t)^2 + (vt * t)^2)
    @test maximum(abs, ϖ_t .- ϖ_ref) < 1e-8
    @test maximum(ϖ_t) - minimum(ϖ_t) > 0.15                 # ~0.18 mas across DR4
    @info "Barnard-like parallax change over the window" Δϖ = maximum(ϖ_t) - minimum(ϖ_t)

    # With the light-time solve on, the light-time terms change the prediction
    # by a small, finite amount on top (second order for these kinematics).
    w_lt = dr4f_simulate(dr4f_absolute_system(obs, DR4F_BARNARD; lighttime=:on)).along_scan
    @test all(isfinite, w_lt)
    @info "Barnard-like light-time effect on the abscissae" max_abs = maximum(abs, w_lt .- w_rv)
end

@testset "rv = 0, slow distant star: absolute frame ≡ linear parallax-only model" begin
    k = (ra=150.0, dec=-30.0, plx=10.0, pmra=20.0, pmdec=-15.0, rv=0.0)
    function linear_model(Δα_mas, Δδ_mas)
        obs = GaiaDR4AstromObs(dr4f_table(); name="GaiaDR4",
            variables=@variables begin
                ra_offset_mas = $Δα_mas
                dec_offset_mas = $Δδ_mas
                pmra = $(k.pmra)
                pmdec = $(k.pmdec)
                ref_epoch = $DR4F_REF
            end)
        return dr4f_simulate(dr4f_parallax_system(obs, k.plx)).along_scan
    end
    absolute_model(ra0, dec0) = dr4f_simulate(dr4f_absolute_system(
        GaiaDR4AstromObs(dr4f_table(; ra0, dec0); name="GaiaDR4"), k)).along_scan

    # (a) Tangent point at the frame's reference position. With rv = 0 a
    # uniformly moving star's gnomonic track about its own reference direction
    # is *exactly* μ t (central projection maps the straight 3D worldline to a
    # straight line, and with no radial motion the projection's denominator is
    # constant), and its parallax changes only at O((μ t)²) ≈ 10⁻¹³ mas. What
    # is left is Float64 rounding of `frame_ra` in degrees, ~10⁻⁷ mas.
    w_lin = linear_model(0.0, 0.0)
    dw_a = maximum(abs, absolute_model(k.ra, k.dec) .- w_lin)
    @test dw_a < 1e-6
    @test maximum(abs, w_lin) > 50                       # not a vacuous comparison

    # (b) A tangent point a few hundred mas off, as DR4's (ra0, dec0) is the
    # source's catalogue position rather than the barycentre's. The offset is
    # the frame's exact gnomonic position about the tangent point, computed
    # here directly with unit vectors. What remains is the rotation between the frame's
    # local axes and the tangent point's, Δα sin δ ≈ 7×10⁻⁷ rad, acting on
    # the ≲ 75 mas of proper motion over the window: ≲ 5×10⁻⁵ mas.
    ra0 = k.ra - 0.25 / 3600 / cosd(k.dec)          # frame is 250 mas east of it
    dec0 = k.dec + 0.40 / 3600                      # …and 400 mas south
    u = [cosd(k.dec) * cosd(k.ra), cosd(k.dec) * sind(k.ra), sind(k.dec)]
    e_α0 = [-sind(ra0), cosd(ra0), 0.0]
    e_δ0 = [-sind(dec0) * cosd(ra0), -sind(dec0) * sind(ra0), cosd(dec0)]
    e_r0 = [cosd(dec0) * cosd(ra0), cosd(dec0) * sind(ra0), sind(dec0)]
    rad2mas = rad2deg(1.0) * 3.6e6
    ξ0 = dot(u, e_α0) / dot(u, e_r0) * rad2mas
    η0 = dot(u, e_δ0) / dot(u, e_r0) * rad2mas
    dw_b = maximum(abs, absolute_model(ra0, dec0) .- linear_model(ξ0, η0))
    @test dw_b < 1e-4
    # The small-angle offsets (Δα cos δ0, Δδ) would instead leave the O(θ²)
    # gnomonic term, ~ ξ0 · Δδ[rad] · tan δ ≈ 3×10⁻⁴ mas — visible, which is
    # why the reference point is projected rather than differenced.
    dw_small = maximum(abs, absolute_model(ra0, dec0) .-
                            linear_model((k.ra - ra0) * cosd(dec0) * 3.6e6,
                                         (k.dec - dec0) * 3.6e6))
    @test dw_small > dw_b
    @info "absolute vs linear model, rv = 0" tangent_at_frame = dw_a tangent_offset_exact = dw_b tangent_offset_small_angle = dw_small
end

@testset "absolute frame: build-time contract" begin
    tbl = dr4f_table(ra0=DR4F_BARNARD.ra, dec0=DR4F_BARNARD.dec)

    # Reference-point variables on the observation: refused, by name.
    for vars in (
        (@variables begin pmra = 1.0; pmdec = 2.0 end),
        (@variables begin ra_offset_mas ~ Normal(0, 10) end),
        (@variables begin ref_epoch = 57936.375 end),
    )
        obs = GaiaDR4AstromObs(tbl; name="GaiaDR4", variables=vars)
        err = try
            dr4f_absolute_system(obs, DR4F_BARNARD)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("absolute", sprint(showerror, err))
        @test occursin("declares", sprint(showerror, err))
    end
    # `astrometric_jitter` is not a reference-point variable and is fine.
    obs_j = GaiaDR4AstromObs(tbl; name="GaiaDR4",
        variables=@variables begin astrometric_jitter = 0.02 end)
    @test isfinite(dr4f_lnlike(dr4f_absolute_system(obs_j, DR4F_BARNARD)))

    # No tangent point: refused at build time, not at the first evaluation.
    obs_nt = GaiaDR4AstromObs(dr4f_table(); name="GaiaDR4")
    err = try
        dr4f_absolute_system(obs_nt, DR4F_BARNARD)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("tangent point", sprint(showerror, err))

    # …and the keyword supplies one for a table without the columns, giving
    # the same prediction as the columns.
    obs_kw = GaiaDR4AstromObs(dr4f_table(); name="GaiaDR4",
        tangent_point=(DR4F_BARNARD.ra, DR4F_BARNARD.dec))
    obs_col = GaiaDR4AstromObs(tbl; name="GaiaDR4")
    @test dr4f_simulate(dr4f_absolute_system(obs_kw, DR4F_BARNARD)).along_scan ==
          dr4f_simulate(dr4f_absolute_system(obs_col, DR4F_BARNARD)).along_scan

    # The frame describes the system barycentre, so `ref` must be it.
    obs_ref = GaiaDR4AstromObs(tbl; target=:A, ref=:A, name="GaiaDR4")
    @test_throws ErrorException dr4f_absolute_system(obs_ref, DR4F_BARNARD)

    # Tangent-point data problems are caught at construction.
    bad = Table(merge(Octofitter.Tables.columntable(tbl),
        (; ra0=[DR4F_BARNARD.ra + (k == 7) * 1e-6 for k in 1:DR4F_N])))
    @test_throws ErrorException GaiaDR4AstromObs(bad; name="GaiaDR4")       # not constant
    onlyra = Table(merge(Octofitter.Tables.columntable(dr4f_table()),
        (; ra0=fill(DR4F_BARNARD.ra, DR4F_N))))
    @test_throws ErrorException GaiaDR4AstromObs(onlyra; name="GaiaDR4")    # half a pair
    @test_throws ErrorException GaiaDR4AstromObs(tbl; name="GaiaDR4",
        tangent_point=(DR4F_BARNARD.ra + 1e-3, DR4F_BARNARD.dec))                     # disagrees
    @test_throws ErrorException GaiaDR4AstromObs(dr4f_table(); name="GaiaDR4",
        tangent_point=(10.0, 95.0))                                         # not a direction
end

@testset "absolute frame: subsets, simulation round trip, gradients" begin
    tbl = dr4f_table(ra0=DR4F_BARNARD.ra, dec0=DR4F_BARNARD.dec)
    obs = GaiaDR4AstromObs(tbl; name="GaiaDR4")

    # The tangent point survives the reconstruction paths.
    @test Octofitter.likeobj_from_epoch_subset(obs, 1:20).tangent_point ==
          (DR4F_BARNARD.ra, DR4F_BARNARD.dec)
    obs_kw = GaiaDR4AstromObs(dr4f_table(); name="GaiaDR4",
        tangent_point=(DR4F_BARNARD.ra, DR4F_BARNARD.dec))
    @test Octofitter.likeobj_from_epoch_subset(obs_kw, 1:20).tangent_point ==
          (DR4F_BARNARD.ra, DR4F_BARNARD.dec)

    # Noiseless simulation from the absolute-frame model is fitted exactly.
    sys = dr4f_absolute_system(obs, DR4F_BARNARD)
    θ = Octofitter.drawfrompriors(sys)
    clean = Octofitter.generate_from_params(sys, θ; add_noise=false)
    @test clean.observations[1].tangent_point == (DR4F_BARNARD.ra, DR4F_BARNARD.dec)
    expected = -DR4F_N * log(2π) / 2 - DR4F_N * log(0.08^2) / 2
    @test dr4f_lnlike(clean) ≈ expected rtol = 1e-10

    # A sampled absolute-frame model differentiates through the frame track.
    simtab = clean.observations[1].table
    obs_fit = GaiaDR4AstromObs(simtab; name="GaiaDR4")
    sysfree = Octofitter.System(name="absfree", bodies=[dr4f_star()],
        observations=[obs_fit], barycentric_lighttime=:off, observing_geometry=:off,
        variables=@variables begin
            plx ~ Normal(DR4F_BARNARD.plx, 1.0)
            ra = $(DR4F_BARNARD.ra)
            dec = $(DR4F_BARNARD.dec)
            pmra ~ Normal(DR4F_BARNARD.pmra, 1.0)
            pmdec ~ Normal(DR4F_BARNARD.pmdec, 1.0)
            rv ~ Normal(DR4F_BARNARD.rv, 1e3)
            ref_epoch = $DR4F_REF
        end)
    model = Octofitter.LogDensityModel(sysfree; verbosity=0)
    θt = model.link(Octofitter.sample_priors(Random.Xoshiro(5), sysfree))
    ℓ, g = model.∇ℓπcallback(θt)
    @test isfinite(ℓ)
    @test all(isfinite, g)
    gfd = FiniteDiff.finite_difference_gradient(model.ℓπcallback, θt)
    @test norm(g .- gfd) / norm(g) < 1e-5
end

@testset "parallax-only models ignore the tangent point" begin
    vars = @variables begin
        ra_offset_mas = 0.5
        dec_offset_mas = -0.3
        pmra = 5.3
        pmdec = -24.2
        ref_epoch = 57388.0
    end
    with = GaiaDR4AstromObs(dr4f_table(ra0=DR4F_BARNARD.ra, dec0=DR4F_BARNARD.dec);
        name="GaiaDR4", variables=vars)
    without = GaiaDR4AstromObs(dr4f_table(); name="GaiaDR4", variables=vars)
    @test dr4f_simulate(dr4f_parallax_system(with, 20.0)).along_scan ==
          dr4f_simulate(dr4f_parallax_system(without, 20.0)).along_scan
end

# ---------------------------------------------------------------------------
# Wide pairs: the excursion is composed with the frame on the sphere
# ---------------------------------------------------------------------------

# A resolved wide pair at high declination: A (1 Msun) and B (0.5 Msun) at
# 300 au, 13.5 pc — ~22″ apart, so each source sits 7–15″ from the
# barycentre. Fixed values: the tests below are about geometry.
const DR4F_WIDE = (ra=150.0, dec=70.0, plx=74.19, pmra=-172.57, pmdec=-381.32, rv=-26.9e3)

function dr4f_wide_bodies(; mass_B=0.5)
    A = Octofitter.Body(name="A", variables=@variables begin
        mass = 1.0
    end)
    B = Octofitter.Body(name="B", about=A, variables=@variables begin
        mass = $mass_B
        a = 300.0; e = 0.2; i = 1.0; ω = 0.5; Ω = 1.0; tp = 50000.0
    end)
    return (A, B)
end

# Spherical-trig inverse gnomonic: the direction [deg] whose gnomonic
# coordinates about (α0, δ0) [deg] are (x, y) [mas]. Written from the
# textbook formulae, independently of `_dr4_project`'s vector form.
function dr4f_ungnomonic(x, y, α0, δ0)
    xr, yr = (x, y) ./ (rad2deg(1.0) * 3.6e6)
    ρ = hypot(xr, yr)
    ρ == 0 && return (α0, δ0)
    c = atan(ρ)
    δ = asind(cos(c) * sind(δ0) + yr * sin(c) * cosd(δ0) / ρ)
    α = α0 + rad2deg(atan(xr * sin(c), ρ * cosd(δ0) * cos(c) - yr * sind(δ0) * sin(c)))
    return (α, δ)
end
# …and the forward projection, in its trigonometric form.
function dr4f_gnomonic(α, δ, α0, δ0)
    D = sind(δ) * sind(δ0) + cosd(δ) * cosd(δ0) * cosd(α - α0)
    s = rad2deg(1.0) * 3.6e6 / D
    return (cosd(δ) * sind(α - α0) * s,
            (sind(δ) * cosd(δ0) - cosd(δ) * sind(δ0) * cosd(α - α0)) * s)
end

@testset "wide pair: excursion composed on the sphere, source's own parallax" begin
    A, B = dr4f_wide_bodies()
    k = DR4F_WIDE
    mksys(obs) = Octofitter.System(name="wide", bodies=[A, B], observations=[obs],
        barycentric_lighttime=:off, observing_geometry=:off,
        variables=@variables begin
            plx = $(k.plx); ra = $(k.ra); dec = $(k.dec)
            pmra = $(k.pmra); pmdec = $(k.pmdec); rv = $(k.rv)
            ref_epoch = $DR4F_REF
        end)

    # Where A is at the reference epoch, from the model itself, so the tangent
    # point sits near the source as DR4's does.
    sys0 = mksys(GaiaDR4AstromObs(dr4f_table(); name="GaiaDR4", target=:A,
        tangent_point=(k.ra, k.dec)))
    nt = Octofitter.make_arr2nt(sys0)(Float64[])
    posys = Octofitter.make_ln_like(sys0).build(nt)
    bary = PlanetOrbits.barycentre(posys)
    sol_ref = orbitsolve(posys, DR4F_REF; Octofitter._solvekw(sys0)...)
    tA = dr4f_ungnomonic(raoff(sol_ref, :A, bary), decoff(sol_ref, :A, bary),
                         frame_ra(sol_ref), frame_dec(sol_ref))
    tA = (tA[1] + 1.0 / 3.6e6, tA[2] - 0.5 / 3.6e6)     # …give or take a mas

    for (name, tgt) in (("A", :A), ("B", :B))
        obs = GaiaDR4AstromObs(dr4f_table(ra0=tA[1], dec0=tA[2]); name="GaiaDR4", target=tgt)
        sim = dr4f_simulate(mksys(obs))
        traj = orbitsolve(posys, DR4F_EPOCHS; Octofitter._solvekw(sys0)...)

        exact = map(eachindex(DR4F_EPOCHS)) do i
            sol = traj[i]
            αβ = dr4f_ungnomonic(raoff(sol, tgt, bary), decoff(sol, tgt, bary),
                                 frame_ra(sol), frame_dec(sol))
            dr4f_gnomonic(αβ..., tA...)
        end
        first_order = map(eachindex(DR4F_EPOCHS)) do i
            sol = traj[i]
            ξB, ηB = dr4f_gnomonic(frame_ra(sol), frame_dec(sol), tA...)
            (ξB + raoff(sol, tgt, bary), ηB + decoff(sol, tgt, bary))
        end
        err_exact = maximum(i -> hypot(sim.ra_offset[i] - exact[i][1],
                                       sim.dec_offset[i] - exact[i][2]), eachindex(exact))
        err_first = maximum(i -> hypot(first_order[i][1] - exact[i][1],
                                       first_order[i][2] - exact[i][2]), eachindex(exact))
        sep = maximum(i -> hypot(raoff(traj[i], tgt, bary), decoff(traj[i], tgt, bary)),
                      eachindex(DR4F_EPOCHS)) / 1000
        @info "wide pair, source $name: composition error" excursion_arcsec = sep exact_mas = err_exact first_order_mas = err_first
        @test err_exact < 1e-5            # mas: rounding only
        @test err_first > 0.3             # mas: the error this replaces is real here

        # The parallax multiplying f_al is the source's own: inverting it with
        # the (separately written) `barycentre_parallax` recovers the
        # barycentre's, and its sign follows the line-of-sight offset.
        ϖ_src = Octofitter._dr4_source_plx(sol_ref, tgt, bary)
        ϖ_bar = Octofitter._dr4_frame_plx(sol_ref)
        dz = posz(sol_ref, tgt, bary)
        @test barycentre_parallax(ϖ_src, dz) ≈ ϖ_bar rtol = 1e-12
        @test sign(ϖ_src - ϖ_bar) == -sign(dz)
        # …and is what the abscissae use.
        ϖ_i = [Octofitter._dr4_source_plx(traj[i], tgt, bary) for i in eachindex(DR4F_EPOCHS)]
        w_expect = @. sim.ra_offset * sind(DR4F_PSI) + sim.dec_offset * cosd(DR4F_PSI) + ϖ_i * DR4F_FAL
        @test maximum(abs, sim.along_scan .- w_expect) < 1e-9
    end
end

# The usage example in docs/src/gaia-dr4-prerelease.md ("Two sources, one
# frame"), with simulated tables: an anchored absolute frame whose anchor
# position is a mas offset from the catalogue, and the two resolved
# components of a wide pair each observed by DR4 against its own tangent
# point.
@testset "two DR4 sources sharing an anchored absolute frame (docs example)" begin
    A, B = dr4f_wide_bodies()
    k = DR4F_WIDE
    # Catalogue positions (and so tangent points) of the two sources: A near
    # the frame centre, B ~22″ away. Made up, as a catalogue's would be off.
    RA_A, DEC_A = k.ra, k.dec
    RA_B, DEC_B = k.ra + 15.0 / 3600 / cosd(k.dec), k.dec + 16.0 / 3600

    obs_A = GaiaDR4AstromObs(dr4f_table(ra0=RA_A, dec0=DEC_A); target=A, ref=Barycentre,
        name="DR4_A")
    obs_B = GaiaDR4AstromObs(dr4f_table(ra0=RA_B, dec0=DEC_B); target=B, ref=Barycentre,
        name="DR4_B")

    sys = Octofitter.System(name="wide pair", bodies=[A, B], observations=[obs_A, obs_B],
        barycentric_lighttime=:off, observing_geometry=:off,
        variables=AnchoredFrame(A; ref_epoch=DR4F_REF, variables=@variables begin
            # A's position as mas offsets from its catalogue position…
            Δra_A  ~ Normal(0, 1000)
            Δdec_A ~ Normal(0, 1000)
            # …and the angles the frame needs, derived from them.
            ra_A  = $RA_A + Δra_A / (3.6e6 * cosd($DEC_A))
            dec_A = $DEC_A + Δdec_A / 3.6e6
            plx_A   ~ Uniform(60, 90)
            pmra_A  ~ Normal(k.pmra, 50)
            pmdec_A ~ Normal(k.pmdec, 50)
            rv_A    ~ Normal(k.rv, 5e3)
        end))

    # A noiseless draw is fitted exactly by both observations…
    θ = Octofitter.drawfrompriors(sys)
    clean = Octofitter.generate_from_params(sys, θ; add_noise=false)
    expected = 2 * (-DR4F_N * log(2π) / 2 - DR4F_N * log(0.08^2) / 2)
    @test Octofitter.make_ln_like(clean)(clean, θ) ≈ expected rtol = 1e-10

    # …and the model differentiates through the anchored frame and both
    # tracks. Checked near a noiseless truth, where ℓ is moderate: at a raw
    # prior draw ℓ ~ -1e11 and central differences are rounding-limited.
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    θv = Octofitter.sample_priors(Random.Xoshiro(7), sys)
    truth = Octofitter.generate_from_params(sys, model.arr2nt(θv); add_noise=false)
    model_t = Octofitter.LogDensityModel(truth; verbosity=0)
    θt = model_t.link(θv .* (1 .+ 1e-6 .* randn(Random.Xoshiro(8), length(θv))))
    ℓ, g = model_t.∇ℓπcallback(θt)
    @test isfinite(ℓ)
    @test all(isfinite, g)
    gfd = FiniteDiff.finite_difference_gradient(model_t.ℓπcallback, θt, Val(:central); relstep=1e-5)
    @info "anchored two-source model: gradient check" ℓ rel_err = norm(g .- gfd) / norm(g)
    @test norm(g .- gfd) / norm(g) < 1e-5
end

# `skytrackplot` draws `_dr4_sky_track` on a fine grid and places each transit
# with it, so at the data epochs it must *be* the likelihood's model.
@testset "the sky-track helper is the likelihood's model (skytrackplot)" begin
    A, B = dr4f_wide_bodies()
    k = DR4F_WIDE
    abs_sys = Octofitter.System(name="wide", bodies=[A, B],
        observations=[GaiaDR4AstromObs(dr4f_table(ra0=k.ra + 2e-4, dec0=k.dec - 1e-4);
                                       name="GaiaDR4", target=:A)],
        barycentric_lighttime=:off, observing_geometry=:off,
        variables=@variables begin
            plx = $(k.plx); ra = $(k.ra); dec = $(k.dec)
            pmra = $(k.pmra); pmdec = $(k.pmdec); rv = $(k.rv)
            ref_epoch = $DR4F_REF
        end)
    plx_sys = Octofitter.System(name="wideplx", bodies=[A, B],
        observations=[GaiaDR4AstromObs(dr4f_table(); name="GaiaDR4", target=:A,
            variables=@variables begin
                ra_offset_mas = 120.0; dec_offset_mas = -40.0
                pmra = $(k.pmra); pmdec = $(k.pmdec); ref_epoch = $DR4F_REF
            end)],
        barycentric_lighttime=:off, observing_geometry=:off,
        variables=@variables begin plx = $(k.plx) end)

    for (label, sys) in (("absolute", abs_sys), ("parallax-only", plx_sys))
        obs = only(sys.observations)
        sim = dr4f_simulate(sys)
        nt = Octofitter.make_arr2nt(sys)(Float64[])
        posys = Octofitter.make_ln_like(sys).build(nt)
        traj = orbitsolve(posys, DR4F_EPOCHS; Octofitter._solvekw(sys)...)
        sols = [traj[i] for i in eachindex(DR4F_EPOCHS)]
        key = Symbol(Octofitter.likelihoodname(obs))
        obsns = hasproperty(nt, :observations) ? nt.observations : (;)
        θ_obs = hasproperty(obsns, key) ? getproperty(obsns, key) : (;)
        trk = Octofitter._dr4_sky_track(obs, posys, sols, DR4F_EPOCHS, nt, θ_obs)
        w = @. trk.ξ * sind(DR4F_PSI) + trk.η * cosd(DR4F_PSI) + trk.ϖ * DR4F_FAL
        @test maximum(abs, trk.ξ .- sim.ra_offset) < 1e-9
        @test maximum(abs, trk.η .- sim.dec_offset) < 1e-9
        @test maximum(abs, w .- sim.along_scan) < 1e-9
        # The exaggeration knob scales the excursion and nothing else: the
        # parallax is untouched. On the flat parallax-only model doubling it
        # adds one more excursion; under the absolute frame the doubled
        # excursion is composed on the sphere like any other (checked with
        # the independent spherical trigonometry above).
        trk2 = Octofitter._dr4_sky_track(obs, posys, sols, DR4F_EPOCHS, nt, θ_obs;
                                         keplerian_mult=2.0)
        @test trk2.ϖ == trk.ϖ
        bary = PlanetOrbits.barycentre(posys)
        if label == "absolute"
            tp = obs.tangent_point
            err = maximum(eachindex(sols)) do i
                sol = sols[i]
                d = dr4f_ungnomonic(2raoff(sol, :A, bary), 2decoff(sol, :A, bary),
                                    frame_ra(sol), frame_dec(sol))
                e = dr4f_gnomonic(d..., tp...)
                hypot(trk2.ξ[i] - e[1], trk2.η[i] - e[2])
            end
            @test err < 1e-5
        else
            exc = [raoff(s, :A, bary) for s in sols]
            @test maximum(abs, (trk2.ξ .- trk.ξ) .- exc) < 1e-9
        end
        @info "sky track vs likelihood, $label" max_dw = maximum(abs, w .- sim.along_scan)
    end
end
