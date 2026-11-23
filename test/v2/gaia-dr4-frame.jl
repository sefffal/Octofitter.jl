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
    # here with unit vectors (independently of the implementation's
    # trigonometric form). What remains is the rotation between the frame's
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
