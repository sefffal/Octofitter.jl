#!/usr/bin/env julia
#
# Gaia DR4 epoch astrometry in the v2 model surface — a worked example, and the
# before/after benchmark for the migration.
#
#   julia --project=. examples/gaia-dr4-v2.jl
#
# The dataset is a synthetic DR4-like scan law (220 transits over ~5.5 yr,
# checked in at test/v2/fixtures/dr4-like-scans.csv) so this runs offline and
# deterministically. Everything else is a real model.
#
# Three things are demonstrated:
#
#   1. The model surface. The host star is a node with a mass and a flux; the
#      companions are nodes with `about=`; the observation names what it
#      observes (`target=Photocentre`) and what it is measured against
#      (`ref=Barycentre`). There is no `basis=`, no `M`, no
#      `θ_at_epoch_to_tperi`, and no per-planet observation attachment.
#
#   2. The cost. The same model, parameter for parameter, against Octofitter
#      v8.2.5 / PlanetOrbits 0.11.4.
#
#   3. The N-body backend. `method=AHL21(h=…, t0=…)` swaps the propagator with
#      no other change: not one likelihood, reference, or variable moves.

using Octofitter
using Distributions
using BenchmarkTools
using DelimitedFiles
using Printf
using PlanetOrbits: KeplerianApprox, AHL21

const SCANS_CSV = joinpath(@__DIR__, "..", "test", "v2", "fixtures", "dr4-like-scans.csv")
const REF_EPOCH = 57388.0        # J2016.0

rows, _ = readdlm(SCANS_CSV, ',', Float64, '\n'; header=true)
scans = Table(
    epoch=rows[:, 1],                    # MJD
    scan_pos_angle=rows[:, 2],           # rad
    parallax_factor_al=rows[:, 3],
    centroid_pos_al=rows[:, 4],          # mas
    centroid_pos_error_al=rows[:, 5],    # mas
    outlier_flag=Int.(rows[:, 6]),
)
@info "Gaia DR4-like scan law" transits = length(scans.epoch) baseline_days = round(
    maximum(scans.epoch) - minimum(scans.epoch), digits=1)

# ---------------------------------------------------------------------------
# 1. The model
# ---------------------------------------------------------------------------

"""
The model as you would actually write it today: the star is a body with a
fitted mass, each companion is a body with its own mass, and the row's
gravitating mass follows from the bodies it binds. No `M`, no manual
`M = M_pri + M_sec*mjup2msol`.
"""
function dr4_model(Np; method=KeplerianApprox(), observing_geometry=true)
    A = Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.1), lower=0.2)     # M⊙
        # The host defines the flux scale; a dark companion leaves the
        # photocentre exactly on it, and a luminous one moves it — the same
        # code path either way.
        flux = 1.0
    end)

    companion(nm) = Body(name=nm, about=A, variables=@variables begin
        mass_jup ~ LogUniform(0.01, 1000)
        mass = mass_jup * mjup                            # M⊙; `mjup` is a constant
        a ~ LogUniform(0.01, 100)                         # AU
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        # The `θ` phase parametrization: sky-plane position angle at an epoch.
        # In v1 this needed `tp = θ_at_epoch_to_tperi(θ, epoch; M=system.M, e,
        # a, i, ω, Ω)` in every model block — and its `P` branch was broken.
        # It is now a constructor group: give `θ` and its `epoch`.
        θ ~ Uniform(0, 2pi)
        epoch = $REF_EPOCH
    end)

    bodies = Np == 1 ? [A, companion("b")] : [A, companion("b"), companion("c")]

    # The observation says what it looks at. For blended Gaia astrometry that
    # is the system's photocentre, measured against the barycentre — which is
    # what the fitted position and proper motion below describe.
    gaia = GaiaDR4AstromObs(scans;
        target=Photocentre, ref=Barycentre, name="GaiaDR4",
        variables=@variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)  # mas
            ra_offset_mas ~ Normal(0, 1000)
            dec_offset_mas ~ Normal(0, 1000)
            pmra ~ Uniform(-1000, 1000)                   # mas/yr
            pmdec ~ Uniform(-1000, 1000)                  # mas/yr
            ref_epoch = $REF_EPOCH
        end)

    return System(name="dr4-like", bodies=bodies, observations=[gaia],
        variables=(@variables begin
            plx ~ Uniform(0.01, 100)                      # mas → parallax frame
        end),
        method=method, observing_geometry=observing_geometry)
end

println("\n", "="^78)
println("The model")
println("="^78)
show(stdout, MIME"text/plain"(), dr4_model(2))

# ---------------------------------------------------------------------------
# 2. Before and after
# ---------------------------------------------------------------------------
#
# For a like-for-like comparison the benchmarked model is the *v1-equivalent*
# one: v1 treated `M` as the total system mass, fixed, and scaled every
# companion's reflex by it. Reproducing that needs the companion masses hoisted
# to the system block (so the host's mass can be written as the remainder) and
# the `M=` compatibility override on each row. Both are things v2 lets you say;
# neither is how you would write a new model.

function v1_equivalent_model(Np; method=KeplerianApprox(), observing_geometry=true)
    sysvars = Np == 1 ?
              (@variables begin
                  plx ~ Uniform(0.01, 100)
                  M_tot = 1.0
                  mass_b ~ LogUniform(0.01, 1000)
                  m_companions = mass_b * mjup
              end) :
              (@variables begin
                  plx ~ Uniform(0.01, 100)
                  M_tot = 1.0
                  mass_b ~ LogUniform(0.01, 1000)
                  mass_c ~ LogUniform(0.01, 1000)
                  m_companions = (mass_b + mass_c) * mjup
              end)

    A = Body(name="A", variables=@variables begin
        mass = system.M_tot - system.m_companions
        flux = 1.0
    end)
    companion(nm, massvar) = Body(name=nm, about=A, variables=@variables begin
        mass = getproperty(system, $massvar) * mjup
        M = system.M_tot          # compatibility: v1's fixed total mass
        a ~ LogUniform(0.01, 100)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        θ ~ Uniform(0, 2pi)
        epoch = $REF_EPOCH
    end)
    bodies = Np == 1 ? [A, companion("b", :mass_b)] :
             [A, companion("b", :mass_b), companion("c", :mass_c)]

    gaia = GaiaDR4AstromObs(scans; target=Photocentre, ref=Barycentre, name="GaiaDR4",
        variables=@variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)
            ra_offset_mas ~ Normal(0, 1000)
            dec_offset_mas ~ Normal(0, 1000)
            pmra ~ Uniform(-1000, 1000)
            pmdec ~ Uniform(-1000, 1000)
            ref_epoch = $REF_EPOCH
        end)

    return System(name="dr4-like", bodies=bodies, observations=[gaia],
        variables=sysvars, method=method, observing_geometry=observing_geometry)
end

# Reference point, natural domain: system priors → node priors → obs priors.
theta(Np) = Np == 1 ?
            [16.0, 10.0, 2.0, 0.2, 0.6, 1.1, 2.4, 3.5, 0.05, 0.5, -0.3, 5.3, -24.2] :
            [16.0, 10.0, 4.0,
             2.0, 0.2, 0.6, 1.1, 2.4, 3.5,
             6.0, 0.1, 1.7, 0.9, 0.4, 1.2,
             0.05, 0.5, -0.3, 5.3, -24.2]

# Recorded from `bench/bench-v1.jl` on Octofitter v8.2.5 / PlanetOrbits 0.11.4,
# same machine, same data.
const V1 = Dict(
    1 => (D=13, ll=-1.81324644778998e7, t_l=45.6, t_g=574.1),
    2 => (D=20, ll=-1.79909269168743e7, t_l=87.3, t_g=581.0),
)

function measure(sys, θ)
    model = Octofitter.LogDensityModel(sys; verbosity=0)
    @assert length(θ) == model.D "expected $(model.D) parameters, got $(length(θ))"
    θt = model.link(θ)
    nt = model.arr2nt(θ)
    ll = Octofitter.make_ln_like(sys)(sys, nt)
    bl = @benchmark $(model.ℓπcallback)($θt)
    bg = @benchmark $(model.∇ℓπcallback)($θt)
    return (; D=model.D, ll,
        t_l=minimum(bl).time / 1e3, t_g=minimum(bg).time / 1e3,
        alloc=minimum(bl).allocs)
end

println("\n", "="^78)
println("Before and after — 220 transits, same data, same parameters")
println("="^78)
@printf("%-3s %-26s %3s %10s %10s %8s %8s\n",
    "Np", "configuration", "D", "ℓπ [µs]", "∇ℓπ [µs]", "ℓπ×", "∇ℓπ×")
for Np in (1, 2)
    v1 = V1[Np]
    @printf("%-3d %-26s %3d %10.1f %10.1f %8s %8s\n",
        Np, "v1 (v8.2.5 / PO 0.11.4)", v1.D, v1.t_l, v1.t_g, "—", "—")
    for (label, og) in (("v2, v1 sky geometry", false), ("v2, full geometry", true))
        θ = theta(Np)
        r = measure(v1_equivalent_model(Np; observing_geometry=og), θ)
        # The `false` configuration selects exactly v1's sky geometry, so the
        # log likelihoods must agree; `true` adds per-body light-travel
        # retardation, per-body depth scaling and line-of-sight projection.
        tag = og ? "" : @sprintf("   Δll/ll = %.1e", abs(r.ll - v1.ll) / abs(v1.ll))
        @printf("%-3d %-26s %3d %10.1f %10.1f %7.1f× %7.1f×%s\n",
            Np, label, r.D, r.t_l, r.t_g, v1.t_l / r.t_l, v1.t_g / r.t_g, tag)
        @assert r.alloc == 0 "the hot path should not allocate; got $(r.alloc)"
    end
end

# ---------------------------------------------------------------------------
# 3. The N-body backend
# ---------------------------------------------------------------------------
#
# `method=` is the only change. The likelihood, the references, the variables
# and the parameter vector are all identical — the propagator seam is below
# every one of them.
#
# h ≲ P_min/20: the inner companion is at a = 2 AU about ~1 M⊙, so P ≈ 1030 d.
# `t0` is the epoch the elements are osculating at, and must be given
# explicitly for a parallax-only system (an absolute frame defaults it to
# `ref_epoch`).

println("\n", "="^78)
println("The same model under the AHL21 N-body propagator")
println("="^78)
@printf("%-3s %-26s %10s %10s %24s\n", "Np", "propagator", "ℓπ [µs]", "∇ℓπ [µs]", "log likelihood")
for Np in (1, 2)
    θ = theta(Np)
    kep = measure(v1_equivalent_model(Np), θ)
    nb = measure(v1_equivalent_model(Np; method=AHL21(h=40.0, t0=REF_EPOCH)), θ)
    @printf("%-3d %-26s %10.1f %10.1f %24.6f\n", Np, "KeplerianApprox", kep.t_l, kep.t_g, kep.ll)
    @printf("%-3d %-26s %10.1f %10.1f %24.6f\n", Np, "AHL21(h=40 d)", nb.t_l, nb.t_g, nb.ll)
    if Np == 1
        @printf("     two bodies: the map is exact, so the two must agree — Δ = %.2e\n",
            abs(nb.ll - kep.ll))
    else
        @printf("     three bodies: the difference (%.3g in log likelihood) is the\n",
            abs(nb.ll - kep.ll))
        println("     companion-companion interaction the Keplerian approximation drops.")
    end
end

println("\nDone.")
