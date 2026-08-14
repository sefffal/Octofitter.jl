# Numerical agreement with Octofitter v1 on the model both versions can express.
#
# The reference values below were produced by `bench/bench-v1.jl` in the
# maintenance workspace, running the *same* 220-transit DR4-like dataset
# (checked in as `fixtures/dr4-like-scans.csv`) through Octofitter v8.2.5 and
# PlanetOrbits v0.11.4. Since v1 is deleted on this branch they are recorded
# fixtures, in the same spirit as PlanetOrbits' own v1 regression fixtures.
#
# The comparison is run with `observing_geometry=false`, which selects exactly
# v1's sky geometry: one shared AU→mas scale per epoch, no per-body
# light-travel retardation, no line-of-sight projection. That is the honest
# like-for-like configuration; the v2 default is *more* accurate, not
# equivalent, and is checked separately below to differ by the expected
# (small) amount rather than by nothing.

const V1_REFERENCE = Dict(
    1 => (D=13, ll=-1.81324644778998e7),
    2 => (D=20, ll=-1.79909269168743e7),
)

const V1_REF_EPOCH = 57388.0

function dr4_scans()
    rows, _ = readdlm(joinpath(@__DIR__, "fixtures", "dr4-like-scans.csv"),
        ',', Float64, '\n'; header=true)
    # The fixture is stored in v1's convention (scan angle in radians),
    # unchanged, because the reference log-likelihoods below were generated
    # from exactly these numbers. `GaiaDR4AstromObs` now ingests the scan
    # angle in degrees — the Gaia archive's own unit — so convert on the way
    # in. This is a unit relabelling only: `sind(rad2deg(ψ)) == sin(ψ)` to
    # within rounding, far inside the `rtol` asserted below.
    return Table(epoch=rows[:, 1], scan_pos_angle=rad2deg.(rows[:, 2]),
        parallax_factor_al=rows[:, 3], centroid_pos_al=rows[:, 4],
        centroid_pos_error_al=rows[:, 5], outlier_flag=Int.(rows[:, 6]))
end

function v1_equivalent_model(Np; observing_geometry=false)
    scans = dr4_scans()
    sysvars = Np == 1 ?
              @variables(begin
                  plx ~ Uniform(0.01, 100)
                  M_tot = 1.0
                  mass_b ~ LogUniform(0.01, 1000)
                  m_companions = mass_b * mjup
              end) :
              @variables(begin
                  plx ~ Uniform(0.01, 100)
                  M_tot = 1.0
                  mass_b ~ LogUniform(0.01, 1000)
                  mass_c ~ LogUniform(0.01, 1000)
                  m_companions = (mass_b + mass_c) * mjup
              end)

    A = Octofitter.Body(name="A", variables=@variables begin
        mass = system.M_tot - system.m_companions
    end)
    companion(nm, mv) = Octofitter.Body(name=nm, about=A, variables=@variables begin
        mass = getproperty(system, $mv) * mjup
        # v1 drove Kepler's third law from the *total* system mass for every
        # companion. That is the compatibility case `M=` exists for.
        M = system.M_tot
        a ~ LogUniform(0.01, 100)
        e ~ Uniform(0, 0.99)
        ω ~ Uniform(0, 2pi)
        i ~ Sine()
        Ω ~ Uniform(0, 2pi)
        θ ~ Uniform(0, 2pi)
        epoch = $V1_REF_EPOCH
    end)
    bodies = Np == 1 ? [A, companion("b", :mass_b)] :
             [A, companion("b", :mass_b), companion("c", :mass_c)]

    obs = GaiaDR4AstromObs(scans; target=A, ref=Barycentre, name="GaiaDR4",
        variables=@variables begin
            astrometric_jitter ~ LogUniform(0.00001, 10)
            ra_offset_mas ~ Normal(0, 1000)
            dec_offset_mas ~ Normal(0, 1000)
            pmra ~ Uniform(-1000, 1000)
            pmdec ~ Uniform(-1000, 1000)
            ref_epoch = $V1_REF_EPOCH
        end)

    return Octofitter.System(name="dr4like", bodies=bodies, observations=[obs],
        variables=sysvars, observing_geometry=observing_geometry)
end

v1_theta(Np) = Np == 1 ?
               [16.0, 10.0, 2.0, 0.2, 0.6, 1.1, 2.4, 3.5, 0.05, 0.5, -0.3, 5.3, -24.2] :
               [16.0, 10.0, 4.0,
                2.0, 0.2, 0.6, 1.1, 2.4, 3.5,
                6.0, 0.1, 1.7, 0.9, 0.4, 1.2,
                0.05, 0.5, -0.3, 5.3, -24.2]

@testset "v1 agreement, Np=$Np" for Np in (1, 2)
    ref = V1_REFERENCE[Np]
    sys = v1_equivalent_model(Np)
    nt = Octofitter.make_arr2nt(sys)(v1_theta(Np))
    @test length(Octofitter._list_priors(sys)) == ref.D
    ll = Octofitter.make_ln_like(sys)(sys, nt)
    # Bit-level: the two implementations do the same arithmetic in the same
    # order, and the epicyclic superposition v1 wrote out by hand is exactly
    # what the A⁻¹ transform produces for an astrocentric row set.
    @test ll ≈ ref.ll rtol = 1e-12
end

@testset "the default geometry differs, and by the expected amount" begin
    # `observing_geometry=true` adds per-body light-travel retardation, the
    # per-body depth scale, and the line-of-sight projection. On a 16 mas / 10
    # Mjup system the photocentre excursion is small, so the difference must be
    # present but tiny — a test that would fail both if the pass were a no-op
    # and if it were wrong by orders of magnitude.
    sys_cheap = v1_equivalent_model(1; observing_geometry=false)
    sys_full = v1_equivalent_model(1; observing_geometry=true)
    θ = v1_theta(1)
    ll_cheap = Octofitter.make_ln_like(sys_cheap)(sys_cheap, Octofitter.make_arr2nt(sys_cheap)(θ))
    ll_full = Octofitter.make_ln_like(sys_full)(sys_full, Octofitter.make_arr2nt(sys_full)(θ))
    @test ll_cheap != ll_full
    @test 1e-9 < abs(ll_full - ll_cheap) / abs(ll_cheap) < 1e-5
end
