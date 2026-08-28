# The Gaia NSS helpers on the v9 model surface: `query_nss`,
# `nss_to_starting_point`, `initialize_from_nss!` and `nss_to_model_chain`.
#
# Everything here is offline. `query_nss` normally hits the ESA archive, so the
# one testset that exercises it writes the response it would have received into
# the on-disk cache (inside a temporary working directory) and lets the helper
# read it back — no network, and nothing written outside the temp dir.
#
# The NSS row itself is synthesised from a *known* orbit rather than copied from
# the catalogue, which is what makes the round trip checkable: the Thiele-Innes
# constants below are exactly those of `nss_reference_system()`, so a correctly
# seeded model must reproduce them.

using Test
using Octofitter
using Octofitter: Octofitter, Body, System
using Octofitter.TypedTables: Table
using Distributions
using PlanetOrbits
using Logging
using Random
using Statistics

# --- the known orbit the synthetic solution describes -------------------------

const NSS_PLX       = 25.0     # mas
const NSS_HOST_MASS = 1.10     # M⊙
const NSS_A         = 1.60     # AU  (relative orbit)
const NSS_E         = 0.35
const NSS_I         = 0.98     # rad
const NSS_Ω         = 2.40     # rad
const NSS_ω         = 1.30     # rad
const NSS_TP        = 57500.0  # MJD

"The PlanetOrbits system the synthetic NSS solution was generated from."
function nss_reference_system()
    A = PlanetOrbits.Body(mass=NSS_HOST_MASS, name=:A)
    b = PlanetOrbits.Body(mass=0.0, name=:b)
    return PlanetOrbits.System((A, b),
        (PlanetOrbits.Orbit(b, about=A; a=NSS_A, e=NSS_E, i=NSS_I,
                            ω=NSS_ω, Ω=NSS_Ω, tp=NSS_TP),);
        plx=NSS_PLX)
end

const NSS_REF_SYS = nss_reference_system()
const NSS_TI      = PlanetOrbits.thieleinnes(NSS_REF_SYS, 1; plx=NSS_PLX)  # mas
const NSS_P_DAYS  = period(NSS_REF_SYS, 1)

# The NSS table stores periastron as days from JD 2457389.0 (the DR3 reference),
# so this is the value whose MJD is `NSS_TP`.
const NSS_T_PERIASTRON = NSS_TP - (2457389.0 - 2400000.5)

const NSS_SOURCE_ID = 4295745059252873600

"The synthetic solution, as `query_nss` would return it."
nss_solution(; errors::Bool=true) = errors ?
    (; source_id = NSS_SOURCE_ID,
       period = NSS_P_DAYS,            period_error = 2.0,
       eccentricity = NSS_E,           eccentricity_error = 0.01,
       t_periastron = NSS_T_PERIASTRON, t_periastron_error = 3.0,
       a_thiele_innes = NSS_TI.A,      a_thiele_innes_error = 0.05,
       b_thiele_innes = NSS_TI.B,      b_thiele_innes_error = 0.05,
       f_thiele_innes = NSS_TI.F,      f_thiele_innes_error = 0.05,
       g_thiele_innes = NSS_TI.G,      g_thiele_innes_error = 0.05,
       parallax = NSS_PLX) :
    (; source_id = NSS_SOURCE_ID,
       period = NSS_P_DAYS,
       eccentricity = NSS_E,
       t_periastron = NSS_T_PERIASTRON,
       a_thiele_innes = NSS_TI.A,
       b_thiele_innes = NSS_TI.B,
       f_thiele_innes = NSS_TI.F,
       g_thiele_innes = NSS_TI.G)

const NSS_EPOCHS = collect(range(56800.0, 58200.0, length=8))

# The v8 helpers put Kepler's third law in *julian* years
# (`a = ∛(M (P/365.25)²)`, and the mirror-image `M = a³/P²` in
# `nss_to_model_chain`), where PlanetOrbits' own third law uses the Kepler year,
# 365.2568983840419 days. The two differ by 1.9e-5, so a seeded semi-major axis
# is high by `NSS_YEAR_SLIP^(2/3)` and the orbit it implies has a period long by
# `NSS_YEAR_SLIP`. That is far inside any NSS error bar and irrelevant for a
# *starting point*, which is all these values are — the formulas are carried
# over unchanged rather than quietly corrected. Pinned exactly here so that the
# day someone does reconcile them, the test says which way and by how much.
const NSS_YEAR_SLIP = PlanetOrbits.kepler_year_to_julian_day_conversion_factor /
                      PlanetOrbits.year2day_julian

# A small Campbell-parameterized v9 model of the same system, with exact
# relative astrometry of the true orbit. Masses, `plx` and the astrometric
# jitter are all fixed, so the model's *free* variables are exactly the six an
# NSS solution maps onto — which is what lets `startingpoints!` and the
# all-fixed path through `initialize!` run without any optimization.
function nss_test_model()
    tr = orbitsolve(NSS_REF_SYS, NSS_EPOCHS)
    astrom = RelAstromObs(
        Table(epoch = NSS_EPOCHS,
              ra    = [raoff(tr[k], :b, :A) for k in eachindex(NSS_EPOCHS)],
              dec   = [decoff(tr[k], :b, :A) for k in eachindex(NSS_EPOCHS)],
              σ_ra  = fill(2.0, length(NSS_EPOCHS)),
              σ_dec = fill(2.0, length(NSS_EPOCHS)));
        target=:b, ref=:A, name="astrom",
        variables=@variables begin
            jitter = 0.0
        end)

    A = Body(name="A", variables=@variables begin
        mass = $NSS_HOST_MASS
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass = 0.0
        e  ~ Uniform(0, 0.9)
        a  ~ LogUniform(0.1, 100)
        i  ~ Sine()
        Ω  ~ Uniform(0, 2pi)
        ω  ~ Uniform(0, 2pi)
        tp ~ Uniform(56000, 59000)
    end)
    sys = System(name="nsstest", bodies=[A, b], observations=[astrom],
        verbosity=0,
        variables=@variables begin
            plx = $NSS_PLX
        end)
    return Octofitter.LogDensityModel(sys, verbosity=0)
end

# Log-likelihood of `model` at the (constrained) parameter vector `θ`.
function nss_ln_like(model, θ)
    nt = model.arr2nt(θ)
    return Octofitter.make_ln_like(model.system, nt)(model.system, nt)
end

# ---------------------------------------------------------------------------

@testset "_ti_to_campbell" begin
    # Face-on prograde orbit
    i, Ω, ω, α = Octofitter._ti_to_campbell(10.0, 0.0, 0.0, 10.0)
    @test α ≈ 10.0
    @test i ≈ 0.0 atol=1e-10

    # Face-on retrograde orbit
    i2, _, _, α2 = Octofitter._ti_to_campbell(10.0, 0.0, 0.0, -10.0)
    @test α2 ≈ 10.0
    @test i2 ≈ π atol=1e-10

    # Round trip against the known orbit. `(ω, Ω)` and `(ω+π, Ω+π)` give
    # identical Thiele-Innes constants, so only the sums and differences of the
    # angles are pinned — the branch is not.
    i3, Ω3, ω3, α3 = Octofitter._ti_to_campbell(NSS_TI.A, NSS_TI.B, NSS_TI.F, NSS_TI.G)
    @test α3 ≈ NSS_A * NSS_PLX rtol=1e-10
    @test i3 ≈ NSS_I rtol=1e-10
    @test sin(ω3 + Ω3) ≈ sin(NSS_ω + NSS_Ω) atol=1e-10
    @test cos(ω3 + Ω3) ≈ cos(NSS_ω + NSS_Ω) atol=1e-10
    @test sin(ω3 - Ω3) ≈ sin(NSS_ω - NSS_Ω) atol=1e-10
    @test cos(ω3 - Ω3) ≈ cos(NSS_ω - NSS_Ω) atol=1e-10

    # It agrees with PlanetOrbits' own inversion up to exactly that branch:
    # `ThieleInnes` folds the node into [0, π), this one does not. Both describe
    # the same orbit. The v8 copy is kept so a v8 seed is reproduced exactly.
    ti = PlanetOrbits.ThieleInnes(A=NSS_TI.A, B=NSS_TI.B, F=NSS_TI.F, G=NSS_TI.G,
                                  plx=NSS_PLX)
    @test ti.a ≈ α3 / NSS_PLX rtol=1e-12
    @test ti.i ≈ i3 rtol=1e-12
    @test mod(ti.Ω - Ω3, π) ≈ 0 atol=1e-10
    @test mod(ti.ω - ω3, π) ≈ 0 atol=1e-10
end

@testset "_nss_get_float" begin
    nt = (; a=1.5, b=missing, c=nothing, d=NaN, e=Inf)
    @test Octofitter._nss_get_float(nt, :a) ≈ 1.5
    @test isnothing(Octofitter._nss_get_float(nt, :b))
    @test isnothing(Octofitter._nss_get_float(nt, :c))
    @test isnothing(Octofitter._nss_get_float(nt, :d))
    @test isnothing(Octofitter._nss_get_float(nt, :e))
    @test isnothing(Octofitter._nss_get_float(nt, :nonexistent))
end

@testset "nss_to_starting_point: Campbell body" begin
    model = nss_test_model()
    seed = with_logger(NullLogger()) do
        Octofitter.nss_to_starting_point(nss_solution(), model; body=:b)
    end

    # v9 nests under `bodies`, not `planets`.
    @test hasproperty(seed, :bodies)
    @test !hasproperty(seed, :planets)
    @test hasproperty(seed.bodies, :b)
    vals = seed.bodies.b

    @test vals.e ≈ NSS_E
    @test vals.tp ≈ NSS_TP
    # `a` comes from the NSS period and the model's own total mass, which is
    # fixed at NSS_HOST_MASS here — so Kepler's third law returns the true
    # value, up to the julian-vs-Kepler year slip described above.
    @test vals.a ≈ NSS_A * NSS_YEAR_SLIP^(2/3) rtol=1e-10
    @test vals.i ≈ NSS_I rtol=1e-10
    # This solution's node lands on the far branch — both are the same orbit.
    @test mod(vals.Ω - NSS_Ω, π) ≈ 0 atol=1e-10
    @test mod(vals.ω - NSS_ω, π) ≈ 0 atol=1e-10

    # A `Body` node may be passed instead of its name.
    node = model.system.nodes[findfirst(n -> n.name === :b, model.system.nodes)]
    seed2 = with_logger(NullLogger()) do
        Octofitter.nss_to_starting_point(nss_solution(), model; body=node)
    end
    @test seed2.bodies.b == vals

    @test_throws r"does not have a body named" Octofitter.nss_to_starting_point(
        nss_solution(), model; body=:nosuchbody)

    # A body with no free variables maps nothing at all.
    empty_seed = with_logger(NullLogger()) do
        Octofitter.nss_to_starting_point(nss_solution(), model; body=:A)
    end
    @test isempty(propertynames(empty_seed))
end

@testset "nss_to_starting_point: Thiele-Innes body" begin
    A = Body(name="A", variables=@variables begin
        mass = $NSS_HOST_MASS
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass = 0.0
        e ~ Uniform(0, 0.9)
        A ~ Normal(0, 1000)   # mas
        B ~ Normal(0, 1000)
        F ~ Normal(0, 1000)
        G ~ Normal(0, 1000)
        ti = PlanetOrbits.ThieleInnes(A=A, B=B, F=F, G=G, plx=system.plx)
        a = ti.a
        i = ti.i
        ω = ti.ω
        Ω = ti.Ω
        tp ~ Uniform(56000, 59000)
    end)
    sys = System(name="nsstest_ti", bodies=[A, b], observations=(),
        verbosity=0,
        variables=@variables begin
            plx = $NSS_PLX
        end)
    model = Octofitter.LogDensityModel(sys, verbosity=0)

    vals = with_logger(NullLogger()) do
        Octofitter.nss_to_starting_point(nss_solution(), model; body=:b)
    end.bodies.b

    # The constants are set directly, and nothing tries to set the derived
    # Campbell elements.
    @test vals.A ≈ NSS_TI.A
    @test vals.B ≈ NSS_TI.B
    @test vals.F ≈ NSS_TI.F
    @test vals.G ≈ NSS_TI.G
    @test vals.e ≈ NSS_E
    @test vals.tp ≈ NSS_TP
    @test !hasproperty(vals, :a)
    @test !hasproperty(vals, :i)
end

@testset "nss_to_starting_point: UniformCircular angles" begin
    A = Body(name="A", variables=@variables begin
        mass = $NSS_HOST_MASS
    end)
    b = Body(name="b", about=A, variables=@variables begin
        mass = 0.0
        e ~ Uniform(0, 0.9)
        P ~ LogUniform(10, 10_000)
        i ~ Sine()
        Ω ~ UniformCircular()
        ω ~ UniformCircular()
        θ ~ UniformCircular()
        epoch = 57000.0
    end)
    sys = System(name="nsstest_uc", bodies=[A, b], observations=(),
        verbosity=0,
        variables=@variables begin
            plx = $NSS_PLX
        end)
    model = Octofitter.LogDensityModel(sys, verbosity=0)

    vals = with_logger(NullLogger()) do
        Octofitter.nss_to_starting_point(nss_solution(), model; body=:b)
    end.bodies.b

    # `Ω ~ UniformCircular()` samples `Ωx`/`Ωy`, so those are what get set.
    @test !hasproperty(vals, :Ω)
    @test !hasproperty(vals, :ω)
    @test hypot(vals.Ωx, vals.Ωy) ≈ 1 atol=1e-12
    @test hypot(vals.ωx, vals.ωy) ≈ 1 atol=1e-12
    @test atan(vals.Ωy, vals.Ωx) ≈ NSS_Ω - π atol=1e-10
    @test atan(vals.ωy, vals.ωx) ≈ NSS_ω - π atol=1e-10

    # A free `P` is mapped in preference to a semi-major axis.
    @test vals.P ≈ NSS_P_DAYS
    @test !hasproperty(vals, :a)
    # `tp` is not a free variable of a θ-phased body, so it is left alone.
    @test !hasproperty(vals, :tp)
end

@testset "a seeded starting point reproduces the NSS solution" begin
    model = nss_test_model()
    seed = with_logger(NullLogger()) do
        Octofitter.nss_to_starting_point(nss_solution(), model; body=:b)
    end
    startingpoints!(model, seed; verbosity=0)

    θ = model.invlink(model.starting_points[1])
    nt = model.arr2nt(θ)

    # The whole point: the seeded model's orbit has the NSS Thiele-Innes
    # constants back, which is the branch-independent statement of "same orbit".
    # The scale carries the year slip (it came in through `a`); the shape,
    # orientation and phase are exact.
    posys = Octofitter.construct_system(model, nt)
    back = PlanetOrbits.thieleinnes(posys, 1; plx=NSS_PLX)
    @test back.A ≈ NSS_TI.A * NSS_YEAR_SLIP^(2/3) rtol=1e-9
    @test back.B ≈ NSS_TI.B * NSS_YEAR_SLIP^(2/3) rtol=1e-9
    @test back.F ≈ NSS_TI.F * NSS_YEAR_SLIP^(2/3) rtol=1e-9
    @test back.G ≈ NSS_TI.G * NSS_YEAR_SLIP^(2/3) rtol=1e-9
    @test period(posys, 1) ≈ NSS_P_DAYS * NSS_YEAR_SLIP rtol=1e-9
    @test eccentricity(posys, 1) ≈ NSS_E rtol=1e-9
    @test periastron(posys, 1) ≈ NSS_TP rtol=1e-9

    # …and the likelihood is finite there. The data were generated noiselessly
    # from the same orbit, so the seed sits essentially at the optimum — within
    # the 1.3e-5 scale slip, which on 2 mas errors is worth ~1e-5 in log
    # likelihood.
    lnl = nss_ln_like(model, θ)
    @test isfinite(lnl)
    @test isfinite(model.ℓπcallback(model.starting_points[1]))

    # Free variables are in declaration order: e, a, i, Ω, ω, tp.
    truth = model.arr2nt([NSS_E, NSS_A, NSS_I, NSS_Ω, NSS_ω, NSS_TP])
    @test lnl ≈ Octofitter.make_ln_like(model.system, truth)(model.system, truth) rtol=1e-5
end

@testset "query_nss and initialize_from_nss!" begin
    sol = nss_solution()
    cols = collect(keys(sol))
    csv = join(string.(cols), ",") * "\n" *
          join((string(getproperty(sol, c)) for c in cols), ",") * "\n"

    mktempdir() do dir
        cd(dir) do
            # `_gaia_cache_path` prefers a `<cwd>/<dir>/<name>` file if one
            # exists, so this stands in for the archive response.
            mkpath("_gaia_nss_dr3")
            write(joinpath("_gaia_nss_dr3", "source-$NSS_SOURCE_ID.csv"), csv)

            parsed = query_nss(gaia_id=NSS_SOURCE_ID)
            @test parsed isa NamedTuple
            @test parsed.period ≈ NSS_P_DAYS
            @test parsed.eccentricity ≈ NSS_E
            @test parsed.a_thiele_innes ≈ NSS_TI.A
            @test parsed.source_id == NSS_SOURCE_ID

            @test_throws r"Unsupported catalog" query_nss(gaia_id=NSS_SOURCE_ID, catalog=:dr5)

            model = nss_test_model()
            # Every free variable is covered by the NSS mapping, so no global
            # optimization runs; pathfinder is off to keep this a unit test.
            with_logger(NullLogger()) do
                initialize_from_nss!(model; gaia_id=NSS_SOURCE_ID, body=:b,
                    use_pathfinder=false, ndraws=10, verbosity=0)
            end
            @test length(model.starting_points) == 10

            nt = model.arr2nt(model.invlink(model.starting_points[1]))
            posys = Octofitter.construct_system(model, nt)
            back = PlanetOrbits.thieleinnes(posys, 1; plx=NSS_PLX)
            @test back.A ≈ NSS_TI.A * NSS_YEAR_SLIP^(2/3) rtol=1e-9
            @test back.G ≈ NSS_TI.G * NSS_YEAR_SLIP^(2/3) rtol=1e-9
            @test period(posys, 1) ≈ NSS_P_DAYS * NSS_YEAR_SLIP rtol=1e-9
            @test isfinite(model.ℓπcallback(model.starting_points[1]))
        end
    end
end

@testset "nss_to_model_chain" begin
    # With no reported errors every prior is effectively a delta, so draw 1 is
    # the published solution and the derived quantities must agree with the
    # orbit the row was generated from.
    m, c = with_logger(NullLogger()) do
        Octofitter.nss_to_model_chain(nss_solution(errors=false); plx=NSS_PLX, N=50)
    end
    @test m isa Octofitter.LogDensityModel
    @test m.system.name === :NSS
    @test m.system.bodynames == [:A, :b]
    @test size(c, 1) == 50

    posys = Octofitter.construct_system(m, c, 1)
    # Here the *mass* is the quantity Kepler's third law produced, so it — and
    # through it the period — is what carries the julian-vs-Kepler year slip.
    # `a` and the Thiele-Innes constants come straight off the row and are exact.
    @test totalmass(posys) ≈ NSS_HOST_MASS / NSS_YEAR_SLIP^2 rtol=1e-5
    @test semimajoraxis(posys, 1) ≈ NSS_A rtol=1e-5
    @test period(posys, 1) ≈ NSS_P_DAYS * NSS_YEAR_SLIP rtol=1e-5
    @test eccentricity(posys, 1) ≈ NSS_E rtol=1e-5
    @test inclination(posys, 1) ≈ NSS_I rtol=1e-5
    back = PlanetOrbits.thieleinnes(posys, 1; plx=NSS_PLX)
    @test back.A ≈ NSS_TI.A rtol=1e-5
    @test back.B ≈ NSS_TI.B rtol=1e-5
    @test back.F ≈ NSS_TI.F rtol=1e-5
    @test back.G ≈ NSS_TI.G rtol=1e-5

    # With errors, the chain is a scatter centred on the solution. Column names
    # follow the usual `<body>_<var>` convention, derived variables included.
    m2, c2 = with_logger(NullLogger()) do
        Octofitter.nss_to_model_chain(nss_solution(); plx=NSS_PLX, N=2000)
    end
    cn = string.(names(c2))
    for col in ("b_e", "b_A", "b_B", "b_F", "b_G", "b_tp", "b_a", "b_i", "A_mass")
        @test col in cn
    end
    @test median(c2["b_e"][:]) ≈ NSS_E atol=0.005
    @test median(c2["b_A"][:]) ≈ NSS_TI.A atol=0.05
    @test std(c2["b_A"][:]) ≈ 0.05 rtol=0.15

    # The parallax comes off the NSS row when it is not given explicitly.
    m3, _ = with_logger(NullLogger()) do
        Octofitter.nss_to_model_chain(nss_solution(); N=10)
    end
    @test m3.system.derived.variables[:plx] ≈ NSS_PLX

    @test_throws r"does not contain Thiele-Innes" Octofitter.nss_to_model_chain(
        (; period=100.0, eccentricity=0.1); plx=NSS_PLX, N=10)
end
