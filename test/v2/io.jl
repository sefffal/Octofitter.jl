# Chain persistence, orbitize!/whereistheplanet interop, and the substellar
# flux tables.
#
# The round trips are the point: every one of these formats is a place where a
# column can quietly land under the wrong name, and nothing downstream errors
# when it does.

using FITSIO: FITSIO, FITS, FITSHeader
using HDF5: HDF5, h5open, attrs, create_dataset
using MCMCChains: MCMCChains, Chains
using ForwardDiff
using Random

# A small but complete v2 model: one host, one companion, priors on everything
# an orbitize! export needs.
function io_test_model()
    A = Octofitter.Body(name="A", variables=@variables begin
        mass ~ truncated(Normal(1.0, 0.05), lower=0.1)
    end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass ~ LogUniform(0.5mjup, 50mjup)
        a ~ LogUniform(5, 50)
        e ~ Uniform(0, 0.5)
        i ~ Sine()
        ω ~ Uniform(0, 2pi)
        Ω ~ Uniform(0, 2pi)
        tp ~ Uniform(50000, 60000)
    end)
    sys = Octofitter.System(name="iotest", bodies=[A, b], observations=(),
        variables=@variables begin
            plx ~ truncated(Normal(50.0, 1.0), lower=1.0)
        end)
    return Octofitter.LogDensityModel(sys, verbosity=0)
end

function io_test_chain(model, n=16)
    rng = Random.Xoshiro(1234)
    nts = [model.arr2nt(collect(model.sample_priors(rng))) for _ in 1:n]
    return Octofitter.result2mcmcchain(nts)
end

@testset "FITS chain round trip" begin
    # Unicode parameter names are the reason for the latexify/restore dance.
    chn = Chains(randn(Random.Xoshiro(1), 32, 5, 1), [:a, :e, :ω, :Ω, :θ])
    fname = tempname() * ".fits"
    Octofitter.savechain(fname, chn)
    loaded = Octofitter.loadchain(fname)

    @test keys(chn) == keys(loaded)
    for nm in keys(chn)
        @test vec(chn[nm]) ≈ vec(loaded[nm])
    end
    # The file says which model surface wrote it.
    @test loaded.info.chainformat == Octofitter.CHAIN_FORMAT_VERSION
    rm(fname)
end

@testset "FITS chain: multiple chains, metadata, sections" begin
    chn = Chains(randn(Random.Xoshiro(2), 12, 3, 4), [:x, :loglike, :logpost],
                 Dict(:internals => [:loglike, :logpost]))
    chn = MCMCChains.setinfo(chn, (;
        model_name=:HD12345,
        sampler="nuts",
        draws=12,
        acceptance=[0.8, 0.81, 0.79, 0.82],   # numeric vector: expanded across cards
    ))
    fname = tempname() * ".fits"
    Octofitter.savechain(fname, chn)
    loaded = Octofitter.loadchain(fname)

    @test size(loaded) == size(chn)
    for nm in keys(chn)
        @test collect(chn[nm]) ≈ collect(loaded[nm])
    end
    @test String(loaded.info.model_name) == "HD12345"
    @test loaded.info.sampler == "nuts"
    @test loaded.info.draws == 12
    @test loaded.info.acceptance ≈ [0.8, 0.81, 0.79, 0.82]

    # v1 dropped the section map, so `tree_depth`-style internals came back
    # masquerading as model parameters.
    @test Set(MCMCChains.names(loaded, :internals)) == Set([:loglike, :logpost])
    @test MCMCChains.names(loaded, :parameters) == [:x]
    rm(fname)
end

@testset "FITS chain: non-1-based iteration numbering" begin
    # A sliced or thinned chain does not start at iteration 1. v1 used the
    # iteration number as an array index, which threw `BoundsError` here.
    full = Chains(randn(Random.Xoshiro(3), 40, 2, 2), [:a, :b])
    sliced = full[21:40]
    fname = tempname() * ".fits"
    Octofitter.savechain(fname, sliced)
    loaded = Octofitter.loadchain(fname)
    @test size(loaded) == size(sliced)
    for nm in keys(sliced)
        @test collect(sliced[nm]) ≈ collect(loaded[nm])
    end
    @test first(MCMCChains.range(loaded)) == first(MCMCChains.range(sliced))
    rm(fname)
end

@testset "FITS chain: info values that are not header-shaped" begin
    # v1 stringified whatever was in `info`, then sliced the string by *byte*
    # at 255 — which throws when a multi-byte character straddles that byte,
    # and which formats the entire object first regardless.
    chn = Chains(randn(Random.Xoshiro(4), 8, 2, 1), [:a, :b])
    chn = MCMCChains.setinfo(chn, (;
        samples_transformed=[randn(Random.Xoshiro(5), 7) for _ in 1:2000],
        unicode_padding="Ω"^400,
    ))
    fname = tempname() * ".fits"
    @test Octofitter.savechain(fname, chn) == fname
    loaded = Octofitter.loadchain(fname)
    @test occursin("Vector", String(loaded.info.samples_transformed))
    @test length(String(loaded.info.unicode_padding)) <= 255
    rm(fname)
end

# Writes the two-HDU layout Octofitter v1 produced: no format stamp, no
# section table. Used to check that `loadchain` recognises it rather than
# handing back a chain whose columns silently will not match a v2 model.
function write_v1_format_chain(fname, coltitles, data)
    FITS(fname, "w") do fits
        h = FITSHeader(["model_name"], Any["HD12345"], ["chain-info"])
        write(fits, zeros(1, 1), header=h)
        write(fits, String.(coltitles), data; hdutype=FITSIO.TableHDU)
    end
    return fname
end

@testset "FITS chain: a v1 file is recognised, not silently accepted" begin
    fname = tempname() * ".fits"
    write_v1_format_chain(fname,
        ["iteration", "chain", "b_a", "b_GPI_jitter"],
        Any[collect(1:6), ones(Int, 6), collect(1.0:6.0), collect(11.0:16.0)])

    loaded = @test_logs (:warn,) match_mode = :any Octofitter.loadchain(fname)
    @test loaded.info.chainformat == 1
    # The samples themselves are intact — the file is not corrupt, it is old.
    @test vec(loaded[:b_a]) ≈ collect(1.0:6.0)
    rm(fname)
end

@testset "checkchain: model ↔ chain agreement" begin
    model = io_test_model()
    chain = io_test_chain(model)
    @test Octofitter.checkchain(model, chain) === chain

    # v1's three-level naming for a planet-owned observation is the failure
    # mode `mcmcchain2result` turns into `missing` rather than an error.
    v1ish = Chains(randn(Random.Xoshiro(6), 4, 2, 1), [:b_GPI_jitter, :b_a])
    A = Octofitter.Body(name="A", variables=@variables begin mass = 1.0 end)
    b = Octofitter.Body(name="b", about=A, variables=@variables begin
        mass = 1mjup
        a = 10.0; e = 0.1; i = 0.5; ω = 0.1; Ω = 0.2; tp = 55000.0
    end)
    obs = Octofitter.RelAstromObs(
        Table(epoch=[57000.0], ra=[100.0], dec=[50.0], σ_ra=[1.0], σ_dec=[1.0]);
        target=b, ref=A, name="GPI",
        variables=@variables begin jitter ~ LogUniform(0.1, 10) end)
    sys2 = Octofitter.System(name="v1names", bodies=[A, b], observations=(obs,),
        variables=@variables begin plx = 50.0 end)
    err = try
        Octofitter.checkchain(sys2, v1ish)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("GPI_jitter", err.msg)
    @test occursin("b_GPI_jitter", err.msg)
    @test occursin("three-level", err.msg)

    # And the non-fatal spelling still returns the chain.
    @test (@test_logs (:warn,) match_mode = :any Octofitter.checkchain(sys2, v1ish; strict=false)) === v1ish
end

@testset "loadchain(fname; model=…) checks as it loads" begin
    model = io_test_model()
    chain = io_test_chain(model)
    fname = tempname() * ".fits"
    Octofitter.savechain(fname, chain)
    @test Octofitter.loadchain(fname; model) isa Chains

    other = Chains(randn(Random.Xoshiro(7), 4, 1, 1), [:not_a_parameter])
    fname2 = tempname() * ".fits"
    Octofitter.savechain(fname2, other)
    @test_throws ErrorException Octofitter.loadchain(fname2; model)
    rm(fname); rm(fname2)
end

# ---------------------------------------------------
# orbitize! interop
# ---------------------------------------------------

const KEPLER_F = PlanetOrbits.kepler_year_to_julian_day_conversion_factor

function write_orbitize_posterior(fname, labels, cols; tau_ref_epoch=58849)
    dat = permutedims(reduce(hcat, cols))   # (nparams, nsamples)
    h5open(fname, "w") do f
        dset = create_dataset(f, "post", Float64, size(dat))
        write(dset, dat)
        attrs(f)["parameter_labels"] = String.(labels)
        attrs(f)["tau_ref_epoch"] = tau_ref_epoch
    end
    return fname
end

@testset "loadhdf5: orbitize! standard basis" begin
    n = 5
    sma = collect(range(9.0, 13.0, length=n))
    ecc = fill(0.2, n)
    inc = fill(0.7, n)
    aop = fill(1.1, n)
    pan = fill(2.2, n)
    tau = collect(range(0.1, 0.9, length=n))
    plx = fill(25.0, n)
    m0 = fill(1.0, n)
    m1 = collect(range(0.005, 0.02, length=n))

    fname = tempname() * ".hdf5"
    write_orbitize_posterior(fname,
        ["sma1", "ecc1", "inc1", "aop1", "pan1", "tau1", "plx", "m0", "m1"],
        [sma, ecc, inc, aop, pan, tau, plx, m0, m1])

    chn = Octofitter.loadhdf5(fname)
    @test vec(chn[:b_a]) ≈ sma
    @test vec(chn[:b_e]) ≈ ecc
    @test vec(chn[:b_i]) ≈ inc
    @test vec(chn[:b_ω]) ≈ aop
    @test vec(chn[:b_Ω]) ≈ pan
    @test vec(chn[:b_τ]) ≈ tau
    @test vec(chn[:plx]) ≈ plx

    # Masses are solar masses, as orbitize! stored them. v1 divided the
    # companion column by `mjup2msol` because a v1 `Planet.mass` was in
    # Jupiter masses; v2 has one mass unit.
    @test vec(chn[:A_mass]) ≈ m0
    @test vec(chn[:b_mass]) ≈ m1

    # τ → tp, against the total mass v1's superposition loop also produced for
    # a single companion (`m0 + m1`). This is the bit-for-bit claim.
    M_v1 = m0 .+ m1
    P_v1 = @. √(sma^3 / M_v1) * KEPLER_F
    @test vec(chn[:b_tp]) ≈ @.(tau * P_v1 + 58849)
    rm(fname)
end

@testset "loadhdf5: mtot-only posterior, and multiple chains" begin
    n = 8
    sma = collect(range(9.0, 13.0, length=n))
    tau = collect(range(0.05, 0.95, length=n))
    mtot = fill(1.2, n)
    fname = tempname() * ".hdf5"
    write_orbitize_posterior(fname,
        ["sma1", "ecc1", "inc1", "aop1", "pan1", "tau1", "plx", "mtot"],
        [sma, fill(0.1, n), fill(0.6, n), fill(1.0, n), fill(2.0, n), tau,
         fill(30.0, n), mtot])

    chn = Octofitter.loadhdf5(fname)
    P = @. √(sma^3 / mtot) * KEPLER_F
    @test vec(chn[:b_tp]) ≈ @.(tau * P + 58849)
    # No per-body masses are recoverable from a total, and the total keeps its
    # own name rather than being renamed onto a body.
    @test !haskey(chn, :A_mass)
    @test haskey(chn, :mtot)

    # v1's multi-chain path built a 3-D array and then `hcat`ed 2-D `tp`
    # columns onto it, which cannot work; `numchains > 1` always errored.
    chn2 = Octofitter.loadhdf5(fname, 2)
    @test size(chn2, 1) == 4
    @test size(chn2, 3) == 2
    @test vec(chn2[:b_a]) ≈ sma
    rm(fname)
end

@testset "savehdf5 → loadhdf5 round trip" begin
    model = io_test_model()
    chain = io_test_chain(model)
    fname = tempname() * ".hdf5"
    Octofitter.savehdf5(fname, model, chain)
    back = Octofitter.loadhdf5(fname)

    # Stored as Float32, so this is a single-precision comparison.
    for (oct, orb) in ((:b_a, :b_a), (:b_e, :b_e), (:b_i, :b_i),
                       (:b_ω, :b_ω), (:b_Ω, :b_Ω), (:plx, :plx))
        @test vec(collect(chain[oct])) ≈ vec(collect(back[orb])) rtol = 1e-6
    end
    @test vec(collect(back[:mtot])) ≈
          vec(collect(chain[:A_mass])) .+ vec(collect(chain[:b_mass])) rtol = 1e-6

    # `tp` only survives modulo the period: orbitize! stores the phase, not the
    # epoch. Same limitation v1 had.
    M = vec(collect(chain[:A_mass])) .+ vec(collect(chain[:b_mass]))
    a_in = vec(collect(chain[:b_a]))
    tp_in = vec(collect(chain[:b_tp]))
    tp_out = vec(collect(back[:b_tp]))
    P = @. √(a_in^3 / M) * KEPLER_F
    ϕ_in = @. mod((tp_in - 58849) / P, 1)
    ϕ_out = @. mod((tp_out - 58849) / P, 1)
    @test ϕ_in ≈ ϕ_out rtol = 1e-5

    # v1 wrote `col_names` but read `parameter_labels`, so its own export only
    # loaded through the fallback guess (with a warning).
    h5open(fname, "r") do f
        @test HDF5.read_attribute(f, "parameter_labels")[1] == "sma1"
    end
    rm(fname)
end

@testset "savehdf5: which bodies an orbit is about" begin
    # A Jacobi companion exports too: orbitize! parametrizes its own
    # multi-planet fits in Jacobi coordinates, so the total of the interior is
    # exactly what its `mtot` means.
    Aa = Octofitter.Body(name="Aa", variables=@variables begin mass = 1.0 end)
    Ab = Octofitter.Body(name="Ab", about=Aa, variables=@variables begin
        mass = 0.5
        a = 1.0; e = 0.0; i = 0.0; ω = 0.0; Ω = 0.0; tp = 55000.0
    end)
    c = Octofitter.Body(name="c", about=(Aa, Ab), variables=@variables begin
        mass = 1mjup
        a = 30.0; e = 0.1; i = 0.4; ω = 0.2; Ω = 0.3; tp = 55000.0
    end)
    sys = Octofitter.System(name="binary", bodies=[Aa, Ab, c], observations=(),
        variables=@variables begin plx = 20.0 end)
    @test Octofitter._interior_of(sys, :c) === (:Aa, :Ab)
    @test Octofitter._interior_of(sys, :Ab) === (:Aa,)
    @test_throws r"no orbit placing body" Octofitter._interior_of(sys, :Aa)
    @test Octofitter._first_companion(sys) === :Ab
end

# ---------------------------------------------------
# whereistheplanet.com
# ---------------------------------------------------

const WITP_DIR = joinpath(get(ENV, "DATADEPS_LOAD_PATH",
                              joinpath(homedir(), ".julia", "datadeps")),
                          "Whereistheplanet")

if isdir(WITP_DIR)
    @testset "whereistheplanet astrometry" begin
        catalog = WITP_DIR
        fname = Octofitter.Whereistheplanet_search("51erib", catalog)
        @test isfile(fname)

        A = Octofitter.Body(name="A", variables=@variables begin mass = 1.75 end)
        b = Octofitter.Body(name="b", about=A, variables=@variables begin
            mass = 3mjup
            a = 12.0; e = 0.4; i = 2.3; ω = 1.0; Ω = 1.0; tp = 55000.0
        end)
        obs = Octofitter.Whereistheplanet_astrom("51erib", catalog; target=b, ref=A)
        # v1 built this vector and never returned it.
        @test obs isa Vector{Octofitter.RelAstromObs}
        @test !isempty(obs)
        # Names must differ or the two cannot go into one System.
        @test allunique(Octofitter.likelihoodname.(obs))
        sys = Octofitter.System(name="51eri", bodies=[A, b], observations=Tuple(obs),
            variables=@variables begin plx = 33.4 end)
        @test length(sys.observations) == length(obs)
    end
else
    @info "Skipping whereistheplanet tests: the `Whereistheplanet` DataDep is not cached."
end

# ---------------------------------------------------
# Substellar flux tables
# ---------------------------------------------------

# These two exercise the parts of the table code that are pure bookkeeping, so
# they run whether or not the grids themselves are cached — which matters
# because the BHAC15 DataDep is not fetched in this environment and would
# otherwise leave that file with no coverage at all.

@testset "mass unit conversion" begin
    # The grids are tabulated in Jupiter masses; v2 body `mass` variables are
    # M⊙, so this is the whole of the v1→v2 calling-convention change.
    @test Octofitter._mass_to_mjup(:Mjup) == 1.0
    @test Octofitter._mass_to_mjup(:Msol) * PlanetOrbits.mjup ≈ 1.0
    @test 12mjup * Octofitter._mass_to_mjup(:Msol) ≈ 12.0
    @test Octofitter._mass_to_mjup(:Mearth) ≈ PlanetOrbits.mearth / PlanetOrbits.mjup
    @test_throws r"mass_unit" Octofitter._mass_to_mjup(:kg)
end

@testset "BHAC15 mass column lookup" begin
    # `M/Ms` survives CSV.jl's `normalizenames` differently across versions, so
    # the column is looked up by any of its plausible spellings rather than
    # hard-coded — a wrong guess would otherwise fail on the last line of a
    # two-minute load.
    for nm in (Symbol("M/Ms"), :M_Ms, :MMs, :M_Ms_)
        tbl = FlexTable(Teff=[3000.0, 3200.0])
        setproperty!(tbl, nm, [0.07, 0.09])
        @test Octofitter._bhac_masscol(tbl) == [0.07, 0.09]
    end
    @test_throws r"mass column" Octofitter._bhac_masscol(FlexTable(Teff=[3000.0]))
end

const SONORA_DIR = joinpath(get(ENV, "DATADEPS_LOAD_PATH",
                                joinpath(homedir(), ".julia", "datadeps")),
                            "SonoraBobcatEvoPhot")

if isdir(SONORA_DIR)
    # Built once: each of these grids a few thousand RBF samples onto a dense
    # rectangular mesh, which takes seconds.
    const SONORA_ABSMAG_L = Octofitter.sonora_photometry_interpolator(
        :Keck_L′; catalog=SONORA_DIR)
    const SONORA_COOLING = Octofitter.sonora_cooling_interpolator(catalog=SONORA_DIR)

    @testset "Sonora interpolators" begin
        m_msol = 12mjup
        @test isfinite(SONORA_ABSMAG_L(1200.0, m_msol))
        @test isfinite(SONORA_COOLING(15.0, m_msol))

        # The default input unit is M⊙, matching every other mass in v2. The
        # v1 calling convention is one keyword away and must agree.
        absmag_jup = Octofitter.sonora_photometry_interpolator(
            :Keck_L′; catalog=SONORA_DIR, mass_unit=:Mjup)
        @test absmag_jup(1200.0, 12.0) ≈ SONORA_ABSMAG_L(1200.0, m_msol)
        cool_jup = Octofitter.sonora_cooling_interpolator(
            catalog=SONORA_DIR, mass_unit=:Mjup)
        @test cool_jup(15.0, 12.0) ≈ SONORA_COOLING(15.0, m_msol)

        @test_throws r"mass_unit" Octofitter.sonora_photometry_interpolator(
            :Keck_L′; catalog=SONORA_DIR, mass_unit=:kg)

        # Out of grid gives NaN, and does so *in the input's own number type* —
        # a bare `Float64` NaN would make the return `Union{Float64,Dual}` and
        # poison inference through the whole `@variables` block.
        @test isnan(SONORA_ABSMAG_L(1e6, m_msol))
        @test SONORA_ABSMAG_L(1e6, m_msol) isa Float64

        # ForwardDiff has to make it through: this sits inside the sampler's
        # hot loop by way of a body's `flux_<band>` variable.
        d = ForwardDiff.derivative(m -> SONORA_ABSMAG_L(1200.0, m), m_msol)
        @test isfinite(d)
        @test d != 0
        dc = ForwardDiff.derivative(m -> SONORA_COOLING(15.0, m), m_msol)
        @test isfinite(dc)
    end

    @testset "Sonora feeds a flux_<band> body variable" begin
        # The mass-photometry bridge, spelled the way the tutorial spells it:
        # mass → temperature → absolute magnitude → *linear* flux. The
        # magnitude→flux step is not optional; photocentres weight bodies by
        # `flux_<band>` linearly.
        cooling = SONORA_COOLING
        absmag = SONORA_ABSMAG_L

        A = Octofitter.Body(name="A", variables=@variables begin
            mass = 1.0
            flux_L = 1.0
        end)
        b = Octofitter.Body(name="b", about=A, variables=@variables begin
            mass ~ LogUniform(5mjup, 30mjup)
            tempK = $cooling(system.age, mass)
            flux_L = 10^(-0.4 * $absmag(tempK, mass))
            a = 15.0; e = 0.1; i = 0.6; ω = 0.2; Ω = 0.4; tp = 55000.0
        end)
        astrom = Octofitter.RelAstromObs(
            Table(epoch=[57000.0], ra=[100.0], dec=[50.0], σ_ra=[1.0], σ_dec=[1.0]);
            target=Photocentre(:L), ref=A, name="phot")
        sys = Octofitter.System(name="massphot", bodies=[A, b], observations=(astrom,),
            variables=@variables begin
                plx = 25.0
                age = 20.0
            end)
        model = Octofitter.LogDensityModel(sys, verbosity=0)
        θ = model.arr2nt(collect(model.sample_priors(Random.Xoshiro(11))))
        @test isfinite(θ.bodies.b.tempK)
        @test θ.bodies.b.flux_L > 0
        @test isfinite(Octofitter.make_ln_like(sys)(sys, θ))
    end
else
    @info "Skipping Sonora tests: the `SonoraBobcatEvoPhot` DataDep is not cached."
end

const BHAC_DIR = joinpath(get(ENV, "DATADEPS_LOAD_PATH",
                              joinpath(homedir(), ".julia", "datadeps")),
                          "BHAC15_GAIA")

if isdir(BHAC_DIR)
    @testset "BHAC15 isochrone interpolator" begin
        fname = only(filter(isfile, readdir(BHAC_DIR, join=true)))
        itp = Octofitter.bhac15_mass_age_interpolator(fname; key=:Teff)
        @test isfinite(itp(15.0, 0.08))
        itpj = Octofitter.bhac15_mass_age_interpolator(fname; key=:Teff, mass_unit=:Mjup)
        @test itpj(15.0, 0.08 / mjup) ≈ itp(15.0, 0.08)
    end
else
    @info "Skipping BHAC15 tests: the `BHAC15_GAIA` DataDep is not cached (it is not " *
          "downloaded on demand here — see the io agent's report)."
end
