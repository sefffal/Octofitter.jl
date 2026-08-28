using Random

# Regression test for issue #129.
#
# `guess_starting_position` used to size its per-thread buffers with
# `Threads.nthreads()` but index them with `Threads.threadid()`. Thread IDs run up
# to `Threads.maxthreadid()`, which exceeds `nthreads()` whenever an interactive
# threadpool is present, so the threaded branch threw a `BoundsError`.
#
# This has to run in a subprocess: the failure only appears when the default
# threadpool is offset by an interactive thread (`--threads=2,1` puts the two
# default threads on IDs 2 and 3), and the test suite itself normally runs
# single-threaded.
@testset "guess_starting_position with an interactive threadpool" begin
    script = """
    using Octofitter, Distributions, Random

    # `guess_starting_position` only takes the threaded branch when the Kepler
    # solver is not itself using threads. This is the default, but be explicit.
    Octofitter._kepsolve_use_threads[] = false

    @assert Threads.maxthreadid() > Threads.nthreads() "test subprocess needs an interactive threadpool"

    b = Planet(
        name="b",
        basis=Visual{KepOrbit},
        observations=[],
        variables=@variables begin
            a ~ Uniform(1, 20)
            e ~ Uniform(0.0, 0.5)
            i ~ Sine()
            ω ~ Uniform(0, 2pi)
            Ω ~ Uniform(0, 2pi)
            θ ~ Uniform(0, 2pi)
            tp = θ_at_epoch_to_tperi(θ, 50000; M=system.M, e, a, i, ω, Ω)
        end
    )
    sys = System(
        name="Issue129TestSys",
        companions=[b],
        observations=[],
        variables=@variables begin
            M ~ truncated(Normal(1.0, 0.1), lower=0.1)
            plx ~ truncated(Normal(50.0, 1.0), lower=0.1)
        end
    )
    model = Octofitter.LogDensityModel(sys; verbosity=0)

    # Cover both fewer draws than chunks and more, to exercise the partitioning.
    for N in (1, 3, 200)
        params, logpost = Octofitter.guess_starting_position(Xoshiro(1), model, N)
        @assert isfinite(logpost) "logpost was not finite for N=\$N"
        # The returned parameters must be the ones that produced the returned
        # logpost -- this is what a torn read across threads would break.
        @assert model.ℓπcallback(model.link(params)) == logpost "returned params do not match returned logpost for N=\$N"
    end
    """
    cmd = `$(Base.julia_cmd()) --threads=2,1 --project=$(Base.active_project()) --startup-file=no -e $script`
    @test success(pipeline(cmd; stdout=stdout, stderr=stderr))
end
