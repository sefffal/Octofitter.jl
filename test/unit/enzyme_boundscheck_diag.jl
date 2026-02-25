# Diagnostic tests to isolate Enzyme + @boundscheck interaction
# Run both ways to compare:
#   julia --project=. test/unit/enzyme_boundscheck_diag.jl
#   julia --project=. --check-bounds=yes test/unit/enzyme_boundscheck_diag.jl

# Stack the test environment so OctofitterRadialVelocity (test-only dep)
# is available alongside Enzyme/Octofitter (main project deps)
let test_env = joinpath(@__DIR__, "..")
    if isfile(joinpath(test_env, "Project.toml")) && !(test_env in LOAD_PATH)
        push!(LOAD_PATH, test_env)
    end
end

using Test
using Enzyme
using Octofitter
using OctofitterRadialVelocity
using Distributions
using TypedTables
using Random
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

cb = Base.JLOptions().check_bounds
mode = cb == 0 ? "auto" : cb == 1 ? "yes" : "no"
@info "Running with --check-bounds=$mode"

@testset "Enzyme + @boundscheck diagnostics" begin

    # Test 1: Can Enzyme differentiate through length()?
    @testset "length() under Enzyme" begin
        function f_length(x::Vector{Float64})
            n = length(x)
            return sum(x) / n
        end
        x = [1.0, 2.0, 3.0]
        dx = zeros(3)
        Enzyme.autodiff(Reverse, f_length, Active, Duplicated(x, dx))
        @test dx ≈ [1/3, 1/3, 1/3]
        @test length(x) == 3
    end

    # Test 2: Does @boundscheck + length() return the right value?
    @testset "@boundscheck length() in regular function" begin
        function f_boundscheck(x::Vector{Float64})
            @boundscheck if length(x) != 3
                error("Expected 3, got $(length(x))")
            end
            @inbounds return x[1] + x[2] + x[3]
        end
        x = [1.0, 2.0, 3.0]
        dx = zeros(3)
        Enzyme.autodiff(Reverse, f_boundscheck, Active, Duplicated(x, dx))
        @test dx ≈ [1.0, 1.0, 1.0]
    end

    # Test 3: Same but with a RuntimeGeneratedFunction (matching arr2nt pattern)
    @testset "@boundscheck in RuntimeGeneratedFunction" begin
        rgf = @RuntimeGeneratedFunction(:(function (arr)
            l = 3
            @boundscheck if length(arr) != l
                error("Expected exactly $l elements in array (got $(length(arr)))")
            end
            @inbounds return arr[1] + arr[2] + arr[3]
        end))

        # Primal should work
        @test rgf([1.0, 2.0, 3.0]) == 6.0

        # Now differentiate
        x = [1.0, 2.0, 3.0]
        dx = zeros(3)
        Enzyme.autodiff(Reverse, rgf, Active, Duplicated(x, dx))
        @test dx ≈ [1.0, 1.0, 1.0]
    end

    # Test 4: Is it the shadow array's length() that's wrong?
    @testset "shadow array length sanity" begin
        observed_lengths = Int[]
        function f_capture_length(x::Vector{Float64})
            push!(observed_lengths, length(x))
            return sum(x)
        end
        x = [1.0, 2.0, 3.0]
        dx = zeros(3)
        empty!(observed_lengths)
        Enzyme.autodiff(Reverse, f_capture_length, Active, Duplicated(x, dx))
        @test all(==(3), observed_lengths)
    end

    # Test 5: Actual Octofitter arr2nt under Enzyme
    # This tests the real code path that fails in CI
    @testset "Octofitter arr2nt under Enzyme" begin
        rv_table = TypedTables.Table(
            epoch=[50000.0, 50100.0, 50200.0],
            rv=[100.0, 110.0, 105.0],
            σ_rv=[10.0, 10.0, 10.0],
        )
        rvlike = MarginalizedStarAbsoluteRVObs(rv_table, name="TestRV")

        b = Planet(
            name="b",
            basis=RadialVelocityOrbit,
            observations=[],
            variables=@variables begin
                e ~ Uniform(0.0, 0.5)
                ω ~ UniformCircular()
                τ ~ Uniform(0, 1)
                P ~ Uniform(0.001, 10)
                a = ∛(P^2 * system.M)
                tp = τ * P * 365.25 + 50000
                mass ~ Uniform(0, 100)
            end
        )

        sys = System(
            name="DiagSys",
            companions=[b],
            observations=[rvlike],
            variables=@variables begin
                M ~ truncated(Normal(1.2, 0.1), lower=0.1)
            end
        )

        # Test 5a: Does arr2nt work on primal path?
        @testset "arr2nt primal" begin
            arr2nt = Octofitter.make_arr2nt(sys)
            θ = collect(Octofitter.sample_priors(sys))
            result = arr2nt(θ)
            @test result isa NamedTuple
            @info "arr2nt primal OK, got $(length(θ))-element -> $(typeof(result))"
        end

        # Test 5b: Does arr2nt work under Enzyme directly?
        @testset "arr2nt under Enzyme.autodiff" begin
            arr2nt = Octofitter.make_arr2nt(sys)
            θ = collect(Octofitter.sample_priors(sys))
            n = length(θ)

            # Wrap in a scalar function so Enzyme can differentiate
            function f_arr2nt(x::Vector{Float64})
                nt = arr2nt(x)
                # Just sum all values to get a scalar
                return Float64(sum(values(nt.planets.b)))
            end

            # Primal check
            val = f_arr2nt(θ)
            @test isfinite(val)
            @info "f_arr2nt primal = $val"

            # Enzyme gradient
            dx = zeros(n)
            try
                Enzyme.autodiff(Reverse, f_arr2nt, Active, Duplicated(θ, dx))
                @test any(!=(0.0), dx)
                @info "Enzyme gradient succeeded: $(dx)"
            catch e
                @test false  # mark failure
                @error "Enzyme autodiff of arr2nt failed" exception=(e, catch_backtrace())
            end
        end

        # Test 5c: Full LogDensityModel gradient (the actual CI failure path)
        @testset "LogDensityModel gradient" begin
            try
                model = Octofitter.LogDensityModel(sys, verbosity=0)
                θ = collect(model.sample_priors(Random.Xoshiro(42)))
                θt = model.link(θ)
                ll = model.ℓπcallback(θt)
                @info "Primal log-density = $ll"

                if isfinite(ll)
                    ll2, grad = model.∇ℓπcallback(θt)
                    @test isfinite(ll2)
                    @test all(isfinite, grad)
                    @info "Gradient succeeded: $(grad)"
                else
                    # Try more starting points
                    for attempt in 1:20
                        θ = collect(model.sample_priors(Random.Xoshiro(42 + attempt)))
                        θt = model.link(θ)
                        ll = model.ℓπcallback(θt)
                        if isfinite(ll)
                            ll2, grad = model.∇ℓπcallback(θt)
                            @test isfinite(ll2)
                            @test all(isfinite, grad)
                            @info "Gradient succeeded on attempt $attempt"
                            break
                        end
                    end
                end
            catch e
                @test false
                @error "LogDensityModel gradient failed" exception=(e, catch_backtrace())
            end
        end
    end
end
