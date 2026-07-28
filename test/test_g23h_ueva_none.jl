#=
`ueva_mode = :none` — opting a source out of the UEVA channel.

Some G23H sources have no sig_AL / sig_att_radec / sig_cal calibration (in the
published catalog these are predominantly the very brightest stars, where the
calibration was not extrapolated). Their σ priors are
`truncated(Normal(sig_X, sig_X_sigma))`, so a NaN in either half makes G23HObs
unconstructable. `:none` drops the UEVA datum, replaces those three sampled
nuisances with fixed constants, and disables the UEVA-driven DR3 deflation.

Uses the same target as test_g23h_simulation.jl so no extra data is needed.
=#

using Octofitter, Distributions, Test

@testset "G23H ueva_mode=:none" begin
    ruwe = G23HObs(; hip_id=18512, freeze_epochs=true, include_rv=false,
                     ueva_mode=:RUWE)
    none = G23HObs(; hip_id=18512, freeze_epochs=true, include_rv=false,
                     ueva_mode=:none)

    @testset "observation table" begin
        @test :ueva_dr3 ∈ ruwe.table.kind
        @test :ueva_dr3 ∉ none.table.kind
        # exactly one row removed, and only that row
        @test length(none.table.kind) == length(ruwe.table.kind) - 1
        @test collect(none.table.kind) ==
              filter(!=(:ueva_dr3), collect(ruwe.table.kind))
    end

    @testset "sigma nuisances become constants" begin
        for k in (:σ_AL, :σ_att, :σ_calib)
            @test k ∈ keys(ruwe.priors.priors)      # sampled under :RUWE
            @test k ∉ keys(none.priors.priors)      # not sampled under :none
            @test k ∈ keys(none.derived.variables)  # fixed instead
        end
    end

    @testset "input validation" begin
        @test_throws Exception G23HObs(; hip_id=18512, ueva_mode=:nonsense)
    end
end
