using Enzyme
Enzyme.API.maxtypeoffset!(4096)

using Octofitter, Distributions, CSV, PlanetOrbits, Profile

df = CSV.read(joinpath(@__DIR__, "..", "docs", "src", "target_1.csv"), FlexTable)
ref_epoch_mjd = 57936.375
gaia_dr4_obs = GaiaDR4AstromObs(df, gaia_id=4373465352415301632,
    variables=@variables begin
        astrometric_jitter ~ LogUniform(0.00001, 10)
        ra_offset_mas ~ Normal(0, 10000); dec_offset_mas ~ Normal(0, 10000)
        pmra ~ Uniform(-1000, 1000); pmdec ~ Uniform(-1000, 1000)
        plx = system.plx; ref_epoch = $ref_epoch_mjd
    end)
orbit_ref_epoch = mean(gaia_dr4_obs.table.epoch)
b = Planet(name="b", basis=Visual{KepOrbit}, observations=[],
    solver=PlanetOrbits.Markley(),
    variables=@variables begin
        a ~ LogUniform(0.01, 100); e ~ Uniform(0, 0.99)
        ω ~ Uniform(0,2pi); i ~ Sine(); Ω ~ Uniform(0,2pi); θ ~ Uniform(0,2pi)
        tp = θ_at_epoch_to_tperi(θ, $orbit_ref_epoch; M=system.M, e, a, i, ω, Ω)
        mass ~ LogUniform(0.01, 1000)
    end)
sys = System(name="target_1", companions=[b], observations=[gaia_dr4_obs],
    variables=@variables begin; M = 1.0; plx ~ Uniform(0.01,100); end)
model = Octofitter.LogDensityModel(sys, verbosity=4)
t = [-1.6495356151026312, -6.491790601708638, 0.029738626170779867,
     0.0832180348920664, 0.010814871781942643, -0.04826697440973689,
     -3.357219015634567, 0.1154498724014276, -0.012059033319817353,
     -0.4030642907378578, -2.391709856749899, 0.21726752142249556, -4.866247384968242]
model.∇ℓπcallback(t)

println("@allocated: $((@allocated model.∇ℓπcallback(t))) bytes")

# Scaling test
function test_scaling()
    df_full = CSV.read(joinpath(@__DIR__, "..", "docs", "src", "target_1.csv"), FlexTable)
    for nrows in [10, 20, 50, 102]
        df_sub = df_full[1:nrows]
        obs = GaiaDR4AstromObs(df_sub, gaia_id=4373465352415301632,
            variables=@variables begin
                astrometric_jitter ~ LogUniform(0.00001, 10)
                ra_offset_mas ~ Normal(0, 10000); dec_offset_mas ~ Normal(0, 10000)
                pmra ~ Uniform(-1000, 1000); pmdec ~ Uniform(-1000, 1000)
                plx = system.plx; ref_epoch = $ref_epoch_mjd
            end)
        ore = mean(obs.table.epoch)
        b2 = Planet(name="b", basis=Visual{KepOrbit}, observations=[],
            solver=PlanetOrbits.Markley(),
            variables=@variables begin
                a ~ LogUniform(0.01, 100); e ~ Uniform(0, 0.99)
                ω ~ Uniform(0,2pi); i ~ Sine(); Ω ~ Uniform(0,2pi); θ ~ Uniform(0,2pi)
                tp = θ_at_epoch_to_tperi(θ, $ore; M=system.M, e, a, i, ω, Ω)
                mass ~ LogUniform(0.01, 1000)
            end)
        s = System(name="t", companions=[b2], observations=[obs],
            variables=@variables begin; M = 1.0; plx ~ Uniform(0.01,100); end)
        m = Octofitter.LogDensityModel(s, verbosity=0)
        m.∇ℓπcallback(t)
        a = @allocated m.∇ℓπcallback(t)
        println("  $nrows rows: $a bytes")
    end
end
println("\n=== Scaling ===")
test_scaling()

# Allocation profiler
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1.0 begin
    for _ in 1:10
        model.∇ℓπcallback(t)
    end
end
results = Profile.Allocs.fetch()
allocs = results.allocs
println("\n=== Profiler: $(length(allocs)÷10) allocs/call ===")

type_counts = Dict{String, Tuple{Int,Int}}()
for alloc in allocs
    tname = string(alloc.type)
    tname = length(tname) > 80 ? tname[1:77] * "..." : tname
    prev = get(type_counts, tname, (0, 0))
    type_counts[tname] = (prev[1] + 1, prev[2] + alloc.size)
end
sorted = sort(collect(type_counts), by=x->x[2][2], rev=true)
for (tname, (count, total_sz)) in sorted
    pc, pb = count÷10, total_sz÷10
    pb < 8 && continue
    println("  $(rpad(tname,80)) cnt/call=$(lpad(pc,3))  bytes/call=$(lpad(pb,6))")
end
