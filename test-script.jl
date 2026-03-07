using Revise

using Octofitter, Distributions, DifferentiationInterface, BenchmarkTools
astrom_dat_1 = Table(;
                    epoch= [50000,  50120, 50240, 50360,50480, 50600, 50720, 50840,], # MJD (days)
                    ra   = [-505.764, -502.57, -498.209, -492.678,-485.977, -478.11, -469.08, -458.896,], # mas
                    dec  = [-66.9298, -37.4722, -7.92755, 21.6356, 51.1472,  80.5359,  109.729,  138.651, ], # mas
                    # Tip! Type this as \sigma + <TAB key>!
                    σ_ra = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, ],  # mas
                    σ_dec = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, ], # mas
                    cor =  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ]
                )
astrom_obs_1 = PlanetRelAstromObs(astrom_dat_1, name="relastrom")
planet_1 = Planet(
                    name="b",
                    basis=Visual{KepOrbit},
                    observations=[astrom_obs_1],
                    variables=@variables begin
                        plx = system.plx
                        M ~ truncated(Normal(1.2, 0.1), lower=0.1)
                        a ~ Uniform(0, 100)
                        e ~ Uniform(0.0, 0.5)
                        i ~ Sine()
                        ω ~ UniformCircular()
                        Ω ~ UniformCircular()
                        θ ~ UniformCircular()
                        tp = θ_at_epoch_to_tperi(θ, 50420; M, e, a, i, ω, Ω)
                        mass ~ LogUniform(0.1, 1000)
                    end
                )
sys = System(
                    name = "test",
                    companions=[planet_1],
                    variables=@variables begin
                        plx ~ Uniform(1, 1000)
                    end
                )
using Enzyme
model = Octofitter.LogDensityModel(sys,autodiff=AutoEnzyme(function_annotation=Enzyme.Const))
 t = [
          -1.8560955014065845
           0.10032287862639093
           0.35221836650994237
           0.22851762104412168
          -0.06018231374718679
           7.377580190515705e-5
          -0.9805758390595423
          -0.9559911808277645
          -0.2181967437562318
           0.07271812577956628
          -0.9778757788661318
          -4.605170581358881
         ]

using Profile

  # Warm up
  model.∇ℓπcallback(t)
  model.∇ℓπcallback(t)

  # Profile allocations
  Profile.Allocs.clear()
  Profile.Allocs.@profile sample_rate=1.0 model.∇ℓπcallback(t)
  results = Profile.Allocs.fetch()

  # Print allocation sites
  for alloc in results.allocs
      println("$(alloc.type) $(alloc.size) bytes")
      for frame in alloc.stacktrace[1:min(3, length(alloc.stacktrace))]
          println("  ", frame)
      end
      println()
  end

