using Octofitter
using Distributions
using CairoMakie

# gaialike = Octofitter.GaiaDR3_v2(;gaia_id=6412595290592307840)
# gaialike = Octofitter.GaiaDivergenceLikelihood_3(;source_id_dr3=2306965202564744064, source_id_dr2=2306965202564506752)
# gaialike = Octofitter.GaiaDivergenceLikelihood_3(;source_id_dr3=2859970279471313664 )
gaialike = Octofitter.GaiaDivergenceLikelihood_3(;
    source_id_dr2=756291174721509376,
    source_id_dr3=756291174721509376,
    scanlaw_table=(
        times = [1752.939291486925, 1780.9299890574835, 1805.1738582301318,
          1947.4970887103518, 2002.2357210388893, 2160.219305272992,
          2179.968745028478, 2274.074257403782, 2274.324422622525,
          2284.5790943458715, 2360.027254399675, 2377.2837472550445,
          2456.551534753258, 2456.8017012214445, 2457.0518657347043,
          2457.302028538316, 2457.552189675778, 2462.305141876399,
          2462.5553033801257, 2462.8054667692613, 2463.0556322547227,
          2504.786485959369, 2539.5082638525732, 2554.261714037261,
          2647.9225728711995, 2664.6757757563023, 2699.9006412338667,
          1753.0132981737922, 1781.0039954466986, 1804.9977018284137,
          1947.571095272144, 2002.3097276801418, 2128.8119324926156,
          2160.2933116484633, 2274.1482652277123, 2274.3984287092217,
          2284.402938355658, 2284.65310160712, 2323.6188693930517,
          2359.8510959981436, 2360.101260944946, 2377.3577542716034,
          2456.375374735121, 2456.6255425302807, 2456.875708411386,
          2457.125872414931, 2457.37603471282, 2462.379147802961,
          2462.629309842445, 2462.8794738507027, 2463.1296399802663,
          2554.0855575107107, 2647.996579910515, 2664.4996192808903,
          2699.974647725961],
        angles = [-101.46008591416796, 159.97589739485045, -126.77752867716161,
          -13.384254386965052, -28.14169154911564, 154.3070339279813,
          -137.86883377094713, 36.73534292462481, 37.99945908961216,
          78.30811815393223, 43.80355321521394, -22.957292362328722,
          -122.43849042222561, -123.38766456474876, -124.33993517584,
          -125.29599071429305, -126.25597256359735, -145.81780532501,
          -146.94447010571182, -148.08333716410172, -149.23468715310455,
          -107.22862650445674, 151.34686022831323, -152.17647703519282,
          18.250756771981205, 78.2819873497146, -23.88729914011067,
          -101.60070471917078, 159.69173492419765, -127.40316006545035,
          -13.101615679369662, -28.23541332400306, -103.58896586134013,
          154.15066122402698, 37.111067110037084, 38.37022837508607,
          77.7862066533938, 78.52559009291367, -19.306679731630958,
          44.63152128601726, 43.45303395836092, -23.129793635393447,
          -121.77189960544872, -122.71899126644449, -123.66899673686468,
          -124.62242228315536, -125.5795074936485, -146.14984716155917,
          -147.28010028016888, -148.42266910267136, -149.57767532854038,
          -153.11917615352792, 18.64190147862433, 78.02616348139668,
          -23.744671384538258
        ],
    )
)
##
out = Octofitter.derive_IAD!(gaialike)
##
@planet b AbsoluteVisual{KepOrbit} begin
    mass  = 0. #~ Uniform(0, 1000)
    e = 0.#        0.2209
    ω = 0.#        4.3257
    a = 1. #~ Uniform(0,10)#       27.4579
    i = 0. #~ Sine()#        1.8528
    Ω = 0. #~ Uniform(0,2pi)#         4.1001
    # θ = 0. #~ UniformCircular()
    tp = 0.0#θ_at_epoch_to_tperi(system,b,Octofitter.meta_gaia_DR2.ref_epoch_mjd) 
    # τ =         0.0334
    # P =       162.4281
    # tp =     62107.4063
end
@system sys begin
    M = 0.7847 # Host mass not important for this example
    rv = 0.0 # system RV not significant for this example
    plx ~ truncated(Normal(30,5),lower=1,upper=50)
    pmra ~ Uniform(-1000, 1000)
    pmdec ~ Uniform(-1000, 1000)


    # It is convenient to put a prior of the catalog value +- 10,000 mas on position
    # gaia_ra_offset_mas ~  Normal(0, 10000000)
    # gaia_dec_offset_mas ~ Normal(0, 10000000)
    # dec = gaialike.dr3.dec + system.gaia_ra_offset_mas/60/60/1000
    # ra = gaialike.dr3.ra + system.gaia_dec_offset_mas/60/60/1000/cosd(gaialike.dr3.dec)
    gaia_ra_offset_deg ~  Normal(0, 0.1)
    gaia_dec_offset_deg ~ Normal(0, 0.1)
    dec = gaialike.dr3.dec + system.gaia_ra_offset_deg
    ra = gaialike.dr3.ra + system.gaia_dec_offset_deg
    # dec = gaialike.dr2.dec + system.gaia_ra_offset_deg
    # ra = gaialike.dr2.ra + system.gaia_dec_offset_deg

    along_scan_uncertainty_mas = 0.4#~ truncated(Normal(0,1),lower=0.0001)
    # along_scan_uncertainty_mas = 21.

    # ref_epoch = 57388.5
    ref_epoch = Octofitter.meta_gaia_DR3.ref_epoch_mjd
    # ref_epoch = Octofitter.meta_gaia_DR2.ref_epoch_mjd

    N ~ MvNormal([1,2,3,4])

end b
using FiniteDiff
model = Octofitter.LogDensityModel(sys,autodiff=:ForwardDiff, verbosity=4)

##
model.starting_points=fill(model.link([
    gaialike.dr3.parallax,
    gaialike.dr3.pmra,
    gaialike.dr3.pmdec,
    0.0,
    0.0,
    # 10.0,
    # 10000
]),60)
model.ℓπcallback(model.starting_points[1])
##
solnt = Octofitter.default_initializer!(model, verbosity=4, nruns=64)

# solnt =  (;logpost=0,model.arr2nt(model.invlink(model.starting_points[1]))...,)
chain = Octofitter.result2mcmcchain([solnt],  Dict(:internals => [:logpost]))

##
chain = octofit(model, max_depth=5)
##
using Pigeons
chain, pt = octofit_pigeons(model,n_rounds=5)
##



using DataFrames
samples_dr3 = DataFrame(
    rand(gaialike.dist_dr3, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
samples_dr3.ra .-= gaialike.dr3.ra
samples_dr3.dec .-= gaialike.dr3.dec
post_dr3 = DataFrame(
    rand(out.P_dr3, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
post_dr3.ra .-= gaialike.dr3.ra
post_dr3.dec .-= gaialike.dr3.dec

samples_dr2 = DataFrame(
    rand(gaialike.dist_dr2, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
samples_dr2.ra .-= gaialike.dr2.ra
samples_dr2.dec .-= gaialike.dr2.dec
post_dr2 = DataFrame(
    rand(out.P_dr2, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
post_dr2.ra .-= gaialike.dr2.ra
post_dr2.dec .-= gaialike.dr2.dec


# θ_systems_from_chain = Octofitter.mcmcchain2result(model, chain)
# i = argmax(chain[:logpost][:])
# i = rand(1:size(chain,1))
# orbits = map(keys(model.system.planets)) do planet_key
#     Octofitter.construct_elements(chain, planet_key, i)
# end
# P = Octofitter.simulate(gaialike, θ_systems_from_chain[i], orbits)


post = (;
    parallax=chain[:plx][:],
    ra=chain[:ra][:],
    dec=chain[:dec][:],
    pmra=chain[:pmra][:],
    pmdec=chain[:pmdec][:],
)
post.ra .-= gaialike.dr3.ra
post.dec .-= gaialike.dr3.dec
# post.ra .-= gaialike.dr2.ra
# post.dec .-= gaialike.dr2.dec


fig = Figure(size=(700,700))

pairplot(
    fig[1:5,1:5],
    PairPlots.Series(samples_dr3,label="DR3",color=:blue)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(post_dr3,label="DR3 Model",color=:blue,linestyle=:dash)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(samples_dr2,label="DR2",color=:red)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5,),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(post_dr2,label="DR2 Model",color=:red,linestyle=:dash)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(post, label="Second Model", color=:black,linestyle=:dash)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
)
fig
##
# sol = Octofitter._refine(model, model.starting_points[1])#get_sample(pt)[end][1:end-1])
sol = Octofitter._refine(model, sol)#get_sample(pt)[end][1:end-1])
@show model.ℓπcallback(sol)
##
solnt =  (;logpost=0,model.arr2nt(model.invlink(sol))...,)
chain = Octofitter.result2mcmcchain([solnt],  Dict(:internals => [:logpost]))
##
let chain = chain[10:end,:,:]
    using LinearAlgebra, StatsBase
    fig = Figure(size=(1080,720))
    j = i = 1
    for prop in (
        (;chain=:ra, gaia=:ra, gaia_err=:ra_error), 
        (;chain=:dec, gaia=:dec, gaia_err=:dec_error),
        (;chain=:plx, gaia=:parallax, gaia_err=:parallax_error), 
        (;chain=:pmra, gaia=:pmra, gaia_err=:pmra_error), 
        (;chain=:pmdec, gaia=:pmdec, gaia_err=:pmdec_error)
    )
        # i, j, ax
        ax = Axis(
            fig[j,i],
            xlabel=string(prop.chain),
        )
        i+=1
        if i > 3
            j+=1
            i = 1
        end
        unc = gaialike.dr3[prop.gaia_err]
        if prop.chain == :ra
            unc /= 60*60*1000 * cosd(gaialike.dr3.dec)
        end
        if prop.chain == :dec
            unc /= 60*60*1000
        end
        if prop.gaia == :zero
            n = Normal(0, unc)
        else
            mu = gaialike.dr3[prop.gaia]
            n = Normal(mu, unc)
        end
        n0,n1=quantile.(n,(1e-4, 1-1e-4))
        nxs = range(n0,n1,length=200)
        h = fit(Histogram, chain[prop.chain][:], nbins=15)
        h = normalize(h, mode=:pdf)
        barplot!(ax, (h.edges[1][1:end-1] .+ h.edges[1][2:end])./2, h.weights, gap=0, color=:red, label="posterior")
        lines!(ax, nxs, pdf.(n,nxs), label="Gaia Catalog", color=:black, linewidth=2)
    end

    ax = Axis(
        fig[j,i],
        xlabel="mass",
    )
    h = fit(Histogram, chain[:b_mass][:], nbins=15)
    h = normalize(h, mode=:pdf)
    barplot!(ax, (h.edges[1][1:end-1] .+ h.edges[1][2:end])./2, h.weights, gap=0, color=:red, label="posterior")

    Legend(fig[i-1,j+1],ax,tellwidth=false)
    fig

end
##

function optimtest(args,gaialike)
    (
        plx,
        ra_off,
        dec_off,
        # rv,
        pmra,
        pmdec,
        # mass,
        unc_mas
    ) = args
    # ra += 1
    # pmra *= 2
    ra = gaialike.dr3.ra + ra_off
    dec = gaialike.dr3.dec + dec_off

    outer_loop_skypath = AbsoluteVisual{KepOrbit}(;plx,ra,dec,rv=0,pmra,pmdec,
        ref_epoch=57388.5,
        M=1,a=1,e=0,i=0,ω=0,Ω=0,tp=0
    )
    mass=zero(typeof(plx))  + 100
    ll = Octofitter.ln_like(gaialike, (;typical_along_scan_uncertainty_mas=abs(unc_mas),planets=(;b=(;mass=abs(mass)))), (outer_loop_skypath,))
    @show ll
    return -ll
end
catalog_best = [
    gaialike.dr3.parallax,
    0.01,#gaialike.dr3.ra,
    0.0,#gaialike.dr3.dec,
    gaialike.dr3.pmra,
    gaialike.dr3.pmdec,
    # 1.0,
    0.35
]
guess = catalog_best
optimtest(guess,gaialike)


##
using Optimization, OptimizationOptimJL
func = OptimizationFunction(optimtest)# AutoForwardDiff())
prob = OptimizationProblem(func, guess, gaialike)
# @time sol = solve(prob, LBFGS(
#     m=4,
#     linesearch=OptimizationOptimJL.Optim.LineSearches.BackTracking(),
#     alphaguess=OptimizationOptimJL.Optim.LineSearches.InitialHagerZhang()
# ), g_tol=1e-9);
@time sol = solve(prob, NelderMead(), g_tol=1e-8, iterations=10_000,)
## Compare the simulated and output models
catalog_path = AbsoluteVisual{KepOrbit}(;
    plx=gaialike.dr3.parallax,
    ref_epoch=Octofitter.meta_gaia_DR2.ref_epoch_mjd,
    ra=gaialike.dr3.ra,
    dec=gaialike.dr3.dec,
    rv=0e3,
    pmra=gaialike.dr3.pmra,
    pmdec=gaialike.dr3.pmdec,
    M=1,a=1,e=0,i=0,ω=0,Ω=0,tp=0
)
x,y,αₘ,δₘ=Octofitter._simulate_skypath_observations(gaialike, catalog_path, 0*Octofitter.mjup2msol)
fig = Figure()
ax = Axis(
    fig[1,1,],
    autolimitaspect=1
)
scatterlines!(ax, x .- gaialike.dr3.ra, y .-gaialike.dr3.dec, label="DR2 catalog values", linewidth=4)

(
    plx,
    ra_off,
    dec_off,
    # rv,
    pmra,
    pmdec,
    # mass,
    unc_mas
) = sol.u
ra = gaialike.dr3.ra + ra_off
dec = gaialike.dr3.dec+ dec_off
args = (
    plx,
    ra,
    dec,
    # rv,
    pmra,
    pmdec,
    # mass,
    unc_mas
)
# args = (
#     gaialike.dr3.parallax,
#     gaialike.dr3.ra,
#     gaialike.dr3.dec,
#     # rv,
#     gaialike.dr3.pmra,
#     gaialike.dr3.pmdec,
#     # mass,
#     unc_mas
# )
x2,y2,resid=Octofitter._simulate_gaia_5_param_model(args,(;
    αₘ,δₘ,gaialike.table,ref_epoch=Octofitter.meta_gaia_DR2.ref_epoch_mjd,along_scan_uncertainty_mas=0,
    full3D=true,
    catalog_rv_m_s=gaialike.dr3.radial_velocity*1e3
))
# std(resid.*60 .*60 .*1000)

scatterlines!(ax, x2 .- gaialike.dr3.ra, y2 .-gaialike.dr3.dec,label="best fit output?")
axislegend(ax)
fig
##
using FiniteDiff
H = FiniteDiff.finite_difference_hessian(u->optimtest(u,gaialike), sol.u)
sqrt_or_inf(num) = num >=0 ? sqrt(num) : Inf
unc = sqrt_or_inf.(diag(Diagonal(inv(H))))
catalog_best = [
    gaialike.dr3.parallax,
    gaialike.dr3.ra,
    gaialike.dr3.dec,
    gaialike.dr3.pmra,
    gaialike.dr3.pmdec,
    # 1.0,
    1.0
]

catalog_error = [
    gaialike.dr3.parallax_error,
    gaialike.dr3.ra_error,
    gaialike.dr3.dec_error,
    gaialike.dr3.pmra_error,
    gaialike.dr3.pmdec_error
]

##
println("                        plx             ra              dec              pmra           pmdec           unc_mas")
println("Catalog values:         ", join(round.(catalog_best[1:6],sigdigits=8),"\t"))
println("Posterior values:       ", join(round.(sol.u[1:6],sigdigits=8),"\t"))
println("Catalog uncertainties:  ", join(round.(catalog_error,sigdigits=8),"\t"))
println("Posterior uncertainties:", join(round.(unc,sigdigits=8),"\t"))
println("")
