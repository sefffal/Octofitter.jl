using Octofitter
using Distributions
using CairoMakie
using FillArrays
using OctofitterRadialVelocity
cd(@__DIR__)
# Just for catalog access for now.
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
astrom_like = PlanetRelAstromLikelihood(
    (epoch=mjd("2016-12-15"), ra=133., dec=-174., σ_ra=07.0, σ_dec=07., cor=0.2),
    (epoch=mjd("2017-03-12"), ra=126., dec=-176., σ_ra=04.0, σ_dec=04., cor=0.3),
    (epoch=mjd("2017-03-13"), ra=127., dec=-172., σ_ra=04.0, σ_dec=04., cor=0.1),
    (epoch=mjd("2018-02-08"), ra=083., dec=-133., σ_ra=10.0, σ_dec=10., cor=0.4),
    (epoch=mjd("2018-11-28"), ra=058., dec=-122., σ_ra=10.0, σ_dec=20., cor=0.3),
    (epoch=mjd("2018-12-15"), ra=056., dec=-104., σ_ra=08.0, σ_dec=08., cor=0.2),
)
scatter(astrom_like.table.ra, astrom_like.table.dec, axis=(;autolimitaspect=1,xreversed=true))
##
# hiplike = HipparcosIADLikelihood(;hip_id=51658)

rvlike = PlanetRelativeRVLikelihood(
    (epoch=mjd("2008-05-01"), rv=1300, σ_rv=150, inst_idx=1),
    (epoch=mjd("2010-02-15"), rv=700, σ_rv=150, inst_idx=1),
    (epoch=mjd("2016-03-01"), rv=-2700, σ_rv=150, inst_idx=1),
)

@planet b AbsoluteVisual{KepOrbit} begin
    # a ~ LogUniform(0.1,400)
    # e ~ Uniform(0,0.999)
    # ω ~ UniformCircular()
    # i ~ Sine()
    # Ω ~ UniformCircular()

    # mass = system.M_sec

    # θ ~ UniformCircular()
    # tp = θ_at_epoch_to_tperi(system,b,57737.0) # epoch of astrometry

    # jitter ~ truncated(Normal(0,400),lower=0)


    a  = 1.0
    e  = 0.0
    ω  = 0.0
    i  = 0.0
    Ω  = 0.0
    tp  = 0.0
    jitter  = 0.0
    mass = 0.0

end rvlike astrom_like # Note the relative astrometry added here!

@system sys begin
    # M = 0.7847 # Host mass not important for this example
    M_pri ~ truncated(Normal(1.61, 0.1), lower=0)
    # M_sec ~ LogUniform(0.5, 1000) # MJup
    M_sec ~ truncated(Normal(300,100),lower=100)
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol
    rv = 0.0 # system RV not significant for this example

    ref_epoch = Octofitter.meta_gaia_DR2.ref_epoch_mjd
    plx  ~ truncated(Normal(30,5),lower=1,upper=50)
    pmra ~ Normal(-137, 10)
    pmdec ~ Normal(2,  10)
    gaia_ra_offset_deg  ~  Normal(0, 1000/1000/60/60)
    gaia_dec_offset_deg  ~ Normal(0, 1000/1000/60/60)
    dec = $gaialike.dr2.dec + system.gaia_ra_offset_deg
    ra = $gaialike.dr2.ra + system.gaia_dec_offset_deg

    # σ_scan ~ truncated(Normal(0,10), lower=0)
    excess_noise_dr2  ~ truncated(Normal(0,10), lower=0)
    excess_noise_dr3  ~ truncated(Normal(0,10), lower=0)

end gaialike b
# end absastromlike b
using FiniteDiff
model = Octofitter.LogDensityModel(sys,autodiff=:FiniteDiff, verbosity=4)

##
Octofitter.default_initializer!(model,ntries=0,verbosity=4)
##
model.starting_points = fill(
    model.link([
        gaialike.dr2.parallax,
        gaialike.dr2.pmra,
        gaialike.dr2.pmdec,
        0.0,0.0,
        1.0, 1.0]
    ),
    100
)
##
using Pigeons
chain,pt = octofit_pigeons(model,n_rounds=8, explorer=SliceSampler())
##
fig = octoplot(model,chain,show_mass=true)
xlims!(fig.content[1], 400,-400)
ylims!(fig.content[1], -400,400)
fig
##
increment_n_rounds!(pt,1)
chain,pt = octofit_pigeons(pt)
##
res = Octofitter.mcmcchain2result(model,chain)
els = Octofitter.construct_elements(chain,:b,:);

using DataFrames
samples_dr2 = DataFrame(
    rand(gaialike.dist_dr2, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
samples_dr2.ra .-= gaialike.dr2.ra
samples_dr2.dec .-= gaialike.dr2.dec
# post_dr2 = DataFrame(
#     stack(getproperty.(res, :P_DR2))',
#     [:parallax,:ra,:dec,:pmra, :pmdec]
# )
# post_dr2.ra .-= gaialike.dr2.ra
# post_dr2.dec .-= gaialike.dr2.dec



samples_dr3 = DataFrame(
    rand(gaialike.dist_dr3, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
samples_dr3.ra .-= gaialike.dr3.ra
samples_dr3.dec .-= gaialike.dr3.dec

# post_dr3 = DataFrame(
#     stack(getproperty.(res, :P_DR3))',
#     [:parallax,:ra,:dec,:pmra, :pmdec]
# )
# post_dr3.ra .-= gaialike.dr3.ra
# post_dr3.dec .-= gaialike.dr3.dec

model_outputs = map(res,els) do res, el
    ll,  (dr3post, dr2post) = Octofitter.simulate(gaialike,res,el)
    return (dr3post, dr2post)
end
post_dr2 = DataFrame(
    stack(getindex.(model_outputs, 2))',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
post_dr2.ra .-= gaialike.dr2.ra
post_dr2.dec .-= gaialike.dr2.dec
post_dr3 = DataFrame(
    stack(getindex.(model_outputs, 1))',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
post_dr3.ra .-= gaialike.dr3.ra
post_dr3.dec .-= gaialike.dr3.dec

post_model = DataFrame(
    parallax=chain[:plx][:],
    ra=chain[:ra][:],
    dec=chain[:dec][:],
    pmra=chain[:pmra][:],
    pmdec=chain[:pmdec][:],
)
post_model.ra .-= gaialike.dr2.ra
post_model.dec .-= gaialike.dr2.dec


fig = Figure(size=(700,700))

pairplot(
    fig[1:5,1:5],
    PairPlots.Series(post_dr2,label="Model: DR2",color=:red,)=>(
        # PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.Scatter(markersize=2),
        # PairPlots.MarginDensity(bandwidth=1.5,linestyle=:dash)
    ),
    PairPlots.Series(post_dr3,label="Model: DR3",color=:blue,)=>(
        # PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.Scatter(markersize=2),
        # PairPlots.MarginDensity(bandwidth=1.5,linestyle=:dash)
    ),
    PairPlots.Series(samples_dr2,label="DR2",color=:red)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5,color=:darkred),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(samples_dr3,label="DR3",color=:blue)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5,color=:darkblue),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    # PairPlots.Series(post_model,label="Model",color=:black,)=>(
    #     # PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
    #     PairPlots.Scatter(markersize=2),
    #     PairPlots.MarginDensity(bandwidth=1.5,linestyle=:solid)
    # )
)
fig
##

fig = Figure()
ax = Axis(
    # fig[1:2,4:5],
    fig[1,1],
    xlabel="[mas]",
    ylabel="[mas]",
    backgroundcolor="#EEE",
    autolimitaspect=1
)
l = 30
xlims!(ax,-l,l)
ylims!(ax,-l,l)

# scatterlines!(ax,gaialike.table.epoch, resids)
N = length(gaialike.table.epoch)
i = rand(1:length(res))
resids = res[i].along_scan_resid
# uncs = fill(res[i].σ_scan,N)
uncs = fill(3,N)
for i in eachindex(resids)
    s,c = sincos(gaialike.table.var" scanAngle[rad]"[i])
    lines!(ax,
        [resids[i] .* c - uncs[i].*c, resids[i] .* c + uncs[i].*c],
        [resids[i] .* s - uncs[i].*s, resids[i] .* s + uncs[i].*s],
        color=:black
        # 
    )
end
p = scatter!(ax,
    vec(resids .* cos.(gaialike.table.var" scanAngle[rad]")),
    vec(resids .* sin.(gaialike.table.var" scanAngle[rad]")),
    color=vec(gaialike.table.epoch),
    colormap=:turbo
)
# Colorbar(fig[1,3], p, label="MJD")
Colorbar(fig[1,2], p, label="MJD")
fig