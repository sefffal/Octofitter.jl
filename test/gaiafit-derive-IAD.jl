using FillArrays
# Just for catalog access for now.
# gaialike = Octofitter.GaiaDivergenceLikelihood_3(;source_id_dr2=756291174721509376,source_id_dr3=756291174721509376)

gaialike = Octofitter.GaiaDivergenceLikelihood_3(;source_id_dr2=4373465352415301632,source_id_dr3=4373465352415301632)

mask_DR2 = vec(Octofitter.meta_gaia_DR2.start_mjd .< gaialike.table.epoch .< Octofitter.meta_gaia_DR2.stop_mjd)
mask_DR3 = vec(Octofitter.meta_gaia_DR3.start_mjd .< gaialike.table.epoch .< Octofitter.meta_gaia_DR3.stop_mjd)
##
dr3like = Octofitter.ParallacticMotionLikelihood_DR3(
    # gaialike.table,
    source_id_dr3=4373465352415301632,
    catalog_parameters=:P_DR3,
    along_scan_residuals=:along_scan_resid_dr3,
    excess_noise=:excess_noise_dr3,
    σ_scan=:σ_scan, 
    ref_epoch=Octofitter.meta_gaia_DR3.ref_epoch_mjd,
)
# ##
# dr2like = Octofitter.ParallacticMotionLikelihood_DR2(
#     # gaialike.table,
#     source_id_dr2=4373465352415301632,
#     catalog_parameters=:P_DR2,
#     along_scan_residuals=:along_scan_resid_dr2,
#     excess_noise=:excess_noise_dr2,
#     σ_scan=:σ_scan, 
#     ref_epoch=Octofitter.meta_gaia_DR2.ref_epoch_mjd,
#     initial_α = dr3like.table.initial_α[mask_DR2],
#     initial_δ = dr3like.table.initial_δ[mask_DR2],
# )
# ##

##
# Given the residuals x and σ_scan, calculate the likelihood P(x{dr2}|)

# Is it better to calculate the loglikihood of the data given the parameters
# Better to go data -> fit DR2, fit DR3. Compare to the catalog values.
#  Add both likelihoods together. 

##
absastromlike = Octofitter.StarAbsoluteAstrometryLikelihood_v1(
    gaialike.table,
    along_scan_residuals=:along_scan_resid_dr3,
    σ_scan=:σ_scan, 
    initial_α = dr3like.table.initial_α,
    initial_δ = dr3like.table.initial_δ,
)
##
# astrom_like = PlanetRelAstromLikelihood(
#     (epoch=mjd("2016-12-15"), ra=133., dec=-174., σ_ra=07.0, σ_dec=07., cor=0.2),
#     (epoch=mjd("2017-03-12"), ra=126., dec=-176., σ_ra=04.0, σ_dec=04., cor=0.3),
#     (epoch=mjd("2017-03-13"), ra=127., dec=-172., σ_ra=04.0, σ_dec=04., cor=0.1),
#     (epoch=mjd("2018-02-08"), ra=083., dec=-133., σ_ra=10.0, σ_dec=10., cor=0.4),
#     (epoch=mjd("2018-11-28"), ra=058., dec=-122., σ_ra=10.0, σ_dec=20., cor=0.3),
#     (epoch=mjd("2018-12-15"), ra=056., dec=-104., σ_ra=08.0, σ_dec=08., cor=0.2),
# )
# scatter(astrom_like.table.ra, astrom_like.table.dec, axis=(;autolimitaspect=1,xreversed=true))
##
# @planet b AbsoluteVisual{KepOrbit} begin
#     mass  = 0. 
#     e = 0.
#     ω = 0.
#     a = 1. 
#     i = 0. 
#     Ω = 0. 
#     tp = 0.0
# end

# hiplike = HipparcosIADLikelihood(;hip_id=51658)

# rvlike = PlanetRelativeRVLikelihood(
#     (epoch=mjd("2008-05-01"), rv=1300, σ_rv=150, inst_idx=1),
#     (epoch=mjd("2010-02-15"), rv=700, σ_rv=150, inst_idx=1),
#     (epoch=mjd("2016-03-01"), rv=-2700, σ_rv=150, inst_idx=1),
# )

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

end #rvlike astrom_like # Note the relative astrometry added here!

@system sys begin
    # M = 0.7847 # Host mass not important for this example
    M_pri ~ truncated(Normal(1.61, 0.1), lower=0)
    M_sec ~ LogUniform(0.5, 1000) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol
    rv = 0.0 # system RV not significant for this example

    ref_epoch = Octofitter.meta_gaia_DR2.ref_epoch_mjd
    plx  ~ truncated(Normal(30,5),lower=1,upper=50)
    pmra ~ Normal(-137, 10)
    pmdec ~ Normal(2,  10)
    gaia_ra_offset_deg  ~  Normal(0, 0.1)
    gaia_dec_offset_deg  ~ Normal(0, 0.1)
    dec = $gaialike.dr2.dec + system.gaia_ra_offset_deg
    ra = $gaialike.dr2.ra + system.gaia_dec_offset_deg
    
    # P_DR2 ~ gaialike.dist_dr2
    # P_DR3 ~ gaialike.dist_dr3

    # along_scan_resid ~ MvNormal(fill(10., count(mask_DR3)))
    σ_scan ~ truncated(Normal(0,10), lower=0)

    # excess_noise_dr2  ~ truncated(Normal(0,10), lower=0)
    # along_scan_resid_dr2 = system.along_scan_resid[$(findfirst(mask_DR2):findlast(mask_DR2))] #~ MvNormal(fill(10.0, count(mask_DR2)))
    excess_noise_dr3  ~ truncated(Normal(0,10), lower=0)
    # along_scan_resid_dr3 = system.along_scan_resid[$(findfirst(mask_DR3):findlast(mask_DR3))] #~ MvNormal(fill(10.0, count(mask_DR2)))
    along_scan_resid_dr3 ~ MvNormal(fill(10., count(mask_DR3)))

# end hiplike dr2like dr3like absastromlike b
end dr3like absastromlike b
# end absastromlike b
model = Octofitter.LogDensityModel(sys,autodiff=:ForwardDiff, verbosity=4,chunk_sizes=[63])

##
@system sysnoabs begin
    # M = 0.7847 # Host mass not important for this example
    M_pri ~ truncated(Normal(1.61, 0.1), lower=0)
    M_sec ~ LogUniform(0.5, 1000) # MJup
    M = system.M_pri + system.M_sec*Octofitter.mjup2msol
    rv = 0.0 # system RV not significant for this example

    ref_epoch = Octofitter.meta_gaia_DR2.ref_epoch_mjd
    plx  ~ truncated(Normal(30,5),lower=1,upper=50)
    pmra ~ Normal(-137, 10)
    pmdec ~ Normal(2,  10)
    gaia_ra_offset_deg  ~  Normal(0, 0.1)
    gaia_dec_offset_deg  ~ Normal(0, 0.1)
    dec = $gaialike.dr2.dec + system.gaia_ra_offset_deg
    ra = $gaialike.dr2.ra + system.gaia_dec_offset_deg
    
    P_DR2 ~ gaialike.dist_dr2
    P_DR3 ~ gaialike.dist_dr3

    along_scan_resid ~ MvNormal(fill(10., count(mask_DR3)))
    σ_scan ~ truncated(Normal(0,10), lower=0)

    excess_noise_dr2  ~ truncated(Normal(0,10), lower=0)
    along_scan_resid_dr2 = system.along_scan_resid[$(findfirst(mask_DR2):findlast(mask_DR2))] #~ MvNormal(fill(10.0, count(mask_DR2)))
    excess_noise_dr3  ~ truncated(Normal(0,10), lower=0)
    along_scan_resid_dr3 = system.along_scan_resid[$(findfirst(mask_DR3):findlast(mask_DR3))] #~ MvNormal(fill(10.0, count(mask_DR2)))

end b
modelnoabs = Octofitter.LogDensityModel(sysnoabs,autodiff=:ForwardDiff, verbosity=4,chunk_sizes=85)
## Use a model with no absolute astrometry to initialize the orbit solution
Octofitter.default_initializer!(modelnoabs,ntries=0,verbosity=4)
model.starting_points = copy(modelnoabs.starting_points)
##
using Pigeons
chain,pt = octofit_pigeons(model,n_rounds=3, explorer=SliceSampler())
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
##
i = rand(1:length(res))
f,a,p,=scatter(gaialike.table.epoch[:],res[i].along_scan_resid)
scatter!(gaialike.table.epoch[mask_DR2],res[i].along_scan_resid_dr2)
f
##
resids = stack(getproperty.(res, :along_scan_resid))'
m = mean(resids,dims=1)[:]
σ = std(resids,dims=1)[:]
f,a,p=errorbars(gaialike.table.epoch[:],m,σ)
scatter!(gaialike.table.epoch[:],m)
f
# ##
# pairplot(
#     stack(getproperty.(res, :along_scan_resid_dr2))',
# )
##
using DataFrames
samples_dr2 = DataFrame(
    rand(gaialike.dist_dr2, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
samples_dr2.ra .-= gaialike.dr2.ra
samples_dr2.dec .-= gaialike.dr2.dec
post_dr2 = DataFrame(
    stack(getproperty.(res, :P_DR2))',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
post_dr2.ra .-= gaialike.dr2.ra
post_dr2.dec .-= gaialike.dr2.dec



samples_dr3 = DataFrame(
    rand(gaialike.dist_dr3, 100000)',
    [:parallax,:ra,:dec,:pmra, :pmdec]
)
samples_dr3.ra .-= gaialike.dr3.ra
samples_dr3.dec .-= gaialike.dr3.dec
post_dr3 = DataFrame(
    stack(getproperty.(res, :P_DR3))',
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
    PairPlots.Series(post_dr2,label="DR2 Model",color=:red,)=>(
        # PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.Scatter(markersize=2),
        PairPlots.MarginDensity(bandwidth=1.5,linestyle=:dash)
    ),
    PairPlots.Series(post_dr3,label="DR3 Model",color=:blue,)=>(
        # PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.Scatter(markersize=2),
        PairPlots.MarginDensity(bandwidth=1.5,linestyle=:dash)
    ),
    PairPlots.Series(samples_dr2,label="DR2",color=:red)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5,color=:darkred),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(samples_dr3,label="DR3",color=:blue)=>(
        PairPlots.Contour(sigmas=1:1,bandwidth=1.5,color=:darkblue),
        PairPlots.MarginDensity(bandwidth=1.5)
    ),
    PairPlots.Series(post_model,label="Model",color=:black,)=>(
        # PairPlots.Contour(sigmas=1:1,bandwidth=1.5),
        PairPlots.Scatter(markersize=2),
        PairPlots.MarginDensity(bandwidth=1.5,linestyle=:solid)
    )
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