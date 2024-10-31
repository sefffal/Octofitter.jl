
using Octofitter
# using OctofitterRadialVelocity
using CairoMakie
# hip_like = Octofitter.HipparcosIADLikelihood(;hip_id=1234,renormalize=true,full_3D_model=false,)
# hip_like_b = Octofitter.HipparcosIADLikelihood(;hip_id=27321,renormalize=true,full_3D_model=true)#27321)1234

# Beta pic: 27321
hip_like_bpic = Octofitter.HipparcosIADLikelihood(;hip_id=27321,renormalize=true,)
# hip_like.table.reject .= false
# hip_like.table.reject[[3,66,84]] .= true

# # Eps Indi: 108870
hip_like_epsind = Octofitter.HipparcosIADLikelihood(;hip_id=108870,renormalize=true,full_3D_model=true)
# hip_like.table.reject .= hip_like.table.sres .<= 0
# hip_like.table.reject[[117;119;120;129;132:137]] .= true
##
# hip_like.table.reject[117] = true
##
# plxs = hip_like.table.plx_vs_time
scatterlines(hip_like.table.Δα✱, hip_like.table.Δδ)
##
hip_like =hip_like_bpic
# hip_like.table.plx_vs_time .= hip_like.hip_sol.plx
using Distributions
@planet b AbsoluteVisual{KepOrbit} begin
# @planet b Visual{KepOrbit} begin
    e = 0. 
    ω = 0. 
    mass = 0.
    a = 1.
    i = 0
    Ω = 0.
    tp = 0.
end

# Conv
@system betapichip begin
    plx ~ Uniform(0,1000)
    M = 1.0#~ truncated(Normal(0.76, 0.04),lower=0)
    pmra ~ Uniform(-10000, 10000)
    pmdec ~  Uniform(-10000, 10000)

    rv = -40.4e3#-40.0e3
    ra_hip_offset_mas ~  Uniform(-10000, 10000)
    dec_hip_offset_mas ~ Uniform(-10000, 10000)
    dec = hip_like.hip_sol.dedeg + (system.dec_hip_offset_mas/60/60/1000)
    ra = hip_like.hip_sol.radeg + system.ra_hip_offset_mas/60/60/1000#/cosd(system.dec)
    ref_epoch = Octofitter.hipparcos_catalog_epoch_mjd

end hip_like b

model = Octofitter.LogDensityModel(betapichip; autodiff=:ForwardDiff, verbosity=4) 
##
Octofitter.default_initializer!(model, initial_samples=5_000,verbosity=4)
##
chain_hip_simp = octofit(model, iterations=2_000, max_depth=6)

## Find best fitting solution
# model.system.observations[1].table.rv_kms .= -400 # -40.4
##
# sol = Octofitter._refine(model, chain_hip_simp.info.samples_transformed[argmax(chain_hip_simp[:logpost][:])])
sol = Octofitter._refine(model, model.starting_points[1])
solnt =  (;logpost=0,model.arr2nt(model.invlink(sol))...,)#ra=hip_like.hip_sol.radeg+100)
# solnt =  (;logpost=0,model.arr2nt(model.invlink(sol))...,plx=10)#ra=hip_like_b.hip_sol.radeg+100)
chn_hip_bf = Octofitter.result2mcmcchain([solnt],  Dict(:internals => [:logpost]))
## Print out residuals
Octofitter.idenfity_rejectable_scans(model.system.observations[1], solnt, [Octofitter.construct_elements(chn_hip_bf,:b,1)])
##
octoplot(model, chn_hip_bf, show_astrom=false, show_astrom_time=false)
# plot_hip(model, chn_hip_bf)
##
using LinearAlgebra, StatsBase
fig = Figure(size=(900,600).*1.2)
j = i = 1
ax = nothing
for prop in (
    (;chain=:ra, hip=:radeg, hip_err=:e_ra), 
    (;chain=:dec, hip=:dedeg, hip_err=:e_de),
    # (;chain=:ra_hip_offset_mas, hip=:zero, hip_err=:e_ra), 
    # (;chain=:dec_hip_offset_mas, hip=:zero, hip_err=:e_de),
    (;chain=:plx, hip=:plx, hip_err=:e_plx), 
    (;chain=:pmra, hip=:pm_ra, hip_err=:e_pmra), 
    (;chain=:pmdec, hip=:pm_de, hip_err=:e_pmde)
)
    global ax = Axis(
        fig[j,i],
        xlabel=string(prop.chain),
    )
    i+=1
    if i > 3
        j+=1
        i = 1
    end
    unc = hip_like.hip_sol[prop.hip_err]
    if prop.chain == :ra
        unc /= 60*60*1000 * cos(hip_like.hip_sol.dedeg)
    end
    if prop.chain == :dec
        unc /= 60*60*1000
    end
    if prop.hip == :zero
        n = Normal(0, unc)
    else
        mu = hip_like.hip_sol[prop.hip]
        n = Normal(mu, unc)
    end
    n0,n1=quantile.(n,(1e-4, 1-1e-4))
    nxs = range(n0,n1,length=200)
    h = fit(Histogram, chain_hip_simp[prop.chain][:], nbins=55)
    h = normalize(h, mode=:pdf)
    barplot!(ax, (h.edges[1][1:end-1] .+ h.edges[1][2:end])./2, h.weights, gap=0, color=:red, label="posterior")
    lines!(ax, nxs, pdf.(n,nxs), label="Hipparcos Catalog", color=:black, linewidth=2)
    vlines!(ax, chn_hip_bf[prop.chain][:],color=:blue,linewidth=3,label="maximum a-posteriori")
end
Legend(fig[i-1,j+1],ax,tellwidth=false)
fig
##
plot_hip(model,chain_hip_simp)
##
function plot_hip(model, chain)

    hip_like = nothing
    for like_obs in model.system.observations
        if like_obs isa HipparcosIADLikelihood
            if !isnothing(hip_like)
                error("more than one HipparcosIADLikelihood present")
            end
            hip_like = like_obs
        end
    end
    if isnothing(hip_like)
        error("No HipparcosIADLikelihood present")
    end

    fig = Figure(
        size=(800,800)
    )

    ax = Axis(fig[1,1:2],
        xlabel="α* [mas]",
        ylabel="δ [mas]",
        autolimitaspect=1
    )


    scatterlines!(hip_like.table.Δα✱, hip_like.table.Δδ,
        label="Hipparcos model",
        linewidth=4,
        alpha=0.5,
        markersize=4,
        color=:grey
    )

    ax2 = Axis(
        fig[2,1:2],
        xlabel="epoch [MJD]",
        ylabel="along-scan residual [mas]"
    )
    errorbars!(ax2, hip_like.table.epoch,  zeros(size(hip_like.table.epoch)), hip_like.table.sres_renorm, color=:black)

    # TODO:
    planet_keys = keys(chain)
    pk = :b
    nts = Octofitter.mcmcchain2result(model, chain,)

    # ii = rand(axes(chain),1)
    _, ii = findmax(chain[:logpost][:])
    for i in ii
        orbC = Octofitter.construct_elements(chain, pk, i)
        sim = Octofitter.simulate(hip_like, nts[i,], [orbC])
        # Model
        scatterlines!(ax,
            sim.α✱_model_with_perturbation[:,1],
            sim.δ_model_with_perturbation[:,1],
            label="Our Model",
            color=:red,
            linewidth=1,
            markersize=4,
            alpha=0.85,
        )

        # Data
        for i in axes(hip_like.table.α✱ₘ,1)
            # @show point1 point2
            lines!(ax, hip_like.table.α✱ₘ[i][1:2], hip_like.table.δₘ[i][1:2],color=:black,alpha=0.1)
        end

        resid = map(eachindex(hip_like.table.epoch)) do i
            # point = α✱ₘ
            point =  [
                sim.α✱_model_with_perturbation[i],
                sim.δ_model_with_perturbation[i]
            ]
            line_point_1 =  [hip_like.table.α✱ₘ[i][1], hip_like.table.δₘ[i][1]]
            line_point_2 =  [hip_like.table.α✱ₘ[i][2], hip_like.table.δₘ[i][2]]
            Octofitter.distance_point_to_line(point, line_point_1, line_point_2)
        end
        scatter!(ax2, hip_like.table.epoch, resid, markersize=4, color=:red, alpha=1)

        for i in 1:length(resid)
            # Need to draw along scan line and residual on the plot
            # start point at data point x and y
            # move at 90 degree angle to the scan direction, for a length of `resid`
            rise = hip_like.table.δₘ[i][2] - hip_like.table.δₘ[i][1]
            run = hip_like.table.α✱ₘ[i][2] - hip_like.table.α✱ₘ[i][1]            
            x0 = sim.α✱_model_with_perturbation[i]
            y0 = sim.δ_model_with_perturbation[i]

            # Two possible directions, check which is closer
            line_point_1 =  [hip_like.table.α✱ₘ[i][1], hip_like.table.δₘ[i][1]]
            line_point_2 =  [hip_like.table.α✱ₘ[i][2], hip_like.table.δₘ[i][2]]
            angle_1 = atan(rise,run) - π/2
            x1_1 = x0 + resid[i]*cos(angle_1)
            y1_1 = y0 + resid[i]*sin(angle_1)
            d1 = Octofitter.distance_point_to_line([x1_1,y1_1], line_point_1, line_point_2)
            angle_2 = atan(rise,run) + π/2
            x1_2 = x0 + resid[i]*cos(angle_2)
            y1_2 = y0 + resid[i]*sin(angle_2)
            d2 = Octofitter.distance_point_to_line([x1_2,y1_2], line_point_1, line_point_2)
            if d1 < d2
                angle = angle_1
                x1 = x1_1 
                y1 = y1_1
            else
                angle = angle_2
                x1 = x1_2
                y1 = y1_2
            end

            # Now we plot the uncertainties along this direction, centred around the intersection point
            unc_x1 = x1 + hip_like.table.sres_renorm[i]*cos(angle)
            unc_y1 = y1 + hip_like.table.sres_renorm[i]*sin(angle)
            unc_x2 = x1 - hip_like.table.sres_renorm[i]*cos(angle)
            unc_y2 = y1 - hip_like.table.sres_renorm[i]*sin(angle)
            lines!(ax, [unc_x1,unc_x2], [unc_y1,unc_y2], color=:blue, alpha=0.25, linewidth=4)


            # Plot residual line (from model point to intersection point [x1,y1])
            lines!(ax, [x0,x1], [y0,y1], color=:blue, linewidth=1)

        end
    end
    Legend(fig[1,3],ax)
    ylims!(ax2, low=0)
    fig
end
plot_hip(model, chn_hip_bf)#chain_hip_simp)


## Incantation to check intermediate values at the Hipparcos solution
sim = Octofitter.simulate(
    hip_like,
    model.arr2nt(
        # [hip_like.hip_sol.plx, hip_like.hip_sol.pm_ra, hip_like.hip_sol.pm_de, 0, 0]
        sol
    ),
    [AbsoluteVisual{KepOrbit}(;
        hip_like.hip_sol.plx,
        hip_like.hip_sol.pm_ra,
        hip_like.hip_sol.pm_de,
        ref_epoch=Octofitter.hipparcos_catalog_epoch_mjd,
        ra=hip_like.hip_sol.radeg,
        dec=hip_like.hip_sol.dedeg,
        rv=-40.1e3,
        pmra=hip_like.hip_sol.pm_ra,
        pmdec=hip_like.hip_sol.pm_de,
        M=1,a=1,e=0,i=0,ω=0,Ω=0,tp=0
    )]
)
##
scatterlines(hip_like.table.Δα✱, hip_like.table.Δδ)
# sim.α✱_model_with_perturbation
##
orb = AbsoluteVisual{KepOrbit}(;
hip_like.hip_sol.plx,
hip_like.hip_sol.pm_ra,
hip_like.hip_sol.pm_de,
ref_epoch=Octofitter.hipparcos_catalog_epoch_mjd,
ra=hip_like.hip_sol.radeg,
dec=hip_like.hip_sol.dedeg,
rv=-40.1e3,
pmra=hip_like.hip_sol.pm_ra,
pmdec=hip_like.hip_sol.pm_de,
M=1,a=1,e=0,i=0,ω=0,Ω=0,tp=0
)
sol = orbitsolve(orb, mjd("1992-10"))
sol.compensated.ra2 .- hip_like.hip_sol.radeg
sol.compensated.dec2 .- hip_like.hip_sol.dedeg