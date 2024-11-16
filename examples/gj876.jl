dat_harps = OctofitterRadialVelocity.HARPS_RVBank_rvs("GJ876")
dat_hires = OctofitterRadialVelocity.HIRES_rvs("GL876")
dat_lick = OctofitterRadialVelocity.Lick_rvs("GL876")

rv_like_harps_03 = MarginalizedStarAbsoluteRVLikelihood(
    dat_harps[dat_harps.epoch .< mjd("2015")],
    instrument_name="HARPS 03",
    jitter=:j_harps_03
)
rv_like_harps_15 = MarginalizedStarAbsoluteRVLikelihood(
    dat_harps[mjd("2015") .< dat_harps.epoch .< mjd("2020")],
    instrument_name="HARPS 15",
    jitter=:j_harps_15
)
# rv_like_harps_20 = MarginalizedStarAbsoluteRVLikelihood(
#     dat_harps[mjd("2020") .< dat_harps.epoch],
#     instrument_name="HARPS 20",
#     jitter=:j_harps_20
# )


rv_like_hires = MarginalizedStarAbsoluteRVLikelihood(
    dat_hires,
    instrument_name="HIRES",
    jitter=:j_hires
)


lick_split = map(eachrow(dat_lick)) do row
    row = row[]
    if row.epoch < mjd("1994-11")
        # 5
        "lick1"
    elseif row.dewar == 13
        # 6
        "lick2"
    elseif row.dewar == 6
        # 7
        "lick3"
    else
        # 8
        "lick4"
    end
end
rv_like_lick1 = MarginalizedStarAbsoluteRVLikelihood(
    dat_lick[lick_split.=="lick1"],
    instrument_name="Lick1",
    jitter=:j_lick1
)
rv_like_lick2 = MarginalizedStarAbsoluteRVLikelihood(
    dat_lick[lick_split.=="lick2"],
    instrument_name="Lick2",
    jitter=:j_lick2
)
rv_like_lick3 = MarginalizedStarAbsoluteRVLikelihood(
    dat_lick[lick_split.=="lick3"],
    instrument_name="Lick3",
    jitter=:j_lick3
)
# rv_like_lick4 = MarginalizedStarAbsoluteRVLikelihood(
#     dat_lick[lick_split.=="lick4"],
#     instrument_name="Lick4",
#     jitter=:j_lick4
# )

fig= Figure()
ax = Axis(fig[1,1])
scatter!(ax, rv_like_harps_03.table.epoch, rv_like_harps_03.table.rv, rv_like_harps_03.table.σ_rv, label="harps_03")
scatter!(ax, rv_like_harps_15.table.epoch, rv_like_harps_15.table.rv, rv_like_harps_15.table.σ_rv, label="harps_15")
scatter!(ax, rv_like_lick1.table.epoch, rv_like_lick1.table.rv, rv_like_lick1.table.σ_rv, label="lick1")
scatter!(ax, rv_like_lick2.table.epoch, rv_like_lick2.table.rv, rv_like_lick2.table.σ_rv, label="lick2")
scatter!(ax, rv_like_lick3.table.epoch, rv_like_lick3.table.rv, rv_like_lick3.table.σ_rv, label="lick3")
# scatter!(ax, rv_like_lick4.table.epoch, rv_like_lick4.table.rv, rv_like_lick4.table.σ_rv, label="lick4")
scatter!(ax, rv_like_hires.table.epoch, rv_like_hires.table.rv, rv_like_hires.table.σ_rv, label="hires")
axislegend(ax)
fig
##



@planet b RadialVelocityOrbit begin
    e = 0.0#~ Uniform(0, 0.7)
    ω = 0.0#~ Uniform(0, 2pi)

    # To match RadVel, we set a prior on Period and calculate semi-major axis from it
    P ~ Uniform((61.1057 - 5)/Octofitter.julian_year, (61.1057 + 5)/Octofitter.julian_year)
    a = cbrt(system.M * b.P^2) # note the equals sign. 

    τ ~ UniformCircular(1.0)
    tp = b.τ*b.P*365.256360417 + 55000 # reference epoch for τ. Choose an MJD date near your data.
    
    # minimum planet mass [jupiter masses]. really m*sin(i)
    mass ~ LogUniform(0.001, 10)
end



@planet c RadialVelocityOrbit begin
    e ~ Uniform(0, 0.7)
    ω ~ Uniform(0, 2pi)

    # To match RadVel, we set a prior on Period and calculate semi-major axis from it
    P ~ Uniform((30 - 5)/Octofitter.julian_year, (30 + 5)/Octofitter.julian_year)
    a = cbrt(system.M * c.P^2) # note the equals sign. 

    τ ~ UniformCircular(1.0)
    tp = c.τ*c.P*365.256360417 + 55000 # reference epoch for τ. Choose an MJD date near your data.
    
    # minimum planet mass [jupiter masses]. really m*sin(i)
    mass ~ LogUniform(0.001, 10)
end


@planet d RadialVelocityOrbit begin
    e = 0 #~ Uniform(0, 0.7)
    ω = 0 #~ Uniform(0, 2pi)

    # To match RadVel, we set a prior on Period and calculate semi-major axis from it
    P ~ Uniform((1.937790 - 1)/Octofitter.julian_year, (1.937790 + 1)/Octofitter.julian_year)
    a = cbrt(system.M * d.P^2) # note the equals sign. 

    τ ~ UniformCircular(1.0)
    tp = d.τ*d.P*365.256360417 + 55000 # reference epoch for τ. Choose an MJD date near your data.
    
    # minimum planet mass [jupiter masses]. really m*sin(i)
    mass ~ LogUniform(0.00001, 1)
end



@planet e RadialVelocityOrbit begin
    e = 0.0 #~ Uniform(0, 0.7)
    ω = 0.0 #~ Uniform(0, 2pi)

    # To match RadVel, we set a prior on Period and calculate semi-major axis from it
    P ~ Uniform((123.55 - 10)/Octofitter.julian_year, (123.55 + 10)/Octofitter.julian_year)
    a = cbrt(system.M * e.P^2) # note the equals sign. 

    τ ~ UniformCircular(1.0)
    tp = e.τ*e.P*365.256360417 + 55000 # reference epoch for τ. Choose an MJD date near your data.
    
    # minimum planet mass [jupiter masses]. really m*sin(i)
    mass ~ LogUniform(0.00001, 1)
end


@system GL876 begin
    # total mass [solar masses]
    M ~ truncated(Normal(0.346,0.007),lower=0.1)
    
    j_harps_03 ~ LogUniform(0.1, 100)
    j_harps_15 ~ LogUniform(0.1, 100)
    j_hires ~ LogUniform(0.1, 100)
    j_lick1 ~ LogUniform(0.1, 100)
    j_lick2 ~ LogUniform(0.1, 100)
    j_lick3 ~ LogUniform(0.1, 100)
    j_lick4 ~ LogUniform(0.1, 100)

end rv_like_harps_03 rv_like_harps_15 rv_like_hires rv_like_lick1 rv_like_lick2 rv_like_lick3 b c d e


model = Octofitter.LogDensityModel(GL876)
##
Octofitter.default_initializer!(model)
##
chain,pt = octofit_pigeons(model,n_rounds=6)

##
increment_n_rounds!(pt,1)
chain,pt=octofit_pigeons(pt)
##
fig = Octofitter.rvpostplot(model,chain,show_summary=true)
ylims!(fig.content[3], -50, 50)
ylims!(fig.content[4], -240, 240)
ylims!(fig.content[7], -110, 110)
ylims!(fig.content[10], -50, 50)
ylims!(fig.content[13], -20, 20)
save("gj876.pdf",fig)
fig
##
Octofitter.savechain("GJ876.fits",chain)
##
Octofitter.rvpostplot_animated(model,chain,show_summary=true) do fig
    ylims!(fig.content[3], -50, 50)
    ylims!(fig.content[4], -240, 240)
    ylims!(fig.content[7], -110, 110)
    ylims!(fig.content[10], -50, 50)
    ylims!(fig.content[13], -20, 20)
end
##
el_b = Octofitter.construct_elements(chain, :b, argmax(chain[:logpost][:]))
o_b = orbit(;
    el_b.a,
    el_b.e,
    el_b.ω,
    el_b.tp,
    el_b.M,
    i=pi/2,
    Ω=0,
    plx=214.0380
)
ep = range(mjd("2025-01"),mjd("2025-05"),length=2000)
s = projectedseparation.(o_b, ep)
lines(
    ep, s,
    axis=(
        xlabel="MJD",
        ylabel="minimum separation [mas]"
    )
)
##
el_e = Octofitter.construct_elements(chain, :e, argmax(chain[:logpost][:]))
o_e = orbit(;
    el_e.a,
    el_e.e,
    el_e.ω,
    el_e.tp,
    el_e.M,
    i=pi/2,
    Ω=0,
    plx=214.0380
)
ep = range(mjd("2025-01"),mjd("2025-05"),length=2000)
s = projectedseparation.(o_e, ep)
lines(
    ep, s,
    axis=(
        xlabel="MJD",
        ylabel="minimum separation [mas]"
    )
)
##
el_e = Octofitter.construct_elements(chain, :e, argmax(chain[:logpost][:]))
o_e = orbit(;
    el_e.a,
    el_e.e,
    el_e.ω,
    el_e.tp,
    el_e.M,
    i=0,
    Ω=0,
    plx=214.0380
)

el_d = Octofitter.construct_elements(chain, :d, argmax(chain[:logpost][:]))
o_d = orbit(;
    el_d.a,
    el_d.e,
    el_d.ω,
    el_d.tp,
    el_d.M,
    i=0,
    Ω=0,
    plx=214.0380
)

el_c = Octofitter.construct_elements(chain, :c, argmax(chain[:logpost][:]))
o_c = orbit(;
    el_c.a,
    el_c.e,
    el_c.ω,
    el_c.tp,
    el_c.M,
    i=0,
    Ω=0,
    plx=214.0380
)


el_b = Octofitter.construct_elements(chain, :b, argmax(chain[:logpost][:]))
o_b = orbit(;
    el_b.a,
    el_b.e,
    el_b.ω,
    el_b.tp,
    el_b.M,
    i=0,
    Ω=0,
    plx=214.0380
)

##
fig = Figure(

)
ax = Axis(
    fig[1,1],
    autolimitaspect=1,
    xlabel="[au]",
    ylabel="[au]",
)

ts = 0:0.1:300
lines!(ax, posx.(o_b, ts), posy.(o_b, ts),)
lines!(ax, posx.(o_c, ts), posy.(o_c, ts),)
lines!(ax, posx.(o_d, ts), posy.(o_d, ts),)
lines!(ax, posx.(o_e, ts), posy.(o_e, ts),)

t = mjd()
scatter!(ax, posx(o_b,t), posy(o_b,t), markersize=15)
scatter!(ax, posx(o_c,t), posy(o_c,t), markersize=15)
scatter!(ax, posx(o_d,t), posy(o_d,t), markersize=15)
scatter!(ax, posx(o_e,t), posy(o_e,t), markersize=15)

fig