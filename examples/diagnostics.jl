
## Auto-correlation
plot(
    autocor(chains[1].planets[1].e, 1:100:1000)
)


## Trace plots
plot()
plot!(chains[1].planets[1].a[1:end], title="d")
plot!(chains[1].planets[2].a[1:end], title="e")
## Trace plots
plot()
plot!(chains[1].planets[1].e[1:end], title="d")
plot!(chains[1].planets[2].e[1:end], title="e")


##
mean(getproperty.(stats[1], :acceptance_rate))
mean(getproperty.(stats[1], :tree_depth))
maximum(getproperty.(stats[1], :tree_depth))
count(getproperty.(stats[1], :numerical_error))
mean(getproperty.(stats[1], :n_steps))
maximum(getproperty.(stats[1], :log_density))

##
cd(@__DIR__)