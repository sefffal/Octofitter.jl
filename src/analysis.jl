"""
    projectpositions(chains.planets[1], mjd("2020-02-02"))

Given the posterior for a particular planet in the model and a modified julian date(s),
return `ra` and `dec` offsets in mas for each sampling in the posterior.
"""
function projectpositions(planet, times)
    N = size(planet, 1) * size(planet, 3)
    ras = zeros(N)
    decs = zeros(N)
    ras = zeros(length(first(planet)) * length(times))
    decs = zeros(length(first(planet)) * length(times))
    @threads for j = 1:length(first(planet))

        a = planet.a[j]
        inc = planet.i[j]
        e = planet.e[j]
        τ = planet.τ[j]
        ω = planet.ω[j]
        Ω = planet.Ω[j]
        μ = planet.μ[j]
        plx = planet.plx[j]

        el = KeplerianElements(; a, i = inc, e, τ, ω, Ω, μ, plx)
        for (k, t) in enumerate(times)
            i = j * length(times) + k - 1
            ra, dec, _ = kep2cart(el, t)
            ras[i] = ra
            decs[i] = dec
        end
    end
    return ras, decs
end
export projectpositions

"""
    sampleorbits(chains.planets[1], 100)

Given the posterior for a particular planet in the model, and the number of orbits to sample,
return a random subset of length N.
"""
function sampleorbits(planet, N)
    return map(rand(eachindex(planet.a), N)) do i
        return KeplerianElements(;
            a = planet.a[i],
            i = planet.i[i],
            e = planet.e[i],
            ω = planet.ω[i],
            Ω = planet.Ω[i],
            μ = planet.μ[i],
            plx = planet.plx[i],
            τ = planet.τ[i],
        )
    end
end
export sampleorbits

# Utility function for access a nested object via a vector or tuple of keys
nestedkey(obj, key::Symbol) = getproperty(obj, key)
nestedkey(obj, key::Symbol, remainingkeys...) = nestedkey(getproperty(obj,key), remainingkeys...)


function plotposterior end

# Optionally depend on Plots. If the user imports it, this code will be run to set up
# our `imshow` function.
using Requires
function init_plots()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin

        function plotposterior(planet, keys=:a, N=500;
            cmap=:turbo,
            rev=true,
            colorbartitle="semi-major axis (au)",
            clims=nothing,
            kwargs...
        )
            # sampled = sampleorbits(planet, N);
            ii = rand(eachindex(planet.a),N)
            elements = map(ii) do j
                KeplerianElements(;
                    a = planet.a[j],
                    i = planet.i[j],
                    e = planet.e[j],
                    ω = planet.ω[j],
                    Ω = planet.Ω[j],
                    μ = planet.μ[j],
                    plx = planet.plx[j],
                    τ = planet.τ[j],
                )
            end

            Plots.plot(;
                size=(700,700),
                dpi=200,
                fmt=:png,
                framestyle=:box,
                grid=:none,
                minorticks=true,
                aspectratio=1,
                fontfamily="Arial",
            # #     xlims=(-2000,2000),
            # #     ylims=(-2000,2000),
            # #     xticks=(-2000:1000:2000, ["-2", "-1", "0", "1", "2",]),
            # #     yticks=(-2000:1000:2000, ["-2", "-1", "0", "1", "2",]),
            #     xlims=(-500,500),
            #     ylims=(-500,500),
            #     xticks=(-500:500:500, ["-.5", "0", ".5",]),
            #     yticks=(-500:500:500, ["-.5", "0", ".5",]),
                margin=8Plots.mm,
                kwargs...
            )
            Plots.xlabel!("Δ right ascension (as)")
            Plots.ylabel!("Δ declination (as)")

            if keys isa Symbol
                colours = getproperty(planet, keys)[ii]
            else
                colours = nestedkey(planet, keys...)[ii]
            end

            if isnothing(clims)
                clims = extrema(colours)
            end
            
            minc = first(clims)
            dc = last(clims) - minc
            clims = (minc, dc+minc)

            cmap= Plots.cgrad(cmap; rev)
            for i in sortperm(colours, rev=false)
                c = colours[i]
                Plots.plot!(elements[i], label="",color=cmap[(c-minc)/dc], lw=0.3, alpha=0.1,)
            end
            Plots.scatter!([0],[0], marker=(:star,:black, 5), label="")
            Plots.scatter!([0],[0]; marker_z=[0], ms=0, clims=clims, color=cmap, label="", colorbartitle)

        end
    end
end