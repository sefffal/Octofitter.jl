
# The stats base sample function makes it easy to get values from Chains
# but converting these values into KeplerianElements along with any static
# paramters takes a bit of work.
# Here we overload the sample function with a method just for this.
# function sample(::Type{KeplerianElements}, chains::Chains, static, N=1)
#     sampled = sample(chains, ceil(Int, N/size(chains,2)))
#     out = KeplerianElements{Float64}[]

#     proto = namedtuple(keys(chains))

#     sizehint!(out, size(sampled,1)*size(sampled,2))
#     for i in 1:size(sampled,1), j in 1:size(sampled,3)
#         nt = proto(Array(sampled[i,:,j]))
#         el = KeplerianElements(merge(nt, static))
#         push!(out, el)
#     end
#     return out[begin:min(N,end)]
# end


function projectchains(chains,times)
    N = size(chains,1)*size(chains,3)
    ras = zeros(N)
    decs = zeros(N)
    ras = zeros(length(first(chains))*length(times))
    decs = zeros(length(first(chains))*length(times))
    i = 0
    for j=1:length(first(chains))
        
        a = chains.planets[1].a[j]
        inc = chains.planets[1].i[j]
        e = chains.planets[1].e[j]
        τ = chains.planets[1].τ[j]
        ω = chains.planets[1].ω[j]
        Ω = chains.planets[1].Ω[j]
        μ = chains.planets[1].μ[j]
        plx = chains.planets[1].plx[j]
        
        el = KeplerianElements(;a, i=inc, e, τ, ω, Ω, μ, plx)
        for t = times
            i+=1
            ra, dec, _ = kep2cart(el, t)
            ras[i] = ra
            decs[i] = dec
        end
    end
    return ras, decs,i
end

function sample_chain(planet,N)
    return map(rand(eachindex(planet.a),N)) do i
        return KeplerianElements(;
            a=planet.a[i],
            i=planet.i[i],
            e=planet.e[i],
            ω=planet.ω[i],
            Ω=planet.Ω[i],
            μ=planet.μ[i],
            plx=planet.plx[i],
            τ=planet.τ[i],
        )
    end
end