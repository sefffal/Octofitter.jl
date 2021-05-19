
# The stats base sample function makes it easy to get values from Chains
# but converting these values into KeplerianElements along with any static
# paramters takes a bit of work.
# Here we overload the sample function with a method just for this.
function sample(::Type{KeplerianElements}, chains::Chains, static, N=1)
    sampled = sample(chains, ceil(Int, N/size(chains,2)))
    out = KeplerianElements{Float64}[]

    proto = namedtuple(keys(chains))

    sizehint!(out, size(sampled,1)*size(sampled,2))
    for i in 1:size(sampled,1), j in 1:size(sampled,3)
        nt = proto(Array(sampled[i,:,j]))
        el = KeplerianElements(merge(nt, static))
        push!(out, el)
    end
    return out[begin:min(N,end)]
end


function chain2elements(chains::Chains, static,)
    out = KeplerianElements{Float64}[]

    proto = namedtuple(keys(chains))

    sizehint!(out, size(chains,1)*size(chains,2))
    for i in 1:size(chains,1), j in 1:size(chains,3)
        nt = proto(Array(chains[i,:,j]))
        el = KeplerianElements(merge(nt, static))
        push!(out, el)
    end
    return out
end