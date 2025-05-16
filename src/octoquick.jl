
octoquick(model::LogDensityModel; kwargs...) = 
    octoquick(nothing, model; kwargs...)
function octoquick(
    rng,
    model::LogDensityModel;
    verbosity::Int=2,
    kwargs...
)
    @nospecialize

    chain = Octofitter.initialize!(rng, model; verbosity, kwargs...)

    
    return chain
end
export octoquick