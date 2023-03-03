using UnPack: @unpack

using Setfield
import Setfield: ConstructionBase

"""
This implements a custom Hamiltonian Monte Carlo integrator
based on the AdvancedHMC Leapfrog type.

We know that many models exhibit pathalogical behaviour as
eccentricity approaches unity.
To ease sampling, this integrator scales down the step 
size as eccentricity increases.
"""

struct EccentricLeapfrog1{T<:AdvancedHMC.AbstractScalarOrVec} <: AdvancedHMC.AbstractLeapfrog{T}
    "Nominal step size."
    ϵ0::T
    "Index of eccentricity parameter"
    e_idx::Int
end

function Base.show(io::IO, l::EccentricLeapfrog1)
    print(
        io,
        "EccentricLeapfrog1(ϵ0=$(round.(l.ϵ0; sigdigits=3)))",
    )
end

AdvancedHMC.step_size(lf::EccentricLeapfrog1) = lf.ϵ0
# AdvancedHMC.step_size(lf::EccentricLeapfrog1, θ) = lf.ϵ0
function AdvancedHMC.step_size(lf::EccentricLeapfrog1, θ)
    e = θ[lf.e_idx] # This is the eccentricity.
    # But if we have used a bijector, it is actually mapped to some real domain.

    # return a stepsize that decreases with higher eccentricity
    lf.ϵ0/(1+10e)
end

AdvancedHMC.nom_step_size(lf::EccentricLeapfrog1) = lf.ϵ0
AdvancedHMC.update_nom_step_size(lf::EccentricLeapfrog1, ϵ0) = EccentricLeapfrog1(ϵ0, lf.e_idx)


function AdvancedHMC.step(
    lf::EccentricLeapfrog1{T},
    h::AdvancedHMC.Hamiltonian,
    z::P,
    n_steps::Int = 1;
    fwd::Bool = n_steps > 0,  # simulate hamiltonian backward when n_steps < 0
    full_trajectory::Val{FullTraj} = Val(false),
) where {T<:AdvancedHMC.AbstractScalarOrVec{<:AbstractFloat},P<:AdvancedHMC.PhasePoint,FullTraj}
    n_steps = abs(n_steps)  # to support `n_steps < 0` cases

    @unpack θ, r = z
    ϵ = fwd ? AdvancedHMC.step_size(lf, θ) : -AdvancedHMC.step_size(lf, θ)


    ϵ = ϵ'

    res = if FullTraj
        Vector{P}(undef, n_steps)
    else
        z
    end

    @unpack value, gradient = z.ℓπ
    for i = 1:n_steps
        # Tempering
        r = AdvancedHMC.temper(lf, r, (i = i, is_half = true), n_steps)
        # Take a half leapfrog step for momentum variable
        r = r - ϵ / 2 .* gradient
        # Take a full leapfrog step for position variable
        ∇r = AdvancedHMC.∂H∂r(h, r)
        θ = θ + ϵ .* ∇r
        # Take a half leapfrog step for momentum variable
        @unpack value, gradient = AdvancedHMC.∂H∂θ(h, θ)
        r = r - ϵ / 2 .* gradient
        # Tempering
        r = AdvancedHMC.temper(lf, r, (i = i, is_half = false), n_steps)
        # Create a new phase point by caching the logdensity and gradient
        z = AdvancedHMC.phasepoint(h, θ, r; ℓπ = AdvancedHMC.DualValue(value, gradient))
        # Update result
        if FullTraj
            res[i] = z
        else
            res = z
        end
        if !isfinite(z)
            # Remove undef
            if FullTraj
                res = res[isassigned.(Ref(res), 1:n_steps)]
            end
            break
        end
    end
    return res
end

# step_size(lf::AbstractLeapfrog) = lf.ϵ


# abstract type AbstractLeapfrog{T} <: AbstractIntegrator end

# step_size(lf::AbstractLeapfrog) = lf.ϵ
# jitter(::Union{AbstractRNG,AbstractVector{<:AbstractRNG}}, lf::AbstractLeapfrog) = lf
# temper(lf::AbstractLeapfrog, r, ::NamedTuple{(:i, :is_half),<:Tuple{Integer,Bool}}, ::Int) =
#     r
# stat(lf::AbstractLeapfrog) = (step_size = step_size(lf), nom_step_size = nom_step_size(lf))

# update_nom_step_size(lf::AbstractLeapfrog, ϵ) = @set lf.ϵ = ϵ
