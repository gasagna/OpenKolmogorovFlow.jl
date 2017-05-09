abstract type AbstractTrajectory{V<:AbstractArray, T} end

struct Trajectory{V, T, F} <: AbstractTrajectory{V, T}
    seeds::V            # a vector of seeds
    ts::Vector{Float64} # a vector of times
    f::F                # the dynamical system that generates the traj
end

# get state at time t
function (o::Type{Trajectory{V, T, F}}){V, T, F}(x::V, t::Real)
    # find nearest previous time
    t₀ = o.ts[findlast(τ -> τ ≥ t)]
    # integrate from that point to t
    propagate!(f, x, t₀, t)
end