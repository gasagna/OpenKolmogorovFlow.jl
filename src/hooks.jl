import Flows: AbstractTimeStepFromHook

export CFLHook

struct CFLHook <: AbstractTimeStepFromHook
       CFL_max::Float64
       Δt_bounds::Tuple{Float64, Float64}
    function CFLHook(CFL_max::Real=1, Δt_min::Real=0, Δt_max::Real=1)
        CFL_max > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(CFL_max), (Float64(Δt_min), Float64(Δt_max)))
    end
end

(hook::CFLHook)(g::ForwardExplicitTerm, A, z) = 
    clamp(hook.CFL_max * g.β[1], hook.Δt_bounds...)
