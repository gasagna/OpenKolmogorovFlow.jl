import Flows: AbstractTimeStepFromHook

export CFLHook

struct CFLHook <: AbstractTimeStepFromHook
       CFL_max::Float64
    function CFLHook(CFL_max::Real=1)
        CFL_max > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(CFL_max))
    end
end

(hook::CFLHook)(g::ForwardExplicitTerm, A, z) = hook.CFL_max * g.β[1]
