export laminarflow, Energy

"""
    Returns the vorticity field of the laminar flow.
"""
function laminarflow(n::Int, Re::Real, kforcing::Int=4)
    Ω = FTField(n)
    Ω[-kforcing, 0] = -Re/kforcing/2
    Ω[ kforcing, 0] = -Re/kforcing/2
    Ω
end

"""
    Kinetic energy density associated to vorticity field `U`
"""
Energy(U::FTField{n}) where n = 
    sum(abs(U[k, 0])^2 for k=0:n, j=0:n)/4π^2