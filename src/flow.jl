export laminarflow, DissipationRate

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
    Energy dissipation rate density associated to vorticity field `U`
"""
function DissipationRate(U::FTField{n}, Re::Real) where n
    d = n>>1
    @inbounds begin
        val  = 2*sum(abs(U[k, j])^2 for k=-d+1:d, j=1:d-1)
        val +=   sum(abs(U[k, 0])^2 for k=-d+1:d)
        val +=   sum(abs(U[k, d])^2 for k=-d+1:d)
        # count properly contribution of extreme cases
        val += -abs(U[d, d])^2 + 2*abs(U[d, d]/2)^2
        val += -abs(U[0, d])^2 + 2*abs(U[0, d]/2)^2
        val += -abs(U[d, 0])^2 + 2*abs(U[d, 0]/2)^2
    end
    val/Re
end