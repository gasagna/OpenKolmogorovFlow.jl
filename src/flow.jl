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
    val = zero(Re)
    for k=-d+1:d
        for j=1:d-1
            @inbounds val += 2*abs(U[k, j])^2
        end
        @inbounds val +=  abs(U[k, 0])^2 + abs(U[k, d])^2
    end
    # count properly contribution of extreme cases
    @inbounds val += -abs(U[d, d])^2 + 2*abs(U[d, d]/2)^2
    @inbounds val += -abs(U[0, d])^2 + 2*abs(U[0, d]/2)^2
    @inbounds val += -abs(U[d, 0])^2 + 2*abs(U[d, 0]/2)^2
    val/Re
end