export laminarflow, dissrate

"""
    Returns the vorticity field of the laminar flow.
"""
function laminarflow(n::Int, Re::Real, kforcing::Int=4)
    kforcing ≤ n>>1 || 
        throw(ArgumentError("forcing wave number too large"))
    Ω = FTField(n)
    Ω[-kforcing, 0] = -Re/kforcing/2
    Ω[ kforcing, 0] = -Re/kforcing/2
    Ω
end

"""
    Energy dissipation rate density associated to vorticity field `U`
"""
function dissrate(U::FTField{n}, Re::Real) where n
    @inbounds begin
        d = n>>1
        val = zero(abs(U[0, 0])^2)
        for k=-d+1:d
            for j=1:d-1
                val += 2*abs(U[k, j])^2
            end
            val +=  abs(U[k, 0])^2 + abs(U[k, d])^2
        end
        # count properly contribution of extreme cases
        val += -abs(U[d, d])^2 + 2*abs(U[d, d]/2)^2
        val += -abs(U[0, d])^2 + 2*abs(U[0, d]/2)^2
        val += -abs(U[d, 0])^2 + 2*abs(U[d, 0]/2)^2
    end
    val/Re
end