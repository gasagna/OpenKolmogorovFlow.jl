export laminarflow, dissrate

# Returns the vorticity field of the laminar flow.
function laminarflow(n::Int, m::Int, Re::Real, kforcing::Int=4)
    0 ≤ kforcing ≤ n  ||
        throw(ArgumentError("forcing wave number must be in [0, n]"))
    Ω = FTField(n, m)
    Ω[WaveNumber(-kforcing, 0)] = -Re/kforcing/2
    Ω[WaveNumber( kforcing, 0)] = -Re/kforcing/2
    return Ω
end

# Energy dissipation rate density associated to vorticity field `Ω`
dissrate(Ω::FTField{n, m}, Re::Real) where {n, m} = norm(Ω)^2/Re