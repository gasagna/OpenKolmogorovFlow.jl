export laminarflow, dissrate

# Returns the vorticity field of the laminar flow.
function laminarflow(n::Int, Re::Real, kforcing::Int=4)
    0 ≤ kforcing ≤ n>>1  ||
        throw(ArgumentError("forcing wave number must be in [0, n/2]"))
    Ω = FTField(n)
    Ω[-kforcing, 0] = -Re/kforcing/2
    Ω[ kforcing, 0] = -Re/kforcing/2
    Ω
end

# Energy dissipation rate density associated to vorticity field `Ω`
dissrate(Ω::FTField{n}, Re::Real) where {n} = (@Σ_jk n abs2(Ω[jk]))/Re