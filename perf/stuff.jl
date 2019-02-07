using ProfileView
using OpenKolmogorovFlow
using Flows

# setup
Re, m, n, Δt = 100, 63, 63, 0.004

# initial condition
ω = Field(m, (x, y)->0.0001*randn()); ω .-= mean(ω); Ω = FFT(ω, n)

# right hand side
f = ForwardEquation(n, m, Re)

# and split
N, L = splitexim(f)

# forward map
ϕ = integrator(N, L, Scheme(:CB3e_3R2R, Ω),  Δt)

# run forward
ϕ(Ω, (0, 1))

@profile ϕ(Ω, (0, 10))
ProfileView.view()
