using Base.Test
using OpenKolmogorovFlow
using IMEXRKCB
using PyPlot; pygui(true)

# parameters
n = 64
Re = 40
kforcing = 4

# initial condition
Ω₀ = FTField(n)
Ω₀[-5:5, -5:5] = 1e-3*randn(11, 11)

# get explicit and implicit parts
L = ImplicitTerm(n, Re)
N = ExplicitTerm(n, kforcing)
scheme = IMEXRK4R3R(IMEXRKCB4,  Ω₀, false)

# define T-time forward map
f = forwmap!(N, L, 200, 0.01, scheme)

# monitor the state excited by forcing
m = Monitor(Ω₀, Ω->Energy(Ω))

# map forward 
@time f(Ω₀, m)

# plot
plot(m.times, m.samples)
show()