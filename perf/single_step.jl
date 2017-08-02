using OpenKolmogorovFlow
using IMEXRKCB

# parameters
const Re = 40
const kforcing = 4
const n = 132 

# get explicit and implicit parts
Ld, Nd = imex(VorticityEquation(n, Re, kforcing; dealias=true))

# initial condition
Ω₀ = FTField(n)

# define scheme
scheme = IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4, false), Ω₀)

# measure time it takes to complete a step
println(minimum([@elapsed step!(scheme, Nd, Ld, 0, 0.1, Ω₀) for i = 1:100]))