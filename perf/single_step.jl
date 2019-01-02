using OpenKolmogorovFlow
using Flows
using FFTW
using BenchmarkTools

# parameters
const Re = 40
const kforcing = 4
const n = 64
const m = 96

# get explicit and implicit parts
sys = Flows.System(splitexim(ForwardEquation(n, m, Re, kforcing, FFTW.EXHAUSTIVE))...)

# initial condition
Ω = FTField(n, m)

# define method
method = CB3R2R3e(Ω, :NORMAL)

# measure time it takes to complete a step
@btime Flows.step!($method, $sys, 0, 0.001, $Ω, nothing)