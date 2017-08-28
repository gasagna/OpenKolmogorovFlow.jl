using OpenKolmogorovFlow
using BenchmarkTools
import MacroTools: postwalk, @capture

n = 100
u = randn(n, n); U = FFT(Field(u))
v = randn(n, n); V = FFT(Field(v))

@btime inner($U, $V)
