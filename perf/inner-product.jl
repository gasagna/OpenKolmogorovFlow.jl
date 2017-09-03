using OpenKolmogorovFlow
using BenchmarkTools

n = 100
u = randn(n, n); U = FFT(Field(u))
v = randn(n, n); V = FFT(Field(v))
W = shifted(V, (1, 2))

@btime inner($U, $V)
@btime inner($U, $W)
