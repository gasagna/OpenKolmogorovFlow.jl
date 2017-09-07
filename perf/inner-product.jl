using OpenKolmogorovFlow
using BenchmarkTools

n = 100
u = randn(n, n); U = FFT(Field(u))
v = randn(n, n); V = FFT(Field(v))
cache = DistanceCache(64)
W = shifted(V, (1, 2))

@btime inner($U, $V)
@btime inner($U, $W)

@btime innerdiff($U, $V)
@btime innerdiff($U, $W)
@btime distance!($U, $V, $cache)