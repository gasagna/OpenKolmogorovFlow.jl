using OpenKolmogorovFlow
using BenchmarkTools

n = 100
u = randn(n, n); U = FFT(Field(u))
v = randn(n, n); V = FFT(Field(v))
cache = DistanceCache(64)
W = shifted(V, (1, 2))

@btime dot($U, $V)
@btime dot($U, $W)

@btime dotdiff($U, $V)
@btime dotdiff($U, $W)
@btime distance!($U, $V, $cache)