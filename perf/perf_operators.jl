using BenchmarkTools
using OpenKolmogorovFlow

c = 3.14
U = FTField(96, 128)
V = FTField(96, 128)

foo(U, V) = (U .*= V; U)

# @btime foo($V, $U)
@btime foo($V.data, $U.data)
@btime OpenKolmogorovFlow.ddx!($V, $U)
# @btime OpenKolmogorovFlow.ddy!($V, $U)
# @btime OpenKolmogorovFlow.invlaplacian!($V, $U)
# @btime OpenKolmogorovFlow.invlaplacian!($V, $U, $c)
# @btime OpenKolmogorovFlow.laplacian!($V, $U)
