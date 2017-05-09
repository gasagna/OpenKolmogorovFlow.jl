using OpenKolmogorovFlow
using BenchmarkTools

n = 128
w = FTField(1.0 + im*randn(n, n>>1+1))
u = FTField(1.0 + im*randn(n, n>>1+1))
A = FTField(1.0 + im*randn(n, n>>1+1))
L = ViscousTerm(n, 1.0)

foo(w, A, u) = w .= A .* u

# @btime foo(w, A, u)
foo(w, A, u)

# @time foo(w, A, u)
# @time foo(w, A, u)
# @time foo(w, A, u)

# @btime foo(w, L, u)
foo(w, L, u)
# @time foo(w, L, u)
# @time foo(w, L, u)
# @time foo(w, L, u)
# @time foo(w, L, u)