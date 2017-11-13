# add method to this function
import IMEXRKCB: ImcA!

# ~~~ THE VISCOUS TERM OF THE GOVERNING EQUATIONS ~~~
struct ImplicitTerm{n, S, D<:DiffOperator{n, S}} <: AbstractFTOperator{n, S}
    Δ::D
    ν::Float64 # inverse of Reynolds number
end

# Outer constructor: TODO: merge with above
function ImplicitTerm(n::Int, Re::Real, ::Type{S}) where {S}
    Δ = DiffOperator(n, :xxyy, S)
    ImplicitTerm{n, S, typeof(Δ)}(Δ, 1/Re)
end

# Read-only data structure
Base.@propagate_inbounds @inline Base.getindex(L::ImplicitTerm, i::Int) = L.Δ[i]*L.ν

# Methods to satisfy the IMEXRKCB interface
Base.A_mul_B!(out::A, L::ImplicitTerm{n}, U::A) where {n, A<:AbstractFTField{n}} =
    (out .= L .* U)

ImcA!(L::ImplicitTerm{n}, c::Real, Y::A, Z::A) where {n, A<:AbstractFTField{n}} =
    (Z .= Y ./ (1 .- c .* L))