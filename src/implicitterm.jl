# add method to this function
import Flows

# ~~~ THE VISCOUS TERM OF THE GOVERNING EQUATIONS ~~~
struct ImplicitTerm
    Re::Float64
end

# Methods to satisfy the Flows interface
Base.A_mul_B!(OUT::FT,
              L::ImplicitTerm, U::FT) where {FT<:AbstractFTField} =
    (laplacian!(OUT, U); OUT .*= 1/L.Re; OUT)

Flows.ImcA!(L::ImplicitTerm,
            c::Real,
            Y::FT,
            OUT::FT) where {FT<:AbstractFTField} = invlaplacian!(OUT, Y, c/L.Re)