# add method to this function
import Flows

import LinearAlgebra: mul!

# ~~~ THE VISCOUS TERM OF THE GOVERNING EQUATIONS ~~~

# allow changin the Reynolds number more easily
mutable struct ImplicitTerm
    Re::Float64
end

# Methods to satisfy the Flows interface
mul!(OUT::FT,
     L::ImplicitTerm, U::FT) where {FT<:AbstractFTField} =
    (laplacian!(OUT, U); OUT .*= 1/L.Re; OUT)

Flows.ImcA!(L::ImplicitTerm,
            c::Real,
            Y::FT,
            OUT::FT) where {FT<:AbstractFTField} = invlaplacian!(OUT, Y, c/L.Re)