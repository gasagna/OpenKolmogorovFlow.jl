# ////// FORCINGS FOR THE DIRECT, TANGENT OR ADJOINT EQUATIONS//////
abstract type AbstractForcing{n} end


# ////// DUMMY FORCING - DOES NOTHING //////
struct DummyForcing{n} <: AbstractForcing{n} end

DummyForcing(n::Int) = DummyForcing{n}()

# call for nonlinear equation
(::DummyForcing{n})(t, 
                    Ω::FT, dΩdt::FT) where {n, FT<:AbstractFTField{n}} = dΩdt

# and for the linear equation
(::DummyForcing{n})(t, Ω::FT, 
                    Ω′::FT, dΩ′dt::FT) where {n, FT<:AbstractFTField{n}} = dΩ′dt