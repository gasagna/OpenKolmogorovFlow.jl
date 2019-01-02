# ////// FORCINGS FOR THE DIRECT, TANGENT OR ADJOINT EQUATIONS//////
abstract type AbstractForcing{n} end


# ////// DUMMY FORCING - DOES NOTHING //////
struct DummyForcing{n} <: AbstractForcing{n} end

DummyForcing(n::Int) = DummyForcing{n}()

(::DummyForcing{n})(t::Real,
                    Ω::FT,
                 dΩdt::FT) where {n, FT<:AbstractFTField{n}} = dΩdt

(::DummyForcing{n})(t::Real,
                    Ω::FT,
                    Λ::FT,
                 dΛdt::FT) where {n, FT<:AbstractFTField{n}} = dΛdt

 (::DummyForcing{n})(t::Real,
                     Ω::FT,
                  dΩdt::FT,
                     Λ::FT,
                  dΛdt::FT) where {n, FT<:AbstractFTField{n}} = dΛdt