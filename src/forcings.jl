# ////// FORCINGS FOR THE DIRECT, TANGENT OR ADJOINT EQUATIONS//////

export DummyForcing, DissRateGradientForcing

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

# ////// FORCING FOR ADJOINT EQUATION FOR DISSIPATION RATE AS THE FUNCTIONAL //////
struct DissRateGradientForcing{n} <: AbstractForcing{n} 
    Re::Float64
end

DissRateGradientForcing(n::Int, Re::Real) = DissRateGradientForcing{n}(Re)

(f::DissRateGradientForcing{n})(t::Real,
                                Ω::FT,
                                Λ::FT,
                             dΛdt::FT) where {n, FT<:AbstractFTField{n}} = 
    (dΛdt .+= 2.0.*Ω./f.Re; dΛdt)