# ////// FORCINGS FOR THE DIRECT, TANGENT OR ADJOINT EQUATIONS//////

export DummyForcing, DissRateGradientForcing, WaveNumberForcing

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

# ////// FORCING AT ONE SPECIFIC WAVENUMBER  //////
struct WaveNumberForcing{n, T} <: AbstractForcing{n} 
      k::Int
      j::Int
    val::T
end

WaveNumberForcing(n::Int, k::Int, j::Int, f::T) where {T<:Number} = 
    WaveNumberForcing{n, T}(k, j, f)

# provide multiple methods for the same function
(f::WaveNumberForcing{n})(t::Real,
                          Ω::FT,
                       dΩdt::FT) where {n, FT<:AbstractFTField{n}} = f(t, Ω, Ω, dΩdt)

(f::WaveNumberForcing{n})(t::Real,
                          Ω::FT,
                          Λ::FT,
                       dΛdt::FT) where {n, FT<:AbstractFTField{n}} = f(t, Ω, Ω, Λ, dΛdt)
                          
function (f::WaveNumberForcing{n})(t::Real,
                                   Ω::FT,
                                dΩdt::FT,
                                   Λ::FT,
                                dΛdt::FT) where {n, FT<:AbstractFTField{n}}
                 dΛdt[WaveNumber( f.k, f.j)] +=      f.val/2
    f.j == 0 && (dΛdt[WaveNumber(-f.k, f.j)] += conj(f.val)/2)
    return dΛdt
end

# ////// FORCING FOR ADJOINT EQUATION FOR DISSIPATION RATE AS THE FUNCTIONAL //////
struct DissRateGradientForcing{n} <: AbstractForcing{n} 
    Re::Float64
end

DissRateGradientForcing(n::Int, Re::Real) = DissRateGradientForcing{n}(Re)

(f::DissRateGradientForcing{n})(t::Real,
                                Ω::FT,
                                Λ::FT,
                             dΛdt::FT) where {n, FT<:AbstractFTField{n}} = 
    (dΛdt .-= 2.0.*Ω./f.Re; dΛdt)