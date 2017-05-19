# add method to this function
import IMEXRKCB: ImcA!

export ImplicitTerm, ExplicitTerm

# The viscous term of the governing equations
struct ImplicitTerm{n}
    dx²::DiffOperator{n, Int}
    dy²::DiffOperator{n, Int}
    ν::Float64 # inverse of Reynolds number
    ImplicitTerm{n}(Re::Real) where n =
        new{n}(DiffOperator(n, :xx), DiffOperator(n, :yy), 1/Re)
end

# Outer constructor
ImplicitTerm(n::Int, Re::Real) = ImplicitTerm{n}(Re)

# Methods to satisfy the IMEXRKCB interface
Base.A_mul_B!(out::FTField{n}, V::ImplicitTerm{n}, U::FTField{n}) where n =
    (out .= (V.dx² .+ V.dy²) .* U .* V.ν)

ImcA!(V::ImplicitTerm{n}, c::Real, y::FTField{n}, z::FTField{n}) where n =
    (z .= y ./ (1 .- c .* (V.dx² .+ V.dy²) .* V.ν))


# The nonlinear term of the governing equations plus the forcing
struct ExplicitTerm{n, T}
kforcing::Int
       U::FTField{n, Complex{T}, Matrix{Complex{T}}}
       V::FTField{n, Complex{T}, Matrix{Complex{T}}}
    ∂Ω∂x::FTField{n, Complex{T}, Matrix{Complex{T}}}
    ∂Ω∂y::FTField{n, Complex{T}, Matrix{Complex{T}}}
      dx::DiffOperator{n, Complex{Int}}
      dy::DiffOperator{n, Complex{Int}}
     dx²::DiffOperator{n, Int}
     dy²::DiffOperator{n, Int}
       u::Field{n, T, Matrix{T}}
       v::Field{n, T, Matrix{T}}
    ∂ω∂x::Field{n, T, Matrix{T}}
    ∂ω∂y::Field{n, T, Matrix{T}}
  InvFFT
 ForwFFT
    function ExplicitTerm{n, T}(kforcing::Int, FFTWflags::UInt32) where {n, T}
        a, b = FTField(n, Complex{T}), FTField(n, Complex{T})
        c, d = FTField(n, Complex{T}), FTField(n, Complex{T})
        f, g, h, i = Field(n, T), Field(n, T), Field(n, T), Field(n, T)
        InvFFT, ForwFFT = InverseFFT(a, FFTWflags), ForwardFFT(f, FFTWflags)
        new{n, T}(kforcing, a, b, c, d, DiffOperator(n, :x), DiffOperator(n, :y), 
                  DiffOperator(n, :xx), DiffOperator(n, :yy), f, g, h, i, 
                  InvFFT, ForwFFT)
    end
end

ExplicitTerm(n::Int, kforcing::Int=4, FFTWflags::UInt32=FFTW.MEASURE, ::Type{T}=Float64) where {T} = 
    ExplicitTerm{n, T}(kforcing, FFTWflags)

function (Eq::ExplicitTerm{n})(t::Real, Ω::FTField{n}, ∂Ω∂t::FTField{n}) where {n}
    # ~~~ PRELIMINARIES ~~~
    # set mean to zero
    Ω[0, 0] = 0.0

    # obtain vorticity derivatives
    Eq.∂Ω∂x .= Eq.dx .* Ω
    Eq.∂Ω∂y .= Eq.dy .* Ω

    # obtain velocity components. Set mean to zero.
    Eq.U .= .- (Eq.dx² .+ Eq.dy²) .\ Eq.∂Ω∂y; Eq.U[0, 0] = 0
    Eq.V .=    (Eq.dx² .+ Eq.dy²) .\ Eq.∂Ω∂x; Eq.V[0, 0] = 0

    # ~~~ NONLINEAR TERM ~~
    # inverse transform to physical space into temporaries
    Eq.InvFFT(   Eq.U, Eq.u)
    Eq.InvFFT(   Eq.V, Eq.v)
    Eq.InvFFT(Eq.∂Ω∂x, Eq.∂ω∂x)
    Eq.InvFFT(Eq.∂Ω∂y, Eq.∂ω∂y)

    # multiply in physical space. Overwrite u
    Eq.u .= .- Eq.u.*Eq.∂ω∂x .- Eq.v.*Eq.∂ω∂y

    # forward transform to Fourier space into destination
    Eq.ForwFFT(Eq.u, ∂Ω∂t)

    # ~~~ FORCING TERM ~~~
    ∂Ω∂t[ Eq.kforcing, 0] -= Eq.kforcing/2
    ∂Ω∂t[-Eq.kforcing, 0] -= Eq.kforcing/2
    return nothing
end
