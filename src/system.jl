# add method to this function
import IMEXRKCB: ImcA!

export ImplicitTerm, ExplicitTerm, VorticityEquation, imex

# ~~~ THE VISCOUS TERM OF THE GOVERNING EQUATIONS ~~~
struct ImplicitTerm{n, T}
    dx²dy²::DiffOperator{n, T}
    ν::Float64                 # inverse of Reynolds number
    function ImplicitTerm{n, T}(Re::Real) where {n, T}
        iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
        new(DiffOperator(n, :xxyy, T), 1/Re)
    end
end

# Outer constructor
ImplicitTerm(n::Int, Re::Real, ::Type{T}) where {T} = ImplicitTerm{n, T}(Re)

# Methods to satisfy the IMEXRKCB interface
Base.A_mul_B!(out::FTField{n}, V::ImplicitTerm{n}, U::FTField{n}) where {n} =
    (ν = V.ν; out .= V.dx²dy² .* U .* ν)

ImcA!(V::ImplicitTerm{n}, c::Real, y::FTField{n}, z::FTField{n}) where {n} =
    (z .= y ./ (1 .- c .* V.dx²dy² .* V.ν))

# ~~~ THE NONLINEAR TERM OF THE GOVERNING EQUATIONS PLUS THE FORCING ~~~
struct ExplicitTerm{n, m, T<:Real, IT<:InverseFFT!, FT<:ForwardFFT!}
       ifft!::IT
        ftt!::FT
    kforcing::Int
           U::FTField{n, Complex{T}, Matrix{Complex{T}}}
           V::FTField{n, Complex{T}, Matrix{Complex{T}}}
        ∂Ω∂x::FTField{n, Complex{T}, Matrix{Complex{T}}}
        ∂Ω∂y::FTField{n, Complex{T}, Matrix{Complex{T}}}
          dx::DiffOperator{n, Complex{T}}
          dy::DiffOperator{n, Complex{T}}
      dx²dy²::DiffOperator{n, T}
           u::Field{m, T, Matrix{T}}
           v::Field{m, T, Matrix{T}}
        ∂ω∂x::Field{m, T, Matrix{T}}
        ∂ω∂y::Field{m, T, Matrix{T}}
    function ExplicitTerm{n, m, T}(kforcing::Int, flags::UInt32) where {n, m, T}
        iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
        iseven(m) || throw(ArgumentError("`m` must be even, got $m"))
        m ≥ n     || throw(ArgumentError("`m` must be bigger than `n`, got `n, m= $n, $m`. Are you sure?"))

        # complex fields have size n
        a, b, c, d = FTField.((n, n, n, n), Complex{T})

        # real fields (might) have (larger) size m
        f, g, h, i = Field.((m, m, m, m), T)

        # transforms
        ifft! = InverseFFT!(Field{m, T}, a,  flags)
         fft! = ForwardFFT!(FTField{n, Complex{T}}, f,  flags)

        new{n, m, T, typeof(ifft!), typeof(fft!)}(ifft!, fft!, kforcing,
            a, b, c, d,
            DiffOperator(n, :x, T),  DiffOperator(n, :y, T), DiffOperator(n, :xxyy, T),
            f, g, h, i)
    end
end

# Outer constructor
ExplicitTerm(n::Int, m::Int, kforcing::Int, ::Type{T}, flags::UInt32) where {T} =
     ExplicitTerm{n, m, T}(kforcing, flags)

function (Eq::ExplicitTerm{n})(t::Real, Ω::FTField{n}, Ω̇::FTField{n}, add::Bool=false) where {n}
    # set mean to zero
    Ω[0, 0] = zero(eltype(Ω))

    # obtain vorticity derivatives
    Eq.∂Ω∂x .= Eq.dx .* Ω
    Eq.∂Ω∂y .= Eq.dy .* Ω

    # obtain velocity components. Set mean to zero.
    Eq.U .= .- (Eq.dx²dy²) .\ Eq.∂Ω∂y; Eq.U[0, 0] = zero(eltype(Eq.U))
    Eq.V .=    (Eq.dx²dy²) .\ Eq.∂Ω∂x; Eq.V[0, 0] = zero(eltype(Eq.V))

    # inverse transform to physical space into temporaries
    Eq.ifft!(Eq.u,    Eq.U)
    Eq.ifft!(Eq.v,    Eq.V)
    Eq.ifft!(Eq.∂ω∂x, Eq.∂Ω∂x)
    Eq.ifft!(Eq.∂ω∂y, Eq.∂Ω∂y)

    # multiply in physical space. Overwrite u
    Eq.u .= .- Eq.u.*Eq.∂ω∂x .- Eq.v.*Eq.∂ω∂y

    # forward transform to Fourier space into destination
    if add == true
        Eq.ftt!(Eq.U, Eq.u); Ω̇ .+= Eq.U
    else
        Eq.ftt!(Ω̇,    Eq.u)
    end

    # ~~~ FORCING TERM ~~~
    Ω̇[ Eq.kforcing, 0] -= Eq.kforcing/2
    Ω̇[-Eq.kforcing, 0] -= Eq.kforcing/2

    return nothing
end


# ~~~ THE GOVERNING EQUATIONS ~~~
struct VorticityEquation{n, m, T<:Real}
    imTerm::ImplicitTerm{n, T}
    exTerm::ExplicitTerm{n, m, T}
    function VorticityEquation{n, m, T}(Re::Real,
                                        kforcing::Int,
                                        flags::UInt32) where {n, m, T}
        iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
        iseven(m) || throw(ArgumentError("`m` must be even, got $m"))
        new(ImplicitTerm(n, Re, T), ExplicitTerm(n, m, kforcing, T, flags))
    end
end

# outer constructor: main entry point
function VorticityEquation(n::Int,
                           Re::Real,
                           kforcing::Int=4;
                           T::Type{<:Real}=Float64,
                           flags::UInt32=FFTW.PATIENT,
                           dealias::Bool=true)
    m = dealias == true ? even_dealias_size(n) : n
    VorticityEquation{n, m, T}(Re, kforcing, flags)
end

# evaluate right hand side of governing equations
function (eq::VorticityEquation{n})(t::Real, Ω::FTField{n}, Ω̇::FTField{n}) where {n}
    A_mul_B!(Ω̇, eq.imTerm, Ω)
    eq.exTerm(t, Ω, Ω̇, true)
end

# obtain two components
imex(eq::VorticityEquation) = (eq.imTerm, eq.exTerm)