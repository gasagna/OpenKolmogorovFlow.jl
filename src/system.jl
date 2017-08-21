# add method to this function
import IMEXRKCB: ImcA!

export ImplicitTerm, ExplicitTerm, VorticityEquation, imex

# ~~~ THE VISCOUS TERM OF THE GOVERNING EQUATIONS ~~~
struct ImplicitTerm{n, S}
    dx²dy²::DiffOperator{n, S}
    ν::Float64                 # inverse of Reynolds number
    function ImplicitTerm{n, S}(Re::Real) where {n, S}
        iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
        new(DiffOperator(n, :xxyy, S), 1/Re)
    end
end

# Outer constructor
ImplicitTerm(n::Int, Re::Real, ::Type{S}) where {S} = ImplicitTerm{n, S}(Re)

# Methods to satisfy the IMEXRKCB interface
Base.A_mul_B!(out::FTField{n}, V::ImplicitTerm{n}, U::FTField{n}) where {n} =
    (ν = V.ν; out .= V.dx²dy² .* U .* ν)

ImcA!(V::ImplicitTerm{n}, c::Real, y::FTField{n}, z::FTField{n}) where {n} =
    (z .= y ./ (1 .- c .* V.dx²dy² .* V.ν))

# ~~~ THE NONLINEAR TERM OF THE GOVERNING EQUATIONS PLUS THE FORCING ~~~
struct ExplicitTerm{n, m, T<:Real, S<:Real, IT<:InverseFFT!, FT<:ForwardFFT!}
       ifft!::IT
        ftt!::FT
    kforcing::Int
           U::FTField{n, Complex{T}, Matrix{Complex{T}}}
           V::FTField{n, Complex{T}, Matrix{Complex{T}}}
        ∂Ω∂x::FTField{n, Complex{T}, Matrix{Complex{T}}}
        ∂Ω∂y::FTField{n, Complex{T}, Matrix{Complex{T}}}
          dx::DiffOperator{n, Complex{S}}
          dy::DiffOperator{n, Complex{S}}
      dx²dy²::DiffOperator{n, S}
           u::Field{m, T, Matrix{T}}
           v::Field{m, T, Matrix{T}}
        ∂ω∂x::Field{m, T, Matrix{T}}
        ∂ω∂y::Field{m, T, Matrix{T}}
    function ExplicitTerm{n, m, T, S}(kforcing::Int, flags::UInt32) where {n, m, T, S}
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

        new{n, m, T, S, typeof(ifft!), typeof(fft!)}(ifft!, fft!, kforcing, 
            a, b, c, d, 
            DiffOperator(n, :x, S),  DiffOperator(n, :y, S), DiffOperator(n, :xxyy, S),
            f, g, h, i)
    end
end

# Outer constructor
ExplicitTerm(n::Int, m::Int, kforcing::Int, ::Type{T}, ::Type{S}, flags::UInt32) where {T, S} =
     ExplicitTerm{n, m, T, S}(kforcing, flags)

function (Eq::ExplicitTerm{n})(t::Real, Ω::FTField{n}, Ω̇::FTField{n}, add::Bool=false) where {n}
    # ~~~ PRELIMINARIES ~~~
    # set mean to zero
    Ω[0, 0] = 0.0

    # obtain vorticity derivatives
    Eq.∂Ω∂x .= Eq.dx .* Ω
    Eq.∂Ω∂y .= Eq.dy .* Ω

    # obtain velocity components. Set mean to zero.
    Eq.U .= .- (Eq.dx²dy²) .\ Eq.∂Ω∂y; Eq.U[0, 0] = 0
    Eq.V .=    (Eq.dx²dy²) .\ Eq.∂Ω∂x; Eq.V[0, 0] = 0

    # ~~~ NONLINEAR TERM ~~
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
struct VorticityEquation{n, m, T<:Real, S<:Real}
    imTerm::ImplicitTerm{n, S}
    exTerm::ExplicitTerm{n, m, T, S}
    function VorticityEquation{n, m, T, S}(Re::Real,
                                           kforcing::Int, 
                                           flags::UInt32) where {n, m, T, S}
        iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
        iseven(m) || throw(ArgumentError("`m` must be even, got $m"))
        new(ImplicitTerm(n, Re, S), ExplicitTerm(n, m, kforcing, T, S, flags))
    end
end

# outer constructor: main entry point
function VorticityEquation(n::Int, 
                           Re::Real, 
                           kforcing::Int=4; 
                           T::Type{<:Real}=Float64,
                           flags::UInt32=FFTW.MEASURE, 
                           dealias::Union{Bool, Int}=true)
    # if dealias is an int we interpret it as the size of the
    # larger grid over which we do interpolation, and assume that
    # the user knows what he/she is doing. If it is a boolean we 
    # select the appropriate value. For optimal performance the user 
    # should select an appropriate pair of grid size, based on tests of 
    # the FFTW compiled library on his/her machine.
    m = dealias isa Int ? dealias : (dealias == true ? even_dealias_size(n) : n)
    # compute eltype of differential operator data if we have a 
    # variational number type as input
    S = T <: VarNum ? T.parameters[1] : T
    VorticityEquation{n, m, T, S}(Re, kforcing, flags)
end

# evaluate right hand side of governing equations
function (eq::VorticityEquation{n})(t::Real, Ω::FTField{n}, Ω̇::FTField{n}) where {n}
    A_mul_B!(Ω̇, eq.imTerm, Ω)
    eq.exTerm(t, Ω, Ω̇, true)
end

# obtain two components
imex(eq::VorticityEquation) = (eq.imTerm, eq.exTerm)