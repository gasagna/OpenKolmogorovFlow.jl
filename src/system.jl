export ForwardEquation, imex, ForwardMode, TangentMode

# solve either forward and variational equation or just the forward problem
abstract type AbstractMode end
abstract type AbstractForwardMode <: AbstractMode end
abstract type AbstractTangentMode <: AbstractMode end

# default tangent and forward modes
struct ForwardMode <: AbstractForwardMode end
struct TangentMode <: AbstractTangentMode end

# ~~~ Explicit Term of Forward Equations ~~~
struct ForwardExplicitTerm{n, FT<:AbstractFTField,     F<:AbstractField,
                              D1<:AbstractFTOperator, D2<:AbstractFTOperator,
                             ITT<:InverseFFT!,       FTT<:ForwardFFT!}
     FTFStore::Vector{FT} # storage
      FStore::Vector{F}
          dx::D1         # operators
          dy::D1
           Δ::D2
       ifft!::ITT        # transforms
        fft!::FTT
    kforcing::Int        # forcing wave number
end

# Outer constructor: TODO: merge with above
function ForwardExplicitTerm(n::Int, m::Int, kforcing::Int, mode::M, ::Type{S}, flags::UInt32) where {S, M<:AbstractMode}
    # stupidity check
    m ≥ n || throw(ArgumentError("`m` must be bigger than `n`, got `n, m= $n, $m`. Are you sure?"))

    # get appropriate constructor and eltype
    CT, C = M <: AbstractTangentMode ? (VariationalFTField, VariationalField) : (FTField,    Field)
    ET, E = M <: AbstractTangentMode ? (Dual{Complex{S}}, Dual{S})        : (Complex{S}, S)

    # complex fields have size `n` but real fields might have larger size `m`
    FTFStore = [CT(n, ET) for i = 1:4]
    FStore  = [ C(m, E)  for i = 1:4]

    # transforms
    ifft! = InverseFFT!(m, FTFStore[1], flags)
     fft! = ForwardFFT!(n,  FStore[1], flags)

    # operators
    dx, dy, Δ = DiffOperator(n, :x, S), DiffOperator(n, :y, S), DiffOperator(n, :xxyy, S)

    # construct object
    ForwardExplicitTerm{n, eltype(FTFStore), eltype(FStore), typeof(dx), typeof(Δ), 
                 typeof(ifft!), typeof(fft!)}(FTFStore, FStore, dx, dy, Δ, ifft!, fft!, kforcing)
end

# Adhere to callable interface for IMEXRKCB
function (Eq::ForwardExplicitTerm{n, FT})(t::Real, Ω::FT, dΩdt::FT, add::Bool=false) where {n, FT<:AbstractFTField{n}}
    # extract aliases
    dΩdx, dΩdy, U, V = Eq.FTFStore[1], Eq.FTFStore[2], Eq.FTFStore[3], Eq.FTFStore[4]
    dωdx, dωdy, u, v = Eq.FStore[1],   Eq.FStore[2],   Eq.FStore[3],   Eq.FStore[4]
    dx, dy, Δ = Eq.dx, Eq.dy, Eq.Δ

    # set mean to zero
    Ω[0, 0] = zero(eltype(Ω))

    # obtain vorticity derivatives
    dΩdx .= dx .* Ω
    dΩdy .= dy .* Ω

    # obtain velocity components. Set mean to zero.
    U .= .- dΩdy ./ Δ; U[0, 0] = zero(eltype(Ω))
    V .=    dΩdx ./ Δ; V[0, 0] = zero(eltype(Ω))

    # inverse transform to physical space into temporaries
    Eq.ifft!(u,    U)
    Eq.ifft!(v,    V)
    Eq.ifft!(dωdx, dΩdx)
    Eq.ifft!(dωdy, dΩdy)

    # multiply in physical space. Overwrite u
    u  .= .- u.*dωdx .- v.*dωdy

    # forward transform to Fourier space into destination
    add == true ? (Eq.fft!(U, u); dΩdt .+= U) : Eq.fft!(dΩdt, u)

    # ~~~ FORCING TERM ~~~
    dΩdt[ Eq.kforcing, 0] -= Eq.kforcing/2
    dΩdt[-Eq.kforcing, 0] -= Eq.kforcing/2

    return nothing
end


# ~~~ SOLVER OBJECT FOR THE GOVERNING EQUATIONS ~~~
struct ForwardEquation{n, IT<:ImplicitTerm{n}, ET<:ForwardExplicitTerm{n}}
    imTerm::IT
    exTerm::ET
end

# outer constructor: main entry point
function ForwardEquation(n::Int,
                         Re::Real,
                         kforcing::Int=4,
                         mode::AbstractMode=ForwardMode();
                         numtype::Type{S}=Float64,
                         flags::UInt32=FFTW.PATIENT,
                         dealias::Bool=true) where {S<:Real}
    iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
    m = dealias == true ? even_dealias_size(n) : n
    imTerm = ImplicitTerm(n, Re, S)
    exTerm = ForwardExplicitTerm(n, m, kforcing, mode, S, flags)
    ForwardEquation{n, typeof(imTerm), typeof(exTerm)}(imTerm, exTerm)
end

# evaluate right hand side of governing equations
(eq::ForwardEquation{n})(t::Real, Ω::FTField{n}, dΩdt::FTField{n}) where {n} =
    (A_mul_B!(dΩdt, eq.imTerm, Ω); eq.exTerm(t, Ω, dΩdt, true))

# obtain two components
imex(eq::ForwardEquation) = (eq.imTerm, eq.exTerm)