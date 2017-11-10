export ForwardEquation, imex

# ~~~ Explicit Term of Forward Equations ~~~
struct ForwardExplicitTerm{n, FT<:FTField,             F<:Field,
                              D1<:AbstractFTOperator, D2<:AbstractFTOperator,
                             ITT<:InverseFFT!,       FTT<:ForwardFFT!}
     FTStore::Vector{FT} # storage
      FStore::Vector{F}
          dx::D1         # operators
          dy::D1
           Δ::D2
       ifft!::ITT        # transforms
        fft!::FTT
    kforcing::Int        # forcing wave number
end

# Outer constructor: TODO: merge with above
function ForwardExplicitTerm(n::Int, m::Int, kforcing::Int, ::Type{S}, flags::UInt32) where {S}
    # stupidity check
    m ≥ n || throw(ArgumentError("`m` must be bigger than `n`, got `n, m= $n, $m`. Are you sure?"))

    # complex fields have size `n` but real fields might have larger size `m`
    FTStore = [FTField(n, Complex{S}) for i = 1:4]
    FStore  = [  Field(m, S)          for i = 1:4]

    # transforms
    ifft! = InverseFFT!(m, FTStore[1], flags)
     fft! = ForwardFFT!(n,  FStore[1], flags)

    # operators
    dx, dy, Δ = DiffOperator(n, :x, S), DiffOperator(n, :y, S), DiffOperator(n, :xxyy, S)

    # construct object
    ForwardExplicitTerm{n, eltype(FTStore), eltype(FStore), typeof(dx), typeof(Δ), 
                 typeof(ifft!), typeof(fft!)}(FTStore, FStore, dx, dy, Δ, ifft!, fft!, kforcing)
end

# Adhere to callable interface for IMEXRKCB
function (Eq::ForwardExplicitTerm{n, FT})(t::Real, Ω::FT, dΩdt::FT, add::Bool=false) where {n, FT<:FTField{n}}
    # extract aliases
    dΩdx, dΩdy, U, V = Eq.FTStore[1], Eq.FTStore[2], Eq.FTStore[3], Eq.FTStore[4]
    dωdx, dωdy, u, v = Eq.FStore[1],  Eq.FStore[2],  Eq.FStore[3],  Eq.FStore[4]
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
                         kforcing::Int=4;
                         numtype::Type{S}=Float64,
                         flags::UInt32=FFTW.PATIENT,
                         dealias::Bool=true) where {S<:Real}
    iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
    m = dealias == true ? even_dealias_size(n) : n
    imTerm = ImplicitTerm(n, Re, S)
    exTerm = ForwardExplicitTerm(n, m, kforcing, S, flags)
    ForwardEquation{n, typeof(imTerm), typeof(exTerm)}(imTerm, exTerm)
end

# evaluate right hand side of governing equations
(eq::ForwardEquation{n})(t::Real, Ω::FTField{n}, dΩdt::FTField{n}) where {n} =
    (A_mul_B!(dΩdt, eq.imTerm, Ω); eq.exTerm(t, Ω, dΩdt, true))

# obtain two components
imex(eq::ForwardEquation) = (eq.imTerm, eq.exTerm)