export AdjointEquation, imex

# ~~~ Explicit Term of Adjoint Equations ~~~
struct AdjointExplicitTerm{n, FT<:AbstractFTField,    F<:AbstractField,
                              D1<:AbstractFTOperator, D2<:AbstractFTOperator,
                             ITT<:InverseFFT!,       FTT<:ForwardFFT!, S, C}
    FTFStore::Vector{FT} # storage
      FStore::Vector{F}
          dx::D1         # operators
          dy::D1
           Δ::D2
       ifft!::ITT        # transforms
        fft!::FTT
         sol::S          # forward solution storage
        cost::C          # cost gradient function
end

# Outer constructor: TODO: merge with above
function AdjointExplicitTerm(n::Int, m::Int, sol, cost, ::Type{S}, flags::UInt32) where {S}
    # stupidity check
    m ≥ n || throw(ArgumentError("`m` must be bigger than `n`, got `n, m= $n, $m`. Are you sure?"))

    # complex fields have size `n` but real fields might have larger size `m`
    FTFStore = [FTField(n, Complex{S}) for i = 1:7]
    FStore   = [  Field(m, S)          for i = 1:6]

    # transforms
    ifft! = InverseFFT!(m, FTFStore[1], flags)
     fft! = ForwardFFT!(n,   FStore[1], flags)

    # operators
    dx, dy, Δ = DiffOperator(n, :x, S), DiffOperator(n, :y, S), DiffOperator(n, :xxyy, S)

    # construct object
    AdjointExplicitTerm{n, 
                        eltype(FTFStore), 
                        eltype(FStore), 
                        typeof(dx), 
                        typeof(Δ), 
                        typeof(ifft!), 
                        typeof(fft!), 
                        typeof(sol), 
                        typeof(cost)}(FTFStore, FStore, dx, dy, Δ, ifft!, fft!, sol, cost)
end

# extract forward solution at current time
(Eq::AdjointExplicitTerm{n, FT})(t::Real, Λ::FT, dΛdt::FT, add::Bool=false) where {n, FT<:FTField{n}} = 
    Eq(t, Eq.sol(Eq.FTFStore[1], t), Λ, dΛdt)

# Adhere to callable interface for IMEXRKCB
function (Eq::AdjointExplicitTerm{n, FT})(t::Real, Ω::FT, Λ::FT, dΛdt::FT, add::Bool=false) where {n, FT<:FTField{n}}
    # extract aliases. We start from one because 1 might be taken for Ω
                U,    V    = Eq.FTFStore[2], Eq.FTFStore[3]
                TMP1, TMP2 = Eq.FTFStore[2], Eq.FTFStore[3]
    dΛdx, dΛdy, dΩdx, dΩdy = Eq.FTFStore[4], Eq.FTFStore[5], Eq.FTFStore[6], Eq.FTFStore[7] 
                u,    v    = Eq.FStore[1],   Eq.FStore[2]
                tmp1, tmp2 = Eq.FStore[1],   Eq.FStore[2]
    dλdx, dλdy, dωdx, dωdy = Eq.FStore[3],   Eq.FStore[4],   Eq.FStore[5],   Eq.FStore[6]
    dx, dy, Δ = Eq.dx, Eq.dy, Eq.Δ

    # set mean to zero
    Λ[0, 0] = zero(eltype(Λ))

    # obtain vorticity derivatives
    dΛdx .= dx .* Λ
    dΛdy .= dy .* Λ
    dΩdx .= dx .* Ω
    dΩdy .= dy .* Ω

    # obtain velocity components. Set mean to zero
    U .= .- dΩdy ./ Δ; U[0, 0] = zero(eltype(U))
    V .=    dΩdx ./ Δ; V[0, 0] = zero(eltype(V))

    # inverse transform to physical space into temporaries
    Eq.ifft!(u,    U)
    Eq.ifft!(v,    V)
    Eq.ifft!(dλdx, dΛdx); Eq.ifft!(dωdx, dΩdx)
    Eq.ifft!(dλdy, dΛdy); Eq.ifft!(dωdy, dΩdy)

    # multiply in physical space
    tmp1 .= u.*dλdx .+ v.*dλdy;       Eq.fft!(TMP1, tmp1)
    tmp2 .= dλdx.*dωdy .- dλdy.*dωdx; Eq.fft!(TMP2, tmp2)
    add == true ? (dΛdt .+= .- TMP1 .- TMP2 ./ Δ) : (dΛdt  .= .- TMP1 .- TMP2 ./ Δ)
    dΛdt[0, 0] = zero(eltype(Λ))
    
    # add forcing from cost gradient
    Eq.cost(dΛdt, Ω)

    return nothing
end


# ~~~ SOLVER OBJECT FOR THE ADJOINT EQUATIONS ~~~
struct AdjointEquation{n, IT<:ImplicitTerm{n}, ET<:AdjointExplicitTerm{n}}
    imTerm::IT
    exTerm::ET
end

# outer constructor: main entry point
function AdjointEquation(n::Int,
                         Re::Real,
                         sol,             # the forward solution
                         cost,            # cost gradient function
                         kforcing::Int=4;
                         numtype::Type{S}=Float64,
                         flags::UInt32=FFTW.PATIENT,
                         dealias::Bool=true) where {S<:Real}
    iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
    m = dealias == true ? even_dealias_size(n) : n
    imTerm = ImplicitTerm(n, -Re, S)
    exTerm = AdjointExplicitTerm(n, m, sol, cost, S, flags)
    AdjointEquation{n, typeof(imTerm), typeof(exTerm)}(imTerm, exTerm)
end

# evaluate right hand side of governing equations
(eq::AdjointEquation{n})(t::Real, Λ::FTField{n}, dΛdt::FTField{n}) where {n} =
    (A_mul_B!(dΛdt, eq.imTerm, Λ); eq.exTerm(t, Λ, dΛdt, true))

# obtain two components
imex(eq::AdjointEquation) = (eq.imTerm, eq.exTerm)





# ~~~ SOLVER OBJECT FOR THE ADJOINT EQUATIONS ~~~
struct AdjointEquation{n, IT<:ImplicitTerm{n}, ET<:LinearisedExplicitTerm{n}}
    imTerm::IT
    exTerm::ET
end

# outer constructor: main entry point
function AdjointEquation(n::Int,
                         Re::Real,
                         sol,             # the forward solution
                         cost,            # cost gradient function
                         kforcing::Int=4;
                         numtype::Type{S}=Float64,
                         flags::UInt32=FFTW.PATIENT,
                         dealias::Bool=true) where {S<:Real}
    iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
    m = dealias == true ? even_dealias_size(n) : n
    imTerm = ImplicitTerm(n, -Re, S)
    exTerm = LinearisedExplicitTerm(n, m, sol, cost, S, flags)
    AdjointEquation{n, typeof(imTerm), typeof(exTerm)}(imTerm, exTerm)
end

# evaluate right hand side of governing equations
(eq::AdjointEquation{n})(t::Real, Λ::FTField{n}, dΛdt::FTField{n}) where {n} =
    (A_mul_B!(dΛdt, eq.imTerm, Λ); eq.exTerm(t, Λ, dΛdt, true))

# obtain two components
imex(eq::AdjointEquation) = (eq.imTerm, eq.exTerm)

# multiply in physical space
u .= u.*dω′dx .+ v.*dω′dy.+ u′.*dωdx .+ v′.*dωdy
