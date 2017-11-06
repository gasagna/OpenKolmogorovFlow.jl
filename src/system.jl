# add method to this function
import IMEXRKCB: ImcA!

export VorticityEquation, 
       imex, 
       ForwardMode,
       TangentMode

# ~~~ DISPATCH BASED ON SOLUTION MODE ~~~
abstract type AbstractMode end
abstract type AbstractAugmentedMode <: AbstractMode end
abstract type AbstractAdjointMode   <: AbstractAugmentedMode end

struct ForwardMode <: AbstractMode end
struct TangentMode <: AbstractAugmentedMode end


# ~~~ THE VISCOUS TERM OF THE GOVERNING EQUATIONS ~~~
struct ImplicitTerm{n, T, D<:DiffOperator{n, T}} <: AbstractFTOperator{n, T}
    Δ::D
    ν::Float64 # inverse of Reynolds number
end

# Outer constructor
function ImplicitTerm(n::Int, Re::Real, ::Type{T}) where {T} 
    Δ = DiffOperator(n, :xxyy, T)
    ImplicitTerm{n, T, typeof(Δ)}(Δ, 1/Re)
end

# Read-only data structure
Base.@propagate_inbounds @inline Base.getindex(L::ImplicitTerm, i::Int) = L.Δ[i]*L.ν

# Methods to satisfy the IMEXRKCB interface
Base.A_mul_B!(out::A, L::ImplicitTerm{n}, U::A) where {n, A<:AbstractFTField{n}} =
    (out .= L .* U)

ImcA!(L::ImplicitTerm{n}, c::Real, Y::A, Z::A) where {n, A<:AbstractFTField{n}} =
    (Z .= Y ./ (1 .- c .* L))

# ~~~ THE NONLINEAR TERM OF THE GOVERNING EQUATIONS PLUS THE FORCING ~~~

# We could be stricter on the definition of type below, but would lead to more 
# complicated code, without too much benefit. It is unlikely that anyone will
# instantiate this directly, so we avoid over constraining input parameters
struct ExplicitTerm{n, M<:AbstractMode,                          
                    FT<:AbstractFTField,    F<:AbstractField,       
                    D1<:AbstractFTOperator, D2<:AbstractFTOperator, 
                    ITT<:InverseFFT!,       FTT<:ForwardFFT!}
           U::FT;     V::FT; dΩdx::FT; dΩdy::FT # storage
           u::F;      v::F;  dωdx::F;  dωdy::F  # storage
          dx::D1;    dy::D1;    Δ::D2           # operators
       ifft!::ITT; fft!::FTT                    # transforms
    kforcing::Int                               # forcing wave number
        mode::M                                 # contains solver specific actions
end

# Outer constructor
function ExplicitTerm(n::Int, m::Int, mode::M, kforcing::Int, ::Type{S}, flags::UInt32) where {S, M<:AbstractMode}
    # stupid check
    m ≥ n || throw(ArgumentError("`m` must be bigger than `n`, got `n, m= $n, $m`. Are you sure?"))

    # get appropriate constructor and eltype
    CT, C = M <: AbstractAugmentedMode ? (AugmentedFTField, AugmentedField) : (FTField,    Field)
    ET, E = M <: AbstractAugmentedMode ? (Dual{Complex{S}}, Dual{S})        : (Complex{S}, S)

    # complex fields have size `n` but real fields might have larger size `m`
    a, b, c, d = CT.((n, n, n, n), ET)
    f, g, h, i =  C.((m, m, m, m), E)

    # transforms
    ifft! = InverseFFT!(m, a, flags)
     fft! = ForwardFFT!(n, f, flags)

    # operators
    dx, dy, Δ = DiffOperator(n, :x, S), DiffOperator(n, :y, S), DiffOperator(n, :xxyy, S)

    # construct object
    ExplicitTerm{n, M, typeof(a), typeof(f), 
                 typeof(dx), typeof(Δ), 
                 typeof(ifft!), typeof(fft!)}(a, b, c, d, 
                                              f, g, h, i, 
                                              dx, dy, Δ, 
                                              ifft!, fft!, kforcing, mode)
end

# This satisfies a callable interface suitable for IMEXRKCB. The fourth
# argument `add` specifies whether `dΩdt` is overwritten or added to.
function (Eq::ExplicitTerm{n, M, FT})(t::Real, Ω::FT, dΩdt::FT, add::Bool=false) where {n, M, FT<:AbstractFTField{n}}
    # set mean to zero
    Ω[0, 0] = zero(eltype(Ω))

    # obtain vorticity derivatives
    Eq.dΩdx .= Eq.dx .* Ω
    Eq.dΩdy .= Eq.dy .* Ω

    # obtain velocity components. Set mean to zero.
    Eq.U .= .- Eq.dΩdy ./ Eq.Δ; Eq.U[0, 0] = zero(eltype(Ω))
    Eq.V .=    Eq.dΩdx ./ Eq.Δ; Eq.V[0, 0] = zero(eltype(Ω))

    # inverse transform to physical space into temporaries
    Eq.ifft!(Eq.u,    Eq.U)
    Eq.ifft!(Eq.v,    Eq.V)
    Eq.ifft!(Eq.dωdx, Eq.dΩdx)
    Eq.ifft!(Eq.dωdy, Eq.dΩdy)

    # multiply in physical space. Overwrite u
    Eq.u  .= .- Eq.u.*Eq.dωdx .- Eq.v.*Eq.dωdy

    # forward transform to Fourier space into destination
    if add == true
        Eq.fft!(Eq.U, Eq.u); dΩdt .+= Eq.U
    else
        Eq.fft!(dΩdt,    Eq.u)
    end

    # ~~~ FORCING TERM ~~~
    dΩdt[ Eq.kforcing, 0] -= Eq.kforcing/2
    dΩdt[-Eq.kforcing, 0] -= Eq.kforcing/2

    return nothing
end


# ~~~ SOLVER OBJECT FOR THE GOVERNING EQUATIONS ~~~

struct VorticityEquation{n, IT<:ImplicitTerm{n}, ET<:ExplicitTerm{n}}
    imTerm::IT
    exTerm::ET
end

# outer constructor: main entry point
function VorticityEquation(n::Int,
                           Re::Real,
                           kforcing::Int=4;
                           mode::AbstractMode=ForwardMode(),
                           numtype::Type{S}=Float64,
                           flags::UInt32=FFTW.PATIENT,
                           dealias::Bool=true) where {S<:Real}
    iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
    m = dealias == true ? even_dealias_size(n) : n
    impl = ImplicitTerm(n, Re, S)
    expl = ExplicitTerm(n, m, mode, kforcing, S, flags)
    VorticityEquation{n, typeof(impl), typeof(expl)}(impl, expl)
end

# evaluate right hand side of governing equations
function (eq::VorticityEquation{n})(t::Real, Ω::FTField{n}, dΩdt::FTField{n}) where {n}
    A_mul_B!(dΩdt, eq.imTerm, Ω)
    eq.exTerm(t, Ω, dΩdt, true)
    return dΩdt
end

# obtain two components
imex(eq::VorticityEquation) = (eq.imTerm, eq.exTerm)