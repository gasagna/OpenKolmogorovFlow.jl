export LinearisedEquation,
       splitexim

# This triggers a different evaluation of the linearised
# operator, but the required software infrastructure is
# the same in both cases.
abstract type AbstractLinearMode end

struct TangentMode <: AbstractLinearMode end
struct AdjointMode <: AbstractLinearMode end


# ~~~ Explicit Term of Linearised Equations ~~~
struct LinearisedExTerm{n,
                        m,
                        M<:AbstractLinearMode,
                        FT<:FTField{n, m},
                        F<:Field{m},
                        ITT<:InverseFFT!,
                        FTT<:ForwardFFT!}
    FTFCache::Vector{FT} # storage
      FCache::Vector{F}
       ifft!::ITT        # transforms
        fft!::FTT
end

# Outer constructor:
function LinearisedExTerm(n::Int,
                          m::Int,
                          mode::M,
                          ::Type{S}, flags) where {S<:AbstractFloat, 
                                                   M<:AbstractLinearMode}
    # stupidity check
    m ≥ n || throw(ArgumentError("`m` must be bigger than `n`"))

    # complex fields have size `n` but real fields might have larger size `m`
    FTFCache = [FTField(n, m, S) for i = 1:8]
    FCache   = [  Field(m, S)    for i = 1:8]

    # plan transforms
    ifft! = InverseFFT!(FTFCache[1], flags)
     fft! = ForwardFFT!(  FCache[1], flags)

    # construct object
    LinearisedExTerm{n, m, M, eltype(FTFCache), eltype(FCache),
                                 typeof(ifft!), typeof(fft!)}(FTFCache,
                                                           FCache, ifft!, fft!)
end

# Callable
function (eq::LinearisedExTerm{n, m, M})(t::Real,
                                         Ω::FTField{n, m},
                                         Λ::FTField{n, m},
                                      dΛdt::FTField{n, m},
                                       add::Bool=false) where {n, m}

    # set mean to zero
    Λ[WaveNumber(0, 0)] = 0
    Ω[WaveNumber(0, 0)] = 0

    if M <: TangentMode
            U,   V  = eq.FTFCache[1], eq.FTFCache[2]
           U′,   V′ = eq.FTFCache[3], eq.FTFCache[4]
         dΩdx, dΩdy = eq.FTFCache[5], eq.FTFCache[6]
         dΛdx, dΛdy = eq.FTFCache[7], eq.FTFCache[8]
            u,    v = eq.FCache[1],   eq.FCache[2]
           u′,   v′ = eq.FCache[3],   eq.FCache[4]
         dωdx, dωdy = eq.FCache[5],   eq.FCache[6]
         dλdx, dλdy = eq.FCache[7],   eq.FCache[8]

        # obtain vorticity derivatives
        ddx!(dΛdx,  Λ)
        ddy!(dΛdy,  Λ)
        ddx!(dΩdx,  Ω)
        ddy!(dΩdy,  Ω)

        # obtain velocity components
        invlaplacian!(U,  dΩdy);  U  .*= -1
        invlaplacian!(V,  dΩdx)
        invlaplacian!(U′, dΛdy); U′ .*= -1
        invlaplacian!(V′, dΛdx)

        # inverse transform to physical space into temporaries
        eq.ifft!(u,      U)   ; eq.ifft!(u′,    U′)
        eq.ifft!(v,      V)   ; eq.ifft!(v′,    V′)
        eq.ifft!(dλdx, dΛdx); eq.ifft!(dωdx, dΩdx)
        eq.ifft!(dλdy, dΛdy); eq.ifft!(dωdy, dΩdy)

        # multiply in physical space, overwriting u, then come back
        u .= .-u.*dλdx .- v.*dλdy .- u′.*dωdx .- v′.*dωdy
        eq.fft!(U, u)

        # add or replace
        add == true ? (dΛdt .+= U) : (dΛdt .= U)
    end

    if M <: AdjointMode
                    U,    V    = Eq.FTFStore[1], Eq.FTFStore[2]
                    TMP1, TMP2 = Eq.FTFStore[1], Eq.FTFStore[2]
        dΛdx, dΛdy, dΩdx, dΩdy = Eq.FTFStore[3], Eq.FTFStore[4], Eq.FTFStore[5], Eq.FTFStore[6] 
                    u,    v    = Eq.FStore[1],   Eq.FStore[2]
                    tmp1, tmp2 = Eq.FStore[1],   Eq.FStore[2]
        dλdx, dλdy, dωdx, dωdy = Eq.FStore[3],   Eq.FStore[4],   Eq.FStore[5],   Eq.FStore[6]

        # set mean to zero
        Λ[WaveNumber(0, 0)] = 0

        # obtain vorticity derivatives
        ddx!(dΛdx, Λ)
        ddy!(dΛdy, Λ)
        ddx!(dΩdx, Ω)
        ddy!(dΩdy, Ω)

        # obtain velocity components. Set mean to zero
        invlaplacian!(U,  dΩdy);  U  .*= -1
        invlaplacian!(V,  dΩdx)

        # inverse transform to physical space into temporaries
        Eq.ifft!(u,    U)
        Eq.ifft!(v,    V)
        Eq.ifft!(dλdx, dΛdx); Eq.ifft!(dωdx, dΩdx)
        Eq.ifft!(dλdy, dΛdy); Eq.ifft!(dωdy, dΩdy)

        # multiply in physical space
        tmp1 .= u.*dλdx .+ v.*dλdy;       Eq.fft!(TMP1, tmp1)
        tmp2 .= dλdx.*dωdy .- dλdy.*dωdx; Eq.fft!(TMP2, tmp2)

        # add or replace
        add == true ? (dΛdt .+= .- TMP1 .- invlaplacian!(U, TMP2)) 
                    : (dΛdt  .= .- TMP1 .- invlaplacian!(U, TMP2))
    end

    # reset mean
    dΛdt[WaveNumber(0, 0)] = 0

    return nothing
end


# ~~~ SOLVER OBJECT FOR THE LINEAR EQUATIONS ~~~
struct LinearisedEquation{n, m,
                          IT<:ImplicitTerm,
                          ET<:LinearisedExTerm{n, m},
                          G<:AbstractForcing{n}}
     imTerm::IT
     exTerm::ET
    forcing::G  # forcing for the linearised equations
end

# outer constructor: main entry point
function LinearisedEquation(n::Int,
                            m::Int,
                           Re::Real,
                         mode::AbstractLinearMode,
                        flags::Int=FFT.EXHAUSTIVE,
                      forcing::AbstractForcing=DummyForcing(n),
                             ::Type{S}=Float64) where {S<:AbstractFloat}
    imTerm = ImplicitTerm(Re)
    exTerm = LinearisedExTerm(n, m, mode, S, flags)
    LinearisedEquation{n, m,
                   typeof(imTerm),
                   typeof(exTerm),
                   typeof(forcing)}(imTerm, exTerm, forcing)
end


# /// SPLIT EXPLICIT AND IMPLICIT PARTS ///
function splitexim(eq::LinearisedEquation{n, m}) where {n, m}
    function wrapper(t::Real,
                     Ω::FTField{n, m},
                  dΩdt::FTField{n, m},
                     Λ::FTField{n, m},
                  dΛdt::FTField{n, m},
                   add::Bool=false)
        eq.exTerm(t, Ω, Λ, dΛdt, add)
        eq.forcing(t, Ω, dΩdt, Λ, dΛdt) # note forcing always adds to dΛdt
        return dΛdt
    end

    # for the adjoint equation we need not to pass the dΩdt term
    function wrapper(t::Real,
                     Ω::FTField{n, m},
                     Λ::FTField{n, m},
                 dΛdt::FTField{n, m},
                   add::Bool=false)
             eq.exTerm(t, Ω, Λ, dΛdt, add)
             eq.forcing(t, Ω, Λ, dΛdt) # note forcing always adds to dΛdt
             return dΛdt
             end

    return wrapper, eq.imTerm
end

# /// EVALUATE RIGHT HAND SIDE OF LINEARISED EQUATION ///
function (eq::LinearisedEquation{n, m})(t::Real,
                                        Ω::FTField{n, m},
                                     dΩdt::FTField{n, m},
                                        Λ::FTField{n, m},
                                     dΛdt::FTField{n, m}) where {n, m}
    A_mul_B!(dΛdt, eq.imTerm, Λ)
    eq.exTerm(t, Ω, Λ, dΛdt, true)
    eq.forcing(t, Ω, dΩdt, Λ, dΛdt) # note forcing always adds to dVdt
    return dΛdt
end

# also provide function withouth dΩdt parameter
function (eq::LinearisedEquation{n, m})(t::Real,
                                        Ω::FTField{n, m},
                                        Λ::FTField{n, m},
                                     dΛdt::FTField{n, m}) where {n, m}
    A_mul_B!(dΛdt, eq.imTerm, Λ)
    eq.exTerm(t, Ω, Λ, dΛdt, true)
    eq.forcing(t, Ω, Λ, dΛdt) # note forcing always adds to dVdt
    return dΛdt
end