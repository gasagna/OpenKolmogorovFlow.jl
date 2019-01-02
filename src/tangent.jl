export TangentEquation,
       splitexim

# ~~~ Explicit Term of Linearised Equations ~~~
struct TangentExplicitTerm{n, m, FT<:AbstractFTField, F<:AbstractField,
                                ITT<:InverseFFT!,   FTT<:ForwardFFT!}
    FTFCache::Vector{FT} # storage
      FCache::Vector{F}
       ifft!::ITT        # transforms
        fft!::FTT
end

# Outer constructor:
function TangentExplicitTerm(n::Int, m::Int, ::Type{S}) where {S<:AbstractFloat}
    # stupidity check
    m ≥ n || throw(ArgumentError("`m` must be bigger than `n`"))

    # complex fields have size `n` but real fields might have larger size `m`
    FTFCache = [FTField(n, m, S) for i = 1:8]
    FCache   = [  Field(m, S)    for i = 1:8]

    # transforms
    ifft! = InverseFFT!(FTFCache[1])
     fft! = ForwardFFT!(  FCache[1])

    # construct object
    TangentExplicitTerm{n, m, eltype(FTFCache), eltype(FCache),
                              typeof(ifft!), typeof(fft!)}(FTFCache,
                                                           FCache, ifft!, fft!)
end

# Callable
function (eq::TangentExplicitTerm{n, m})(t::Real,
                                         Ω::FTField{n, m},
                                  NOT_USED::FTField{n, m},
                                        Ω′::FTField{n, m},
                                     dΩ′dt::FTField{n, m},
                                       add::Bool=false) where {n, m}
    # extract aliases. We start from one because 1 might be taken for Ω
        U,    V  = eq.FTFCache[1], eq.FTFCache[2]
       U′,    V′ = eq.FTFCache[3], eq.FTFCache[4]
     dΩdx,  dΩdy = eq.FTFCache[5], eq.FTFCache[6]
    dΩ′dx, dΩ′dy = eq.FTFCache[7], eq.FTFCache[8]
        u,     v = eq.FCache[1],   eq.FCache[2]
       u′,    v′ = eq.FCache[3],   eq.FCache[4]
     dωdx,  dωdy = eq.FCache[5],   eq.FCache[6]
    dω′dx, dω′dy = eq.FCache[7],   eq.FCache[8]

    # obtain vorticity derivatives
    ddx!(dΩ′dx, Ω′)
    ddy!(dΩ′dy, Ω′)
    ddx!(dΩdx,  Ω)
    ddy!(dΩdy,  Ω)

    # obtain velocity components
    invlaplacian!(U,  dΩdy);  U  .*= -1
    invlaplacian!(V,  dΩdx)
    invlaplacian!(U′, dΩ′dy); U′ .*= -1
    invlaplacian!(V′, dΩ′dx)

    # inverse transform to physical space into temporaries
    eq.ifft!(u,      U)   ; eq.ifft!(u′,    U′)
    eq.ifft!(v,      V)   ; eq.ifft!(v′,    V′)
    eq.ifft!(dω′dx, dΩ′dx); eq.ifft!(dωdx, dΩdx)
    eq.ifft!(dω′dy, dΩ′dy); eq.ifft!(dωdy, dΩdy)

    # multiply in physical space, overwriting u, then come back
    u .= .-u.*dω′dx .- v.*dω′dy .- u′.*dωdx .- v′.*dωdy
    eq.fft!(U, u)

    add == true ? (dΩ′dt .+= U) : (dΩ′dt .= U)

    return nothing
end


# ~~~ SOLVER OBJECT FOR THE LINEAR EQUATIONS ~~~
struct TangentEquation{n, m,
                       IT<:ImplicitTerm,
                       ET<:TangentExplicitTerm{n, m},
                       G<:AbstractForcing{n}}
     imTerm::IT
     exTerm::ET
    forcing::G  # forcing for the linearised equations
end

# outer constructor: main entry point
function TangentEquation(n::Int,
                         m::Int,
                         Re::Real,
                         forcing::AbstractForcing=DummyForcing(n),
                         ::Type{S}=Float64) where {S<:AbstractFloat}
    imTerm = ImplicitTerm(Re)
    exTerm = TangentExplicitTerm(n, m, S)
    TangentEquation{n, m,
                   typeof(imTerm),
                   typeof(exTerm),
                   typeof(forcing)}(imTerm, exTerm, forcing)
end


# /// SPLIT EXPLICIT AND IMPLICIT PARTS ///
function splitexim(eq::TangentEquation{n, m}) where {n, m}
    function wrapper(t::Real,
                     Ω::FTField{n, m},
                  dΩdt::FTField{n, m},
                    Ω′::FTField{n, m},
                 dΩ′dt::FTField{n, m},
                   add::Bool=false)
        eq.exTerm(t, Ω, dΩdt, Ω′, dΩ′dt, add)
        eq.forcing(t, Ω, dΩdt, Ω′, dΩ′dt) # note forcing always adds to dVdt
        return dΩ′dt
    end
    return wrapper, eq.imTerm
end

# /// EVALUATE RIGHT HAND SIDE OF LINEARISED EQUATION ///
function (eq::TangentEquation{n, m})(t::Real,
                                     Ω::FTField{n, m},
                                  dΩdt::FTField{n, m},
                                    Ω′::FTField{n, m},
                                 dΩ′dt::FTField{n, m}) where {n, m}
    A_mul_B!(dΩ′dt, eq.imTerm, Ω′)
    eq.exTerm(t, Ω, dΩdt, Ω′, dΩ′dt, true)
    eq.forcing(t, Ω, dΩdt, Ω′, dΩ′dt) # note forcing always adds to dVdt
    return dΩ′dt
end