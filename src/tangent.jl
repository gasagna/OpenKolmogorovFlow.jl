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

# Outer constructor: TODO: merge with above
function TangentExplicitTerm(n::Int, m::Int, ::Type{S}) where {S<:Real}
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
function (Eq::TangentExplicitTerm{n, m, FT})(t::Real, Ω::FT, Ω′::FT, dΩ′dt::FT, add::Bool=false) where {n, m, FT<:FTField{n, m}}
    # extract aliases. We start from one because 1 might be taken for Ω
        U,    V  = Eq.FTFCache[1], Eq.FTFCache[2]
       U′,    V′ = Eq.FTFCache[3], Eq.FTFCache[4]
     dΩdx,  dΩdy = Eq.FTFCache[5], Eq.FTFCache[6]
    dΩ′dx, dΩ′dy = Eq.FTFCache[7], Eq.FTFCache[8]
        u,     v = Eq.FCache[1],   Eq.FCache[2]
       u′,    v′ = Eq.FCache[3],   Eq.FCache[4]
     dωdx,  dωdy = Eq.FCache[5],   Eq.FCache[6]
    dω′dx, dω′dy = Eq.FCache[7],   Eq.FCache[8]

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
    Eq.ifft!(u,      U)   ; Eq.ifft!(u′,    U′)
    Eq.ifft!(v,      V)   ; Eq.ifft!(v′,    V′)
    Eq.ifft!(dω′dx, dΩ′dx); Eq.ifft!(dωdx, dΩdx)
    Eq.ifft!(dω′dy, dΩ′dy); Eq.ifft!(dωdy, dΩdy)

    # multiply in physical space, overwriting u, then come back
    u .= .-u.*dω′dx .- v.*dω′dy .- u′.*dωdx .- v′.*dωdy
    Eq.fft!(U, u)

    add == true ? (dΩ′dt .+= U) : (dΩ′dt .= U)

    return nothing
end


# ~~~ SOLVER OBJECT FOR THE LINEAR EQUATIONS ~~~
struct TangentEquation{n, m, 
                       IT<:ImplicitTerm,
                       ET<:TangentExplicitTerm{n, m},
                       G<:AbstractForcing{n},
                       M<:Flows.AbstractMonitor,
                       FT<:AbstractFTField{n, m}}
     imTerm::IT
     exTerm::ET
    forcing::G  # forcing for the linearised equations
        mon::M  # monitor for the forward solution
        TMP::FT # temporary field for the interpolation
          χ::Float64
end

# outer constructor: main entry point
function TangentEquation(n::Int,
                         m::Int,
                         Re::Real,
                         mon::Flows.AbstractMonitor{T, X},
                         forcing::AbstractForcing=DummyForcing(n),
                         χ::Real=0,
                         ::Type{S}=Float64) where {T, X, S<:Real}
    X <: FTField{n, m} || error("invalid monitor object")
    imTerm = ImplicitTerm(Re)
    exTerm = TangentExplicitTerm(n, m, S)
    TMP = FTField(n, m, S)
    TangentEquation{n, m,
                   typeof(imTerm),
                   typeof(exTerm),
                   typeof(forcing),
                   typeof(mon),
                   typeof(TMP)}(imTerm, exTerm, forcing, mon, TMP, χ)
end

# obtain two components
function splitexim(eq::TangentEquation{n, m}) where {n, m}
    function wrapper(t::Real, Ω′::AbstractFTField{n, m}, dΩ′dt::AbstractFTField{n, m})
        # interpolate U and evaluate nonlinear interaction term and forcing
        eq.mon(eq.TMP, t, Val{0}())
        eq.exTerm(t, eq.TMP, Ω′, dΩ′dt, false)
        eq.forcing(t, eq.TMP, Ω′, dΩ′dt)

        # interpolate dUdt if needed
        if eq.χ != 0
            eq.mon(eq.TMP, t, Val{1}())
            dΩ′dt .+= eq.χ .* eq.TMP
        end

        return dΩ′dt
    end
    return wrapper, eq.imTerm
end



