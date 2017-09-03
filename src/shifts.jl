export Shift, shifted

# Small struct to describe shifts
struct Shift{S<:Real, M<:Integer} 
    s::S # continuous x shift
    m::M # discrete   y shift
end

# construct a shifted view of U
struct Shifted{n, T, Sx, Sy, F<:AbstractFTField{n, T}} <: AbstractFTField{n, T}
         U::F  # original field
    jcache::Sx # cache values of exp(im*j*s) for fast calculations
    kcache::Sy #
end

# outer constructor
function shifted(U::FTField{n, T}, Δ::Shift) where {n, T}
    jc, kc = jcache(n, Δ.s), kcache(n, Δ.m)
    Shifted{n, T, typeof(jc), typeof(kc), typeof(U)}(U, jc, kc)
end

# Construct cache of exp(im*j*s) 
jcache(n::Int, s::Real) = exp.(.-im.*js(n).*s)
kcache(n::Int, m::Int)  = exp.(.-im.*ks(n).*m.*π./4)

# Read only data type 
Base.@propagate_inbounds @inline Base.getindex(S::Shifted{n}, i::Int) where {n} = 
    S[i, ItoKJ(i, n)...]
Base.@propagate_inbounds @inline Base.getindex(S::Shifted{n}, k::Int, j::Int) where {n} = 
    S.U[k, j]*S.jcache[JtoI(j)]*S.kcache[KtoI(k, n)]
Base.@propagate_inbounds @inline Base.getindex(S::Shifted{n}, i::Int, k::Int, j::Int) where {n} = 
    S.U[i]*S.jcache[JtoI(j)]*S.kcache[KtoI(k, n)]

# Apply shift `Δ` to field `U`, in place
function Base.shift!(U::FTField{n}, Δ::Shift) where {n}
    jc, kc = jcache(n, Δ.s), kcache(n, Δ.m)
    for jj = 1:n>>1+1, kk = 1:n
        @inbounds U.data[kk, jj] *= jc[jj]*kc[kk]
    end
    U
end