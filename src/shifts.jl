export Shift, shifted

# Small struct to describe the x and y shifts of a field
const Shift = Tuple{T, Int} where {T<:Real}

# construct a shifted view of U
struct Shifted{n, T, F<:AbstractFTField{n, T}} <: AbstractFTField{n, T}
         U::F         # original field
    jcache::Vector{T} # cache values of exp(im*j*s) for fast calculations
    kcache::Vector{T} #
end

# outer constructor
function shifted(U::FTField{n, T}, Δ::Shift) where {n, T}
    jc, kc = jcache(n, Δ[1]), kcache(n, Δ[2])
    Shifted{n, T, typeof(U)}(U, jc, kc)
end

# Construct cache of exp(im*j*s) 
jcache(n::Int, s::Real) = [cis(-j*s) for j = 0:n>>1]
kcache(n::Int,  m::Int) = [cis(-ItoK(i, n)*m*π/4) for i = 1:n]

# Read only data type 
Base.@propagate_inbounds @inline Base.getindex(S::Shifted{n, T}, i::Int) where {n, T} = 
    S[i, ItoKJ(i, n)...]
Base.@propagate_inbounds @inline Base.getindex(S::Shifted{n, T}, k::Int, j::Int) where {n, T} = 
    S.U[k, j]*S.jcache[JtoI(j)]*S.kcache[KtoI(k, n)]
Base.@propagate_inbounds @inline Base.getindex(S::Shifted{n, T}, i::Int, k::Int, j::Int) where {n, T} = 
    S.U[i]*S.jcache[JtoI(j)]*S.kcache[KtoI(k, n)]

# Apply shift `Δ` to field `U`, in place
function Base.shift!(U::FTField{n}, Δ::Shift) where {n}
    jc, kc = jcache(n, Δ[1]), kcache(n, Δ[2])
    for jj = 1:n>>1+1, kk = 1:n
        @inbounds U.data[kk, jj] *= jc[jj]*kc[kk]
    end
    U
end