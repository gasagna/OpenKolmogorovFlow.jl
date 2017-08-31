export Shift, shifted

# Small struct to describe shifts
struct Shift{S<:Real, M<:Integer} 
    s::S # continuous x shift
    m::M # discrete   y shift
end

# construct a shifted view of U
struct Shifted{n, T, Sx, Sy, F<:AbstractFTField{n, T}} <: AbstractFTField{n, T}
         U::F  # original field
    xcache::Sx #
    ycache::Sy #
end

# outer constructor
function shifted(U::FTField{n, T}, Δ::Shift) where {n, T}
    cx, cy = _xcache(Δ.s, n), _ycache(Δ.m, n)
    Shifted{n, T, typeof(cx), typeof(cy), typeof(U)}(U, cx, cy)
end

# build cache with complex exponentials for fast shifts
_xcache(s, n) = ntuple(jj->cis((jj-1)*s/2π), n>>1+1)
_ycache(m, n) = ntuple(jj->cis((jj-1)*m/2π), n)

# Read only data type 
Base.@propagate_inbounds @inline Base.getindex(S::Shifted{n}, i) where {n} = 
    S[i, ItoKJ(i, n)...]
Base.@propagate_inbounds @inline Base.getindex(S::Shifted, k, j) = 
    S.U[k, j]*S.xcache[j+1]*S.ycache[k+1]
Base.@propagate_inbounds @inline Base.getindex(S::Shifted, i, k, j) = 
    S.U[i]*S.xcache[j+1]*S.ycache[k+1]