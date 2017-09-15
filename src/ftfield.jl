export AbstractFTField, FTField, growto!, shrinkto!

# ~~~ ABSTRACT TYPE FOR STRUCTURED, ARRAY-LIKE OBJECTS FROM FFT ~~~
abstract type AbstractFTArray{n, T}    <: AbstractMatrix{T}  end
abstract type AbstractFTField{n, T}    <: AbstractFTArray{n, T} end
abstract type AbstractFTOperator{n, T} <: AbstractFTArray{n, T} end

Base.size(::AbstractFTArray{n})      where {n} = (n, n>>1+1)
Base.eltype(::AbstractFTArray{n, T}) where {n, T} = T
Base.IndexStyle(::Type{<:AbstractFTArray}) = IndexLinear()
fieldsize(::Type{<:AbstractFTArray{n}}) where {n} = n
fieldsize(U::AbstractFTArray{n}) where {n} = n

function Base.checkbounds(U::AbstractFTArray{n}, k::Int, j::Int) where {n}
    d = n>>1
    -d+1 ≤ j ≤ d || throw(BoundsError(U, (k, j)))
    -d+1 ≤ k ≤ d || throw(BoundsError(U, (k, j)))
    return nothing
end


# ~~~ MAIN TYPE ~~~
struct FTField{n, T<:Complex, M<:AbstractMatrix{T}} <: AbstractFTField{n, T}
    data::M
    FTField{n}(data::M) where {n, T, M<:AbstractMatrix{T}} = 
        new{n, T, M}(_checksize(data, n))
end

# extra constructors
FTField(n::Int) = FTField(n, Complex128)
FTField(n::Int, ::Type{T}) where {T} = FTField(zeros(T, n, n>>1+1))
FTField(data::AbstractMatrix) = FTField{size(data, 1)}(data)

# size checker
function _checksize(data::AbstractMatrix{<:Complex}, n::Int)
    M, N = size(data)
    iseven(n)   || throw(ArgumentError("`n` must be even, got $n"))
    M == n      || throw(ArgumentError("wrong row number, got $M but wanted $n"))
    N == n>>1+1 || throw(ArgumentError("wrong column number, got $N but wanted $n"))
    return data
end

# ~~~ array interface ~~~
@inline function Base.getindex(U::FTField{n}, k::Int, j::Int) where {n}
    @boundscheck checkbounds(U, k, j)
    @inbounds ret = U.data[KJtoI(k, j, n)]
    rectify(ret, j)
end

@inline function Base.setindex!(U::FTField{n}, val::Number, k::Int, j::Int) where {n}
    @boundscheck checkbounds(U, k, j)
    @inbounds U.data[KJtoI(k, j, n)] = rectify(val, j)
    val
end

@inline function Base.getindex(U::FTField, i::Int)
    P = U.data
    @boundscheck checkbounds(P, i)
    @inbounds ret = P[i]
    ret
end

@inline function Base.setindex!(U::FTField, val::Number, i::Int)
    P = U.data
    @boundscheck checkbounds(P, i)
    @inbounds P[i] = val
    val
end

# accessors functions
@inline Base.parent(U::FTField) = U.data

# allow constructing similar fields. Used by IMEXRKCB to allocate storage.
Base.similar(U::FTField{n, T}, m::Int=n) where {n, T} = FTField(m, T)

# ~~~ Copy one field to the other, e.g. for zero padding or truncation ~~~

# same size is just a copy
growto!(W::FTField{n}, U::FTField{n}) where {n} = W .= U

# different size requires special handling of boundary terms
function growto!(W::AbstractFTField{m}, U::AbstractFTField{n}) where {m, n}
    m >= n || throw(ArgumentError("output `W` should be larger or equal than input `U`"))
    @inbounds begin
        dU = n>>1
        W .= 0
        for j = 0:dU-1, k = -dU+1:dU
            W[k, j] = U[k, j]
        end
        # move the waves at last column to waves at positive k frequency only
        for k = 0:dU
            W[k, dU] = U[k, dU]
        end
        # make sure we preserve the appropriate weight for extreme frequencies
        W[ dU, 0] *= 0.5
        W[-dU, 0]  = W[ dU, 0]
        W[0,  dU] *= 0.5
        W[dU, dU] *= 0.5
    end
    W
end

# same size is a copy
shrinkto!(W::FTField{n}, U::FTField{n}) where {n} = W .= U

function shrinkto!(W::FTField{n}, U::FTField{m}) where {n, m}
    m >= n || throw(ArgumentError("input `U` should be larger or equal than output `W`"))
    @inbounds begin
        dW = n>>1
        for j = 0:dW, k = (-dW+1):dW
            W[k, j] = U[k, j]
        end
        # make last column conjugate symmetric
        for k = (-dW+1):dW
            val = U[k, dW] + U[-k, dW]
            W[ k, dW] =      0.5*val
            W[-k, dW] = conj(0.5*val)
        end
        # multiply by two extreme frequencies
        W[0,  dW] = 2*real(U[0,  dW])
        W[dW,  0] = 2*real(U[dW,  0])
        W[dW, dW] = 2*real(U[dW, dW])
    end
    W
end