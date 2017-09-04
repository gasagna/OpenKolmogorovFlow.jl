export AbstractFTField, FTField, growto!, shrinkto!, fieldsize

# AbstractType for Fourier Transformed fields
abstract type AbstractFTField{n, T} <: AbstractMatrix{T} end
Base.indices(::AbstractFTField{n}) where {n} = (d = n>>1; (-d:d, -d:d))
Base.linearindices(U::AbstractFTField{n}) where {n} = 1:(n*(n>>1+1))
Base.IndexStyle(::Type{<:AbstractFTField}) = IndexLinear()

# Main type
struct FTField{n, T<:Complex, M<:AbstractMatrix{T}} <: AbstractFTField{n, T}
    data::M
    function FTField{n, T, M}(data::M) where {n, T, M}
        FTField_checksize(data, n)
        new(data)
    end
end

fieldsize(::Type{FTField{n}}) where {n} = n
fieldsize(::FTField{n})       where {n} = n

FTField(n::Int) = FTField(n, Complex{Float64})
FTField(n::Int, T::Type) = FTField(zeros(T, n, n>>1+1))
FTField(data::AbstractMatrix) = FTField{size(data, 1), eltype(data), typeof(data)}(data)

function FTField_checksize(data::AbstractMatrix, n::Int)
    M, N = size(data)
    iseven(n)   || throw(ArgumentError("`n` must be even, got $n"))
    M == n      || throw(ArgumentError("wrong row number, got $M"))
    N == n>>1+1 || throw(ArgumentError("wrong column number, got $N"))
end

# ~~~ array interface ~~~
@inline function Base.getindex(U::FTField{n}, k::Int, j::Int) where n
    @boundscheck checkbounds(U, k, j)
    @inbounds ret = U.data[KJtoI(k, j, n)]
    rectify(ret, j)
end

@inline function Base.setindex!(U::FTField{n}, val::Number, k::Int, j::Int) where n
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

Base.unsafe_get(U::FTField) = U.data

# allow constructing similar fields. Used by IMEXRKCB to allocate storage.
# TODO. allow, different, type and shape
Base.similar(U::FTField) = FTField(zeros(U.data))

# get/set perturbation for variational analysis. There is no
# need for checks of FTField have same `n`
VariationalNumbers.set_pert!(A::FTField{n, Complex{VarNum{T}}}, 
                             p::FTField{n, Complex{T}}) where {n, T} = 
    (unsafe_set_pert!(A.data, p.data))

VariationalNumbers.get_pert!(A::FTField{n, Complex{VarNum{T}}}, 
                             p::FTField{n, Complex{T}}) where {n, T} = 
    (unsafe_get_pert!(A.data, p.data))

# ~~~ Copy one field to the other, e.g. for zero padding or truncation ~~~

# same size is just a copy
growto!(W::FTField{n}, U::FTField{n}) where {n} = W .= U

# different size requires special handling of boundary terms
function growto!(W::FTField{m}, U::FTField{n}) where {m, n}
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
        # halve energy of the last modes
        for k = (-dW+1):dW
            W[k, dW] *= 0.5
        end
    end
    W
end