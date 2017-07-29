# TODO
# ~ check that setindex with k and j does not invalidate the fft

export FTField, growto!, shrinkto!

struct FTField{n, T<:Complex, M<:AbstractMatrix{T}}
    data::M
    function FTField{n, T, M}(data::M) where {n, T, M}
        checksize(data, n)
        new(data)
    end
end

FTField(n::Int) = FTField(n, Complex{Float64})
FTField(n::Int, ::Type{T}) where T = FTField(zeros(T, n, n>>1+1))
FTField(data::AbstractMatrix) = FTField{size(data, 1), eltype(data), typeof(data)}(data)

function checksize(data::AbstractMatrix, n::Int)
    M, N = size(data)
    iseven(n)   || throw(ArgumentError("`n` must be even, got $n"))
    M == n      || throw(ArgumentError("wrong row number, got $M"))
    N == n>>1+1 || throw(ArgumentError("wrong column number, got $N"))
end

# ~~~ array interface ~~~
@inline function Base.getindex(U::FTField{n}, k::Int, j::Int) where n
    P = U.data
    @boundscheck checkbounds(U, k, j)
    @inbounds ret = P[KJtoI(k, j, n)]
    rectify(ret, j)
end

@inline function Base.setindex!(U::FTField{n}, val::Number, k::Int, j::Int) where n
    P = U.data
    @boundscheck checkbounds(U, k, j)
    @inbounds P[KJtoI(k, j, n)] = rectify(val, j)
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

# `indices` is used for bounds checking
Base.indices(::FTField{n}) where n = (-n>>1:n>>1, -n>>1:n>>1)
Base.linearindices(U::FTField) = eachindex(U.data)
Base.IndexStyle(::Type{<:FTField}) = IndexLinear()
Base.unsafe_get(U::FTField) = U.data

# allow constructing similar fields. Used by IMEXRKCB to allocate storage.
Base.similar(U::FTField{n}, T::Type, shape::Tuple{Range, Range}) where n =
    FTField(similar(U.data))

# ~~~ Copy one field to the other, e.g. for zero padding or truncation ~~~
function growto!(W::FTField{m}, U::FTField{n}) where {m, n}
    m >= n || throw(ArgumentError("output `W` should be larger than input `U`"))
    @inbounds begin
        dU = n>>1
        W .= 0
        for j = 0:dU, k = -dU:dU
            W[k, j] = U[k, j]
        end
    end
    W
end

function shrinkto!(W::FTField{n}, U::FTField{m}) where {n, m}
    m >= n || throw(ArgumentError("input `U` should be larger than output `W`"))
    @inbounds begin
        dW = n>>1
        for j = 0:dW, k = -dW:dW
            W[k, j] = U[k, j]
        end
        # enforce symmetries in fft
        W[dW,  0] = real(W[dW,  0])
        W[0,  dW] = real(W[0,  dW])
        W[dW, dW] = real(W[dW, dW])
    end
    W
end