export Field, fieldsize, make_grid

abstract type AbstractField{n, T} <: AbstractMatrix{T} end
Base.indices(f::AbstractField{n}) where {n} = (0:n-1, 0:n-1)
Base.linearindices(f::AbstractField{n}) where {n} =  1:n^2
Base.IndexStyle(::Type{<:AbstractField}) = Base.IndexLinear()
fieldsize(::Type{<:AbstractField{n}}) where {n} = n
fieldsize(u::AbstractField{n}) where {n} = n

struct Field{n, T<:Real, M<:AbstractMatrix{T}} <: AbstractField{n, T}
    data::M
    function Field{n, T, M}(data::M) where {n, T, M}
        Field_checksize(data, n)
        new{n, T, M}(data)
    end
end

Field(data::AbstractMatrix) =
    Field{size(data, 1), eltype(data), typeof(data)}(data)
Field(n::Int, ::Type{T}=Float64) where {T} = Field(zeros(T, n, n))

function Field_checksize(data::AbstractMatrix, n::Int)
    M, N = size(data)
    N == M    || throw(ArgumentError("input matrix must be square: got $M×$N"))
    iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
    M == n    || throw(ArgumentError("wrong row number, got $M, should be $n"))
    N == n    || throw(ArgumentError("wrong column number, got $N, should be $n"))
end

# accessors functions
@inline Base.parent(U::Field) = U.data

# ~~~ array interface ~~~
@inline function reindex(n, i, j)
    ii, jj = i%n + 1, j%n + 1
    ii ≤ 0 && (ii += n)
    jj ≤ 0 && (jj += n)
    ii, jj
end

@inline function Base.getindex(f::Field{n}, i::Int, j::Int) where {n}
    P = f.data
    ii, jj = reindex(n, i, j)
    @inbounds ret = P[ii, jj]
    return ret
end

@inline function Base.setindex!(f::Field{n}, val::Number, i::Int, j::Int) where {n}
    P = f.data
    ii, jj = reindex(n, i, j)
    @inbounds P[ii, jj] = val
    return val
end

@inline function Base.getindex(f::Field, i::Int)
    P = f.data
    @boundscheck checkbounds(P, i)
    @inbounds ret = P[i]
    return ret
end

@inline function Base.setindex!(f::Field, val::Number, i::Int)
    P = f.data
    @boundscheck checkbounds(P, i)
    @inbounds P[i] = val
    return val
end

Base.similar(u::Field{n, T}, m::Int=n) where {n, T} = Field(m, T)

# ~~~ GRID FUNCTIONALITY ~~~
function make_grid(n)
    Δ = 2π/n
    x = reshape(collect(0:Δ:2π-Δ), 1, n)
    y = reshape(collect(0:Δ:2π-Δ), n, 1)
    return x, y
end