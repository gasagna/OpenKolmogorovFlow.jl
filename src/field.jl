export Field, fieldsize, make_grid

struct Field{n, T<:Real, M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    data::M
    function Field{n, T, M}(data::M) where {n, T, M}
        checksize(data, n)
        new{n, T, M}(data)
    end
end

Field(data::AbstractMatrix) = 
    Field{size(data, 1), eltype(data), typeof(data)}(data)
Field(n::Int, ::Type{T}=Float64) where {T} = Field(zeros(T, n, n))

function checksize(data::AbstractMatrix{<:Real}, n::Int)
    M, N = size(data)
    N == M    || throw(ArgumentError("input matrix must be square: got $M×$N"))
    iseven(n) || throw(ArgumentError("`n` must be even, got $n"))
    M == n    || throw(ArgumentError("wrong row number, got $M, should be $n"))
    N == n    || throw(ArgumentError("wrong column number, got $N, should be $n"))
end

fieldsize(::Type{Field{n}}) where {n} = n
fieldsize(::Field{n})       where {n} = n

# ~~~ array interface ~~~

@inline function Base.getindex(f::Field, i::Int, j::Int)
    P = f.data
    @boundscheck checkbounds(f, i, j)
    @inbounds ret = P[i+1, j+1]
    return ret
end

@inline function Base.setindex!(f::Field, val::Number, i::Int, j::Int)
    P = f.data
    @boundscheck checkbounds(f, i, j)
    @inbounds P[i+1, j+1] = val
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

Base.indices(f::Field{n})       where n = (0:n-1, 0:n-1)
Base.linearindices(f::Field{n}) where n =  1:n^2
Base.IndexStyle(::Type{Field}) = Base.IndexLinear()
Base.similar(u::Field) = Field(similar(u.data))

# ~~~ GRID FUNCTIONALITY ~~~
function make_grid(n)
    Δ = 2π/n
    x = reshape(collect(0:Δ:2π-Δ), 1, n)
    y = reshape(collect(0:Δ:2π-Δ), n, 1)
    return x, y
end