export Field, make_grid

abstract type AbstractField{m, T} <: AbstractMatrix{T} end

Base.size(f::AbstractField{m}) where {m} = (2m+2, 2m+2)
Base.IndexStyle(::Type{<:AbstractField}) = Base.IndexLinear()

struct Field{m, T<:Real, M<:AbstractMatrix{T}} <: AbstractField{m, T}
    data::M
       α::Float64
    function Field(data::M, α::Real=1) where {T<:Real, M<:AbstractMatrix{T}}
        _checksize(data)
        new{size(data, 1)>>1 - 1, T, M}(data, Float64(α))
    end
end


# OUTER CONSTRUCTORS
Field(m::Int, ::Type{T}=Float64, α::Real=1) where {T} = Field(zeros(T, 2m+2, 2m+2))

# provide function
Field(m::Int, fun::Base.Callable, α::Real=1) = Field(fun.(make_grid(m)...))

function _checksize(data::AbstractMatrix{<:Real})
    M, N = size(data)
    N == M    || throw(ArgumentError("input matrix must be square: got $M×$N"))
    iseven(M) || throw(ArgumentError("size must be even"))
    return nothing
end

# ~~~ array interface ~~~
@inline function Base.getindex(f::Field{n}, i::Int, j::Int) where {n}
    @boundscheck checkbounds(f, i, j)
    @inbounds ret = f.data[i, j]
    return ret
end

@inline function Base.setindex!(f::Field{n}, val::Number, i::Int, j::Int) where {n}
    @boundscheck checkbounds(f, i, j)
    @inbounds f.data[i, j] = val
    return val
end

# Linear indexing
@inline function Base.getindex(f::Field, i::Int)
    @boundscheck checkbounds(f, i)
    @inbounds ret = f.data[i]
    return ret
end

@inline function Base.setindex!(f::Field, val::Number, i::Int)
    @boundscheck checkbounds(f, i)
    @inbounds f.data[i] = val
    return val
end

# accessors functions
Base.parent(U::Field) = U.data

Base.similar(u::Field{m, T}) where {m, T} = Field(m, T)
Base.copy(u::Field{m, T}) where {m, T} = (v = similar(u); v .= u; v)

# ~~~ GRID FUNCTIONALITY ~~~
function make_grid(m::Int, α::Real=1)
    x = range(0, stop=2π,   length=2m+3)[1:(2m+2)]
    y = range(0, stop=2π/α, length=2m+3)[1:(2m+2)]
    return reshape(x, 1, 2m+2), reshape(y, 2m+2, 1)
end

make_grid(u::Field{m}) where {m} = make_grid(m, u.α)