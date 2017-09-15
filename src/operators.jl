struct DiffOperator{n, T, M<:AbstractMatrix{T}} <: AbstractFTOperator{n, T}
    data::M
end

# read only
Base.@propagate_inbounds @inline Base.getindex(d::DiffOperator, i::Int) = d.data[i]

function DiffOperator(n::Int, i::Symbol, ::Type{T}=Float64) where {T<:Union{AbstractFloat, Integer}}
    i == :x    && return DiffOperator{n, Complex{T}, Matrix{Complex{T}}}(     Complex{T}[im*j for k=      1:n,    j=0:n>>1])
    i == :y    && return DiffOperator{n, Complex{T}, Matrix{Complex{T}}}(vcat(Complex{T}[im*k for k=      0:n>>1, j=0:n>>1],
                                                          Complex{T}[im*k for k=-n>>1+1:-1,   j=0:n>>1]))
    i == :xx   && return DiffOperator{n, T, Matrix{T}}(                       T[-j*j for k=      1:n,    j=0:n>>1])
    i == :yy   && return DiffOperator{n, T, Matrix{T}}(                  vcat(T[-k*k for k=      0:n>>1, j=0:n>>1],
                                                                   T[-k*k for k=-n>>1+1:-1,   j=0:n>>1]))
    i == :xxyy && return DiffOperator{n, T, Matrix{T}}(                       T[-j*j for k=      1:n,    j=0:n>>1] +
                                                              vcat(T[-k*k for k=      0:n>>1, j=0:n>>1],
                                                                   T[-k*k for k=-n>>1+1:-1,   j=0:n>>1]))
    throw(ArgumentError("differentiation not understood, got $i"))
end