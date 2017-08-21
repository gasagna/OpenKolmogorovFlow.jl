export DiffOperator

struct DiffOperator{n, T}
    data::Matrix{T}
end

function DiffOperator(n::Int, i::Symbol, ::Type{T}) where {T<:Union{AbstractFloat, Integer}} 
    i == :x    && return DiffOperator{n, Complex{T}}(     Complex{T}[im*j for k=      1:n,    j=0:n>>1])
    i == :y    && return DiffOperator{n, Complex{T}}(vcat(Complex{T}[im*k for k=      0:n>>1, j=0:n>>1], 
                                                          Complex{T}[im*k for k=-n>>1+1:-1,   j=0:n>>1]))
    i == :xx   && return DiffOperator{n, T}(                       T[-j*j for k=      1:n,    j=0:n>>1])
    i == :yy   && return DiffOperator{n, T}(                  vcat(T[-k*k for k=      0:n>>1, j=0:n>>1], 
                                                                   T[-k*k for k=-n>>1+1:-1,   j=0:n>>1]))
    i == :xxyy && return DiffOperator{n, T}(                       T[-j*j for k=      1:n,    j=0:n>>1] +
                                                              vcat(T[-k*k for k=      0:n>>1, j=0:n>>1], 
                                                                   T[-k*k for k=-n>>1+1:-1,   j=0:n>>1]))
    throw(ArgumentError("differentiation not understood, got $i"))
end