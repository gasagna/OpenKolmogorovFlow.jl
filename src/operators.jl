export DiffOperator

struct DiffOperator{n, T<:Number}
    data::Matrix{T}
end

function DiffOperator(n::Int, i::Symbol)
    i == :x  && return DiffOperator{n, Complex{Int}}(     [im*j for k=      0:0,    j=0:n>>1])
    i == :y  && return DiffOperator{n, Complex{Int}}(vcat([im*k for k=      0:n>>1, j=0:0], 
                                                          [im*k for k=-n>>1+1:-1,   j=0:0]))
    i == :xx && return DiffOperator{n, Int}(              [-j*j for k=0:0,          j=0:n>>1])
    i == :yy && return DiffOperator{n, Int}(         vcat([-k*k for k=      0:n>>1, j=0:0], 
                                                          [-k*k for k=-n>>1+1:-1,   j=0:0]))
    throw(ArgumentError("differentiation not understood, got $i"))
end