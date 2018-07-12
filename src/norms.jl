import LinearAlgebra: dot, norm

# Definitions of inner products, norms and distances
export dotdiff

# Inner product between two vorticity fields (Ω₁, Ω₂)
function dot(Ω₁::AbstractFTField{n, m, T}, Ω₂::AbstractFTField{n, m, T}) where {n, m, T}
    s = zero(T)
    @inbounds for j = 1:n, k=-n:n
        s += 2*real(Ω₁[k, j]*conj(Ω₂[k, j]))
    end
    @inbounds for k=-n:n
        s += real(Ω₁[k, 0]*conj(Ω₂[k, 0]))
    end
    return s
end

# Inner product of the difference of two vorticity fields (Ω₁-Ω₂, Ω₁-Ω₂)
function dotdiff(Ω₁::AbstractFTField{n, m, T}, Ω₂::AbstractFTField{n, m, T}) where {n, m, T}
    s = zero(T)
    @inbounds for j = 1:n, k=-n:n
        s += 2*abs2(Ω₁[k, j]-Ω₂[k, j])
    end
    @inbounds for k=-n:n
        s += abs2(Ω₁[k, 0]-Ω₂[k, 0])
    end
    return s
end

# The two norm of a field ||Ω|| = sqrt(inner(Ω, Ω))
norm(Ω::AbstractFTField, n::Int=2) =
    (n==2 || throw(ArgumentError("only the 2-norm is defined"));
    sqrt(dot(Ω, Ω)))