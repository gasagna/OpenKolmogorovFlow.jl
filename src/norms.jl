import LinearAlgebra: dot, norm

# Definitions of inner products, norms and distances
export dotdiff

# Inner product between two vorticity fields (Ω₁, Ω₂) with possible non-zero mean
function dot(Ω₁::AbstractFTField{n, m, T}, Ω₂::AbstractFTField{n, m, T}) where {n, m, T}
    # we split the loop because different part carry different weight
    s = zero(T)
    @inbounds for _j = 1:n
        @loop_k n m (s += 2 * real(Ω₁[_k, _j] * conj(Ω₂[_k, _j])))
    end
    @loop_k n m (s += real(Ω₁[_k, 0] * conj(Ω₂[_k, 0])))
    # do not count the mean flow twice
    s -= 0.5*real(Ω₁[0, 0] * conj(Ω₂[0, 0]))
    return s
end

# Inner product of the difference of two vorticity fields (Ω₁-Ω₂, Ω₁-Ω₂)
function dotdiff(Ω₁::AbstractFTField{n, m, T}, Ω₂::AbstractFTField{n, m, T}) where {n, m, T}
    s = zero(T)
    @inbounds for _j = 1:n
        @loop_k n m (s += 2*abs2(Ω₁[_k, _j] - Ω₂[_k, _j]))
    end
    @loop_k n m (s += abs2(Ω₁[_k, 0] - Ω₂[_k, 0]))
    s -= 0.5*abs2(Ω₁[0, 0] - Ω₂[0, 0])
    return s
end

# The two norm of a field ||Ω|| = sqrt(inner(Ω, Ω))
norm(Ω::AbstractFTField, n::Int=2) =
    (n==2 || throw(ArgumentError("only the 2-norm is defined"));
    sqrt(dot(Ω, Ω)))