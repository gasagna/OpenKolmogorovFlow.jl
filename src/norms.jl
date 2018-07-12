import LinearAlgebra: dot, norm

# Definitions of inner products, norms and distances
export dotdiff

# Inner product between two vorticity fields (Ω₁, Ω₂)
function dot(Ω₁::AbstractFTField{n, m}, Ω₂::AbstractFTField{n, m}) where {n, m}
    s = real(Ω₁[_reindex(0, 0, m)...]*conj(Ω₂[_reindex(0, 0, m)...]))
    @inbounds for j = 1:n, k=-n:n
        kk, jj = _reindex(k, j, m)
        s += 2*real(Ω₁[kk, jj]*conj(Ω₂[kk, jj]))
    end
    @inbounds for k=-n:n
        kk, jj = _reindex(k, 0, m)
        s += real(Ω₁[kk, jj]*conj(Ω₂[kk, jj]))
    end
    return s
end


# Inner product of the difference of two vorticity fields (Ω₁-Ω₂, Ω₁-Ω₂)
function dotdiff(Ω₁::AbstractFTField{n, m}, Ω₂::AbstractFTField{n, m}) where {n, m}
    s = abs2(Ω₁[_reindex(0, 0, m)...]-Ω₂[_reindex(0, 0, m)...])
    @inbounds for j = 1:n, k=-n:n
        kk, jj = _reindex(k, j, m)
        s += 2*abs2(Ω₁[kk, jj]-Ω₂[kk, jj])
    end
    @inbounds for k=-n:n
        kk, jj = _reindex(k, 0, m)
        s += abs2(Ω₁[kk, jj]-Ω₂[kk, jj])
    end
    return s
end

# The two norm of a field ||Ω|| = sqrt(inner(Ω, Ω))
norm(Ω::AbstractFTField, n::Int=2) =
    (n==2 || throw(ArgumentError("only the 2-norm is defined"));
    sqrt(dot(Ω, Ω)))