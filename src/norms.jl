# Definitions of inner products, norms and distances
export dotdiff

# Inner product between two vorticity fields (Ω₁, Ω₂). When the fields are 
# augmented by perturbations this is what gets computed
# 
#     dot(Ω₁+η₁, Ω₂+η₂) = dot(Ω₁, Ω₂) + ε*[dot(Ω₁, η₂) + dot(η₁, Ω₂)]
#
#
Base.dot(Ω₁::AbstractFTField{n}, Ω₂::AbstractFTField{n}) where {n} =
    4π^2 * @Σ_jk n real(Ω₁[jk]*conj(Ω₂[jk]))

# Inner product of the difference of two vorticity fields (Ω₁-Ω₂, Ω₁-Ω₂)
dotdiff(Ω₁::AbstractFTField{n}, Ω₂::AbstractFTField{n}) where {n} =
    4π^2 * @Σ_jk n abs2(Ω₁[jk]-Ω₂[jk])

# The two norm of a field ||Ω|| = sqrt(inner(Ω, Ω))
Base.norm(Ω::AbstractFTField, n::Int=2) =
    (n==2 || throw(ArgumentError("only the 2-norm is defined"));
    sqrt(dot(Ω, Ω)))