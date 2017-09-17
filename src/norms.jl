# Definitions of inner products, norms and distances
export inner, innerdiff

# Inner product between two vorticity fields (Ω₁, Ω₂). When the fields are 
# augmented by perturbations this is what gets computed
# 
#     inner(Ω₁+η₁, Ω₂+η₂) = inner(Ω₁, Ω₂) + ε*[inner(Ω₁, η₂) + inner(η₁, Ω₂)]
#
#
inner(Ω₁::AbstractFTField{n}, Ω₂::AbstractFTField{n}) where {n} =
    4π^2 * @Σ_jk n real(Ω₁[jk]*conj(Ω₂[jk]))

# Inner product of the difference of two vorticity fields (Ω₁-Ω₂, Ω₁-Ω₂)
innerdiff(Ω₁::AbstractFTField{n}, Ω₂::AbstractFTField{n}) where {n} =
    4π^2 * @Σ_jk n abs2(Ω₁[jk]-Ω₂[jk])

# The two norm of a field ||Ω|| = sqrt(inner(Ω, Ω))
Base.norm(Ω::AbstractFTField, n::Int=2) =
    (n==2 || throw(ArgumentError("only the 2-norm is defined"));
    sqrt(inner(Ω, Ω)))