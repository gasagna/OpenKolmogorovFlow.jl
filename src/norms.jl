# Definitions of dot products norms and distances

export inner, innerdiff

# This enables \cdot notation
Base.dot(Ω₁::FTField{n}, Ω₂::FTField{n}) where {n} =
    inner(Ω₁, Ω₂)

# Inner product between two vorticity fields
function inner(Ω₁::FTField{n, Complex{T}}, Ω₂::FTField{n, Complex{T}})::T where {n, T}
    @inbounds begin
        # initialise
        Σ1, Σ2 = zero(Complex{T}), zero(Complex{T})
        # loop
        d = n>>1
        # count j = 0 with weight 1, skip the (0, 0) mode, assumed zero
        @simd for k = 2:n
            Σ1 += Ω₁[k]*conj(Ω₂[k])
        end
        # count all other with weight 2
        @simd for k = (n+1):(n*d+n)
            Σ2 += Ω₁[k]*conj(Ω₂[k])
        end
        Σ = 2*Σ2 + Σ1
        # remove some elements counted with wrong weight
        Σ -= Ω₁[d, 0]*conj(Ω₂[d, 0])*0.5
        Σ -= Ω₁[0, d]*conj(Ω₂[0, d])*1.5
        Σ -= Ω₁[d, d]*conj(Ω₂[d, d])*1.5
    end
    # remember ∫∫eⁱ⁰dxdy = 4π^2
    real(Σ*4π^2)
end

# Norm of a field
Base.norm(Ω::FTField, n::Int=2) = 
    (n==2 || throw(ArgumentError("only the 2-norm is defined"));
    sqrt(inner(Ω, Ω)))

# Inner product of the difference (u-v, u-v)
function innerdiff(Ω₁::FTField{n, Complex{T}}, Ω₂::FTField{n, Complex{T}})::T where {n, T}
    @inbounds begin
        Σ1, Σ2 = zero(T), zero(T)
        d = n>>1
        @simd for k = 2:n
            Σ1 += abs2(Ω₁[k]-Ω₂[k])
        end
        @simd for k = (n+1):(n*d+n)
            Σ2 += abs2(Ω₁[k]-Ω₂[k])
        end
        Σ = 2*Σ2 + Σ1
        Σ -= abs2(Ω₁[d, 0] -Ω₂[d, 0])*0.5
        Σ -= abs2(Ω₁[0, d] -Ω₂[0, d])*1.5
        Σ -= abs2(Ω₁[d, d] -Ω₂[d, d])*1.5
    end
    Σ*4π^2
end