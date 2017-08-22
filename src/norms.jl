# Definitions of dot products norms and distances

function Base.dot(Ω₁::FTField{n, Complex{T}}, Ω₂::FTField{n, Complex{T}})::T where {n, T}
    @inbounds begin
        # initialise, skipping the (0, 0) mode, assumed zero
        Σ1, Σ2 = zero(Complex{T}), zero(Complex{T})
        # loop
        d = n>>1
        # count j = 0 with weight 1
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