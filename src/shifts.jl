export shift!, xshift!, yshift!

function xshift!(U::FTField{n}, s::Real) where {n}
    @inbounds for j = 0:n
        # precompute this, since it's expensive
        val = cis(j*s)
        @simd for k = -n:n
            U[WaveNumber(k, j)] *= val
        end
    end
    return U
end

function yshift!(U::FTField{n}, m::Int) where {n}
    @inbounds for k = -n:n
        # precompute this, since it's expensive. Note also this
        # assumes that kf = 4
        val = cis(k*m*π/2)
        @simd for j = 0:n
            U[WaveNumber(k, j)] *= val
        end
    end
    return U
end

shift!(U::FTField, s::Real, m::Int) = yshift!(xshift!(U, s), m)