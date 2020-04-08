export ddx!, ddy!, invlaplacian!, laplacian!

function ddx!(OUT::FTField{n,m}, U::FTField{n,m}) where {n, m}
    α = U.α
    @loop_jk n m OUT[_k, _j] = im * j * α * U[_k, _j]
    return OUT
end

function ddy!(OUT::FTField{n,m}, U::FTField{n,m}) where {n, m}
    @loop_jk n m OUT[_k, _j] = im * k * U[_k, _j]
    return OUT
end

function invlaplacian!(OUT::FTField{n,m}, U::FTField{n,m}) where {n, m}
    α = U.α
    @loop_jk n m OUT[_k, _j] = - U[_k, _j] / ((α*j)^2 + k^2)
    @inbounds OUT[WaveNumber(0, 0)] = 0
    return OUT
end

function invlaplacian!(OUT::FTField{n,m}, U::FTField{n,m}, c::Real) where {n, m}
    α = U.α
    @loop_jk n m OUT[_k, _j] = U[_k, _j] / (1 + c * ((α*j)^2 + k^2))
    @inbounds OUT[WaveNumber(0, 0)] = 0
    return OUT
end

function laplacian!(OUT::FTField{n,m}, U::FTField{n,m}) where {n, m}
    α = U.α
    @loop_jk n m OUT[_k, _j] = -U[_k, _j] * ((α*j)^2 + k^2)
    return OUT
end