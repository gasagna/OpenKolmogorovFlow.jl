export ddx!, ddy!, invlaplacian!, laplacian!

function ddx!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m}
    @inbounds for j = 0:n, k = -n:n
        OUT[k, j] = im*j*U[k, j]
    end
    return OUT
end

function ddy!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m}
    @inbounds for j = 0:n, k = -n:n
        OUT[k, j] = im*k*U[k, j]
    end
    return OUT
end

function invlaplacian!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m}
    @inbounds for j = 0:n, k = -n:n
        OUT[k, j] = -U[k, j]/(j^2 + k^2)
    end
    @inbounds OUT[0, 0] = 0
    return OUT
end

function invlaplacian!(OUT::FTField{n, m}, U::FTField{n, m}, c::Real) where {n, m}
    @inbounds for j = 0:n, k = -n:n
        OUT[k, j] = U[k, j]/(1 + c*(j^2 + k^2))
    end
    @inbounds OUT[0, 0] = 0
    return OUT
end

function laplacian!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m}
    @inbounds for j = 0:n, k = -n:n
        OUT[k, j] = -U[k, j]*(j^2 + k^2)
    end
    return OUT
end