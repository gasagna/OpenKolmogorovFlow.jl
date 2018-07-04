function ddx!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m}
    for j = 0:n, k = -n:n
        kk, jj = _reindex(k, j, m)
        @inbounds OUT[kk, jj] = im*j*U[kk, jj]
    end
    return OUT
end

function ddy!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m}
    for j = 0:n, k = -n:n
        kk, jj = _reindex(k, j, m)
        @inbounds OUT[kk, jj] = im*k*U[kk, jj]
    end
    return OUT
end

function invlaplacian!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m}
    for j = 0:n, k = -n:n
        kk, jj = _reindex(k, j, m)
        @inbounds OUT[kk, jj] = -U[kk, jj]/(j^2 + k^2)
    end
    @inbounds OUT[0, 0] = 0
    return OUT
end

function invlaplacian!(OUT::FTField{n, m}, U::FTField{n, m}, c::Real) where {n, m}
    for j = 0:n, k = -n:n
        kk, jj = _reindex(k, j, m)
        @inbounds OUT[kk, jj] = U[kk, jj]/(1 + c*(j^2 + k^2))
    end
    @inbounds OUT[0, 0] = 0
    return OUT
end

function laplacian!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m}
    for j = 0:n, k = -n:n
        kk, jj = _reindex(k, j, m)
        @inbounds OUT[kk, jj] = -U[kk, jj]*(j^2 + k^2)
    end
    return OUT
end