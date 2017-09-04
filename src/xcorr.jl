export XCorrCache, distance!

struct XCorrCache{m, T, rT<:Field{m, T}, RT<:FTField{m, Complex{T}}, IP}
    r::rT       # the cross correlation field
    R::RT       # the cross correlation transform
    irfft!::IP  # inverse Fourier transform object
    function XCorrCache(m::Int, ::Type{T}=Float64, flags::UInt32=FFTW.PATIENT) where {T}
        # we want to index exact location for vertical integer shifts
        rem(m, 4) == 0 || error("cache size / 4 must be an even number")
        r, R = Field(m, T), FTField(m, Complex{T})
        irfft! = InverseFFT!(typeof(r), R, flags)
        new{m, T, typeof(r), typeof(R), typeof(irfft!)}(r, R, irfft!)
    end
end

function distance!(U::FTField{n}, V::FTField{n}, cache::XCorrCache{nc, T}) where {n, nc, T}
    evalconjprod!(U, V, cache)              # evaluate product
    cache.irfft!(cache.r, cache.R)          # inverse transform
    rp, s, m_ = locatepeak(cache.r)         # find peak
    norm(U)^2 + norm(V)^2 - 8*rp*π^2, (s/nc*2π, m_) # recover shift in dimensional terms
end

# Same size, use broadcasting
evalconjprod!(U::FTField{n}, V::FTField{n}, cache::XCorrCache{n}) where {n} =
    (cache.R .= conj.(U).*V; nothing)

# Different size, use indexing
function evalconjprod!(U::FTField{n}, V::FTField{n}, cache::XCorrCache{nc}) where {nc, n}
    d = nc>>1
    for j = 0:d, k = -d+1:d
        @inbounds cache.R[k, j] = conj(U[k, j])*V[k, j]
    end
end

# Find peak of correlation function for y shifts m = 0, 1, 2, 3
function locatepeak(r::Field{nc}) where {nc}
    rmax, mmax, smax = r[0, 0], 0, 0
    for m = (0, 1, 2, 3), j = 0:nc-1
        @inbounds val = r[m*div(nc, 4), j]
        if val > rmax
            # m is a y-shift by π/4
            rmax, mmax, smax = val, 2*m, j
        end
    end
    # the maximum correlation, and the s and m shifts
    rmax, smax, mmax
end