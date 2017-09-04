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
    evalconjprod!(U, V, cache)                 # evaluate product
    cache.irfft!(cache.r, cache.R)             # inverse transform
    s, m_ = peakdetect(cache.r)                # find peak
    innerdiff(shifted(U, (s, m_)), V), (s, m_) # calculate distance
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

# Find peak of correlation function for y shifts m = 0, 1, 2, 3,
# using a quadratic interpolation of the data points
function peakdetect(r::Field{nc}) where {nc}
    rmax, mmax, smax = r[0, 0], 0, 0
    @inbounds begin
        for m2 = (0, 1, 2, 3), j = 0:nc-1
            # see https://www.dsprelated.com/freebooks/sasp/
            #            Quadratic_Interpolation_Spectral_Peaks.html
            row = m2*div(nc, 4)
            α, β, γ  = r[row, j-1], r[row, j], r[row, j+1]
            peakloc = clamp(0.5*(α-γ)/(α - 2β + γ), -0.5, 0.5)
            peakval = β - 0.25*(α-γ)*peakloc
            if peakval > rmax
                # m is a y-shift by π/4
                rmax, mmax, smax = peakval, 2*m2, j+peakloc
            end
        end
    end
    # return the `s` and `m` shifts in proper units
    smax/nc*2π, mmax
end