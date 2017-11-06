export DistanceCache, distance!

# Storage type to assist the computation of distance between two vorticity field
struct DistanceCache{m, T, S<:Field{m, T}, D<:FTField{m, Complex{T}}}
            r::S                    # the cross correlation field
       r_line::Vector{T}            # cross correlation along Δy = m*π/4
            R::D                    # the cross correlation transform
       R_line::Vector{Complex{T}}   # Fourier series or r_line
         plan                       # 2d plan for FFT
    plan_line                       # 1d plan for FFT
    function DistanceCache(m::Int, ::Type{T}=Float64, flags::UInt32=FFTW.PATIENT) where {T}
        # we want to index exact location for vertical integer shifts
        rem(m, 4) == 0 || error("cache size / 4 must be an even number")
        r, R = Field(m, T), FTField(m, Complex{T})
        r_line, R_line = zeros(T, m), zeros(Complex{T}, m)
             plan = plan_brfft(R.data,  m, [2, 1], flags=flags) # inverse transform
        plan_line = plan_rfft(r_line,              flags=flags) # forward transform
        new{m, T, typeof(r), typeof(R)}(r, r_line, R, R_line, plan, plan_line)
    end
end

"""
    distance!(Ω₁, Ω₂, cache)

Calculate the minimum distance (L2 norm) between vorticity fields `Ω₁` and `Ω₂`,
among all possible continuous shift `s ∈ [0, 2π)` along the `x` direction and
for discrete vertical shifts `m*π/4, m ∈ {0, 2, 4, 6}`. It returns the minimum
distance `d` and the shifts `s` and `m`, such that

julia> dotdiff(shift!(Ω₁, (s, m)), Ω₂) ≈ d

where the equality holds approximately because different algorithms are used.
This function exploits the correlation theorem to compute the optimal shift
more efficiently.

We currently use the algorithm that is fastest and allocates less memory, but
can be significantly less precise to to cancellation errors. This is not a big
issue, because the main purpose of this function is to detect near recurrences
and we do not care if the distance is not calculated to 16 decimal digits, but
only, say, 11.

The object `cache` is created as

julia> cache = DistanceCache(nc)

where `nc` is the size of the grid where the correlation function is evaluated
on. This can be the same size as the fields `Ω₁` and `Ω₂` above, or lower, for
faster, but maybe less accurate calculations. The number `nc/4` must be even.
"""
function distance!(Ω₁::FTField{n}, Ω₂::FTField{n}, cache::DistanceCache{nc, T}) where {n, nc, T}
    nc ≤ n || throw(ArgumentError("cache size must be lower or equal that field size"))
    _evalconjprod!(Ω₁, Ω₂, cache)                                # evaluate product
    unsafe_execute!(cache.plan, cache.R.data, cache.r.data)      # inverse transform
    s_opt, m_opt, r_opt = _peakdetect(cache)               # find peak
    return norm(Ω₁)^2 + norm(Ω₂)^2 - 8*π^2*r_opt, (s_opt, m_opt) # apply definition
end


# ~~~~ Below is the private interface, mainly helper functions ~~~

# Same size, use broadcasting
_evalconjprod!(Ω₁::FTField{n}, Ω₂::FTField{n}, cache::DistanceCache{n}) where {n} =
    (cache.R .= conj.(Ω₁).*Ω₂; nothing)

# Different size, use indexing
function _evalconjprod!(Ω₁::FTField{n}, Ω₂::FTField{n}, cache::DistanceCache{nc}) where {nc, n}
    d = nc>>1
    @inbounds begin
        for j = 0:d, k = -d+1:d
            cache.R[k, j] = conj(Ω₁[k, j])*Ω₂[k, j]
        end
        # satisfy symmetries
        for k = 0:d
            cache.R[k, d] = conj(cache.R[-k, d])
        end
        cache.R[0, d] = 2*real(cache.R[0, d])
        cache.R[d, 0] = 2*real(cache.R[d, 0])
        cache.R[d, d] = 2*real(cache.R[d, d])
    end
end

# Copy value of correlation function at `Δy = m_opt*π/4` into a temporary
# buffer and compute its Fourier transform, so that we can use it for
# interpolation and for finding the correlation peak.
function _fillline!(cache::DistanceCache{nc}, m_opt::Int) where {nc}
    I = div(m_opt, 2)*div(nc, 4)
    @inbounds @simd for j = 0:nc-1 # copy line
        cache.r_line[j+1] = cache.r[I, j]
    end
    unsafe_execute!(cache.plan_line, cache.r_line, cache.R_line)
    cache.R_line .*= 1/nc
    return nothing
end

# Find peak of correlation function for y shifts m = 0, 1, 2 and 3.
function _peakdetect(cache::DistanceCache{nc}) where {nc}
    # scan correlation data to find maximum
    r_opt, m_opt, s_opt_int = cache.r[0, 0], 0, 0
    @inbounds begin
        for m = (0, 2, 4, 6), s = 0:nc-1
            p = cache.r[div(m, 2)*div(nc, 4), s]
            if p > r_opt
                r_opt, m_opt, s_opt_int = p, m, s
            end
        end
    end
    # now get line of the optimum and compute FFT
    _fillline!(cache, m_opt)
    # use newton method to find correlation peak accurately
    s_opt, r_max = _peroptim(cache.R_line, s_opt_int/nc*2π)
    # return x shift, y shift and correlation maximum
    return s_opt, m_opt, r_max
end

# Interpolate periodic function at `x` from its Fourier series `R`. Returns
# value and first and second derivatives at `x`. Only works for even sized
# input arrays (cover the use cases in this package). We assume the grid
# is over the domain [0, 2π].
function _perinterp(R::Vector{Complex{T}}, x::Real) where {T}
    lX = length(R)
    lX % 2 == 0 || throw(ArgumentError("array length must be even"))
    @inbounds begin
        val, grad, hess  = zero(T), zero(T), zero(T)
        for kk = 2:lX-1
            coskx, sinkx = cos((kk-1)*x), sin((kk-1)*x)
            rexk, imxk = reim(R[kk])
            val  +=  (rexk*coskx - imxk*sinkx)
            grad -=  (rexk*sinkx + imxk*coskx)*(kk-1)
            hess -=  (rexk*coskx - imxk*sinkx)*(kk-1)^2
        end
        val  = 2*val  + real(R[lX])*cos((lX-1)*x) + real(R[1])
        grad = 2*grad - real(R[lX])*sin((lX-1)*x)*(lX-1)
        hess = 2*hess - real(R[lX])*cos((lX-1)*x)*(lX-1)^2
    end
    val, grad, hess
end

# Use Newton method to find an extremum of the periodic function defined
# by the Fourier Series `R` near `x`. We assume that this is going to
# converge in less than 10 iteration, otherwise we raise an error.
function _peroptim(R::Vector{<:Complex}, x::Real, stol::Real=1e-12)
    for i = 1:10 # this must end
        val, grad, hess = _perinterp(R, x)
        Δ = grad/hess
        x = x - Δ
        abs(Δ) < stol && return (x, _perinterp(R, x)[1])
    end
    throw(error("iterations did not converge"))
end