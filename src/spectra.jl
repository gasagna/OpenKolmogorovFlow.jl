export radial_mean

# Returns number of bins in which wave vectors fall
_nbins(n) = (d = n>>1; round(Int, sqrt(d^2 + d^2)))

# Calculate mean of function `f`, accepting a single complex number argument,
# acting on the amplitudes of the Fourier coefficients, over the wave vectors
# that have norm between `n-1/2` and `n+1/2`. The `out` is overwritten
# such that the element `n` contains this information.
function radial_mean!(f, Ω::FTField{n}, out::Vector, counts::Vector{Int}) where n
    # check input vector has correct length
    length(out) == length(counts) || error("wrong input length")

    # maximum value of component of wave number vector (±d, ±d)
    d = n >> 1

    # zero arrays
    out    .= 0
    counts .= 0

    # loop over wave numbers
    for j = -d:d, k = -d:d

        # integer wave number
        m = round(Int, sqrt(j^2+k^2))

        # increment value and count except for `(0, 0)`
        # mode, for which `m=0`.
        if m > 0
               out[m] += f(Ω[j, k])
            counts[m] += 1
        end
    end

    # average
    return out ./= counts
end

# Allocating version. Use this version, unless you already have buffers
# of memory for the two work vector allocated.
radial_mean(f, Ω::FTField{n}) where {n} =
    radial_mean!(f, Ω, zeros(_nbins(n)), zeros(Int, _nbins(n)))