export radial_mean

# Calculate mean of function `f`, accepting a single complex number argument,
# acting on the amplitudes of the Fourier coefficients. The vector `out` is 
# overwritten such that the element `n` contains this information.
function radial_mean!(f, Ω::FTField{n, m}, out::Vector, counts::Vector{Int}) where {n, m}
    # check input vector has correct length
    length(out) == length(counts) || error("wrong input length")

    # zero arrays
    out    .= 0
    counts .= 0

    # loop over wave numbers
    for j = 1:n, k = -n:n
        if k != 0
            # integer wave number
            nint = round(Int, sqrt(j^2+k^2))

            # increment value and count except for `(0, 0)`
            # mode, for which `nint=0`.
            if nint > 0 && nint ≤ n
                   out[nint] += f(Ω[WaveNumber(k, j)])
                counts[nint] += 1
            end
        end
    end

    # average
    return out ./= counts
end

# Allocating version. Use this version, unless you already have buffers
# of memory for the two work vector allocated.
radial_mean(f, Ω::FTField{n, m}) where {n, m} =
    radial_mean!(f, Ω, zeros(n), zeros(Int, n))