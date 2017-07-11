export radial_mean

# Returns number of bins in which wave vectors fall
_nbins(n) = (d = n>>1; round(Int, sqrt(d^2 + d^2)))

# Calculate mean of function `f` acting on the amplitudes of the Fourier
# coefficient, and average over values for all wave vectors that have norm 
# between `n-1/2` and `n+1/2`. This version does not allocate any memory. 
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
        
        # integer wave vector length
        m = round(Int, sqrt(j^2+k^2))
        
        # increment except for (0, 0) mode
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