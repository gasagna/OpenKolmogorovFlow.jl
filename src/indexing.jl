export WaveNumber, @loop_jk, @loop_k

# Transform the wave number k into a row position
# @inline _reindex(k::Int, j::Int, m::Int) = (ifelse(k≥0, k+1, 2m+3+k), j+1)
@inline _reindex(k::Int, j::Int, m::Int) = (ifelse(k≥0, k, 2m+2+k), j)

# allow indexing a specific wavenumber
struct WaveNumber
    k::Int 
    j::Int 
end

# various macros used to simplify indexing
macro loop_k(n, m, expr)
    quote
        for $(esc(:k)) = 0:$(esc(n))
            $(esc(:_k)) = $(esc(:k))
            $(esc(expr))
        end
        for $(esc(:k)) = -$(esc(n)):-1
            $(esc(:_k)) = 2 * $(esc(m)) + 2 + $(esc(:k))
            $(esc(expr))
        end
    end
end

macro loop_jk(n, m, expr)
    quote
        @inbounds for $(esc(:j)) = 0:$(esc(n))
            $(esc(:_j)) = $(esc(:j))
            for $(esc(:k)) = 0:$(esc(n))
                $(esc(:_k)) = $(esc(:k))
                $(esc(expr))
            end
            for $(esc(:k)) = -$(esc(n)):-1
                $(esc(:_k)) = 2 * $(esc(m)) + 2 + $(esc(:k))
                $(esc(expr))
            end
        end
    end
end


# # Location 'i' of wave number vector  (`k`, `j`), for a grid with `2m+2` points.
# @inline KJtoI(k::Int, j::Int, m::Int) =
#      j≥0 ? j*(2m+2) + ifelse(k≥0,  k+1,  k+2m+3) : -j*m + ifelse(k≤0, -k+1, -k+m+1)

# # Return the 'k' or 'j' wave numbers corresponding to storage location
# # `i`, for a grid with `n` points
# # @inline ItoJ(i::Int, n::Int) = div(i-1, n)
# # @inline ItoK(i::Int, n::Int) = (k = rem(i-1, n); ifelse(k > n>>1, k-n, k))
# # @inline ItoKJ(i, n) = (ItoK(i, n), ItoJ(i, n))

# # Return `i`-th column or row where wave numbers `j` or `k` are stored
# @inline JtoI(j::Int) = abs(j) + 1
# @inline KtoI(k::Int, m::Int) = k ≥ 0 ? k + 1 : 2m+2 + k + 1

# # Return conjugate of `val` if `j` is negative. Other symmetries are taken
# # care of by the `KJtoI` function
# @inline rectify(val::Number, j::Int) = j < 0 ? conj(val) : val

# # Tuples of j and k wave numbers
# # js(n::Int) = ntuple(i->i-1,   n>>1+1)
# # ks(n::Int) = ntuple(i->ItoK(i, n), n)