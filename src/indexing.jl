# Return the storage location 'i' corresponding to wave numbers `k`, `j`,
# for a grid with `n` points. Since the data in the storage has some 
# symmetries (i.e. for j=0 and j=n/2 the entries for k and -k are conjugate)
# we point to the exact location indicated by the input arguments, so that
# looping over `k` and `j` spans the entire underlying data
@inline KJtoI(k::Int, j::Int, n::Int) =
     j≥0 ? j*n + ifelse(k≥0,  k+1,  k+n+1) : -j*n + ifelse(k≤0, -k+1, -k+n+1)

# Return the 'k' or 'j' wave numbers corresponding to storage location 
# `i`, for a grid with `n` points
@inline ItoJ(i::Int, n::Int) = div(i-1, n)
@inline ItoK(i::Int, n::Int) = (k=rem(i, n)-1; ifelse(k>n>>1, k-n, k))

# Return conjugate of `val` if `j` is negative. Other symmetries are taken
# care of by the `KJtoI` function
@inline rectify(val::Number, j::Int) = ifelse(j<0, conj(val), val)