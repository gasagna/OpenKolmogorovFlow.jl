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
@inline ItoK(i::Int, n::Int) = (k = rem(i-1, n); ifelse(k > n>>1, k-n, k))
@inline ItoKJ(i, n) = (ItoK(i, n), ItoJ(i, n))

# Return `i`-th column or row where wave numbers `j` or `k` are stored
@inline JtoI(j::Int) = abs(j) + 1
@inline KtoI(k::Int, n::Int) = k ≥ 0 ? k + 1 : n + k + 1

# Return conjugate of `val` if `j` is negative. Other symmetries are taken
# care of by the `KJtoI` function
@inline rectify(val::Number, j::Int) = j < 0 ? conj(val) : val

# Tuples of j and k wave numbers
js(n::Int) = ntuple(i->i-1,   n>>1+1)
ks(n::Int) = ntuple(i->ItoK(i, n), n)