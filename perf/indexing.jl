using OpenKolmogorovFlow
using BenchmarkTools

const n = 128
const U = FTField(n)
const V = FTField(n)

function baz(dest::FTField{n}, src::FTField{n}) where n
    d = n>>1
    for j = 0:d, k=-d:d
        @inbounds dest[k, j] = src[k, j]
    end
    dest
end

# function bar(dest, src)
#     @simd for i in linearindices(dest)
#         @inbounds dest[i] = src[i]
#     end
#     dest
# end

# @btime bar($U, $V)
@btime baz($U, $V)
# @btime bar($U.data, $V.data)