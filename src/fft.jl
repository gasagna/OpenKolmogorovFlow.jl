import FFTW: unsafe_execute!, plan_rfft, plan_brfft

export FFT, IFFT, ForwardFFT!, InverseFFT!, even_dealias_size

# ~~~ UTILS ~~~

# Return next even number.
_next_even(n::Int) = ifelse(iseven(n), n, n+1)

# Return even `m`, the minimum size of a `Field{m}` that
# avoids aliasing on a `FTField{n}` of size `n`. The returned
# size is not necessarily good for performance, and one should
# rather select a size that is optimal for his/her own hardware
even_dealias_size(n::Int) = _next_even(3n>>1 + 1)

# ~~~ ALLOCATING VERSIONS - Always Aliased ~~~
# We need the copy on IFFT because irfft does not preserve input
 FFT(u::Field{n, T}) where {n, T} =
    ForwardFFT!(FTField{n, Complex{T}}, similar(u))(FTField(n, Complex{T}), u)
IFFT(U::FTField{n, Complex{T}}, m::Int=n) where {n, T} =
    InverseFFT!(Field{m, T}, similar(U))(Field(m, T), copy(U))


# ~~~ NON ALLOCATING VERSION ~~~

# ~~~ FORWARD TRANSFORM ~~~
struct ForwardFFT!{m, F, P}
    tmp::F
    plan::P
end

# used to call the inner constructor
function _ForwardFFT!(v::Field{m}, tmp::Union{FTField{m}, Void}, flags::UInt32) where {m}
    plan = plan_rfft(v.data, [2, 1], flags=flags)
    ForwardFFT!{m, typeof(tmp), typeof(plan)}(tmp, plan)
end

# For details about the following two methods see below, for InverseFFT

# different size - de-aliased calculations - create temporary
ForwardFFT!(outType::Type{<:FTField{n, Complex{T}}}, u::Field{m, T}, flags::UInt32=FFTW.MEASURE) where {m, n, T} =
    _ForwardFFT!(u, FTField(m, Complex{T}), flags)

# same size - aliased calculations - no temporary needed
ForwardFFT!(outType::Type{<:FTField{m, Complex{T}}}, u::Field{m, T}, flags::UInt32=FFTW.MEASURE) where {m, T} =
    _ForwardFFT!(u, nothing, flags)

# same size - no worries
(f::ForwardFFT!{m})(U::FTField{m}, u::Field{m}) where {m} =
    (unsafe_execute!(f.plan, u.data, U.data); U .*= 1/m^2; U)

# different size - use temporary and shrink
(f::ForwardFFT!{m})(U::FTField, u::Field{m}) where {m} =
    shrinkto!(U, f(f.tmp, u))


# ~~~ INVERSE TRANSFORM ~~~
struct InverseFFT!{m, F, P}
    tmp::F
    plan::P
end

# Used to call the inner constructor. These three arguments are: `V` - a `FTField{m}`
# to create the plan on. `m` will be the size of the grid in physical space. `tmp` -
# this can be an `FTField{m}`, or Void. The first case is for when we do de aliased
# calculations and we need to prolong an `FTField{n}` , `n < m` to a larger grid, before
# performing the transform. It is `nothing` for when we do not need de-aliasing. `flags`
# are the usual FFTW flags.
function _InverseFFT!(V::FTField{m}, tmp::Union{FTField{m}, Void}, flags::UInt32)  where {m}
    plan = plan_brfft(V.data, m, [2, 1], flags=flags)
    InverseFFT!{m, typeof(tmp), typeof(plan)}(tmp, plan)
end

# The following two methods are those used to create the FFTW plans. The input
# arguments are
#    - outType : the type of the output data that the transform produces
#    - U       : an instance of the input data that the transform accepts
#    - flags   : FFTW flags. See FFTW documentation. Default is MEASURE.
#
# When an InverseFFT object has been created, it will adhere to a callable
# interface. For instance:
#
# julia> f = InverseFFT!(Field{m, T}, U::FTField{m, Complex{T}})
# julia> f(out, U)
#
# Note that the fieldsize of outType and U is used to determine whether
# the transform requires dealiasing or not.
#
# different size - dealised calculations - needs padding, so create a temporary
InverseFFT!(outType::Type{<:Field{m, T}}, U::FTField{n, Complex{T}}, flags::UInt32=FFTW.MEASURE) where {m, n, T} =
    (tmp = FTField(m, Complex{T}); _InverseFFT!(tmp, tmp, flags))

# same size - aliased calculations - do not create temporary
InverseFFT!(outType::Type{<:Field{m, T}}, U::FTField{m, Complex{T}}, flags::UInt32=FFTW.MEASURE) where {m, T} =
    _InverseFFT!(U, nothing, flags)

# same size - no worries
@inline (i::InverseFFT!{m})(u::Field{m}, U::FTField{m}) where {m} =
    (unsafe_execute!(i.plan, U.data, u.data); u)

# different size - grow `n` up to `m` then transform. Up to user
# to make sure that `U` has size that makes sense for dealiasing.
@inline (i::InverseFFT!{m})(u::Field{m}, U::FTField) where {m} =
    i(u, growto!(i.tmp, U))