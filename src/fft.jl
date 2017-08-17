export FFT, IFFT, ForwardFFT!, InverseFFT!, even_dealias_size

# ~~~ UTILS ~~~

# Return next even number. 
_next_even(n::Int) = ifelse(iseven(n), n, n+1)

# Return even `m`, the minimum size of a `Field{m}` that 
# avoids aliasing on a `FTField{n}` of size `n`.
even_dealias_size(n::Int) = _next_even(3n>>1 + 1)

# ~~~ ALLOCATING VERSIONS - Always Aliased ~~~ 
# We need the copy on IFFT because irfft does not preserve input
 FFT(u::Field{n})   where {n} = ForwardFFT!(FTField{n}, similar(u))(FTField(n),    u)
IFFT(U::FTField{n}) where {n} = InverseFFT!(Field{n},   similar(U))(Field(n), copy(U))


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

# same size - aliased calculations - no temporary needed
ForwardFFT!(::Type{<:FTField{m}}, u::Field{m}, flags::UInt32=FFTW.MEASURE) where {m} =
    _ForwardFFT!(u, nothing, flags)

# different size - de-aliased calculations - create temporary
ForwardFFT!(::Type{<:FTField}, u::Field{m}, flags::UInt32=FFTW.MEASURE) where {m} =
    _ForwardFFT!(u, FTField(m), flags)

# same size - no worries
(f::ForwardFFT!{m})(U::FTField{m}, u::Field{m}) where {m} = 
    (Base.FFTW.unsafe_execute!(f.plan, u.data, U.data); U .*= 1/m^2; U)

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

# different size - dealised calculations - needs padding
InverseFFT!(::Type{<:Field{m}}, U::FTField, flags::UInt32=FFTW.MEASURE) where {m} =
    (tmp = FTField(m); _InverseFFT!(tmp, tmp, flags))

# same size - aliased calculations - do not create temporary
InverseFFT!(::Type{<:Field{m}}, U::FTField{m}, flags::UInt32=FFTW.MEASURE) where {m} =
    _InverseFFT!(U, nothing, flags)

# same size - no worries
@inline (i::InverseFFT!{m})(u::Field{m}, U::FTField{m}) where {m} =
    (Base.FFTW.unsafe_execute!(i.plan, U.data, u.data); u)

# different size - grow `n` up to `m` then transform. Up to user
# to make sure that `U` has size that makes sense for dealiasing.
@inline (i::InverseFFT!{m})(u::Field{m}, U::FTField) where {m} = 
    i(u, growto!(i.tmp, U))