import FFTW: unsafe_execute!, plan_rfft, plan_brfft

export FFT, IFFT, ForwardFFT!, InverseFFT!, even_dealias_size

# ~~~ UTILS ~~~

# Construct an instance of the type that is the output of 
# application of FFT or inverse FFT on the input type `u`
_oftransftype(u::Field{m, T}, n::Int=m)                           where {m, T} = FTField(n, Complex{T})
_oftransftype(u::FTField{m, Complex{T}}, n::Int=m)                where {m, T} = Field(n, T)
_oftransftype(u::VariationalField{m, Dual{T}}, n::Int=m)            where {m, T} = VariationalFTField(n, Dual{Complex{T}})
_oftransftype(u::VariationalFTField{m, Dual{Complex{T}}}, n::Int=m) where {m, T} = VariationalField(n, Dual{T})

# Return next even number.
_next_even(n::Int) = ifelse(iseven(n), n, n+1)

# Return even `m`, the minimum size of a `Field{m}` 
# that avoids aliasing on a `FTField{n}` of size `n`. 
even_dealias_size(n::Int) = _next_even(3n>>1 + 1)


# ~~~ ALLOCATING VERSIONS ~~~

# We need copies because the plan creation destroys the input
FFT(u::AbstractField{m}) where {m} = 
     ForwardFFT!(m, similar(u))(_oftransftype(u), deepcopy(u))

IFFT(U::AbstractFTField{n}, m::Int=n) where {n} = 
    InverseFFT!(m, similar(U))(_oftransftype(U, m), deepcopy(U))


# ~~~ NON ALLOCATING VERSION ~~~

# forward transform
struct ForwardFFT!{m, F, P}
    tmp::F
    plan::P
    ForwardFFT!{m}(tmp::F, plan::P) where {m, F, P} = new{m, F, P}(tmp, plan)
end

ForwardFFT!(n::Int, u::AbstractField{m}, flags::UInt32=FFTW.PATIENT) where {m} =
    ForwardFFT!{m}(n == m ? nothing : _oftransftype(u), 
                   plan_rfft(parent(state(u)), [2, 1], flags=flags))

# same size - no worries
(f::ForwardFFT!{m})(U::FTField{m}, u::Field{m}) where {m} =
    (unsafe_execute!(f.plan, parent(u), parent(U)); U .*= 1/m^2; U)

(f::ForwardFFT!{m})(U::VariationalFTField{m}, u::VariationalField{m}) where {m} =
    (unsafe_execute!(f.plan, parent(state(u)), parent(state(U)));
     unsafe_execute!(f.plan, parent(prime(u)), parent(prime(U))); U .*= 1/m^2; U)

# different size - use temporary and shrink
(f::ForwardFFT!{m})(U::AbstractFTField, u::AbstractField{m}) where {m} =
    shrinkto!(U, f(f.tmp, u))


# inverse transform
struct InverseFFT!{m, F, P}
    tmp::F # Union{Void, Field}
    plan::P
    InverseFFT!{m}(tmp::F, plan::P) where {m, F, P} = new{m, F, P}(tmp, plan)
end

function InverseFFT!(m::Int, U::AbstractFTField{n}, flags::UInt32=FFTW.PATIENT) where {n}
    tmp   = m == n ? nothing : similar(U, m)
    input = m == n ? U       : tmp
    plan = plan_brfft(parent(state(input)), m, [2, 1], flags=flags)
    InverseFFT!{m}(tmp, plan)
end

# same size - no worries
@inline (i::InverseFFT!{m})(u::Field{m}, U::FTField{m}) where {m} =
    (unsafe_execute!(i.plan, parent(state(U)), parent(state(u))); return u)

# Keep this in mind: FFT(U+ε*u) = FFT(U) + ε*FFT(u)
@inline (i::InverseFFT!{m})(u::VariationalField{m}, U::VariationalFTField{m}) where {m} =
    (unsafe_execute!(i.plan, parent(state(U)), parent(state(u)));
     unsafe_execute!(i.plan, parent(prime(U)), parent(prime(u))); return u)

# different size - grow `n` (unspecified) up to `m` then transform. Up to 
# the user to make sure that `U` has size that makes sense for dealiasing.
@inline (i::InverseFFT!{m})(u::AbstractField{m}, U::AbstractFTField) where {m} =
    i(u, growto!(i.tmp, U))