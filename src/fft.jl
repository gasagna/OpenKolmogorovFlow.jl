import FFTW: unsafe_execute!, plan_rfft, plan_brfft

export FFT, IFFT, ForwardFFT!, InverseFFT!, up_dealias_size, down_dealias_size

# ~~~ UTILS ~~~
# set to zero Fourier coefficients
function _apply_mask(U::AbstractFTField{n, m, T}) where {n, m, T}
    @inbounds begin
        # middle block
        for jj = 1:n+1
            @simd for kk = (n+2):2m+2-n
                U.data[kk, jj] = zero(Complex{T})
            end
        end
        # vertical block
        for jj = n+2:m+2
            @simd for kk = 1:2m+2
                U.data[kk, jj] = zero(Complex{T})
            end
        end
        # zero mean
        U.data[1, 1] = 0
    end
    return U
end

# ensure field satisfies symmetry at all times
function _apply_symmetry(U::AbstractFTField{n, m, T}) where {n, m, T}
    @inbounds @simd for k = 1:m
        pos = U.data[k+1, 1]
        neg = U.data[2*(m+1)-k+1, 1]
        _re = 0.5 * (real(pos) + real(neg))
        _im = 0.5 * (imag(pos) - imag(neg))
        U.data[k+1, 1] = _re + im * _im 
        U.data[2*(m+1)-k+1, 1] = _re - im * _im 
    end
    return U
end

# Return the smallest `m` that avoids aliasing on a `FTField{n, ⋅}`
up_dealias_size(n::Int) = n + n>>1

# Return the largest `n`that avoids aliasing on a `Field{m}`
down_dealias_size(m::Int) = findlast(n->(up_dealias_size(n) ≤ m), 1:m)

# ~~~ NON ALLOCATING VERSION ~~~

# forward transform
struct ForwardFFT!{m, P}
    plan::P
    function ForwardFFT!(u::AbstractField{m}, flags=FFTW.EXHAUSTIVE) where {m}
        plan = plan_rfft(parent(u), [2, 1], flags=flags)
        new{m, typeof(plan)}(plan)
    end
end

# callable interface
(f::ForwardFFT!{m})(U::FTField{n, m}, u::Field{m}) where {n, m} =
    (unsafe_execute!(f.plan, parent(u), parent(U)); 
        U .*= 1/(2m+2)^2; _apply_symmetry(_apply_mask(U)))


# inverse transform
struct InverseFFT!{m, P}
    plan::P
    function InverseFFT!(U::AbstractFTField{n, m}, flags=FFTW.EXHAUSTIVE) where {n, m}
        plan = plan_brfft(parent(U), 2m+2, [2, 1], flags=flags)
        new{m, typeof(plan)}(plan)
    end
end

# callable interface
(i::InverseFFT!{m})(u::Field{m}, U::FTField{n, m}) where {n, m} =
    (unsafe_execute!(i.plan, parent(_apply_mask(U)), parent(u)); u)


# ~~~ ALLOCATING VERSIONS ~~~

# We need copies because the plan destroys the input
function FFT(u::AbstractField{m, T}, n::Int) where {m, T}
    v = copy(u)
    fun = ForwardFFT!(v, FFTW.ESTIMATE); v .= u
    return fun(FTField(n, m, T), v)
 end

function IFFT(U::AbstractFTField{n, m, T}) where {n, m, T}
    V = copy(U)
    fun = InverseFFT!(V, FFTW.ESTIMATE); V .= U
    return fun(Field(m, T), V)
end
