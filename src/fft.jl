export FT, IFT

# ~~~ ALLOCATING VERSION ~~~ 
 FT(u::Field{n})   where n = FTField( rfft(u.data,    [2, 1])/n^2)
IFT(U::FTField{n}) where n =   Field(brfft(U.data, n, [2, 1]))


# ~~~ NON ALLOCATING VERSION ~~~

# ~ Forward transform
struct ForwardFFT{n, P}
    plan::P
end

function ForwardFFT(u::Field{n}, FFTWflags::UInt32=FFTW.MEASURE) where {n}
    p = plan_rfft(u.data, [2, 1], flags=FFTWflags)
    ForwardFFT{n, typeof(p)}(p)
end

function (f::ForwardFFT{n})(U::FTField{n}, u::Field{n}) where {n}
    A_mul_B!(U.data, f.plan, u.data)
    U .*= 1/n^2
end


# ~ Inverse Transform
struct InverseFFT{n, P}
    plan::P
end

function InverseFFT(U::FTField{n}, FFTWflags::UInt32=FFTW.MEASURE) where {n}
    p = plan_brfft(U.data, n, [2, 1], flags=FFTWflags)
    InverseFFT{n, typeof(p)}(p)
end

function (i::InverseFFT{n})(u::Field{n}, U::FTField{n}) where {n}
    A_mul_B!(u.data, i.plan, U.data)
end