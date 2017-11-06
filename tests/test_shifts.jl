using Base.Test
using OpenKolmogorovFlow

@testset "norm is shift invariant                " begin
    n = 4
    s, m = 0.123, 3
    data = rand(n, n)
    U = FFT(Field(data))
    V = shifted(U, (s, m))
    @test norm(U) ≈ norm(V)
end

@testset "shift cos right by π/2 get a sine      " begin
    n = 40
    x, y = make_grid(n)
    data = cos.(x .+ 0.*y)
    Ucos = FFT(Field(data))

    data = sin.(x .+ 0.*y)
    Usin = FFT(Field(data))

    s, m = π/2, 0
    @test abs(dotdiff(Usin, shifted(Ucos, (s, m)))) < 1e-20
end

@testset "do it up, by π/4                       " begin
    n = 40
    x, y = make_grid(n)
    data = cos.(0.*x .+ 2.*y)
    Ucos = FFT(Field(data))

    data = sin.(0.*x .+ 2.*y)
    Usin = FFT(Field(data))

    s, m = 0, 1
    @test abs(dotdiff(Usin, shifted(Ucos, (s, m)))) < 1e-20
end

@testset "shift by exact fraction of grid size   " begin
    n = 4
    s, m = π/2, 2

    data = Float64[1 2 3 4;
                   5 6 7 8;
                   9 1 1 2;
                   3 4 5 6]

    U  = FFT(Field(data))
    shift!(U, (s, m))
    @test IFFT(U).data ≈ Float64[6 3 4 5;
                                 4 1 2 3;
                                 8 5 6 7;
                                 2 9 1 1]
end