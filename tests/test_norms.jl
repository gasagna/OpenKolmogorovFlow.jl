using OpenKolmogorovFlow
using Base.Test

@testset "dot product                          " begin
    n, m = 5, 10
    x, y = make_grid(m)

    # tolerance on integrals
    TOL = 1e-14

    # init rng
    srand(0)

    # all cosine waves fit properly on the grid
    for j = -n:n, k =-n:n
        a, b = rand(), rand()
        u = a.*cos.(j.*x.+k.*y); U = FFT(Field(u), n)
        v = b.*cos.(j.*x.+k.*y); V = FFT(Field(v), n)
        j == k == 0 || @test abs(dot(U, V) - a*b/2) < TOL
    end

    # orthogonal fields have zero dot product
    u = sin.(1.*x.+1.*y); U = FFT(Field(u), n)
    v = sin.(1.*x.+2.*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 0) < TOL

    u = cos.(5.*x.+1.*y) .+ sin.(4.*x.+1.*y); U = FFT(Field(u), n)
    v = cos.(1.*x.+2.*y) .+ sin.(1.*x.+2.*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 0) < TOL

    # count only what is not orthogonal
    u = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); U = FFT(Field(u), n)
    v = cos.(1.*x.+0.*y) .+ sin.(1.*x.+2.*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 0.5) < TOL

    u = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); U = FFT(Field(u), n)
    v = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 2*0.5) < TOL

    u = cos.(1.*x.+5.*y) .+ sin.(1.*x.+5.*y); U = FFT(Field(u), n)
    v = cos.(1.*x.+5.*y) .+ sin.(1.*x.+5.*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 2*0.5) < TOL

    u = 0.4*cos.(5.*x.+1.*y) + 0.4*sin.(5.*x.+2.*y); U = FFT(Field(u), n)
    v = 0.3*cos.(5.*x.+1.*y) + 0.3*sin.(5.*x.+2.*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 0.24*0.5) < TOL
end

@testset "dot product performance              " begin
    m, n = 49, 49
    U = FFT(Field(m, (x, y)->rand()), n)
    V = FFT(Field(m, (x, y)->rand()), n)
    @test minimum([@elapsed dot(U, V) for i = 1:100000]) < 8*10.0^(-6)
end

@testset "norm                                   " begin
    n, m = 10, 10
    x, y = make_grid(m)

    # tolerance on integrals
    TOL = 1e-14

    for (val, u) in [(sqrt(2)*π, sin.(1.*x.+1.*y)),
                     (2*π,       cos.(1.*x.+2.*y) .+ sin.(1.*x.+2.*y)),
                     (2*π,       cos.(0.*x.+5.*y) .+ sin.(1.*x.+5.*y)),
                     (2*π,       cos.(0.*x.+5.*y) .+ sin.(0.*x.+5.*y))]
        U = FFT(Field(u), n)
        @test abs(norm(U) - val/2π) < TOL
        @test_throws ArgumentError norm(U, 1)
    end
end

@testset "diff                                   " begin
    n, m = 5, 10
    x, y = make_grid(m)

    # tolerance on integrals
    TOL = 1e-14

    u = sin.(1.*x.+1.*y); U = FFT(Field(u), n)
    v = sin.(1.*x.+1.*y); V = FFT(Field(v), n)
    @test dotdiff(U, V) == 0

    u = 2*sin.(1.*x.+1.*y); U = FFT(Field(u), n)
    v =   sin.(1.*x.+1.*y); V = FFT(Field(v), n)
    @test abs(dotdiff(U, V) - 0.5) < TOL
end