using OpenKolmogorovFlow
using Base.Test

@testset "inner product                          " begin
    x, y = make_grid(10)

    # tolerance on integrals
    TOL = 1e-14
    
    # init rng
    srand(0)

    # all cosine waves fit properly on the grid
    for j = -5:5, k =-5:5
        a, b = rand(), rand()
        u = a.*cos.(j.*x.+k.*y); U = FFT(Field(u))
        v = b.*cos.(j.*x.+k.*y); V = FFT(Field(v))
        j == k == 0 || @test abs(inner(U, V) - a*b*2π^2) < TOL
    end

    # these does not fit on the grid, so it is zero
    u = sin.(0.*x.+5.*y); U = FFT(Field(u))
    v = sin.(0.*x.+5.*y); V = FFT(Field(v))
    @test abs(inner(U, V) - 0) < TOL

    u = sin.(5.*x.+0.*y); U = FFT(Field(u))
    v = sin.(5.*x.+0.*y); V = FFT(Field(v))
    @test abs(inner(U, V) - 0) < TOL

    u = sin.(5.*x.+5.*y); U = FFT(Field(u))
    v = sin.(5.*x.+5.*y); V = FFT(Field(v))
    @test abs(inner(U, V) - 0) < TOL

    # orthogonal fields have zero inner product
    u = sin.(1.*x.+1.*y); U = FFT(Field(u))
    v = sin.(1.*x.+2.*y); V = FFT(Field(v))
    @test abs(inner(U, V) -   0) < TOL

    u = cos.(5.*x.+1.*y) .+ sin.(4.*x.+1.*y); U = FFT(Field(u))
    v = cos.(1.*x.+2.*y) .+ sin.(1.*x.+2.*y); V = FFT(Field(v))
    @test abs(inner(U, V) -   0) < TOL

    # count only what is not orthogonal
    u = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); U = FFT(Field(u))
    v = cos.(1.*x.+0.*y) .+ sin.(1.*x.+2.*y); V = FFT(Field(v))
    @test abs(inner(U, V) - 2π^2) < TOL

    u = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); U = FFT(Field(u))
    v = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); V = FFT(Field(v))
    @test abs(inner(U, V) - 2*2π^2) < TOL

    u = cos.(1.*x.+5.*y) .+ sin.(1.*x.+5.*y); U = FFT(Field(u))
    v = cos.(1.*x.+5.*y) .+ sin.(1.*x.+5.*y); V = FFT(Field(v))
    @test abs(inner(U, V) - 2*2π^2) < TOL

    u = 0.4*cos.(5.*x.+1.*y) + 0.4*sin.(5.*x.+2.*y); U = FFT(Field(u))
    v = 0.3*cos.(5.*x.+1.*y) + 0.3*sin.(5.*x.+2.*y); V = FFT(Field(v))
    @test abs(inner(U, V) - 0.24*2π^2) < TOL
end

@testset "inner product performance              " begin
    n = 100
    u = randn(n, n); U = FFT(Field(u))
    v = randn(n, n); V = FFT(Field(v))
    @test minimum([@elapsed inner(U, V) for i = 1:100000]) < 4*10.0^(-6)
end

@testset "norm                                   " begin
    x, y = make_grid(10)

    # tolerance on integrals
    TOL = 1e-14
 
    for (val, u) in [(sqrt(2)*π, sin.(1.*x.+1.*y)), 
                     (2*π,       cos.(1.*x.+2.*y) .+ sin.(1.*x.+2.*y)),
                     (2*π,       cos.(0.*x.+5.*y) .+ sin.(1.*x.+5.*y)),
                     (sqrt(2)*π, cos.(0.*x.+5.*y) .+ sin.(0.*x.+5.*y))]
        U = FFT(Field(u))
        @test abs(norm(U) - val) < TOL
        @test_throws ArgumentError norm(U, 1)
    end
end

@testset "diff                                   " begin
    x, y = make_grid(10)

    # tolerance on integrals
    TOL = 1e-14
 
    u = sin.(1.*x.+1.*y); U = FFT(Field(u))
    v = sin.(1.*x.+1.*y); V = FFT(Field(v))
    @test innerdiff(U, V) == 0

    u = 2*sin.(1.*x.+1.*y); U = FFT(Field(u))
    v =   sin.(1.*x.+1.*y); V = FFT(Field(v))
    @test innerdiff(U, V) == 2π^2

    # corner case
    u = 2*sin.(5.*x.+1.*y); U = FFT(Field(u))
    v =   sin.(5.*x.+1.*y); V = FFT(Field(v))
    @test innerdiff(U, V) == 2π^2
end