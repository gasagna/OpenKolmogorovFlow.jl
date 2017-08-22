using OpenKolmogorovFlow
using Base.Test

@testset "dot product                            " begin
    x, y = make_grid(10)

    # tolerance on integrals
    TOL = 1e-14
 
    u = sin.(1.*x.+1.*y); U = FFT(Field(u))
    v = sin.(1.*x.+2.*y); V = FFT(Field(v))
    @test abs(dot(U, V) -   0) < TOL

    u = cos.(1.*x.+1.*y) .+ sin.(1.*x.+1.*y); U = FFT(Field(u))
    v = cos.(1.*x.+2.*y) .+ sin.(1.*x.+2.*y); V = FFT(Field(v))
    @test abs(dot(U, V) -   0) < TOL

    u = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); U = FFT(Field(u))
    v = cos.(1.*x.+0.*y) .+ sin.(1.*x.+2.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 2π^2) < TOL

    u = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); U = FFT(Field(u))
    v = cos.(1.*x.+0.*y) .+ sin.(1.*x.+1.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 2*2π^2) < TOL

    u = cos.(0.*x.+5.*y); U = FFT(Field(u))
    v = cos.(0.*x.+5.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 2π^2) < TOL

    u = cos.(1.*x.+5.*y) .+ sin.(1.*x.+5.*y); U = FFT(Field(u))
    v = cos.(1.*x.+5.*y) .+ sin.(1.*x.+5.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 2*2π^2) < TOL

    u = cos.(1.*x.+5.*y); U = FFT(Field(u))
    v = cos.(1.*x.+5.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 2π^2) < TOL

    u = 0.2*cos.(0.*x.+5.*y); U = FFT(Field(u))
    v = 0.4*cos.(0.*x.+5.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 0.08*2π^2) < TOL

    u = 0.2*cos.(4.*x.+4.*y); U = FFT(Field(u))
    v = 0.4*cos.(4.*x.+4.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 0.08*2π^2) < TOL

    u = 0.4*cos.(5.*x.+5.*y); U = FFT(Field(u))
    v = 0.2*cos.(5.*x.+5.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 0.08*2π^2) < TOL

    u = 0.4*cos.(5.*x.+0.*y); U = FFT(Field(u))
    v = 0.2*cos.(5.*x.+0.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 0.08*2π^2) < TOL

    # this does not fit on the grid, so it is zero
    u = sin.(0.*x.+5.*y); U = FFT(Field(u))
    v = sin.(0.*x.+5.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 0) < TOL

    u = sin.(5.*x.+0.*y); U = FFT(Field(u))
    v = sin.(5.*x.+0.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 0) < TOL

    u = sin.(5.*x.+5.*y); U = FFT(Field(u))
    v = sin.(5.*x.+5.*y); V = FFT(Field(v))
    @test abs(dot(U, V) - 0) < TOL
end

@testset "dot product performance                " begin
    u = randn(100, 100); U = FFT(Field(u))
    v = randn(100, 100); V = FFT(Field(v))
    @test minimum([@elapsed dot(U, V) for i = 1:100000]) < 4*10.0^(-6)
end