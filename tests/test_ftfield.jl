using OpenKolmogorovFlow
using Base.Test

@testset "input                                  " begin
    @test_throws TypeError     FTField(rand(5, 5))
    @test_throws ArgumentError FTField(5)
    @test_throws ArgumentError FTField(im*rand(5, 5))
    @test_throws ArgumentError FTField(im*rand(5, 3))
    @test_throws ArgumentError FTField(im*rand(4, 4))
    @test_throws ArgumentError FTField(im*rand(4, 2))
end

@testset "symmetries                             " begin
    for n = [4, 24]
        U = FFT(Field(randn(n, n)))
        d = n>>1
        for k = (-d+1):(d-1), j = (-d+1):(d-1)
            # it is not exactly zero, because of FFTW does
            # not ensure exact conjugate symmetry
            @test abs(U[j, k] - conj(U[-j, -k])) < 1e-16
        end
    end
end

@testset "cartesian indexing                     " begin
    # construct data with appropriate symmetries
    data = [1+2im  9+10im 17+0im
            3+4im 11+12im 19+20im
            5+0im 13+14im 21+0im
            3-4im 15+16im 19-20im]
    u = FTField(data)

    # inbounds
    @test u[ 0, -1] ==  9-10im
    @test u[ 1, -1] == 15-16im
    @test u[ 2, -1] == 13-14im
    @test u[-1, -1] == 11-12im

    @test u[ 0,  0] ==  1+2im
    @test u[ 1,  0] ==  3+4im
    @test u[ 2,  0] ==  5+0im
    @test u[-1,  0] ==  3-4im

    @test u[ 0,  1] ==  9+10im
    @test u[ 1,  1] == 11+12im
    @test u[ 2,  1] == 13+14im
    @test u[-1,  1] == 15+16im

    @test u[ 0,  2] == 17+0im
    @test u[ 1,  2] == 19+20im
    @test u[ 2,  2] == 21+0im
    @test u[-1,  2] == 19-20im

    # out of bounds
    @test_throws BoundsError u[ 0, -3]
    @test_throws BoundsError u[ 0,  3]
    @test_throws BoundsError u[ 0, -2]
    @test_throws BoundsError u[ 1, -2]
    @test_throws BoundsError u[-2, -2]
    @test_throws BoundsError u[-2,  0]
    @test_throws BoundsError u[-3,  0]
    @test_throws BoundsError u[ 3,  0]
end

@testset "linear indexing                        " begin
    data = rfft(randn(6, 6), [2, 1])
    u = FTField(data)
    for i in eachindex(data)
        @test u[i] == data[i]
    end
end

@testset "similar                                " begin
    for n = [2, 4, 8]
        for T in [Float64, Float32]
            U = FTField(n, Complex{T})
            V = similar(U)
            @test typeof(V) == FTField{n, Complex{T}, Matrix{Complex{T}}}
        end
    end
    # different size
    U = FTField(6, Complex{Float64})
    V = similar(U, 8)
    @test typeof(V) == FTField{8, Complex{Float64}, Matrix{Complex{Float64}}}
end

@testset "transform                              " begin
    n = 4
    x, y = make_grid(n)

    # ~~~ Along k ~~~
    # for k < N/2 the transform has half the value
    f = cos.(0.*x .+ y) .+ 2.*sin.(0.*x .+ y)
    U = FFT(Field(f))
    @test U[ 1, 0] ≈ 0.5*(1 - 2*im)
    @test U[-1, 0] ≈ 0.5*(1 + 2*im)

    # for k = N/2 we count the value in full. The coefficient
    # must be real, hence the sine part gets killed.
    f = cos.(0.*x .+ 2.*y) .+ 2.*sin.(0.*x .+ 2.*y)
    U = FFT(Field(f))
    @test U[ 2, 0] ≈ 1

    # for k = 0 we count the value in full. The coefficient
    # must be real
    f = 2.0.*cos.(0.*x .+ 0.*y)
    U = FFT(Field(f))
    @test U[ 0, 0] ≈ 2.0

    # for k = -N/2 the result should be equal to that of k = +N/2
    fp = cos.(0.*x .+ 2.*y); Up = FFT(Field(f))
    fm = cos.(0.*x .- 2.*y); Um = FFT(Field(f))
    @test Up == Um

    # ~~~ Along j ~~~
    # for j < N/2 the transform has half the value
    f = cos.(1.*x .+ 0.*y) .+ 2.*sin.(1.*x .+ 0.*y)
    U = FFT(Field(f))
    @test U[0,  1] ≈ 0.5*(1 - 2*im)
    @test U[0, -1] ≈ 0.5*(1 + 2*im)

    # for j = N/2 we count the value in full. The coefficient
    # must be real, hence the sine part gets killed.
    f = cos.(2.*x .+ 0.*y) .+ 2.*sin.(2.*x .+ 0.*y)
    U = FFT(Field(f))
    @test U[0,  2] ≈ 1

    # for j = -N/2 the result should be equal to that of j = +N/2
    fp = cos.(+ 2.*x .+ 0.*y); Up = FFT(Field(f))
    fm = cos.(- 2.*x .+ 0.*y); Um = FFT(Field(f))
    @test Up == Um

    # ~~~ Along both j and k ~~~
    f = cos.(1.*x .+ 1.*y) .+ 2.*sin.(1.*x .+ 1.*y)
    U = FFT(Field(f))
    @test U[ 1,  1] ≈ 0.5*(1 - 2im)
    @test U[-1, -1] ≈ 0.5*(1 + 2im)

    f = cos.(1.*x .+ 2.*y) .+ 3.*sin.(1.*x .+ 2.*y)
    U = FFT(Field(f))
    @test U[ 2,  1] ≈ 0.5*(1 - 3im)

    f = cos.(2.*x .+ 1.*y) .+ 4.*sin.(2.*x .+ 1.*y)
    U = FFT(Field(f))
    @test U[ 1,  2] ≈ 0.5*(1 - 4im)

    f = 2.*cos.(2.*x .- 1.*y) .+ 1.*sin.(2.*x .- 1.*y)
    U = FFT(Field(f))
    @test U[-1,  2] ≈ 0.5*(2 - 1im)

    # the corner wave number is real
    f = cos.(2.*x .+ 2.*y) .+ 4.*sin.(2.*x .+ 2.*y)
    U = FFT(Field(f))
    @test U[ 2,  2] ≈ 1
end

@testset "grow/shrink-to!                        " begin
    @testset "growto!                            " begin
        # growing a field should not change its energy
        U = FFT(Field(randn(6, 6)))
        for m = [6, 12, 24, 48]
            @test norm(U) ≈ norm(growto!(FTField(m, Complex{Float64}), U))
        end
    end
    @testset "shrinkto!                          " begin
        data = Complex{Float64}[1+0im 5-7im 5+4im 16+0im
                                3+4im 3+3im 4-2im 19+20im
                                4-1im 7+3im 2+3im 22+2im
                                6+0im 5+8im 6-3im 12+0im
                                4+1im 1+3im 4+0im 22-2im
                                3-4im 5-3im 4-4im 19-20im]

        # define field
        v = FTField(data)

        # same size
        w = FTField(6)
        shrinkto!(w, v)
        @test w.data == Complex{Float64}[1+0im 5-7im 5+4im 16+0im
                                         3+4im 3+3im 4-2im 19+20im
                                         4-1im 7+3im 2+3im 22+2im
                                         6+0im 5+8im 6-3im 12+0im
                                         4+1im 1+3im 4+0im 22-2im
                                         3-4im 5-3im 4-4im 19-20im]

        # smaller
        w = FTField(4)
        shrinkto!(w, v)
        @test w.data == Complex{Float64}[1+0*im 5-7im 10+0*im
                                         3+4*im 3+3im  4-3*im
                                         8+0*im 7+3im  4+0*im
                                         3-4*im 5-3im  4+3*im]

        # smaller
        w = FTField(2)
        shrinkto!(w, v)
        @test w.data == Complex{Float64}[1+0*im 10+0*im
                                         6+0*im  6+0*im]
    end
end