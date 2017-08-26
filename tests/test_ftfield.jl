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
    for n = [2, 24]
        d = n>>1
        U = FTField(rfft(randn(n, n), [2, 1]))
        for k = -d:d, j=-d:d
            @test abs(U[j, k] - conj(U[-j, -k])) < 1e-13
        end
    end
end

@testset "cartesian                              " begin
    # construct data with appropriate symmetries
    data = [1+2im  9+10im 17+0im
            3+4im 11+12im 19+20im
            5+0im 13+14im 21+0im
            3-4im 15+16im 19-20im]
    u = FTField(data)

    # inbounds
    @test u[ 0, -2] == 17-0im
    @test u[ 1, -2] == 19+20im
    @test u[ 2, -2] == 21-0im
    @test u[-2, -2] == 21-0im
    @test u[-1, -2] == 19-20im

    @test u[ 0, -1] ==  9-10im
    @test u[ 1, -1] == 15-16im
    @test u[ 2, -1] == 13-14im
    @test u[-2, -1] == 13-14im
    @test u[-1, -1] == 11-12im

    @test u[ 0,  0] ==  1+2im
    @test u[ 1,  0] ==  3+4im
    @test u[ 2,  0] ==  5+0im
    @test u[-2,  0] ==  5+0im
    @test u[-1,  0] ==  3-4im
    
    @test u[ 0,  1] ==  9+10im
    @test u[ 1,  1] == 11+12im
    @test u[ 2,  1] == 13+14im
    @test u[-2,  1] == 13+14im
    @test u[-1,  1] == 15+16im

    @test u[ 0,  2] == 17+0im
    @test u[ 1,  2] == 19+20im
    @test u[ 2,  2] == 21+0im
    @test u[-2,  2] == 21+0im
    @test u[-1,  2] == 19-20im

    # out of bounds
    for k=-4:4, j=-4:4
        if max(abs(k), abs(j)) > 2
            @test_throws BoundsError u[k, j]
        end
    end
end

@testset "symmetry                               " begin
    for n = [2, 24]
        d = n>>1
        U = FTField(rfft(randn(n, n), [2, 1]))
        for k = -d:d, j=-d:d
            @test abs(U[j, k] - conj(U[-j, -k])) < 1e-13
        end
    end
end

@testset "linear                                 " begin
    data = rfft(randn(6, 6), [2, 1])
    u = FTField(data)
    for i in 1:length(data)
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
end

# this is meant to be a documentation of the properties of the FFT 
# for a transform of 2D data over a grid with even number of points.
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
    @test U[-2, 0] ≈ 1

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
    @test U[0, -2] ≈ 1

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
    @test U[-2, -1] ≈ 0.5*(1 + 3im)

    f = cos.(2.*x .+ 1.*y) .+ 4.*sin.(2.*x .+ 1.*y)
    U = FFT(Field(f))
    @test U[ 1,  2] ≈ 0.5*(1 - 4im)
    @test U[-1, -2] ≈ 0.5*(1 + 4im)

    f = 2.*cos.(2.*x .- 1.*y) .+ 1.*sin.(2.*x .- 1.*y)
    U = FFT(Field(f))
    @test U[-1,  2] ≈ 0.5*(2 - 1im)
    @test U[ 1, -2] ≈ 0.5*(2 + 1im)

    # the corner wave number is real
    f = cos.(2.*x .+ 2.*y) .+ 4.*sin.(2.*x .+ 2.*y)
    U = FFT(Field(f))
    @test U[ 2,  2] ≈ 1
    @test U[-2, -2] ≈ 1
end

@testset "grow/shrink-to!                        " begin
    @testset "growto!                                " begin
        # growing a field should not change its energy
        srand(0)
        n = 6
        u = Field(randn(n, n))
        U = FFT(u)
        for m = [6, 8, 10, 12]
            @test norm(U) == norm(growto!(FTField(m, Complex{Float64}), U))
        end

        data = Complex{Float64}[1+0im  9+10im 16+0im
                                3+4im 11+12im 19+20im
                                6+0im 13+14im 22+0im
                                3-4im 15+16im 19-20im]
        u = FTField(data)

        # smaller size
        v = FTField(2)
        @test_throws ArgumentError growto!(v, u)
        
        # same size
        v = FTField(4)
        growto!(v, u)
        @test v == u

        # larger size
        v = FTField(6)
        growto!(v, u)

        # note how the frequencies (2, 0), (0, 2) and (2, 2) are halved
        @test v.data == Complex{Float64}[  1+0im  9+10im   8+0im   0+0im
                                           3+4im 11+12im  19+20im  0+0im
                                           3+0im 13+14im  11+0im   0+0im
                                           0+0im  0+0im    0+0im   0+0im
                                           3+0im  0+0im    0+0im   0+0im
                                           3-4im 15+16im  19-20im  0+0im]
    end
    @testset "shrinkto!                          " begin
        data = Complex{Float64}[1+0im 5-7im 5+4im 16+0im
                                3+4im 3+3im 4-3im 19+20im
                                4-1im 7+3im 0+3im 22+2im
                                6+0im 5+8im 6-3im 12+0im
                                4+1im 1+3im 4+0im 22-2im
                                3-4im 5-3im 4-4im 19-20im]       

        # define field                                
        v = FTField(data)                                

        # same size                                 
        w = FTField(6)
        shrinkto!(w, v)
        @test w.data == Complex{Float64}[1+0im 5-7im 5+4im 16+0im
                                         3+4im 3+3im 4-3im 19+20im
                                         4-1im 7+3im 0+3im 22+2im
                                         6+0im 5+8im 6-3im 12+0im
                                         4+1im 1+3im 4+0im 22-2im
                                         3-4im 5-3im 4-4im 19-20im]                                  

        # smaller                                         
        w = FTField(4)
        shrinkto!(w, v)
        # TODO: preserve conjugate symmetry of modes of last column
        expected = Complex{Float64}[1+0im          5-7im  2sqrt(5^2+4^2)
                                    3+4im          3+3im  4-3im 
                                    2sqrt(4^2+1^2) 7+3im  sqrt(3^2+4^2)
                                    3-4im          5-3im  4-4im]
        @test sum(abs, w.data - expected) < 1e-16

        # smaller                         
        w = FTField(2)
        shrinkto!(w, v)
        @test w.data == Complex{Float64}[1+0im          2sqrt(5^2+7^2)
                                         2sqrt(3^2+4^2) sqrt(3^2+3^2+5^2+3^2)]
    end
end

@testset "fieldsize                              " begin
    for n in 2:2:10
        @test fieldsize(FTField(n)) == n
        @test fieldsize(FTField{n}) == n
    end
end