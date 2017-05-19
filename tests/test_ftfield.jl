using OpenKolmogorovFlow
using Base.Test

@testset "input        " begin
    @test_throws TypeError     FTField(   rand(5, 5))
    @test_throws ArgumentError FTField(5)
    @test_throws ArgumentError FTField(im*rand(5, 5))
    @test_throws ArgumentError FTField(im*rand(5, 3))
    @test_throws ArgumentError FTField(im*rand(4, 4))
    @test_throws ArgumentError FTField(im*rand(4, 2))
end

@testset "symmetries   " begin
    for n = [2, 24]
        d = n>>1
        U = FTField(rfft(randn(n, n), [2, 1]))
        for k = -d:d, j=-d:d
            @test abs(U[j, k] - conj(U[-j, -k])) < 1e-13
        end
    end
end

@testset "cartesian    " begin
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

@testset "symmetry     " begin
    for n = [2, 24]
        d = n>>1
        U = FTField(rfft(randn(n, n), [2, 1]))
        for k = -d:d, j=-d:d
            @test abs(U[j, k] - conj(U[-j, -k])) < 1e-13
        end
    end
end

@testset "linear       " begin
    data = rfft(randn(6, 6), [2, 1])
    u = FTField(data)
    for i in 1:length(data)
        @test u[i] == data[i]
    end
end

@testset "similar      " begin
    for n = [2, 4, 8]
        for T in [Float64, Float32]
            U = FTField(n, Complex{T})
            V = similar(U)
            @test typeof(V) == FTField{n, Complex{T}, Matrix{Complex{T}}}
        end
    end
end

@testset "transform    " begin
    N = 4
    x = linspace(0, 2π, N+1)[1:end-1]'
    y = linspace(0, 2π, N+1)[1:end-1]

    # ~~~ Along j ~~~
    # for j < N/2 we count the value twice
    c = randn() + im*randn()
    û = FT(Field(fun(x, y, c, 1, 0)))
    @test û[0,  1] ≈      c/2
    @test û[0, -1] ≈ conj(c/2)

    # for j = N/2 we count the value in full
    c = randn() + im*randn()
    û = FT(Field(fun(x, y, c, 2, 0)))
    @test û[0,  2] ≈ real(c)
    @test û[0, -2] ≈ real(c)

    # for j = -N/2 the result should be equal to that of j = +N/2
    c = randn() + im*randn()
    ûp = FT(Field(round.(fun(x, y, c, -2, 0), 3)))
    ûm = FT(Field(round.(fun(x, y, c, +2, 0), 3)))
    @test ûp == ûm

    # ~~~ Along k ~~~
    # for k < N/2 we count the value twice
    c = randn() + im*randn()
    û = FT(Field(fun(x, y, c, 0, 1)))
    @test û[ 1, 0] ≈      c/2
    @test û[-1, 0] ≈ conj(c/2)

    # for k = N/2 we count the value in full
    c = randn() + im*randn()
    û = FT(Field(fun(x, y, c, 0, 2)))
    @test û[ 2, 0] ≈ real(c)
    @test û[-2, 0] ≈ real(c)

    # for j = -N/2 the result should be equal to that of j = +N/2
    c = randn() + im*randn()
    ûp = FT(Field(round.(fun(x, y, c, 0, -2), 3)))
    ûm = FT(Field(round.(fun(x, y, c, 0, +2), 3)))
    @test ûp == ûm

    # ~~~ Along both j and k ~~~
    c = randn() + im*randn()
    û = FT(Field(fun(x, y, c, 2, 2)))
    @test û[ 2,  2] ≈ real(c)
    @test û[-2, -2] ≈ real(c)

    c = randn() + im*randn()
    û = FT(Field(fun(x, y, c, 1, 2)))
    @test û[ 2,  1] ≈      c/2
    @test û[-2, -1] ≈ conj(c/2)

    c = randn() + im*randn()
    û = FT(Field(fun(x, y, c, 2, 1)))
    @test û[ 1,  2] ≈      c/2
    @test û[-1, -2] ≈ conj(c/2)

    c = randn() + im*randn()
    û = FT(Field(fun(x, y, c, -1, 2)))
    @test û[ 2, -1] ≈      c/2
    @test û[-2,  1] ≈ conj(c/2)
end