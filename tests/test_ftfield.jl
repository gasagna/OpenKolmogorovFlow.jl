using OpenKolmogorovFlow
using Base.Test

@testset "Constructors                           " begin
    # provide wave numbers
    @test_throws ArgumentError FTField(3, 2)
    
    U = FTField(1, 2, Float64)
    @test typeof(U) == FTField{1, 2, Float64, Matrix{Complex{Float64}}}

    # provide data
    @test_throws ArgumentError FTField(im*rand(5, 5), 1)
    @test_throws ArgumentError FTField(im*rand(6, 3), 1)
    @test_throws ArgumentError FTField(im*rand(6, 4), 3)
    
    U = FTField(im*rand(6, 4), 1)
    @test typeof(U) == FTField{1, 2, Float64, Matrix{Complex{Float64}}}
end


@testset "indexing                               " begin
    @testset "getindex                               " begin
        # test data
        data = [1+2im  9+10im 0+0im
                3+4im 11+12im 0+0im
                0+0im  0+0im  0+0im
                3-4im 15+16im 0+0im]
        U = FTField(data, 1)

        # index over underlying data
        for i = 1:4, j = 1:3
            @test U[i, j] == parent(U)[i, j]
        end
        @test_throws BoundsError U[0, 0]

        # test wavenumbers
        @test U[( 0, 0)] ==  1+2im
        @test U[( 1, 0)] ==  3+4im
        @test U[(-1, 0)] ==  3-4im
        @test U[(-1, 1)] == 15+16im
        @test U[( 1, 1)] == 11+12im
        @test U[( 0, 1)] ==  9+10im
        
        @test_throws BoundsError U[(0, -1)]
        @test_throws BoundsError U[(1, -1)]
        @test_throws BoundsError U[(4,  0)]

        # linear indexing
        @test U[1] == 1+2im
        @test U[2] == 3+4im
        @test U[4] == 3-4im
    end
    @testset "setindex!                              " begin
        # test data
        data = [1+2im  9+10im 0+0im
                3+4im 11+12im 0+0im
                0+0im  0+0im  0+0im
                3-4im 15+16im 0+0im]
        U = FTField(data, 1)

        # index over underlying data
        for i = 1:4, j = 1:3
            U[i, j] = i+j + (i-j)*im
            @test parent(U)[i, j] == i+j + (i-j)*im
        end
        @test_throws BoundsError U[0, 0] = 0

        # test wavenumbers
        U[( 0, 0)] =  (1+2im)*2;   @test U[( 0, 0)] == ( 1+2im)*2
        U[( 1, 0)] =  (3+4im)*2;   @test U[( 1, 0)] == ( 3+4im)*2
        U[(-1, 0)] =  (3-4im)*2;   @test U[(-1, 0)] == ( 3-4im)*2
        U[(-1, 1)] = (15+16im)*2;  @test U[(-1, 1)] == (15+16im)*2
        U[( 1, 1)] = (11+12im)*2;  @test U[( 1, 1)] == (11+12im)*2
        U[( 0, 1)] =  (9+10im)*2;  @test U[( 0, 1)] == ( 9+10im)*2
        
        @test_throws BoundsError U[(0, -1)] == 1
        @test_throws BoundsError U[(1, -1)] == 1
        @test_throws BoundsError U[(4,  0)] == 1

        # linear indexing
        U[1] = 2*(1+2im); @test U[1] == 2*(1+2im)
        U[2] = 2*(3+4im); @test U[2] == 2*(3+4im)
        U[4] = 2*(3-4im); @test U[4] == 2*(3-4im)
    end
end

@testset "similar/copy                           " begin
    U = FTField(1, 2); U[(0, 0)] = 1+2im
    V = similar(U)
    W = copy(U)
    @test typeof(V) == FTField{1, 2, Float64, Matrix{Complex{Float64}}}
    @test W[(0, 0)] == 1+2im
end