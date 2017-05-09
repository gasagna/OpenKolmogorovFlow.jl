using OpenKolmogorovFlow
using Base.Test

@testset "input type   " begin
    @test_throws ArgumentError Field(5)
    @test_throws ArgumentError Field(randn(4, 5))
    @test_throws ArgumentError Field(randn(4, 6))
end

@testset "indexing     " begin
    # getindex
    u = Field([0  4  8 12
               1  5  9 13
               2  6 10 14
               3  7 11 15])
    # lower left corner of domain, at 
    # (x, y) = (0, 0)
    @test u[0, 0] == 0 
    @test u[2, 1] == 6
    @test_throws BoundsError u[-1, 0]
    @test_throws BoundsError u[0, -1]
    @test_throws BoundsError u[4, 0]
    @test_throws BoundsError u[0, 4]

    # setindex
    u[0, 0] = 1
    @test u[0, 0] == 1
    u[0, 0] = 0
    @test u[0, 0] == 0

    # linear indexing
    for i = 0:15
        @test u[i+1] == i
    end

    for i = 1:16
        u[i] = 2i
        @test u[i] == 2i
    end
end