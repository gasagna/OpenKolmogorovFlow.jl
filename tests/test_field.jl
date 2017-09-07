using OpenKolmogorovFlow
using Base.Test

@testset "input type                             " begin
    @test_throws ArgumentError Field(5)
    @test_throws ArgumentError Field(randn(4, 5))
    @test_throws ArgumentError Field(randn(4, 6))
end

@testset "indexing                               " begin
    # getindex
    u = Field([0  4  8 12
               1  5  9 13
               2  6 10 14
               3  7 11 15])
    # lower left corner of domain, at
    # (x, y) = (0, 0)
    @test u[ 0,  0] ==  0
    @test u[ 2,  1] ==  6
    @test u[-1,  0] ==  3
    @test u[ 0, -1] == 12
    @test u[ 4,  0] ==  0
    @test u[ 0,  4] ==  0
    @test u[ 0,  5] ==  4
    @test u[ 5,  0] ==  1

    # linear indexing
    for i = 0:15
        @test u[i+1] == i
    end

    for i = 1:16
        u[i] = 2i
        @test u[i] == 2i
    end

    # setindex
    for i = -20:20, j = -20:20
        c = rand(1:1000)
        u[i, j] = c
        @test u[i, j] == c
    end
end

@testset "fieldsize                              " begin
    for n in 2:2:10
        @test fieldsize(Field(n)) == n
        @test fieldsize(Field{n}) == n
    end
end