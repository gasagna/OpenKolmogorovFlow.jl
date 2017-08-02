using Base.Test
using OpenKolmogorovFlow

@testset "ItoK                                   " begin
    # these are the wave numbers for k
    data_k = [0  0  0
              1  1  1
              2  2  2
             -1 -1 -1]
    for i in 1:12
        @test data_k[i] == OpenKolmogorovFlow.ItoK(i, 4)
    end
end

@testset "ItoJ                                   " begin
    # these are the wave numbers for j
    data_j = [0  1  2
              0  1  2
              0  1  2
              0  1  2]
    for i in 1:12
        @test data_j[i] == OpenKolmogorovFlow.ItoJ(i, 4)
    end
end

@testset "KJtoI                                  " begin
    @test OpenKolmogorovFlow.KJtoI( 0, -2, 4) ==  9
    @test OpenKolmogorovFlow.KJtoI( 1, -2, 4) == 12
    @test OpenKolmogorovFlow.KJtoI( 2, -2, 4) == 11
    @test OpenKolmogorovFlow.KJtoI(-2, -2, 4) == 11
    @test OpenKolmogorovFlow.KJtoI(-1, -2, 4) == 10

    @test OpenKolmogorovFlow.KJtoI( 0, -1, 4) ==  5
    @test OpenKolmogorovFlow.KJtoI( 1, -1, 4) ==  8
    @test OpenKolmogorovFlow.KJtoI( 2, -1, 4) ==  7
    @test OpenKolmogorovFlow.KJtoI(-2, -1, 4) ==  7
    @test OpenKolmogorovFlow.KJtoI(-1, -1, 4) ==  6

    @test OpenKolmogorovFlow.KJtoI( 0,  0, 4) ==  1
    @test OpenKolmogorovFlow.KJtoI( 1,  0, 4) ==  2
    @test OpenKolmogorovFlow.KJtoI( 2,  0, 4) ==  3
    @test OpenKolmogorovFlow.KJtoI(-2,  0, 4) ==  3
    @test OpenKolmogorovFlow.KJtoI(-1,  0, 4) ==  4

    @test OpenKolmogorovFlow.KJtoI( 0,  1, 4) ==  5
    @test OpenKolmogorovFlow.KJtoI( 1,  1, 4) ==  6
    @test OpenKolmogorovFlow.KJtoI( 2,  1, 4) ==  7
    @test OpenKolmogorovFlow.KJtoI(-2,  1, 4) ==  7
    @test OpenKolmogorovFlow.KJtoI(-1,  1, 4) ==  8

    @test OpenKolmogorovFlow.KJtoI( 0,  2, 4) ==  9
    @test OpenKolmogorovFlow.KJtoI( 1,  2, 4) == 10
    @test OpenKolmogorovFlow.KJtoI( 2,  2, 4) == 11
    @test OpenKolmogorovFlow.KJtoI(-2,  2, 4) == 11
    @test OpenKolmogorovFlow.KJtoI(-1,  2, 4) == 12
end

@testset "rectify                                " begin
    @test OpenKolmogorovFlow.rectify(1,   1)     == 1
    @test OpenKolmogorovFlow.rectify(1,   0)     == 1
    @test OpenKolmogorovFlow.rectify(1,  -1)     == 1
    @test OpenKolmogorovFlow.rectify(1 + im,  1) == 1 + im
    @test OpenKolmogorovFlow.rectify(1 + im,  0) == 1 + im
    @test OpenKolmogorovFlow.rectify(1 + im, -1) == 1 - im
end