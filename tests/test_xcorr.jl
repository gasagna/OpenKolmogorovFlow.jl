using OpenKolmogorovFlow
using Base.Test

@testset "locatepeak                             " begin
    u = Field([0 1 2 3;
               4 5 6 7;
               8 9 0 1;
               2 3 4 5])
    @test OpenKolmogorovFlow.locatepeak(u) == (9, 1, 2)

    u = Field([0 1 2 3 0  1 2 3;
               4 5 6 7 4  5 6 7;
               8 9 0 1 8 10 0 1;
               8 9 0 1 8  9 0 1;
               8 9 0 1 8  9 0 1;
               8 9 0 1 8  9 0 1;
               8 9 0 1 8  9 0 1;
               2 3 4 5 2  3 4 5])
    @test OpenKolmogorovFlow.locatepeak(u) == (10, 5, 1)

    u = Field([0 1 2 3 0  1 2 3;
               4 5 6 7 4 10 6 7;
               8 9 0 1 8  5 0 1;
               8 1 0 1 8  3 0 1;
               8 3 0 1 8  5 0 1;
               8 5 0 1 8  2 0 1;
               8 8 0 1 8  1 0 1;
               2 3 4 5 2  3 4 5])
    @test OpenKolmogorovFlow.locatepeak(u) == (9, 1, 1)
end

@testset "distance                               " begin
    @testset "same field                         " begin
        n = 64
        x, y = make_grid(n)
        u = Field(cos.(x .+ y)) # peak must be at grid point,
        U, V = FFT(u), FFT(u)   # else we get sampling artefacts
        # same cache size
        cache = XCorrCache(n)
        @test distance!(U, V, cache) == (0.0, (0.0, 0))
        # reduced cache size
        cache = XCorrCache(48)
        @test distance!(U, V, cache) == (0.0, (0.0, 0))
    end
    @testset "known shift 1                      " begin
        n = 64
        x, y = make_grid(n)
        u = Field(cos.(x .+ 0.*y)) # peak must be at grid point,
        v = Field(sin.(x .+ 0.*y)) # else we get sampling artefacts
        U, V = FFT(u), FFT(v)
        # same cache size
        cache = XCorrCache(n)
        @test distance!(U, V, cache) == (0.0, (π/2, 0))
        # reduced cache size
        cache = XCorrCache(48)
        @test distance!(U, V, cache) == (0.0, (π/2, 0))
    end
    @testset "known shift 2                      " begin
        n = 64
        x, y = make_grid(n)
        u = Field(cos.(0.*x .+ y)) # peak must be at grid point,
        v = Field(sin.(0.*x .+ y)) # else we get sampling artefacts
        U, V = FFT(u), FFT(v)
        # same cache size
        cache = XCorrCache(n)
        @test distance!(U, V, cache) == (0.0, (0.0, 1))
        # reduced cache size
        cache = XCorrCache(48)
        @test distance!(U, V, cache) == (0.0, (0.0, 1))
    end
end