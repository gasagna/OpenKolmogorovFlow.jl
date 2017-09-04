using OpenKolmogorovFlow
using IMEXRKCB
using Base.Test


@testset "locatepeak                             " begin
    u = Field([0 1 2 3;
               4 5 6 7;
               8 9 0 1;
               2 3 4 5])
    @test OpenKolmogorovFlow.locatepeak(u) == (9, 1, 4)

    u = Field([0 1 2 3 0  1 2 3;
               4 5 6 7 4  5 6 7;
               8 9 0 1 8 10 0 1;
               8 9 0 1 8  9 0 1;
               8 9 0 1 8  9 0 1;
               8 9 0 1 8  9 0 1;
               8 9 0 1 8  9 0 1;
               2 3 4 5 2  3 4 5])
    @test OpenKolmogorovFlow.locatepeak(u) == (10, 5, 2)

    u = Field([0 1 2 3 0  1 2 3;
               4 5 6 7 4 10 6 7;
               8 9 0 1 8  5 0 1;
               8 1 0 1 8  3 0 1;
               8 3 0 1 8  5 0 1;
               8 5 0 1 8  2 0 1;
               8 8 0 1 8  1 0 1;
               2 3 4 5 2  3 4 5])
    @test OpenKolmogorovFlow.locatepeak(u) == (9, 1, 2)
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
        @test distance!(U, V, cache) == (0.0, (0.0, 2))
        # reduced cache size
        cache = XCorrCache(48)
        @test distance!(U, V, cache) == (0.0, (0.0, 2))
    end
end


@testset "distance on fields                     " begin
    # setup
    Re, n, Δt = 40, 64, 0.015

    # initial condition
    Ω = laminarflow(n, Re)
    for j = 1:5, k=1:5
        Ω[k, j] = 0.1*(randn() + im*randn())
    end

    # integration schemes
    RK = IMEXRKScheme(IMEXRK3R2R(IMEXRKCB3e, false), Ω)

    # Get system
    L, N = imex(VorticityEquation(n, Re; dealias=true))

    # forward map
    f  = integrator(N, L, RK,  Δt)

    # run forward to get to steady state
    f(Ω, 100)

    # distance cache
    cache = XCorrCache(64)

    for j = 0:70, m = 0:4
        # shift exactly by a multiple of the grid size
        Δ = (j*2π/64, 2*m)
        Ωs = shift!(deepcopy(Ω), Δ)

        # calculate distance
        d, (s, m_) = distance!(Ω, Ωs, cache)
        @test d < 1e-10
        @test abs(s  - Δ[1] % 2π) < 1e-8
        @test m_ == 2*m % 8
    end
end