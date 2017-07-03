using OpenKolmogorovFlow
using Base.Test

@testset "dissipation                            " begin
    # test on laminar flow
    Re = rand()
    kf = 4
    Ω = laminarflow(10, Re, kf)
    @test DissipationRate(Ω, Re) ≈ Re/2/kf^2

    # test on some other flows
    N, d = 4, 2
    x = linspace(0, 2π, N+1)[1:end-1]'
    y = linspace(0, 2π, N+1)[1:end-1]

    # phase shift does not change dissipation
    Ωa = FT(Field(fun(x, y, 1,  1, 0)))
    Ωb = FT(Field(fun(x, y, im, 1, 0)))
    @test DissipationRate(Ωa, 1.0) == DissipationRate(Ωb, 1.0)
   
    # for any wave dissipation is 1/2/Re
    for j = -d+1:d, k=-d+1:d
        if !(j == 0 && k == 0)
            Re = randn()
            Ω = FT(Field(fun(x, y, 1, j, k)))
            @test DissipationRate(Ω, Re) ≈ 1/2/Re
        end
    end

    # test no allocations 
    @test (@allocated DissipationRate(Ω, Re)) == 16
    @test (@allocated DissipationRate(Ω, Re)) == 16
    @test (@allocated DissipationRate(Ω, Re)) == 16
    @test (@allocated DissipationRate(Ω, Re)) == 16
end