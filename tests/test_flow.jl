using OpenKolmogorovFlow
using Base.Test

@testset "dissipation                            " begin
    # test on laminar flow
    Re = rand()
    kf = 4
    @test_throws ArgumentError laminarflow(10, Re, 6)
    Ω = laminarflow(10, Re, kf)
    @test dissrate(Ω, Re) ≈ Re/2/kf^2

    # test on some other flows
    n = 4
    x, y = make_grid(n)

    # phase shift does not change dissipation
    Ωa = FFT(Field(cos.(x.+y)))
    Ωb = FFT(Field(sin.(x.+y)))
    @test dissrate(Ωa, 1.0) == dissrate(Ωb, 1.0)
   
    # for any wave dissipation is 1/2/Re
    d = n>>1 + 1
    for j = -d+1:d, k=-d+1:d
        if !(j == 0 && k == 0)
            Re = randn()
            Ω = FFT(Field(cos.(j.*x.+k.*y)))
            @test dissrate(Ω, Re) ≈ 1/2/Re
        end
    end

    # test no allocations 
    @test (@allocated dissrate(Ω, 1.0)) == 16
    
    # warm up with int
    dissrate(Ω, 1)
    @test (@allocated dissrate(Ω, 1)) == 16
end