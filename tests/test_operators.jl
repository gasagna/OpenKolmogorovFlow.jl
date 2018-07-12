using Base.Test
using OpenKolmogorovFlow

@testset "derivatives                            " begin

    # create data
    data = Complex{Float64}[1+2im  9+10im 0+0im
                            3+4im 11+12im 0+0im
                            0+0im  0+0im  0+0im
                            3-4im 15+16im 0+0im]
    U = FTField(data, 1)

    # do work
    out = ddx!(similar(U), U)
    @test out == [ 0+0im  9im-10 0+0im
                   0+0im 11im-12 0+0im
                   0+0im  0im-0  0+0im
                   0+0im 15im-16 0+0im]

    out = ddy!(similar(U), U)
    @test out == [  0im-0   0im-0   0im+0
                    3im-4  11im-12  0im+0
                    0im-0   0im-0   0im+0
                   -3im-4 -15im+16  0im+0]

    out = invlaplacian!(similar(U), U)
    @test out == [  0+0im   -9-10im 0im+0
                   -3-4im -5.5-6im  0im+0
                    0+0im    0+0im  0im+0
                   -3+4im -7.5-8im  0im+0]
    
    out = laplacian!(similar(U), U)
    @test out == [  0+0im   -9-10im  0im+0
                   -3-4im  -22-24im  0im+0
                    0+0im     0+0im  0im+0
                   -3+4im  -30-32im  0im+0]
end
