using Base.Test
using OpenKolmogorovFlow

@testset "derivatives                            " begin
    
    # create data
    data = [1+2im  9+10im 17+0im
            3+4im 11+12im 19+20im
            5+0im 13+14im 21+0im
            3-4im 15+16im 19-20im]
    U = FTField(data)

    # do work
    ∂U∂x   = DiffOperator(4, :x)  .* U
    ∂U∂y   = DiffOperator(4, :y)  .* U
    ∂²U∂x² = DiffOperator(4, :xx) .* U
    ∂²U∂y² = DiffOperator(4, :yy) .* U

    # Note that these field do not represent valid
    # rfft data, because some of the symmetries have
    # been lost, i.e. some terms that should be zero 
    # are not zero. These are neglected by FFTW.            
    @test ∂U∂x.data   == [  0+0im   9im-10  34im-0
                            0+0im  11im-12  38im-40
                            0+0im  13im-14  42im-0
                            0+0im  15im-16  38im+40]
    @test ∂U∂y.data   == [  0im-0   0im-0    0im+0
                            3im-4  11im-12  19im-20
                           10im-0  26im-28  42im+0
                           -3im-4 -15im+16 -19im-20]

    @test ∂²U∂x².data == [  0+0im  -9-10im -68-0im
                            0+0im -11-12im -76-80im
                            0+0im -13-14im -84-0im
                            0+0im -15-16im -76+80im]

    @test ∂²U∂y².data == [  0+0im   0+00im   0+0im
                           -3-4im -11-12im -19-20im
                          -20+0im -52-56im -84+0im
                           -3+4im -15-16im -19+20im]
end

@testset "analytic                               " begin

    # make grid
    n = 10
    x = linspace(0, 2π, n+1)[1:n]'
    y = linspace(0, 2π, n+1)[1:n]
    
    # select α, β up to n/2-1, to avoid aliasing
    l = n>>1-1

    # make grid
    x, y = make_grid(n)
    
    for α in rand(-l:l, 10), β in rand(-l:l, 10)
        # construct a random field
        cc, cs = randn(), randn()

        u = Field(cc.*cos.(α.*x .+ β.*y) .+ cs.*sin.(α.*x .+ β.*y))

        # transform to Fourier space
        U = FFT(u)

        # calculate derivatives. This also tests broadcast for 
        # operators
        Ux  = DiffOperator(n, :x)  .* U
        Uy  = DiffOperator(n, :y)  .* U
        Uxx = DiffOperator(n, :xx) .* U
        Uyy = DiffOperator(n, :yy) .* U

        # check against analytic value
        @test IFFT(Ux).data  ≈    -α.*cc.*sin.(α.*x .+ β.*y) .+ α.*cs.*cos.(α.*x .+ β.*y)
        @test IFFT(Uy).data  ≈    -β.*cc.*sin.(α.*x .+ β.*y) .+ β.*cs.*cos.(α.*x .+ β.*y)
        @test IFFT(Uxx).data ≈  -α^2.*cc.*cos.(α.*x .+ β.*y) .- α^2.*cs.*sin.(α.*x .+ β.*y)
        @test IFFT(Uyy).data ≈  -β^2.*cc.*cos.(α.*x .+ β.*y) .- β^2.*cs.*sin.(α.*x .+ β.*y)
    end
end