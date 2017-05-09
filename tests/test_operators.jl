using Base.Test
using OpenKolmogorovFlow

@testset "derivatives  " begin
    
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

@testset "analytic     " begin

    # make grid
    n = 10
    x = linspace(0, 2π, n+1)[1:n]'
    y = linspace(0, 2π, n+1)[1:n]
    
    # a complex exponential and its derivatives
    f(x, y, c::Number, α::Integer, β::Integer) = 
        0.5.*real(c.*exp.(im*(α.*x .+ β.*y)) .+ conj(c).*exp.(-im*(α.*x .+ β.*y)))

    fx(x, y, c::Number, α::Integer, β::Integer) = 
        0.5.*real(im.*α.*c.*exp.(im*(α.*x .+ β.*y)) .- im.*α.*conj(c).*exp.(-im*(α.*x .+ β.*y)))

    fy(x, y, c::Number, α::Integer, β::Integer) =     
        0.5.*real(im.*β.*c.*exp.(im*(α.*x .+ β.*y)) .- im.*β.*conj(c).*exp.(-im*(α.*x .+ β.*y)))        

    fxx(x, y, c::Number, α::Integer, β::Integer) = -α^2*f(x, y, c, α, β)
    fyy(x, y, c::Number, α::Integer, β::Integer) = -β^2*f(x, y, c, α, β)

    # select α, β up to n/2-1, to avoid aliasing
    l = n>>1-1
    for α in rand(-l:l, 20), β in rand(-l:l, 20)
        # construct a random field
        c = rand() + im*rand()
        u = Field(f(x, y, c, α, β))

        # transform to Fourier space
        U = FT(u)

        # calculate derivatives. This also tests broadcast for 
        # operators
        Ux  = DiffOperator(n, :x)  .* U
        Uy  = DiffOperator(n, :y)  .* U
        Uxx = DiffOperator(n, :xx) .* U
        Uyy = DiffOperator(n, :yy) .* U

        # check against analytic value
        @test IFT(Ux).data  ≈  fx(x, y, c, α, β)
        @test IFT(Uy).data  ≈  fy(x, y, c, α, β)
        @test IFT(Uxx).data ≈ fxx(x, y, c, α, β)
        @test IFT(Uyy).data ≈ fyy(x, y, c, α, β)
    end
end