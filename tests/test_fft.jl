using OpenKolmogorovFlow
using Base.Test

macro display(ex)
    quote
        display($(esc(ex)))
        println()
    end
end

@testset "allocating - some tests                " begin 
    # make grid
    n = 4
    x, y = make_grid(n)
    
    # constant field
    U = FTField(n)
    U[0, 0] = 1
    u = IFFT(U)
    @test all(u.data .== 1.0)

    # sin(y)
    U = FTField(n)
    U[ 1, 0] =      -0.5im
    U[-1, 0] = conj(-0.5im)
    u = IFFT(U)
    @test maximum(abs, u.data - sin.(0.*x .+ 1.*y)) < 1e-15

    # cos(2y)
    U = FTField(n)
    U[2, 0] = 1
    u = IFFT(U)
    @test maximum(abs, u.data - cos.(0.*x .+ 2.*y)) < 1e-15

    # cos(2x)
    U = FTField(n)
    U[0, 2] = 1
    u = IFFT(U)
    @test maximum(abs, u.data - cos.(2.*x .+ 0.*y)) < 1e-15

    # cos(2x + 2y)
    U = FTField(n)
    U[2, 2] = 1
    u = IFFT(U)
    @test maximum(abs, u.data - cos.(2.*x .+ 2.*y)) < 1e-15

    # cos(x + y)
    U = FTField(n)
    U[1, 1] = 0.5
    u = IFFT(U)
    @test maximum(abs, u.data - cos.(1.*x .+ 1.*y)) < 1e-15
end

@testset "allocating - forw/inv                  " begin
    n = 2
    u = Field(randn(n, n))
    # FIXME norm
    @test maximum(abs, IFFT(FFT(u)).data - u.data) < 1e-15

    U = FFT(Field(randn(n, n)))
    @test maximum(abs, FFT(IFFT(U)).data - U.data) < 1e-15
end

@testset "utils                                  " begin
    # the reasoning is as follows. on a n × n grid the largest wave is
    # n/2. When you square this wave you get a wave at a frequency n.
    # Now, you need a grid large enough such that this wave does not
    # alias to waves equal to n/2 or lower. If it aliases to n/2 + 1,
    # for instance, you do not care, because you care about having alias
    # free in the original n × n grid. So you find an a grid of size m, 
    # such that n - m/2 > n/2, because on such a grid, the wave a frequency
    # n would alias to n - m/2. It is guaranteed that n/m > 1.5, the exact
    # values just depend on picking even grid sizes that are convenient
    # for speed.
    @test even_dealias_size(4) == 8
    @test even_dealias_size(6) == 10
    @test even_dealias_size(8) == 14
    @test even_dealias_size(10) == 16

    for n = 12:2:100
        @test even_dealias_size(n)/n > 1.5
    end
end

@testset "advanced interface                     " begin
    # define cos(2y)
    U = FTField(4)   
    U[2, 0] = 1

    # define fields for dealiased and aliased calculations
    u_aliased   = Field(4)
    u_dealiased = Field(even_dealias_size(4))

    # define de-aliased plans
    ip_dealiased! = InverseFFT!(typeof(u_dealiased), U)

    # these do not destroy input at construction
    @test U[2, 0] == 1

    # but the aliased might
    ip_aliased! = InverseFFT!(typeof(u_aliased), U)

    # now execute plans and go to physical space
    ip_aliased!(u_aliased, U)

    # U will be different, as inverse transforms
    # modify input
    @test U[1, 0] != 0

    # now square. It should be cos(2y)^2 = 0.5*(1 + cos(4y))
    u_aliased .*= u_aliased

    # and transform 
    U_aliased = FFT(u_aliased)

    # of course we cannot hold the 4y wave and this is aliased to the zero frequency
    @test U_aliased[0, 0] != 0.5

    # try with dealiased calculations
    U .= 0
    U[2, 0] = 1
    ip_dealiased!(u_dealiased, U)
    u_dealiased .*= u_dealiased


    # Fourier transform, then shrink to a 4 x 4 grid
    U_dealiased = FFT(u_dealiased)
    U_shrink = shrinkto!(FTField(4), U_dealiased)
    # note the mean is correctly calculated, but clearly we
    # miss the frequency at wave number 4
    @test U_shrink.data == [0.5+0.0im  0.0+0.0im  0.0+0.0im
                            0.0+0.0im  0.0+0.0im  0.0+0.0im
                            0.0+0.0im  0.0+0.0im  0.0+0.0im
                            0.0+0.0im  0.0+0.0im  0.0+0.0im] 
end