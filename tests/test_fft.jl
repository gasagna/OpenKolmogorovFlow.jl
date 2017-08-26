using OpenKolmogorovFlow
using Base.Test
using VariationalNumbers

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

@testset "allocating - forw/inv - interpolate    " begin
    n = 10
    x_coarse, y_coarse = make_grid(n)
    for k = -5:5
        u_coarse = Field(cos.(k.*x_coarse.+0.*y_coarse))
        U = FFT(u_coarse)
    
        # go backwards on a finer grid
        u_fine = IFFT(U, 2n)

        # this should match 
        x_fine, y_fine = make_grid(2*n)
        @test maximum(abs, u_fine.data - cos.(k.*x_fine.+0.*y_fine)) < 3e-15
    end
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

@testset "advanced interface                     " begin
    for n = [4, 6, 8, 64]
        # define cos(n/2y)
        U = FTField(n)   
        U[n>>1, 0] = 1

        # define fields for dealiased  calculations
        u_dealiased = Field(even_dealias_size(n))

        # define de-aliased plans
        ip_dealiased! = InverseFFT!(typeof(u_dealiased), U)

        # execute plan
        ip_dealiased!(u_dealiased, U)

        # square
        u_dealiased .*= u_dealiased

        # Fourier transform, then shrink to a n x n grid
        U² = shrinkto!(FTField(n), FFT(u_dealiased))

        # note the mean is correctly calculated, but clearly we
        # miss the frequency at wave number 4
        @test U²[0, 0] ≈ 0.5

        # test all others frequencies are zero
        U²[0, 0] = 0
        @test maximum(abs, U².data) < 1e-16
    end
end

# try squaring a lower frequency
@testset "lower frequency                        " begin
    n = 6
    U = FTField(n)   
    
    # sin(2y)^2 = 0.5*(1 - cos(4y))
    U[ 2, 0] =  0.5*im
    U[-2, 0] = -0.5*im

    # define fields for dealiased  calculations
    u_dealiased = Field(even_dealias_size(n))

    # define de-aliased plans
    ip_dealiased! = InverseFFT!(typeof(u_dealiased), U)

    # execute plan
    ip_dealiased!(u_dealiased, U)

    # square
    u_dealiased .*= u_dealiased

    # Fourier transform, then shrink to a n x n grid
    U² = shrinkto!(FTField(n), FFT(u_dealiased))

    # note the mean is correctly calculated, but clearly we
    # miss the frequency at wave number 4
    @test U²[0, 0] ≈ 0.5
    U²[0, 0] = 0
    @test maximum(abs, U².data) < 1e-16
end

@testset "product of two waves                   " begin
    n = 6
    U = FTField(n)   
    V = FTField(n)   
    
    # product of two waves
    A = 1+5im  
    B = 2-3im  
    U[ 1,  1] = A  
    V[ 2,  3] = B 

    # define fields for dealiased  calculations
    u_dealiased = Field(even_dealias_size(n))
    v_dealiased = Field(even_dealias_size(n))

    # define de-aliased plans
    ip_dealiased! = InverseFFT!(typeof(u_dealiased), U)
    fp_dealiased! = ForwardFFT!(typeof(U), u_dealiased)

    # execute plan
    ip_dealiased!(u_dealiased, U)
    ip_dealiased!(v_dealiased, V)

    # product
    u_dealiased .*= v_dealiased

    # Fourier transform, then shrink to a n x n grid
    UV = fp_dealiased!(FTField(n), u_dealiased)

    # sum of frequencies is out of bounds
    # UV[3, 4] == A*B

    # difference of frequencies
    @test UV[1, 2] == conj(A)*B

    # ensure there is nothing else
    UV[1, 2] = 0
    @test maximum(abs, UV.data) < 1e-14
end

@testset "ForwardFFT!                            " begin
    @testset "aliased calculations               " begin
        # define field
        U = FTField(4)
        u = Field(4)
        x, y = make_grid(4)
        u.data .= cos.(x .+ y) .+ sin.(x .+ y)

        # define aliased transform, same size
        f! = ForwardFFT!(FTField{4, Complex{Float64}}, u)

        # apply to same field size
        f!(U, u)

        # test value
        @test U[1, 1] ≈ 0.5*(1-im)

        # ensure there is nothing else
        U[1, 1] = 0
        @test maximum(abs, U.data) < 1e-14
    end
    @testset "dealiased calculations             " begin
        # define field
        U = FTField(4)
        u = Field(even_dealias_size(4))
        x, y = make_grid(even_dealias_size(4))
        u.data .= cos.(x .+ y) .+ sin.(x .+ y)

        # define dealiased transform, go to smaller size
        f! = ForwardFFT!(FTField{4, Complex{Float64}}, u)
        
        # apply to same field size
        f!(U, u)

        # test value
        @test U[1, 1] ≈ 0.5*(1-im)

        # ensure there is nothing else
        U[1, 1] = 0
        @test maximum(abs, U.data) < 1e-14
    end
end

@testset "with VariationalNumbers                " begin
    # test allocating versions, as these are based on the non-allocating ones
    # construct random data
    u_ = randn(4, 4) + δ*randn(4, 4)
    U_ = rfft(real.(u_), [2, 1])/16 + δ*rfft(pert.(u_), [2, 1])/16

    # create field and transform
    U = FTField(U_)
    u = IFFT(U)
    # FIXME: make Field one-based
    @test maximum(abs, pert.(u.data) - pert.(u_)) < 1e-15
    @test maximum(abs, real.(u.data) - real.(u_)) < 1e-15

    V = FFT(u)
    @test maximum(abs, real.(real.(V.data)) - real.(real.(U_))) < 1e-15
    @test maximum(abs, real.(imag.(V.data)) - real.(imag.(U_))) < 1e-15
    @test maximum(abs, pert.(real.(V.data)) - pert.(real.(U_))) < 1e-15
    @test maximum(abs, pert.(imag.(V.data)) - pert.(imag.(U_))) < 1e-15
end