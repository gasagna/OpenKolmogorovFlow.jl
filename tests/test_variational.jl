using OpenKolmogorovFlow
using Base.Test
using DualNumbers

@testset "operators                              " begin
    n = 10
    x, y = make_grid(n)
    @testset "laplacian                          " begin
        ux, up = Field(cos.(2x .+ 3y)), Field(sin.(3x .+ 4y))
        Ux, Up = FFT(ux), FFT(up); U = VariationalFTField(Ux, Up); V = similar(U)
        Δ = OpenKolmogorovFlow.ImplicitTerm(n, 1, Float64)

        V  .= Δ .* U
        @test IFFT(state(V)).data ≈ .-4*cos.(2x .+ 3y) .-  9*cos.(2x .+ 3y)
        @test IFFT(prime(V)).data ≈ .-9*sin.(3x .+ 4y) .- 16*sin.(3x .+ 4y)
    end
    @testset "dx/dy                                 " begin
        ux, up = Field(cos.(2x .+ 3y) + sin.(x .+ y)), Field(sin.(3x .+ 4y) .+ cos.(3x .+ 4y))
        Ux, Up = FFT(ux), FFT(up); U = VariationalFTField(Ux, Up); V = similar(U)
        ∂x = OpenKolmogorovFlow.DiffOperator(n, :x)
        ∂y = OpenKolmogorovFlow.DiffOperator(n, :y)

        V  .= ∂x .* U
        @test IFFT(state(V)).data ≈ .-2*sin.(2x .+ 3y) .+   cos.(x  .+ y)
        @test IFFT(prime(V)).data ≈   3*cos.(3x .+ 4y) .- 3*sin.(3x .+ 4y)
        V  .= ∂y .* U
        @test IFFT(state(V)).data ≈ .-3*sin.(2x .+ 3y) .+   cos.(x  .+ y)
        @test IFFT(prime(V)).data ≈   4*cos.(3x .+ 4y) .- 4*sin.(3x .+ 4y)
    end
end

@testset "eltypes                                " begin
    @test eltype(VariationalFTField(4))        == Dual{Complex{Float64}}
    @test eltype(FTField(4))                 == Complex{Float64}

    @test eltype(VariationalField(4))          == Dual{Float64}
    @test eltype(Field(4))                   == Float64

    @test eltype(VariationalFTField(4, Dual{Complex{Int64}})) == Dual{Complex{Int64}}
    @test eltype(FTField(4, Complex{Int64}))                == Complex{Int64}

    @test eltype(VariationalField(4, Dual{Int64}))   == Dual{Int64}
    @test eltype(Field(4, Int64))                  == Int64
end 

@testset "broadcast                              " begin
    n = 10
    x, y = make_grid(n)
    ux, up = Field(cos.(2x .+ 3y)), Field(sin.(3x .+ 4y))
    Ux, Up = FFT(ux), FFT(up); U = VariationalFTField(Ux, Up); 
    V = VariationalFTField(n)

    # set to zero
    V .= zero(eltype(V))

    # assert both parts have been zeroed
    @test norm(V) == 0 + 0*ε
end

@testset "norms/inners                           " begin
    n = 10
    x, y = make_grid(n)

    # state and prime are orthogonal, so result is real
    ux, up = Field(cos.(2x .+ 3y)), Field(sin.(3x .+ 4y))
    U = VariationalFTField(FFT(ux), FFT(up)) 
    @test dot(U, U) == 2*π^2 + 0*ε

    # state and prime are not orthogonal, so result contains a perturbation
    ux, up = Field(cos.(2x .+ 3y)), Field(cos.(2x .+ 3y))
    U = VariationalFTField(FFT(ux), FFT(up)) 
    @test dot(U, U) == 2*π^2 + 4*π^2*ε
end

@testset "transforms                             " begin
    n = 10
    x, y = make_grid(n)

    # define field
    u = VariationalField(Field(cos.(2x .+ 3y)), Field(sin.(3x .+ 4y)))

    # transform
    U = FFT(u)

    d = abs(U[3, 2] - (0.5 + ε*0.0) + im*(0.0 + ε*0.0))
    @test max(value(d), epsilon(d)) < 3e-16
    d = abs(U[4, 3] - (0.0 + ε*0.0) + im*(0.0 - ε*0.5))
    @test max(value(d), epsilon(d)) < 3e-16
end