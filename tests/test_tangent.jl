using VariationalNumbers
using OpenKolmogorovFlow
using Base.Test

@testset "field creation and broadcast           " begin
    U = FTField(4, VarNum{Float64})
    u = Field(4, VarNum{Float64})
    @test eltype(U) == Complex{VarNum{Float64}}
    @test eltype(u) == VarNum{Float64}
end

@testset "operators work on real and pert parts  " begin
    # take example operator, the laplacian
    Δ = DiffOperator(4, :xxyy, Float64)
    
    # define a field with known values
    U = FTField(4, VarNum{Float64}) + 1 + 3*im + 2*δ + 4*δ*im
    
    # apply laplacian
    V = Δ.*U

    # test a few entries
    @test V[0, 0] == -(0^2+0^2)*((1+3*im) + (2+4*im)*δ)
    @test V[1, 1] == -(1^2+1^2)*((1+3*im) + (2+4*im)*δ)
    @test V[2, 0] == -(2^2+0^2)*((1+3*im) + (2+4*im)*δ)

    # take another example
    ∂ₓ = DiffOperator(4, :x, Float64)
    
    # define a field with known values
    U = FTField(4, VarNum{Float64}) + 1 + 3*im + 2*δ + 4*δ*im
    
    # apply laplacian
    V = ∂ₓ.*U

    # test a few entries
    @test V[0, 0] == 0*im*((1+3*im) + (2+4*im)*δ)
    @test V[0, 1] == 1*im*((1+3*im) + (2+4*im)*δ)
    @test V[0, 2] == 2*im*((1+3*im) + (2+4*im)*δ)
end

@testset "FTField perturbation                   " begin
    @testset "set perturbation                       " begin
        U = FTField(4, VarNum{Float64})
        p = FTField(randn(4, 3) + im*randn(4, 3))
        set_pert!(U, p)
        @test pert.(real.(U.data)) == real.(p.data)
        @test pert.(imag.(U.data)) == imag.(p.data)
        @test all(real.(real.(U.data)) .== 0)
        @test all(real.(imag.(U.data)) .== 0)
    end
    @testset "get perturbation                       " begin
        U_ = (randn(4, 3) + δ*randn(4, 3)) + im*(randn(4, 3) + δ*randn(4, 3))
        U = FTField(U_)
        p = FTField(4, Float64)
        get_pert!(U, p)
        @test real.(p.data) == pert.(real.(U_))
        @test imag.(p.data) == pert.(imag.(U_))
    end
end