using VariationalNumbers
using OpenKolmogorovFlow
using Base.Test

@testset "field creation and broadcast           " begin
    U = FTField(4, VarNum{Float64})
    u = Field(4, VarNum{Float64})
    @test eltype(U) == Complex{VarNum{Float64}}
    @test eltype(u) == VarNum{Float64}
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