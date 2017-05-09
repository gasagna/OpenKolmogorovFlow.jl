using Base.Test
using OpenKolmogorovFlow

@testset "broadcast    " begin
    for n = [2, 4, 10]
        u = FTField(randn(n, n>>1+1) + im*randn(n, n>>1+1))
        v = FTField(randn(n, n>>1+1) + im*randn(n, n>>1+1))

        # use FTFields and scalar
        v .= u .+ conj.(u)
        for i in eachindex(v)
            @test imag(v[i]) == 0.0
        end

        v .= .- u .* u .* 2im
        for i in eachindex(v)
            @test v[i] == -u[i]*u[i]*2im
        end

        v .= 0.0 + 1im
        for i in eachindex(v)
            @test v[i] == 0.0 + 1im
        end

        v .-= 0.0 + 1im
        for i in eachindex(v)
            @test v[i] == 0.0 + 0.0*im
        end

        v .= 0.0
        v .+= u ./ u .+ u.^2
        for i in eachindex(v)
            @test v[i] == u[i]/u[i] + u[i]^2
        end
    end
end