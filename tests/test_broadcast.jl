using Base.Test
using OpenKolmogorovFlow

@testset "broadcast                              " begin
    for n = [2, 4, 10]
        u = FTField(randn(n, n>>1+1) + im*randn(n, n>>1+1))
        v = FTField(randn(n, n>>1+1) + im*randn(n, n>>1+1))

        # return v and not the underlying array
        @test typeof(v .*= 1) <: FTField

        # use FTFields and scalar
        v .= u .+ conj.(u)
        for i in linearindices(v)
            @test imag(v[i]) == 0.0
        end

        v .= .- u .* u .* 2im
        for i in linearindices(v)
            @test v[i] == -u[i]*u[i]*2im
        end

        v .= 0.0 + 1im
        for i in linearindices(v)
            @test v[i] == 0.0 + 1im
        end

        v .-= 0.0 + 1im
        for i in linearindices(v)
            @test v[i] == 0.0 + 0.0*im
        end

        v .= 0.0
        v .+= u ./ u .+ u.^2
        for i in linearindices(v)
            @test v[i] == u[i]/u[i] + u[i]^2
        end
    end
end

@testset "broadcast allocation                   " begin
    n = 10
    u = FTField(randn(n, n>>1+1) + im*randn(n, n>>1+1))
    v = FTField(randn(n, n>>1+1) + im*randn(n, n>>1+1))
    c = 1.0

    # use FTFields and scalar
    foo(v, u, c) = (v .= u .+ conj.(u) .* c)

    foo(v, u, c)
    @test (@allocated foo(v, u, c)) == 0 
end