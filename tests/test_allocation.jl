using Base.Test
using OpenKolmogorovFlow
using IMEXRKCB
using BenchmarkTools

# Test that no memory is allocated in calls
@testset "allocation                             " begin
    # example dimension
    n = 100

    # take a Re and a forcing wave number
    Re = 1.2345678
    kforcing = 4

    # initial condition
    Ω₀ = FTField(n)
    
    # get explicit and implicit parts
    L, N = KolmogorovFlowSystem(n, Re, kforcing)

    # forward integration should have not allocation. Expected elapsed time
    # is checked.
    for (scheme, elaps_exp) in [(IMEXRK3R2R(IMEXRKCB3c, false, Ω₀), 0.00155),
                                (IMEXRK3R2R(IMEXRKCB3e, false, Ω₀), 0.00155),
                                (IMEXRK4R3R(IMEXRKCB4,  false, Ω₀), 0.00250)]

        # count allocation and time
        @test (@belapsed  IMEXRKCB.step!($scheme, $N, $L, 0.0, 0.1, $Ω₀)) < elaps_exp
        @test (@allocated IMEXRKCB.step!( scheme,  N,  L, 0.0, 0.1,  Ω₀)) == 0

        # forward integration
        f = integrator(N, L, scheme, 0.01)
        # warm up
        f(Ω₀, 1)
        @test (@allocated f(Ω₀, 1)) == 0
    end
end