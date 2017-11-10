using Base.Test
using OpenKolmogorovFlow
using IMEXRKCB

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
    L, N = imex(ForwardEquation(n, Re, kforcing; dealias=false))

    for (scheme, elaps_exp) in [(IMEXRKScheme(IMEXRK3R2R(IMEXRKCB3c, false), Ω₀), 0.150),
                                (IMEXRKScheme(IMEXRK3R2R(IMEXRKCB3e, false), Ω₀), 0.155),
                                (IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4,  false), Ω₀), 0.320)]
        # get integrator
        f = integrator(N, L, scheme, 0.01)

        # test elapsed time and allocation
        @test minimum([@elapsed f(Ω₀, (0, 1)) for i = 1:10]) < elaps_exp
        @test (@allocated f(Ω₀, (0, 1))) == 0
    end
end