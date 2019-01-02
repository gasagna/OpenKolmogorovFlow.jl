# Test that no memory is allocated in calls
@testset "allocation                             " begin
    # example dimension
    m = 49
    n = down_dealias_size(m)

    # take a Re and a forcing wave number
    Re = 1.2345678
    kforcing = 4

    # initial condition
    Ω₀ = FTField(n, m)

    # get explicit and implicit parts
    N, L = splitexim(ForwardEquation(n, m, Re, kforcing))

    for (scheme, elaps_exp) in [(Scheme(:CB3e_3R2R, Ω₀), 0.13),
                                (Scheme(:CB3c_3R2R, Ω₀), 0.13),
                                (Scheme(:CB4_4R3R,  Ω₀), 0.22)]
        # get integrator
        f = integrator(N, L, scheme, TimeStepConstant(0.01))

        # test allocations
        min_time = minimum([@elapsed f(Ω₀, (0, 1)) for i = 1:10])
        @test min_time < elaps_exp
        @test (@allocated f(Ω₀, (0, 1))) == 0   
    end
end