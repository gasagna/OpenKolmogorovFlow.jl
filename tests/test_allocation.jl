# Test that no memory is allocated in calls
@testset "allocation                             " begin
    # example dimension
    m = 49
    n = down_dealias_size(m)

    # take a Re and a forcing wave number
    Re = 1.2345678
    kforcing = 4

    # initial condition
    Ω = FTField(n, m)

    # get explicit and implicit parts
    f = ForwardEquation(n, m, Re, kforcing)

    for (method, elaps_exp) in [(CB3R2R3e(Ω, :NORMAL), 0.115),
                                (CB3R2R3c(Ω, :NORMAL), 0.115)]

        # define integrator
        ϕ = flow(splitexim(f)..., method, TimeStepConstant(0.01))

        # test allocations
        min_time = minimum([@elapsed ϕ(Ω, (0, 1)) for i = 1:10])
        @test min_time < elaps_exp
        @test (@allocated ϕ(Ω, (0, 1))) == 0   
    end
end