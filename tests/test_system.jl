using OpenKolmogorovFlow
using IMEXRKCB
using Base.Test

# Test that solution converges to the
# laminar flow for small Reynolds numbers
@testset "laminar flow                           " begin
    # example dimension
    n = 10

    # take a Re and a forcing wave number
    Re = 1.2345678
    kforcing = 4

    # initial condition
    Ω₀ = FTField(n)

    for dealias in [true, false]
        # get explicit and implicit parts
        L, N = imex(ForwardEquation(n, Re, kforcing; dealias=dealias))

        # for each integration scheme
        for impl in [IMEXRK3R2R(IMEXRKCB3c, false),
                     IMEXRK3R2R(IMEXRKCB3e, false),
                     IMEXRK4R3R(IMEXRKCB4,  false)]

            # define integrator
            f = integrator(N, L, IMEXRKScheme(impl, Ω₀), 0.005)

            # start from some non zero initial condition
            Ω₀ .= 0.01; Ω₀[0, 0] = 0

            # monitor the state excited by forcing
            m = Monitor(Ω₀, Ω->Ω)

            # map forward
            f(Ω₀, (0, 50),  m)

            # test final value is that predicted by explicit equation
            Δ = m.xs[end] .- laminarflow(n, Re, kforcing)
            @test maximum(abs, Δ) < 1e-15
        end
    end
end