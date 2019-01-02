using OpenKolmogorovFlow
using IMEXRKCB
  

@testset "adjoint allocations                    " begin
    
    # setup
    const Δt = 0.015    # time step
    const n  = 64       # problem size
    const Re = 20.0     # reynolds number
    const kforcing = 4  # forcing wave number
    const T = 5         # time horizon
    const dealias = true

    impl = IMEXRK3R2R(IMEXRKCB3c, false)

    # ~~~ FORWARD PROBLEM ~~~
    # random initial conditions
    srand(0)

    # state
    Ω = FFT(Field(0.1*randn(n, n))); Ω[0, 0] = 0

    # define forward equation
    L_forw, N_forw = imex(ForwardEquation(n, Re, kforcing, ForwardMode(); dealias=false))

    # define forward integrator
    f_forw = integrator(N_forw, L_forw, IMEXRKScheme(impl, Ω), 0.01)

    # storage for forward solution only
    forw_sol = Monitor(Ω, Ω->copy(state(Ω)))

    # now store forward solution
    f_forw(Ω, (0, T), forw_sol)

    # ~~~ ADJOINT ~~~
    cost(dΛdt::FTField, Ω::FTField) = nothing

    # define backward equations
    L_back, N_back = imex(AdjointEquation(n, Re, forw_sol, cost, kforcing; dealias=true))

    # define backward integrator
    f_back = integrator(N_back, L_back, IMEXRKScheme(impl, Ω), 0.01)

    # march backwards the final state
    f_back(Ω, (T, 0))

    # wrap in a function
    foo(f_back, Ω, T) = @allocated f_back(Ω, (T, 0))
    @test foo(f_back, Ω, T) == 0
end