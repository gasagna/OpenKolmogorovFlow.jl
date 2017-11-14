using OpenKolmogorovFlow
using IMEXRKCB
using Base.Test

# We use the adjoint code to verify the identity 
#
# (Λ_T, L Ω′_0) = (L⁺Λ_T, Ω′_0)
# 
# where the subscripts _T and _0 denote the final and initial 
# states of the adjoint variables Λ and perturbation Ω′. The 
# operator L is the linearisation of the forward map, while L⁺ 
# is its adjoint, from integration by parts in space and time.
# The adjoint operator L⁺ runs backward in time. Here we check
# that the identity is verified up to a certain precision. 
# The identity cannot be satisfied exactly, because we use 
# interpolation of the state variables to rune the adjoint, 
# whereas the variational equations are run along with the 
# forward equations.
@testset "adjoint                           " begin
    
    # setup
    const Δt = 0.015    # time step
    const n  = 64       # problem size
    const Re = 20.0     # reynolds number
    const kforcing = 4  # forcing wave number
    const T = 5         # time horizon

    for dealias in [false, true]
        for impl in [IMEXRK3R2R(IMEXRKCB3c, false),
                     IMEXRK3R2R(IMEXRKCB3e, false),
                     IMEXRK4R3R(IMEXRKCB4,  false)]

            # ~~~ FORWARD PROBLEM ~~~
            # random initial conditions
            srand(0)

            # state
            Ωs = FFT(Field(0.1*randn(n, n))); Ωs[0, 0] = 0

            # perturbation in a single direction
            Ωp = FTField(n); Ωp[1, 2] = 1.0

            # augmented state
            Ω = VariationalFTField(Ωs, copy(Ωp));

            # define forward and variational equation
            L_forw, N_forw = imex(ForwardEquation(n, Re, kforcing, TangentMode(); dealias=false))

            # define forward integrator
            f_forw = integrator(N_forw, L_forw, IMEXRKScheme(impl, Ω), 0.01)

            # storage for forward solution only
            forw_sol = Monitor(Ω, Ω->copy(state(Ω)));

            # now store forward solution
            f_forw(Ω, (0, T), forw_sol);

            # calculate 
            LHS = inner(state(Ω), prime(Ω)) 

            # ~~~ ADJOINT ~~~
            cost(dΛdt::FTField, Ω::FTField) = nothing

            # define backward equations
            L_back, N_back = imex(AdjointEquation(n, Re, forw_sol, cost, kforcing; dealias=true))

            # define backward integrator
            f_back = integrator(N_back, L_back, IMEXRKScheme(impl, Ωs), 0.01);

            # monitor the adjoint state
            mon = Monitor(Ωs, copy);

            # march backwards the final state
            f_back(Ωs, (T, 0), mon);

            # calculate 
            RHS = inner(Ωs, Ωp)

            Δ = abs(LHS - RHS)/abs(RHS)
            @test Δ < 1e-6
        end
    end
end