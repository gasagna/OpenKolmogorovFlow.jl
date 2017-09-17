using DualNumbers
using OpenKolmogorovFlow
using IMEXRKCB
using Base.Test

@testset "nonlinear simulation vs Var. Eqs.      " begin
    # seed random number generator
    srand(0)

    # take a simple setup
    n = 32
    Re = 20
    kforcing = 4

    for dealias in [false, true]
        # for each integration scheme
        for impl in [IMEXRK3R2R(IMEXRKCB3c, false),
                     IMEXRK3R2R(IMEXRKCB3e, false),
                     IMEXRK4R3R(IMEXRKCB4,  false)]

            # get explicit and implicit parts
            L, N = imex(VorticityEquation(n, Re, kforcing; dealias=dealias))

            # define integrator
            f = integrator(N, L, IMEXRKScheme(impl, FTField(n, Complex{Float64})), 0.005)

            # start from a random non zero initial conditions
            ω = Field(randn(n, n)); ω .-= mean(ω)
            Ω0_reference = FFT(ω)

            # Run for a while to land on attractor. Otherwise, we are just
            # looking at the viscous decay of the high frequency modes.
            # After this Ω0 will be the initial condition for the simulation
            # of the variational equations.
            f(Ω0_reference, 100)

            # Now integrate in time the nonlinear equations
            Ωf1 = f(copy(Ω0_reference), 10)

            # now perturb Ω0 along a particular direction
            Ω0_reference[1, 1] += 1e-6

            # and integrate in time, for the same time span as before
            Ωf2 = f(copy(Ω0_reference), 10)

            # reset perturbation, so we can reuse the initial condition later
            Ω0_reference[1, 1] -= 1e-6

            # now take the difference between the two fields from the
            # nonlinear simulations and use that as reference. This
            # difference represents the evolution of the perturbation over
            # the reference trajectory.
            diff_nl = (Ωf2-Ωf1)/1e-6

            # Now define integration setup for variational equations
            L, N = imex(VorticityEquation(n, Re, kforcing; mode=TangentMode(), dealias=dealias))
            f = integrator(N, L, IMEXRKScheme(impl, AugmentedFTField(n)), 0.005)

            # define initial condition, using the previous one
            Ω0 = AugmentedFTField(Ω0_reference, FTField(n))

            # perturb mode (1, 1), by the same amount as in the nonlinear
            # simulation. The perturbation size does not matter, because
            # the equations are linearised. We do this because is is easier
            # to check the two. Note that here, we do not perturb the initial
            # condition of the nonlinear simulation. We set the initial condition
            # of the variational equations
            Ω0[1, 1] += 0 + 1e-6*ε

            # integrate in time
            Ωf1_var = f(copy(Ω0), 10)

            # check that amplitude of mode (1, 1) of perturbation is similar
            # to the difference of two perturbed nonlinear simulations. We cannot
            # set too low a tolerance here, because nonlinear effects will kick in
            # at some point reducing the validity of the linearised equations.
            # this checks that we are basically not doing rubbish.
            a = diff_nl[1, 1]
            b = Ωf1_var[1, 1]
            @test abs(a-epsilon(b)) < 1e-6
        end
    end
end