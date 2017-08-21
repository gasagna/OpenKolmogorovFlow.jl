using VariationalNumbers
using OpenKolmogorovFlow
using IMEXRKCB
using Base.Test

@testset "field creation and broadcast           " begin
    U = FTField(4, Complex{VarNum{Float64}})
    u = Field(4, VarNum{Float64})
    @test eltype(U) == Complex{VarNum{Float64}}
    @test eltype(u) == VarNum{Float64}
end

@testset "operators work on real and pert parts  " begin
    # take example operator, the laplacian
    Δ = DiffOperator(4, :xxyy, Float64)
    
    # define a field with known values
    U = FTField(4, Complex{VarNum{Float64}}) + 1 + 3*im + 2*δ + 4*δ*im
    
    # apply laplacian
    V = Δ.*U

    # test a few entries
    @test V[0, 0] == -(0^2+0^2)*((1+3*im) + (2+4*im)*δ)
    @test V[1, 1] == -(1^2+1^2)*((1+3*im) + (2+4*im)*δ)
    @test V[2, 0] == -(2^2+0^2)*((1+3*im) + (2+4*im)*δ)

    # take another example
    ∂ₓ = DiffOperator(4, :x, Float64)
    
    # define a field with known values
    U = FTField(4, Complex{VarNum{Float64}}) + 1 + 3*im + 2*δ + 4*δ*im
    
    # apply laplacian
    V = ∂ₓ.*U

    # test a few entries
    @test V[0, 0] == 0*im*((1+3*im) + (2+4*im)*δ)
    @test V[0, 1] == 1*im*((1+3*im) + (2+4*im)*δ)
    @test V[0, 2] == 2*im*((1+3*im) + (2+4*im)*δ)
end

@testset "FTField perturbation                   " begin
    @testset "set perturbation                       " begin
        U = FTField(4, Complex{VarNum{Float64}})
        p = FTField(randn(4, 3) + im*randn(4, 3))
        set_pert!(U, p)
        @test pert.(real.(U.data)) == real.(p.data)
        @test pert.(imag.(U.data)) == imag.(p.data)
        @test all(real.(real.(U.data)) .== 0)
        @test all(real.(imag.(U.data)) .== 0)
    end
    @testset "get perturbation                       " begin
        U_ = (randn(4, 3) + δ*randn(4, 3)) + im*(randn(4, 3) + δ*randn(4, 3))
        U = FTField(U_)
        p = FTField(4, Complex{Float64})
        get_pert!(U, p)
        @test real.(p.data) == pert.(real.(U_))
        @test imag.(p.data) == pert.(imag.(U_))
    end
end

@testset "DNS is not affected by Var. Eqs.       " begin
    # seed random number generator
    srand(0)

    # take a simple setup
    n = 64
    Re = 1
    kforcing = 4
    
    for dealias in [false, true]
        # for each integration scheme
        for impl in [IMEXRK3R2R(IMEXRKCB3c, false),
                     IMEXRK3R2R(IMEXRKCB3e, false),
                     IMEXRK4R3R(IMEXRKCB4,  false)]

            # get explicit and implicit parts
            L, N = imex(VorticityEquation(n, Re, kforcing; T=Float64, dealias=dealias))

            # define integrator
            f = integrator(N, L, IMEXRKScheme(impl, FTField(n, Complex{Float64})), 0.005)

            # start from a random non zero initial conditions
            ω = Field(randn(n, n)); ω .-= mean(ω) 
            Ω0 = FFT(ω)

            # integrate in time
            Ωf1 = f(copy(Ω0), 1)

            # Now define integration for variational equations
            L, N = imex(VorticityEquation(n, Re, kforcing; T=VarNum{Float64}, dealias=dealias))
            f = integrator(N, L, IMEXRKScheme(impl, FTField(n, Complex{VarNum{Float64}})), 0.005)

            # define initial condition
            Ω0 = FFT(ω) + 0*δ

            # integrate in time
            Ωf1_var = f(copy(Ω0), 1)

            # the real part, (not the perturbation), should be the same, 
            # regardless of the fact that we are integrating the variational 
            # equations as well difference arise in the differences in FFTW
            # algorithms for data that is not layed out in memory in the same way
            a = real.(Ωf1.data)
            b = real.(real.(Ωf1_var.data))
            @test maximum(abs, a-b) < 1e-16
            a = imag.(Ωf1.data)
            b = real.(imag.(Ωf1_var.data))
            @test maximum(abs, a-b) < 1e-16
        end
    end
end

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
            L, N = imex(VorticityEquation(n, Re, kforcing; T=Float64, dealias=dealias))

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
            L, N = imex(VorticityEquation(n, Re, kforcing; T=VarNum{Float64}, dealias=dealias))
            f = integrator(N, L, IMEXRKScheme(impl, FTField(n, Complex{VarNum{Float64}})), 0.005)

            # define initial condition, using the previous one
            Ω0 = Ω0_reference + 0*δ

            # perturb mode (1, 1), by the same amount as in the nonlinear 
            # simulation. The perturbation size does not matter, because 
            # the equations are linearised. We do this because is is easier 
            # to check the two. Note that here, we do not perturb the initial
            # condition of the nonlinear simulation. We set the initial condition
            # of the variational equations
            Ω0[1, 1] += δ

            # integrate in time
            Ωf1_var = f(copy(Ω0), 10)

            # check that amplitude of mode (1, 1) of perturbation is similar
            # to the difference of two perturbed nonlinear simulations. We cannot
            # set too low a tolerance here, because nonlinear effects will kick in
            # at some point reducing the validity of the linearised equations.
            # this checks that we are basically not doing rubbish.
            a = real(diff_nl[1, 1])
            b = pert(real(Ωf1_var[1, 1]))
            @test abs(a-b) < 1e-6

            a = imag(diff_nl[1, 1])
            b = pert(imag(Ωf1_var[1, 1]))
            @test abs(a-b) < 1e-6
        end
    end
end