using OpenKolmogorovFlow
using IMEXRKCB
using BenchmarkTools
using VariationalNumbers

# parameters
Re = 40
kforcing = 4

for n = [32, 64, 128, 256]

    # ~~~ RUN NONLINEAR EQUATIONS ~~~
    # get explicit and implicit parts
    L, N = imex(VorticityEquation(n, Re, kforcing; T=Float64, flags=FFTW.PATIENT, dealias=false))

    # initial condition
    Ω₀ = FTField(n)

    # define scheme
    scheme = IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4, false), Ω₀)

    # run step
    t_nl = @belapsed step!($scheme, $N, $L, 0, 0.1, $Ω₀)


    # ~~~ RUN NONLINEAR+VARIATIONAL EQUATIONS ~~~
    # get explicit and implicit parts
    L, N = imex(VorticityEquation(n, Re, kforcing; T=VarNum{Float64}, flags=FFTW.PATIENT, dealias=false))

    # initial condition
    Ω₀ = FTField(n) + δ

    # define scheme
    scheme = IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4, false), Ω₀)

    # run step
    t_var = @belapsed step!($scheme, $N, $L, 0, 0.1, $Ω₀)

    # print output: nonlinear simulation time, nonlinear+variational time and slowdown
    @printf "%4d %10.2f ms %10.2f ms %10.2f\n" n 10^3*t_nl 10^3*t_var t_var/t_nl
end

# ~~~ COMMENT ~~~ 
# It seems like if running the variational equations + nonlinear equations would 
# cost 3.5 times more than running the nonlinear equations alone, depending on the
# grid size. It seems a bit too large, because:
# - FFTs should be twice as expensive
# - all linear operations should be twice as expensive
# - the evaluation of the nonlinear term should be three times more expensive
# Hence, on average running the var. equations + the nonlinear simulation 
# should be between two and three times more expensive, depending on how
# much importance is the evaluation of the nonlinear term. The larger than
# expected slowdown is probably due to the lesser performance of FFTW on 
# data with split arrays, since the FFTW take a larger fraction of the 
# total running cost. This would need profiling and improvements. We should
# also run profiling of a C implementation of the split array format.
# However, the code has not changed a lot and the logic for solving the 
# variational equations is quite simple.