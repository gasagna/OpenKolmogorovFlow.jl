using OpenKolmogorovFlow
using IMEXRKCB
using BenchmarkTools

function foo()
    # parameters
    Re = 40
    kforcing = 4
    n = 132

    # get explicit and implicit parts
    Ld, Nd = imex(VorticityEquation(n, Re, kforcing, Float64, FFTW.MEASURE, true))

    # initial condition
    Ω₀ = FTField(n)

    # define scheme
    scheme = IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4, false), Ω₀)

    # run step
    for i = 1:100
        step!(scheme, Nd, Ld, 0, 0.1, Ω₀)
    end
end

foo()
@profile foo()
Profile.print()