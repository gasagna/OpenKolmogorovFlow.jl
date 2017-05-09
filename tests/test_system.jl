using Base.Test
using OpenKolmogorovFlow
using IMEXRKCB

# Test that solution converges to the 
# laminar flow for small Reynolds numbers
@testset "laminar flow " begin
    # example dimension
    n = 10

    # take a Re and a forcing wavenumber
    Re = 1.2345678
    kforcing = 4

    # initial condition
    Ω₀ = FTField(n)
    
    # get explicit and implicit parts
    L = ImplicitTerm(n, Re)
    N = ExplicitTerm(n, kforcing)

    # for each integration scheme
    for scheme in [IMEXRK3R2R(IMEXRKCB3c, Ω₀, false),
                   IMEXRK3R2R(IMEXRKCB3e, Ω₀, false),
                   IMEXRK4R3R(IMEXRKCB4,  Ω₀, false)]
    
        # define T-time forward map
        f = forwmap!(N, L, 50, 0.005, scheme)

        # start from some non zero initial condition
        Ω₀ .= 0.01; Ω₀[0, 0] = 0

        # monitor the state excited by forcing
        m = Monitor(Ω₀, Ω->Ω)
    
        # map forward 
        f(Ω₀, m)
        
        # test final value is that predicted by explicit equation
        Δ = m.samples[end] .- laminarflow(n, Re, kforcing)
        @test maximum(abs, Δ) < 1e-15
    end
end