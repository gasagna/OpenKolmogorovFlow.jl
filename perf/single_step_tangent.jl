using OpenKolmogorovFlow
using IMEXRKCB
using BenchmarkTools

# parameters
const Re = 40
const kforcing = 4
const n = 100
const m = even_dealias_size(n)

Lf, Nf = imex(VorticityEquation(n, Re, kforcing; mode=ForwardMode(), dealias=true))
Ω₀f = FTField(n); dΩdtf = similar(Ω₀f)
schemef = IMEXRKScheme(IMEXRK3R2R(IMEXRKCB3e, false), Ω₀f)
@btime step!($schemef, $Nf, $Lf, 0, 0.1, $Ω₀f)

Lt, Nt = imex(VorticityEquation(n, Re, kforcing; mode=TangentMode(), dealias=true))
Ω₀t = AugmentedFTField(n); dΩdtt = similar(Ω₀t)
schemet = IMEXRKScheme(IMEXRK3R2R(IMEXRKCB3e, false), Ω₀t)
@btime step!($schemet, $Nt, $Lt, 0, 0.1, $Ω₀t)