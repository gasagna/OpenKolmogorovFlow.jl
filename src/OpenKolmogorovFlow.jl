__precompile__(false)
module OpenKolmogorovFlow

using VariationalNumbers

include("operators.jl")
include("indexing.jl")
include("ftfield.jl")
include("field.jl")
include("fft.jl")
include("broadcast.jl")
include("system.jl")
include("flow.jl")
include("spectra.jl")

end