__precompile__(false)
module OpenKolmogorovFlow

using VariationalNumbers

include("operators.jl")
include("indexing.jl")
include("utils.jl")
include("ftfield.jl")
include("field.jl")
include("fft.jl")
include("norms.jl")
include("broadcast.jl")
include("system.jl")
include("flow.jl")
include("spectra.jl")
include("shifts.jl")
include("distance.jl")

end