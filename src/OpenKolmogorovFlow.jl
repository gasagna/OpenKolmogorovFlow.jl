__precompile__(false)
module OpenKolmogorovFlow

include("indexing.jl")
include("utils.jl")
include("ftfield.jl")
include("field.jl")
include("augmented.jl")
include("operators.jl")
include("fft.jl")
include("norms.jl")
include("linear.jl")
include("system.jl")
include("adjoint.jl")
include("flow.jl")
include("spectra.jl")
include("shifts.jl")
include("distance.jl")

end