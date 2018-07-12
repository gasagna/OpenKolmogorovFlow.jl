module OpenKolmogorovFlow

include("indexing.jl")
# include("utils.jl")
include("ftfield.jl")
include("field.jl")
include("operators.jl")
include("fft.jl")
include("norms.jl")
include("implicitterm.jl")
include("forcings.jl")
include("system.jl")
include("tangent.jl")
include("hooks.jl")
include("flow.jl")
# include("adjoint.jl")
include("spectra.jl")
# include("shifts.jl")
# include("distance.jl")

end