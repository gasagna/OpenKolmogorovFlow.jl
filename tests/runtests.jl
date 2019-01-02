using Test
using OpenKolmogorovFlow
import Random: seed!
import LinearAlgebra: dot, norm
using Flows

# include("test_adjoint_identity.jl")
include("test_allocation.jl")
include("test_broadcast.jl")
# include("test_distance.jl")
# include("test_fft.jl")        # FIXME
include("test_field.jl")
include("test_flow.jl")
include("test_ftfield.jl")
include("test_indexing.jl")
include("test_nonlinear.jl")
include("test_norms.jl")
include("test_operators.jl")
# include("test_shifts.jl")
# include("test_spectra.jl")
include("test_tangent.jl")
