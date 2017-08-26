# a useful macro for printing arrays
macro display(ex)
    quote
        display($(esc(ex)))
        println()
    end
end

include("test_allocation.jl")
include("test_broadcast.jl")
include("test_fft.jl")
include("test_field.jl")
include("test_flow.jl")
include("test_ftfield.jl")
include("test_indexing.jl")
include("test_norms.jl")
include("test_operators.jl")
include("test_spectra.jl")
include("test_system.jl")
# include("test_tangent.jl")