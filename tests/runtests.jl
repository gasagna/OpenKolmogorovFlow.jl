# an useful macro for printing arrays
macro display(ex)
    quote
        display($(esc(ex)))
        println()
    end
end

include("test_broadcast.jl")
include("test_field.jl")
include("test_ftfield.jl")
include("test_indexing.jl")
include("test_operators.jl")
# include("test_system.jl")
# include("test_flow.jl")