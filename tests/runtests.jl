# a complex exponential and its derivatives
fun(x, y, c::Number, j::Integer, k::Integer) = 
    0.5.*real(c.*exp.(im*(j.*x .+ k.*y)) .+ conj(c).*exp.(-im*(j.*x .+ k.*y)))

funx(x, y, c::Number, j::Integer, k::Integer) = 
    0.5.*real(im.*j.*c.*exp.(im*(j.*x .+ k.*y)) .- im.*j.*conj(c).*exp.(-im*(j.*x .+ k.*y)))

funy(x, y, c::Number, j::Integer, k::Integer) =     
    0.5.*real(im.*k.*c.*exp.(im*(j.*x .+ k.*y)) .- im.*k.*conj(c).*exp.(-im*(j.*x .+ k.*y)))        

funxx(x, y, c::Number, j::Integer, k::Integer) = -j^2*fun(x, y, c, j, k)
funyy(x, y, c::Number, j::Integer, k::Integer) = -k^2*fun(x, y, c, j, k)

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
include("test_system.jl")
include("test_flow.jl")