# a complex exponential and its derivatives
fun(x, y, c::Number, α::Integer, β::Integer) = 
    0.5.*real(c.*exp.(im*(α.*x .+ β.*y)) .+ conj(c).*exp.(-im*(α.*x .+ β.*y)))

funx(x, y, c::Number, α::Integer, β::Integer) = 
    0.5.*real(im.*α.*c.*exp.(im*(α.*x .+ β.*y)) .- im.*α.*conj(c).*exp.(-im*(α.*x .+ β.*y)))

funy(x, y, c::Number, α::Integer, β::Integer) =     
    0.5.*real(im.*β.*c.*exp.(im*(α.*x .+ β.*y)) .- im.*β.*conj(c).*exp.(-im*(α.*x .+ β.*y)))        

funxx(x, y, c::Number, α::Integer, β::Integer) = -α^2*fun(x, y, c, α, β)
funyy(x, y, c::Number, α::Integer, β::Integer) = -β^2*fun(x, y, c, α, β)


include("test_broadcast.jl")
include("test_field.jl")
include("test_ftfield.jl")
include("test_indexing.jl")
include("test_operators.jl")
include("test_system.jl")