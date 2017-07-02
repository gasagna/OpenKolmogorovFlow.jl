# DEFINE STUFF NECESSARY TO USE . NOTATION
import Base.Broadcast: _containertype, 
                       promote_containertype, 
                       broadcast_c!, 
                       broadcast_c

# when you are broadcasting over FTField and DiffOperator always go down to FTField
_containertype(::Type{<:FTField}) = FTField
_containertype(::Type{<:DiffOperator}) = FTField

# These are need for scalars and arrays
promote_containertype(::Type{FTField}, ::Type{Any})   = FTField
promote_containertype(::Type{Any}, ::Type{FTField})   = FTField
promote_containertype(::Type{FTField}, ::Type{Array}) = FTField
promote_containertype(::Type{Array}, ::Type{FTField}) = FTField

# Extract underlying matrices in the fields and operators and then call broadcast 
# for the arrays. See this link
# https://discourse.julialang.org/t/allocation-using-broadcasting-with-custom-type/4592/3
# for the reason why we do not have a single broadcast function using splatting.
for nargs = 1:20
    args  = [Symbol("src", i) for i = 1:nargs]
    calls = [:(unsafe_get($(args[i]))) for i = 1:nargs]
    @eval @generated function broadcast_c!(f, ::Type{FTField}, ::Type{FTField}, dest, $(args...))
            :(broadcast!(f, unsafe_get(dest), $($calls...)); return dest)
          end    
end    

# Allocating version. This is called when you need to create a new object e.g.
# u = v .* 3. This creates a new object and returns it.
@generated function broadcast_c(f, ::Type{FTField}, src::Vararg{Any, N}) where N
    args = [:(unsafe_get(src[$k])) for k = 1:N]
    quote
        FTField(broadcast(f, $(args...)))
    end
end

# extract the underlying matrix representation 
@inline Base.unsafe_get(U::Union{FTField, DiffOperator}) = U.data
